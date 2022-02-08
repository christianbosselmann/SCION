# SCN
# phenotype pre-processing

# load packages
library(librarian)
librarian::shelf(tidyverse,
                 data.table,
                 httr, 
                 jsonlite,
                 matrixcalc,
                 Matrix,
                 ontologyIndex,
                 ontologySimilarity,
                 klic,
                 CEGO,
                 jaccard)

# params
sim_method <- "jaccard" # choose similarity measure: jaccard, lin or resnik
psd_method <- "clip" # choose spectrum method to find nearest psd matrix: clip, shift or flip

# get data
raw_data <- read_csv("data/clean_tbl.csv") %>%
  select(id, pheno)

# inner join disease list (manual) to get OMIM IDs
list_disease <- read_delim("pheno/disease_list.csv", delim = ";")

raw_data$omim <- str_replace_all(raw_data$pheno, setNames(list_disease$omim, list_disease$term))

# remove # (for API call later)
raw_data$omim <- str_replace_all(raw_data$omim, "#", "")

# handle the case of multiple OMIM IDs per variant by separating across rows
raw_data <- raw_data %>%
  separate_rows(omim, sep = "[^[:alnum:].]+", convert = TRUE)

# define API request
base <- "https://hpo.jax.org/api/hpo/disease/"

# loop through all OMIM IDs and get HPO terms from HPO API
list_terms <- list()

for (i in 1:nrow(raw_data)) {
  disease <- paste("OMIM:", raw_data$omim[i], sep = "") # diseaseID
  
  # call
  get_terms <- GET(paste(base, disease, sep = ""))
  
  # get content
  content <- as.data.frame(
    fromJSON(
      content(get_terms, "text"),
      flatten = TRUE)
  )
  
  # extract HPO terms from disease association and store as list by rowindex
  term_vector <- rbindlist(content$catTermsMap.terms)$ontologyId
  list_terms[[i]] <- term_vector
}

# map HPO terms back from rowindex to variant id
set_terms <- list() # set_terms is a list of HPO terms by variant ID

for (i in raw_data$id) {
  ind <- which(raw_data$id == i)
  set_terms[[i]] <- as_vector(do.call(c, list(list_terms[ind])))
}

# remove null elements from list (which are introduced by indexing)
set_terms <- set_terms[lengths(set_terms) != 0]

# load ontology
hpo <- get_ontology("pheno/hp.obo.txt", propagate_relationships = "is_a", extract_tags = "minimal")

# insert for Jaccard similarity coefficient
if (sim_method == "jaccard"){
  # propagate terms to parents
  list_prop <- lapply(set_terms, propagate_relations, ontology = hpo, relations = "parents")
  
  # reshape into a wide dataframe
  df_prop <- list_prop %>% 
    map_dfr(~ .x %>% as_tibble(), .id = "name")
  
  df_prop <- reshape2::dcast(df_prop, as.numeric(name) ~ value)
  
  # remove id column, redundant to rowname
  df_prop <- df_prop[,-1]
  
  # make sure that the matrix is binary (artifact from propagation)
  df_prop[df_prop > 0] <- 1
  
  # calculate pairwise jaccard coefficient
  mat_pheno <- jaccard.test.pairwise(as.matrix(df_prop), compute.qvalue = FALSE)
  
  # force to a square symmetric matrix and set diagonal to 1 (identity matrix)
  mat_pheno <- forceSymmetric(mat_pheno$statistics)
  diag(mat_pheno) <- 1
  
  # fix class, and we are done!
  mat_pheno <- as.matrix(mat_pheno)
  
  # if we are not doing a Jaccard coefficient, carry on:
} else {
  
  # get similarity matrix
  mat_pheno <- get_sim_grid(ontology = hpo, 
                            information_content = descendants_IC(hpo),
                            term_sim_method = sim_method, # or resnik
                            term_sets = set_terms)
  
  # check if phenotypic similarity matrix, mat_pheno, is positive semi-definite
  # which is generally not the case for Lin and Resnik measures
  # this correction has some implications, cf https://authors.library.caltech.edu/69864/1/p145-chen.pdf
  if (is.positive.semi.definite(mat_pheno) == FALSE) {
    print("Phenotypic similarity matrix is not PSD, applying correction.")
    
    if (psd_method == "clip") {
      mat_pheno <- nearPD(mat_pheno, base.matrix = TRUE)
      mat_pheno <- mat_pheno$mat
    }
    
    if (psd_method == "shift") {mat_pheno <- spectrumShift(mat_pheno, coeff = 1.2)}
    
    if (psd_method == "flip") {mat_pheno <- correctionKernelMatrix(mat_pheno, method = "flip", repair = TRUE, tol = 1e-08)
    mat_pheno <- mat_pheno$mat
    }
  }
}

# export for use in kernel matrix precomputation and MTMKL model
write.table(mat_pheno, "mat/hpomatrix.csv")
