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

# set seed
set.seed(42)

# get helper functions
source("func.R")

#' @params sim_method choose similarity measure: jaccard, lin, or resnik
#' @params psd_method choose spectrum method to find nearest psd matrix: clip, shift or flip
#' @params pheno_sim if TRUE, sparse and noisy phenotypes are simulated for each similarity method
#' @return phenotypic similarity matrix as R object
sim_method <- "jaccard" 
psd_method <- "shift"
pheno_sim <- FALSE

# get set of hpo terms for each OMIM ID (helper_omim.R)
load("pheno/set_terms.RData")

# load ontology
ont_hpo <- get_ontology("pheno/hp.obo.txt", 
                        propagate_relationships = "is_a", 
                        extract_tags = "minimal")

# sparse/noisy phenotype simulation
if (pheno_sim == TRUE) {
term_noise <- 2 ^ seq(-2, 2, by = 1) # vector of noise to introduce to set terms
term_len <- unlist(lapply(set_terms, length)) # vector of number of terms for each variant
term_list <- setNames(vector('list', length(term_noise)), paste0('P', term_noise)) # named list to store results as list of lists
term_bag <- ont_hpo$id # vector of all possible hpo terms

# iterate through noise levels and pick appropriate method to remove or add terms
for (i in 1:length(term_noise)) {
  
  # sparse phenotype
  if (term_noise[[i]] < 1){
    term_list[[i]] <- list()
    for (j in 1:length(set_terms)){
      term_list[[i]][[j]] <- sample(set_terms[[j]], term_len[[j]] * term_noise[[i]])
    }
  } 
  
  # standard phenotype
  if (term_noise[[i]] == 1){
    term_list[[i]] <- set_terms
  }
  
  # noisy phenotype
  if(term_noise[[i]] > 1) {
    for (j in 1:length(set_terms)){
      term_target <- term_len[[j]] * term_noise[[i]] - term_len[[j]]  # how many hpo terms should be added to this term set
      term_tmp <- unname(sample(term_bag, term_target)) # grab the terms from term_bag
      term_list[[i]][[j]] <- append(term_tmp, set_terms[[j]])
    }
  }
}
}

# Jaccard similarity coefficient
if (sim_method == "jaccard"){
  # propagate terms to parents
  list_prop <- lapply(set_terms, propagate_relations, ontology = ont_hpo, relations = "parents")
  
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
  mat_pheno <- get_sim_grid(ontology = ont_hpo, 
                            information_content = descendants_IC(ont_hpo),
                            term_sim_method = sim_method,
                            term_sets = set_terms)
  
  # check if phenotypic similarity matrix, mat_pheno, is positive semi-definite
  # which is generally not the case for Lin and Resnik measures
  # this correction has some implications, cf https://authors.library.caltech.edu/69864/1/p145-chen.pdf
  if (is.positive.semi.definite(mat_pheno) == FALSE) {
    print("Phenotypic similarity matrix is not PSD, applying correction.")
    
    if (psd_method == "clip") { # seems to perform worse than the other methods
      mat_pheno <- nearPD(mat_pheno, 
                          base.matrix = TRUE, 
                          doSym = TRUE,
                          ensureSymmetry = TRUE,
                          corr = FALSE) # this object also stores the minimum eigenvalue (for correction diagnostics)
      mat_pheno <- mat_pheno$mat
      mat_pheno <- round(mat_pheno, 10) # this is necessary to avoid floating point errors with isSymmetric and is.positive.semi.definite
    }
    
    if (psd_method == "shift") {
      mat_pheno <- spectrumShift(mat_pheno, coeff = 1.2)
      mat_pheno <- round(mat_pheno, 10)
    }
    
    if (psd_method == "flip") {
      mat_pheno <- correctionKernelMatrix(mat_pheno, method = "flip", repair = TRUE, tol = 1e-08)
      mat_pheno <- mat_pheno$mat
      mat_pheno <- round(mat_pheno, 10)
    }
  }
}

# export for use in kernel matrix precomputation and MTMKL model
hpo <- mat_pheno
save(hpo, file = "mat/hpomatrix.RData")
