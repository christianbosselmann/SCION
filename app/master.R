#' SCN
#' Functional variant prediction for voltage-gated sodium channels

# packages
library("librarian")
librarian::shelf(tidyverse,
                 tidymodels,
                 data.table,
                 openxlsx,
                 yardstick,
                 caret,
                 bestNormalize,
                 kernlab,
                 e1071,
                 data.table,
                 magic,
                 matrixcalc,
                 Matrix,
                 ontologyIndex,
                 ontologySimilarity,
                 klic,
                 CEGO,
                 jaccard,
                 quiet = TRUE)

# set seed
set.seed(42)

# get helper fn
source("func.R")

# get preprocessed training data set
train <- read_csv("training_data.csv") # aka "data_prep.csv"
y <- as.factor(train$y)

# get similarity matrix
sim_matrices <- read_csv("similaritymatrix.csv")
rownames(sim_matrices) <- colnames(sim_matrices)

sim_match <- sim_matrices %>%
  rownames_to_column() %>%
  pivot_longer(cols = -c(1)) %>%
  rename(k = rowname, l = name)

# get pre-processing recipe
load("recipe.rds")

# get feature lookup tables
aa_feats <- read_tsv("aa_feats.tsv")
load("str_feats.rda")
cid_raw <- read_csv("cid.csv")

# collect data features of test observations
n_test <- nrow(df_in)
n_train <- nrow(train)
n_all <- n_test + n_train

# fix aa label of test observations
df_in$aa1 <- str_split(df_in$aa1, " ", 2) %>%
  lapply(., function(x) x[1]) %>%
  unlist()

df_in$aa2 <- str_split(df_in$aa2, " ", 2) %>%
  lapply(., function(x) x[1]) %>%
  unlist()

# merge in features
test <- df_in %>%
  merge(aa_feats, by = c("aa1", "aa2")) %>%
  merge(str_feats, by = c("gene", "pos"))

# split, pivot and merge canonical alignment id (cid) and annotation features
cid_align <- cid_raw[,1:10] # split into domain annotation and actual alignment
cid_feats <- cid_raw[,10:15] # this is currently a magic number

cid_align <- cid_align %>%
  pivot_longer(!cid, names_to = "gene", values_to = "pos") %>%
  group_by(gene) %>%
  arrange(pos) %>%
  distinct(pos, .keep_all = TRUE)

test <- merge(test, cid_align)

test <- test %>%
  merge(cid_feats, by = c("cid"))

# drop unneccessary features
test <- test %>% 
  select(-one_of("pos", 
                 "aa1", 
                 "aa2"))

# fix data types
test <- hablar::retype(test)

# apply prep recipe
test <- bake(dat_rec, test)

# generate instance-level kernel matrix for sequence/structure features (base/union)
t_train <- train$gene %>%
  as.character()
train <- train %>%
  select(-c("y", "gene"))

t_test <- test$gene %>%
  as.character()
test <- test %>%
  select(-gene)

Kb <- as.matrix(
  kernlab::kernelMatrix(
    kernlab::rbfdot(sigma = 0.01),
    as.matrix(rbind(train, test))
  ))

# calculate task similarity kernel matrix Kt
Kt <- matrix(data = 0, nrow = n_all, ncol = n_all,
             dimnames = list(c(t_train, t_test), c(t_train, t_test)))

for (i in 1:n_all) {
  for (j in 1:n_all) { 
    Kt[i,j] <- sim_match$value[sim_match$k == c(t_train, t_test)[[i]] & 
                                 sim_match$l == c(t_train, t_test)[[j]]] 
  }
}

# calculate MTL kernel Km
Km <- as.matrix(Kb) * as.matrix(Kt)

# load phenotypic data from user interface
test_terms <- str_split(df_hpo, " ", 2) %>%
  lapply(., function(x) x[1]) %>%
  unlist()

# assert that we extracted valid non-null HPO terms before proceeding
if(all(test_terms %in% ont_hpo$id) == FALSE | is.null(test_terms) == TRUE){
  flag_mkl <- FALSE # proceed with MTL
}else{
  flag_mkl <- TRUE # proceed with MTMKL
}

if(flag_mkl == TRUE) {
  
  # get set of hpo terms for each training observation
  load("set_terms.RData")
  
  # append test terms
  set_terms <- append(set_terms, list(test_terms))

  # propagate terms to parents
  list_prop <- lapply(set_terms, propagate_relations, ontology = ont_hpo, relations = "parents")

  # reshape into a wide dataframe
  df_prop <- list_prop %>%
    map_dfr(~ .x %>% as_tibble(), .id = "name")

  df_prop <- reshape2::dcast(df_prop, as.numeric(name) ~ value, length)

  # remove id column, redundant to rowname
  df_prop <- df_prop[,-1]

  # make sure that the matrix is binary (artifact from propagation)
  df_prop[df_prop > 0] <- 1

  # calculate pairwise jaccard coefficient
  mat_pheno <- jaccard.test.pairwise(as.matrix(df_prop), compute.qvalue = FALSE)

  # force to a square symmetric matrix and set diagonal to 1 (identity matrix)
  mat_pheno <- forceSymmetric(mat_pheno$statistics)
  diag(mat_pheno) <- 1

  # fix type
  mat_pheno <- as.matrix(mat_pheno)
  
  # scale and center kernel matrices
  Kb <- kernelPreparation(Kb) # instance-level sequence and structure
  Km <- kernelPreparation(mat_pheno) # instance-level phenotype
  Kt <- kernelPreparation(Kt) # instance-level task similarity
  Km <- Kb+Km+Kt
  
  if(flag_exp == TRUE){
    print("this is where the fun begins!")
  }
}

# set up train and test indices
train_indices <- 1:length(t_train)
test_indices <- (length(t_train)+1):(length(t_train)+length(t_test))

# training SVM
model <- e1071::svm(
  x = Km[c(-test_indices), c(-test_indices), drop = FALSE],
  y = y,
  cost = 4,
  probability = TRUE
)

# predicting test set
test <- as.kernelMatrix(Km[test_indices, -test_indices, drop = FALSE])
prediction <- predict(model, test) %>%
  as_tibble() %>%
  cbind(attr(predict(model, test, probability = TRUE), "probabilities")) %>%
  cbind(attr(predict(model, test, decision.value = TRUE), "decision.values")) 

# print output
out <- cbind(df_in, prediction)
# write.xlsx(out, "output.xlsx") # DEBUG

# verbose output for app
verb_name <- paste(df_in$gene, paste("p.", paste(df_in$aa1, df_in$pos, df_in$aa2, sep = ""), sep = ""), sep = " ")
verb_out <- paste(verb_name, "is predicted to be:", out$value, sep = " ")
