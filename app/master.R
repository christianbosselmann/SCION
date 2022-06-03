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

# get helper fn
source("app/func.R")

# get preprocessed training data set
train <- read_csv("app/training_data.csv") # aka "data_prep.csv"
y <- as.factor(train$y)

# get similarity matrix
sim_matrices <- read_csv("app/similaritymatrix.csv")
rownames(sim_matrices) <- colnames(sim_matrices)

sim_match <- sim_matrices %>%
  rownames_to_column() %>%
  pivot_longer(cols = -c(1)) %>%
  rename(k = rowname, l = name)

# get pre-processing recipe
load("app/recipe.rds")

# get feature lookup tables
aa_feats <- read_tsv("app/aa_feats.tsv")
load("app/str_feats.rda")
cid_raw <- read_csv("app/cid.csv")

# get test data, collect data features, then merge in features
df_in <- openxlsx::read.xlsx("app/input.xlsx") # DEBUG
n_test <- nrow(df_in)
n_train <- nrow(train)
n_all <- n_test + n_train

test <- df_in %>%
  merge(aa_feats, by = c("aa1", "aa2")) %>%
  merge(str_feats, by = c("gene", "pos"))

# split, pivot and merge canonical alignment id (cid) and annotation features
cid_align <- cid_raw[,1:10] # split into domain annotation and actual alignment
cid_feats <- cid_raw[,10:15] # sorry, this is currently a magic number

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

Kb <- as.data.frame(
  kernelMatrix(
    rbfdot(sigma = 0.01),
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

# get set of hpo terms for each training observation
load("app/set_terms.RData")

# load ontology
ont_hpo <- get_ontology("app/hp.obo.txt", 
                        propagate_relationships = "is_a", 
                        extract_tags = "minimal")

### TODO pheno kernel deploy

# # propagate terms to parents
# list_prop <- lapply(set_terms, propagate_relations, ontology = ont_hpo, relations = "parents")
# 
# # reshape into a wide dataframe
# df_prop <- list_prop %>% 
#   map_dfr(~ .x %>% as_tibble(), .id = "name")
# 
# df_prop <- reshape2::dcast(df_prop, as.numeric(name) ~ value, length)
# 
# # remove id column, redundant to rowname
# df_prop <- df_prop[,-1]
# 
# # make sure that the matrix is binary (artifact from propagation)
# df_prop[df_prop > 0] <- 1
# 
# # calculate pairwise jaccard coefficient
# mat_pheno <- jaccard.test.pairwise(as.matrix(df_prop), compute.qvalue = FALSE)
# 
# # force to a square symmetric matrix and set diagonal to 1 (identity matrix)
# mat_pheno <- forceSymmetric(mat_pheno$statistics)
# diag(mat_pheno) <- 1
# 
# # fix class, and we are done!
# mat_pheno <- as.matrix(mat_pheno)
# 
# 
# 
# 
# 
# 
# # set up class weights
# weights <- c("GOF" = length(train$y) / sum(train$y == "GOF"),
#              "LOF" = length(train$y) / sum(train$y == "LOF"), 
#              "Neutral" = length(train$y) / sum(train$y == "Neutral"))
# 
# # set up train and test indices
# train_indices <- 1:length(t_train)
# test_indices <- (length(t_train)+1):(length(t_train)+length(t_test))
# 
# # training SVM
# model <- e1071::svm(
#   x = Km[c(-test_indices), c(-test_indices), drop = FALSE],
#   y = y,
#   cost = 1,
#   probability = TRUE,
#   class.weights = weights
# )
# 
# # predicting test set
# test <- as.kernelMatrix(Km[test_indices, -test_indices, drop = FALSE])
# prediction <- predict(model, test) %>%
#   as_tibble() %>%
#   cbind(attr(predict(model, test, probability = TRUE), "probabilities")) %>%
#   cbind(attr(predict(model, test, decision.value = T), "decision.values")) %>%
#   rename(prediction = value, 
#          prob.LOF = LOF, prob.Neutral = Neutral, prob.GOF = GOF, 
#          `conf.LOF_Neutral` = `LOF/Neutral`, `conf.LOF_GOF` = `LOF/GOF`, `conf.Neutral_GOF` = `Neutral/GOF`)
# 
# # print output
# out <- cbind(input, prediction)
# write.xlsx(out, "output.xlsx")
# 
# # verbose output for app
# verb_out <- paste("This variant is predicted to be:", out$prediction, sep = " ")
