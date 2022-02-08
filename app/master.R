# PREFEKT
# PREdicting the Functional Effects of Kv muTations

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
                 quiet = TRUE)

# set parameters
seed <- 42 # random seed
cost_vec <- 1 # C-SVM hyperparameter
sigma_vec <- 0.01 # rbf kernel hyperparameter
weights <- "uniform" # SVM weights for imbalanced classes. "uniform" for uniform weight distribution, "inverse" for weights inversely proportional to class frequency.

# get preprocessed training data set
train <- read_csv("training_data.csv")
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

# read feature tables
aa_feats <- read_tsv("features/aa_feats.tsv")
load("features/str_feats.rda")
cid_raw <- read_csv("features/cid.csv")

# get test data
# input <- read.xlsx("input.xlsx")
input <- df_in # for app

# merge aa
test <- input %>%
  merge(aa_feats, by = c("aa1", "aa2")) 

# unlist and merge str
str_feats <- do.call(rbind, Map(cbind, str_feats, id = names(str_feats))) %>%
  dplyr::rename(gene = "id")

test <- test %>%
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

# drop redundant features
test <- test %>% 
  select(-one_of("pos", 
                 "aa1", 
                 "aa2"))

# fix data types
test <- hablar::retype(test)

# apply prep recipe
test <- bake(dat_rec, test)

# calculate rbf kernel matrix Kb
df_train <- train %>%
  select(-y) %>%
  select(-gene)

t_train <- train %>%
  select(gene) %>%
  as_vector() # training set task vector

df_test <- test %>%
  select(-gene)

t_test <- test %>%
  select(gene) %>%
  mutate_if(is.factor, as.character) %>% # this is necessary to prevent a bug related to factor output
  as_vector() # test set task vector

df_all <- rbind(data.frame(id = "train", df_train),
                data.frame(id = "test", df_test))

t_all <- c(t_train, t_test)

Kb <- as.data.frame(
  kernelMatrix(
    rbfdot(sigma = sigma_vec), 
    as.matrix(df_all[,-1])
  ))

# calculate task similarity kernel matrix Kt
Kt <- matrix(data = 0, nrow = length(t_all), ncol = length(t_all),
             dimnames = list(as.vector(t_all), as.vector(t_all)))

for(i in 1:nrow(sim_match)){
  Kt[rownames(Kt) == sim_match[[i,2]], colnames(Kt) == sim_match[[i,1]]] <- sim_match[[i,3]]
}

# calculate MTL kernel Km
Km <- as.matrix(Kb) * as.matrix(Kt)

# set up class weights
if (weights == "uniform") {
  weights <- c("GOF" = 1,
               "LOF" = 1)
}

if (weights == "inverse") {
  weights <- c("GOF" = length(data$y) / sum(data$y == "GOF"),
               "LOF" = length(data$y) / sum(data$y == "LOF"))
}

# set up train and test indices
train_indices <- 1:length(t_train)
test_indices <- (length(t_train)+1):(length(t_train)+length(t_test))

# training SVM
model <- e1071::svm(
  x = Km[c(-test_indices), c(-test_indices), drop = FALSE],
  y = y,
  cost = cost_vec,
  probability = TRUE,
  class.weights = weights
)

# predicting test set
test <- as.kernelMatrix(Km[test_indices, -test_indices, drop = FALSE])
prediction <- predict(model, test) %>%
  as_tibble() %>%
  cbind(attr(predict(model, test, probability = TRUE), "probabilities")) %>%
  cbind(attr(predict(model, test, decision.value = T), "decision.values"))

# print output
out <- cbind(input, prediction)
# write.xlsx(out, "output.xlsx")

# verbose output for app
verb_out1 <- paste("This variant is predicted to be:", out$value, sep = " ")
verb_out2 <- paste("Probability of gain-of-function:", signif(out$GOF,3), sep = " ")
verb_out3 <- paste("Probability of loss-of-function:", signif(out$LOF,3), sep = " ")
verb_out4 <- paste("Confidence (distance from hyperplane):", signif(out$`GOF/LOF`,3), sep = " ")
