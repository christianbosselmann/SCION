# SCN
# Functional variant prediction for voltage-gated sodium channels
# helper_kernel.R
# this script generates MTL kernel matrices (Km) for later use in MTL-SVM.
# @params sigma for RBF kernel matrix
# @params alpha is assumed from the available list of similarity matrices
# @return Kt task-level similarity matrix
# @return Kf instance-level RBF kernel matrix
# @return Km MTL kernel matrix (Kt*Kf)

# packages
library("librarian")
librarian::shelf(tidyverse,
                 data.table,
                 Matrix,
                 matrixcalc,
                 kernlab,
                 quiet = TRUE)

# set magic values
sigma <- c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10)

# get data
data <- read_csv("data/dat_prep.csv")

# get similarity matrices, maintain order of alpha vector, fix names
path <- list.files(path = "mat", pattern = "^similaritymatrix_", full.names = TRUE)
path <- str_sort(path, numeric = TRUE)
names <- sub("^[^_]*_", "", gsub(path, pattern = ".csv$", replacement = ""))

sim_matrices <- lapply(path, read_csv) 
names(sim_matrices) <- names

# fixing rownames
for (m in 1:length(sim_matrices)) {
  rownames(sim_matrices[[m]]) <- colnames(sim_matrices[[m]])
}

# reshape to long format, yielding the pairwise similarity for all combinations of tasks
# loop through all values of alpha, the baseline similarity, then store as list
sim_match <- vector(mode = "list", length = length(path))

for (a in 1:length(path)) {
  sim_match[[a]] <- sim_matrices[[a]] %>%
    rownames_to_column() %>%
    pivot_longer(cols = -c(1)) %>%
    dplyr::rename(k = rowname, l = name)
}

names(sim_match) <- names(sim_matrices)

# get vector of task identity. store for later, used for task-specific model assessment later on..
t <- data$gene
save(t, file = "t_vec.rda")

### Kf
# constructing Kf, the instance-level kernel matrix (i.e. similarity between observations)
f <- data %>%
  select(-y) %>%
  select(-gene) %>%
  as.matrix()

# prepare list to apply kernelMatrix on
f <- rep(list(f), length(sigma)) 

# iterate over the list and get the rbf kernel matrix for each value of sigma
Kf <- lapply(1:length(f),
             function(x){
               as.matrix(kernelMatrix(rbfdot(sigma = sigma[[x]]), x = f[[x]]))
             })

# export
saveRDS(Kf, "mat/kernelmatrices_instance.rds")

### Kt
# preallocating Kt, the list of task similarity matrices
Kt <- lapply(1:a, matrix, data = NA, nrow = length(t), ncol = length(t))

# lapply through the list of empty matrices and assign values by col/rownames and index from appropriate sim_match list (inner for loop)
Kt <- lapply(seq_along(Kt),
             function(index){
               m <- Kt[[index]]
               colnames(m) <- t
               rownames(m) <- t
               
               d <- sim_match[[index]] 
               
               for(i in 1:nrow(d)){
                 m[rownames(m) == d[[i,2]], colnames(m) == d[[i,1]]] <- d[[i,3]]
               }
               
               colnames(m) <- NULL
               rownames(m) <- NULL
               Kt[[index]] <- m
             })

# export (for MKL later)
saveRDS(Kt[[1]], "mat/taskmatrix.rds")

### Km
# preallocate an empty nested list for each alpha and sigma
Km <- lapply(Km <- vector(mode = 'list', length(sim_match)),
             function(x) x <- vector(mode = 'list', length(sigma)))

# element-wise multiplication of Kf and Kt to obtain the list of MTL kernel matrices Km
Km <- lapply(1:length(Kt), function(a)
  lapply(1:length(Kf), function(s){
    Kf[[s]] * Kt[[a]]
  }))

# give the kernel matrix list a descriptive name
Km <- setNames(lapply(Km, setNames, paste("s", sigma, sep = "")), names)

# save Km (nested list of kernel matrices by alpha*sigma)
saveRDS(Km, "mat/kernelmatrices_mtl.rds")

### Kd
# Dirac kernel, vgl. DOI 10.1093/bioinformatics/btm611
# intuitively the same as the 'single' method from prefeKt, where one SVM
# is trained for each task
Kd <- Kt

Kd <- lapply(seq_along(Kd),
             function(index){
               m <- Kd[[index]]
               colnames(m) <- t
               rownames(m) <- t
               
               # get outer product to check if rowname == colname for each matrix element
               out <- outer(row.names(m), colnames(m), `==`)
               dimnames(out) <- dimnames(m)
               m <- out*1
               
               colnames(m) <- NULL
               rownames(m) <- NULL
               Kt[[index]] <- m
             })

Kd <- lapply(1:length(Kd), function(a)
  lapply(1:length(Kf), function(s){
    Kf[[s]] * Kd[[a]]
  }))

saveRDS(Kd, "mat/kernelmatrices_dirac.rds")

### Ku
# uniform kernel, where each task is weighted the same
# intuitively the same as the 'base' method from prefeKt, where one SVM
# is trained for all tasks
Ku <- lapply(1:a, matrix, data = 1, nrow = length(t), ncol = length(t))

Ku <- lapply(1:length(Ku), function(a)
  lapply(1:length(Kf), function(s){
    Kf[[s]] * Ku[[a]]
  }))

saveRDS(Ku, "mat/kernelmatrices_union.rds")
