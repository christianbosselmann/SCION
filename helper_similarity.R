# SCN
# Functional variant prediction for voltage-gated sodium channels
# helper_similarity.R
# this script leverages taxonomy-based information on task-relatedness by 
# using the msa-distance matrix to obtaining similarity matrices for the MTL kernel
# one similarity matrix for each value of alpha (see below)
# cf. DOI 10.1007/978-3-642-12683-3_34
# requires the distance matrix provided by pipeline_msa.R

# set desired values of alpha
alpha <- c(1, 2, 3, 5, 10, 100) # baseline similarity

# packages
library("librarian")
librarian::shelf(tidyverse,
                 data.table,
                 openxlsx,
                 Matrix,
                 matrixcalc,
                 ggcorrplot,
                 quiet = TRUE)

# set seed
set.seed(42)

# read data
d <- read_csv("mat/distancematrix.csv")
rownames(d) <- colnames(d)

# remember that dist.alignment = sqrt(1 - identity). Let's get the pairwise distance.
d = d^2

# convert to similarity (gamma) by the transformation gamma = alpha − d / dmax
# where alpha ≥ 1 is a hyperparameter to control baseline similarity between tasks
# and dmax is the maximal distance. alpha is set in the master script.
for (a in alpha) {
  gamma <- a - (d / max(d))
  gamma <- as.matrix(gamma)
  
  # matrices need to be positive semi-definite to yield a valid kernel. 
  if (is.positive.semi.definite(as.matrix(gamma)) == FALSE) {
    gamma <- nearPD(as.matrix(gamma), 
                    keepDiag = TRUE)
  }
  
  # fix names
  gamma <- data.frame(as.matrix(gamma$mat))
  rownames(gamma) <- rownames(d)
  colnames(gamma) <- colnames(d)
  
  # final check
  if (is.positive.semi.definite(as.matrix(gamma)) == FALSE) {
    stop("Error, similarity matrix is not PDM. Kernel will be invalid.")
  }
  
  # export
  write_csv(gamma, paste('mat/similaritymatrix_a', a, '.csv', sep = ''))
  
  # assign alpha value as matrix id
  assign(paste("gamma", a, sep = "_a"), gamma)   
}

# visualization
p <- ggcorrplot(cor(gamma_a1))
plot(p)

tiff("fig/similarity.tiff", units="in", width=8, height=8, res=300)
plot(p)
dev.off()
