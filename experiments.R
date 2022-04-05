#' SCN
#' Functional variant prediction for voltage-gated sodium channels
#' setup for experiment control

#' @params seed random seed
#' @params k inner and outer folds for nested cross validation
#' @params cost_vec vector of cost values to tune over
#' @params class_weight SVM weights for imbalanced classes: "uniform" for equal weight distribution, "inverse" for weights inversely proportional to class frequency
#' @params kernel choice of kernel method: "mkl" for MTMKL-SVM, "mtl" for MTL-SVM, "dirac" for Dirac kernel SVM, "union" for union SVM
#' @params mkl_method choice of MKL method for kernel weights, "semkl" for SEMKL, "simple" for simpleMKL, "uniform" for no kernel weights, "group" for unregularized task-level weights
#' @params mkl_cost penalty for MKL kernel prioritization (only applies to SEMKL, SimpleMKL)
#' @params sim_method choose similarity measure: jaccard, lin, or resnik
#' @params psd_method choose spectrum method to find nearest psd matrix: clip, shift or flip
#' @params pheno_sim if TRUE, sparse and noisy phenotypes are simulated for each similarity method
#' @return phenotypic similarity matrix as R object
#' @return saved timestamped objects of parameters, metrics and raw predictions

#' experiment 1: standard model tuning and assessment
seed <- 42
k <- 10
cost_vec <- 2 ^ seq(-5, 5, by = 1)
class_weight <- "uniform"
kernel <- "mkl"
mkl_method <- "uniform"
mkl_cost <- 1
source("model.R")

#' experiment 2: compare different kernel methods
seed <- 42
k <- 10
cost_vec <- 2 ^ seq(-5, 5, by = 1)
class_weight <- "uniform"
mkl_method <- "uniform"
mkl_cost <- 1
iter <- c("mkl", "mtl", "union", "dirac")
for (i in iter) {
  kernel <- i
  source("model.R")
}

#' experiment 3: compare different mkl methods
seed <- 42
k <- 10
cost_vec <- 2 ^ seq(-5, 5, by = 1)
class_weight <- "uniform"
kernel <- "mkl"
mkl_cost <- 1
iter <- c("semkl", "simple", "uniform")
for (i in iter) {
  mkl_method <- i
  source("model.R")
}

#' experiment 4: compare different similarity measures
seed <- 42
k <- 10
cost_vec <- 2 ^ seq(-5, 5, by = 1)
class_weight <- "uniform"
kernel <- "mkl"
mkl_method <- "uniform"
mkl_cost <- 1
psd_method <- "shift"
pheno_sim <- FALSE
iter <- c("jaccard", "lin", "resnik")
for (i in iter) {
  sim_method <- i
  source("helper_pheno.R")
  source("model.R")
}

#' experiment 5: compare differet psd corrections across similarity measures
seed <- 42
k <- 10
cost_vec <- 2 ^ seq(-5, 5, by = 1)
class_weight <- "uniform"
kernel <- "mkl"
mkl_method <- "uniform"
mkl_cost <- 1
pheno_sim <- FALSE
loop_sim <- c("jaccard", "lin", "resnik")
loop_psd <- c("clip", "flip", "shift")

for (sim in loop_sim){
  sim_method <- sim
  for (psd in loop_psd){
    psd_method <- psd
    print(paste(sim, psd, sep = " ")) # print to console to keep order
    source("helper_pheno.R")
    source("model.R")
  }
}

#' experiment 6: simulate sparse/noisy phenotypes and compare similarity measures
seed <- 42
k <- 10
cost_vec <- 2 ^ seq(-5, 5, by = 1)
class_weight <- "uniform"
kernel <- "mkl"
mkl_method <- "uniform"
mkl_cost <- 1
psd_method <- "shift"
pheno_sim <- TRUE
loop_sim <- c("jaccard", "lin", "resnik")
loop_noise <- 2 ^ seq(-2, 2, by = 1) 

for (sim in loop_sim){
  sim_method <- sim
  for (nse in loop_noise){
    term_noise <- nse
    print(paste(sim, nse, sep = ""))
    source("helper_pheno.R")
    source("model.R")
  }
}

#' experiment 7: unregularized group-level MKL
#' TODO bugfixing. block diagonal matrix ok? reason for low performance?
seed <- 42
k <- 10
cost_vec <- 2 ^ seq(-5, 5, by = 1)
class_weight <- "uniform"
kernel <- "mkl"
mkl_method <- "group"
mkl_cost <- 1
source("model.R")

