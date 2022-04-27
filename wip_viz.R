# toy example to visualize data

# pkg
librarian::shelf(tidyverse,
                 ggplot2,
                 pROC,
                 yardstick,
                 gridExtra,
                 quiet = TRUE)

# read data: raw predictions from model script
paths <- list.files(path = "out/experiment 2", 
                    pattern = "preds", 
                    full.names = TRUE)

files <- lapply(paths, read_csv)

# fix names
names(files) <- basename(paths) %>%
  sub("([^.]+)\\.[[:alnum:]]+$", "\\1", .)

# manual labels
labels <- c("Dirac", "MKL", "MTL", "Union") 

### ROC curves
roc <- list()

for (i in 1:length(files)){
  auc_label <- list()
  auc_label[[i]] <- files[[i]] %>% 
    group_by(fold) %>%
    roc_auc(truth = factor(truth, levels = c("GOF", "LOF")), GOF) %>%
    summarize(mean = round(mean(.estimate), 3), sd = round(sd(.estimate), 3))
  
  roc[[i]] <- files[[i]] %>%
    roc_curve(truth = factor(truth, levels = c("GOF", "LOF")), GOF) %>%
    ggplot(aes(x = 1-specificity, y = sensitivity)) +
    geom_path() +
    geom_abline(lty = 3) +
    coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1), expand = FALSE, clip = "on") +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    labs(title = labels[[i]],
         x = "False Positive Rate",
         y = "True Positive Rate") +
    annotate("text", x = 0.750, y = 0.125, label = paste("AUC:", auc_label[[i]]$mean, sep = " "))
}

n <- length(roc)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(roc, ncol = nCol))

# PRC curves
prc <- list()

for (i in 1:length(files)){
  auc_label <- list()
  auc_label[[i]] <- files[[i]] %>% 
    group_by(fold) %>%
    pr_auc(truth = factor(truth, levels = c("GOF", "LOF")), GOF) %>%
    summarize(mean = round(mean(.estimate), 3), sd = round(sd(.estimate), 3))
  
  prc[[i]] <- files[[i]] %>%
    pr_curve(truth = factor(truth, levels = c("GOF", "LOF")), GOF) %>%
    ggplot(aes(x = recall, y = precision)) +
    geom_path() +
    coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1), expand = FALSE, clip = "on") +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    labs(title = labels[[i]],
         x = "Recall",
         y = "Precision") +
    annotate("text", x = 0.750, y = 0.125, label = paste("AUC:", auc_label[[i]]$mean, sep = " "))
}

n <- length(prc)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(prc, ncol = nCol))
