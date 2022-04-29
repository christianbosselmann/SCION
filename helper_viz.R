# data visualization functions

# pkgs
librarian::shelf(tidyverse,
                 ggplot2,
                 pROC,
                 yardstick,
                 gridExtra,
                 data.table,
                 quiet = TRUE)

# read data: raw predictions from model script
paths <- list.files(path = "out/experiment 2", 
                    pattern = "preds", 
                    full.names = TRUE)

files <- lapply(paths, read_csv)

# fix names
names(files) <- basename(paths) %>%
  sub("([^.]+)\\.[[:alnum:]]+$", "\\1", .)

# manual: control order
files <- list(files$`preds dirac`, files$`preds union`, files$`preds mtl`, files$`preds mkl`)

# manual: labels
labels <- c("Dirac", "Union", "MTL", "MKL")
load("t_vec.rda")
gene_labels <- c("SCN1A", "SCN2A", "SCN3A", "SCN4A", "SCN5A", "SCN8A", "SCN9A", "SCN10A", "SCN11A")

### ROC curves
roc <- list()

for (i in 1:length(files)){
  auc_label <- list()
  auc_label[[i]] <- files[[i]] %>% 
    group_by(fold) %>%
    roc_auc(truth = factor(truth, levels = c("GOF", "LOF")), GOF) %>%
    summarize(mean = format(round(mean(.estimate), digits = 3), nsmall = 3), 
              sd = format(round(sd(.estimate), digits = 3), nsmall = 3))
  
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
    annotate("text", x = 0.625, y = 0.125, 
             label = paste("AUC:", auc_label[[i]]$mean, sep = " "),
             size = 4)
}

n <- length(roc)
nCol <- floor(sqrt(n))
roc <- do.call("grid.arrange", c(roc, ncol = nCol))

# PRC curves
prc <- list()

for (i in 1:length(files)){
  auc_label <- list()
  auc_label[[i]] <- files[[i]] %>% 
    group_by(fold) %>%
    pr_auc(truth = factor(truth, levels = c("GOF", "LOF")), GOF) %>%
    summarize(mean = format(round(mean(.estimate), digits = 3), nsmall = 3), 
              sd = format(round(mean(.estimate), digits = 3), nsmall = 3))
  
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
    annotate("text", x = 0.625, y = 0.125, 
             label = paste("AUC:", auc_label[[i]]$mean, sep = " "),
             size = 4)
}

n <- length(prc)
nCol <- floor(sqrt(n))
prc <- do.call("grid.arrange", c(prc, ncol = nCol))

### task-wise box plot
names(files) <- labels
files %>%
  rbindlist(idcol = "model") %>%
  group_by(gene, fold, model) %>%
  accuracy(factor(truth, levels = c("GOF", "LOF")), factor(pred, levels = c("GOF", "LOF"))) %>%
  ggplot(aes(x = gene, y = .estimate)) +
  geom_boxplot(aes(fill = factor(model, levels = c("Dirac", "Union", "MTL", "MKL")))) +
  theme_bw() +
  theme(legend.title=element_blank())

### task-wise bar chart 
files %>%
  rbindlist(idcol = "model") %>%
  group_by(gene, fold, model) %>%
  accuracy(factor(truth, levels = c("GOF", "LOF")), 
      factor(pred, levels = c("GOF", "LOF"))) %>%
  replace(is.na(.), 0) %>%
  group_by(gene, model) %>%
  summarize(mean = mean(.estimate)) %>%
  ggplot(aes(x = factor(gene, levels = gene_labels), 
             y = mean, 
             fill = factor(model, levels = c("Dirac", "Union", "MTL", "MKL")))) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = "Gene",
       y = "Accuracy") +
  scale_fill_brewer(palette = "Pastel1") +
  coord_fixed(ratio = 5, ylim = c(0,1), expand = FALSE, clip = "on")
