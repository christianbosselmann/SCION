# data visualization functions

# pkgs
librarian::shelf(tidyverse,
                 ggplot2,
                 pROC,
                 yardstick,
                 gridExtra,
                 ape,
                 data.table,
                 dendextend,
                 RColorBrewer,
                 egg,
                 rstatix,
                 ggpubr,
                 quiet = TRUE)

# read data: raw predictions from model script
paths <- list.files(path = "out/experiment 1", 
                    pattern = "preds", 
                    full.names = TRUE)

files <- lapply(paths, read_csv)

# fix names
names(files) <- basename(paths) %>%
  sub("([^.]+)\\.[[:alnum:]]+$", "\\1", .)

# manual: labels
labels <- c("Dirac", "Union", "MTL", "MKL")
load("t_vec.rda")
gene_labels <- c("SCN1A", "SCN2A", "SCN3A", "SCN4A", "SCN5A", "SCN8A", "SCN9A", "SCN10A", "SCN11A")

# manual: control order
# files <- list(files$`preds dirac`, files$`preds union`, files$`preds mtl`, files$`preds mkl`)

### ROC curves
roc <- list()

for (i in 1:length(files)){
  auc_label <- list()
  auc_label[[i]] <- files[[i]] %>% 
    group_by(fold) %>%
    roc_auc(truth = factor(truth, levels = c("LOF", "GOF")), LOF) %>%
    summarize(mean = format(round(mean(.estimate), digits = 3), nsmall = 3), 
              sd = format(round(sd(.estimate), digits = 3), nsmall = 3))
  
  roc[[i]] <- files[[i]] %>%
    roc_curve(truth = factor(truth, levels = c("LOF", "GOF")), LOF) %>%
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

# stacked ROC curves
roc_stacked <- ggplot() +
  geom_path(data = files[[1]] %>% roc_curve(truth = factor(truth, levels = c("LOF", "GOF")), LOF),
            aes(x = 1-specificity, y = sensitivity),
            col = brewer.pal(4, "Pastel1")[[1]],
            lwd = 1.1) +
  geom_path(data = files[[2]] %>% roc_curve(truth = factor(truth, levels = c("LOF", "GOF")), LOF),
            aes(x = 1-specificity, y = sensitivity),
            col = brewer.pal(4, "Pastel1")[[2]],
            lwd = 1.1) +
  geom_path(data = files[[3]] %>% roc_curve(truth = factor(truth, levels = c("LOF", "GOF")), LOF),
            aes(x = 1-specificity, y = sensitivity),
            col = brewer.pal(4, "Pastel1")[[3]],
            lwd = 1.1) +
  geom_path(data = files[[4]] %>% roc_curve(truth = factor(truth, levels = c("LOF", "GOF")), LOF),
            aes(x = 1-specificity, y = sensitivity),
            col = brewer.pal(4, "Pastel1")[[4]],
            lwd = 1.1) +
  geom_abline(lty = 3) +
  coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1), expand = FALSE, clip = "on") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(title = "",
       x = "False Positive Rate",
       y = "True Positive Rate") 

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
# names(files) <- labels
# files %>%
#   rbindlist(idcol = "model") %>%
#   group_by(gene, fold, model) %>%
#   accuracy(factor(truth, levels = c("GOF", "LOF")), factor(pred, levels = c("GOF", "LOF"))) %>%
#   ggplot(aes(x = gene, y = .estimate)) +
#   geom_boxplot(aes(fill = factor(model, levels = c("Dirac", "Union", "MTL", "MKL")))) +
#   theme_bw() +
#   theme(legend.title=element_blank())

### task-wise bar chart 
names(files) <- labels
task_bar <- files %>%
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
  theme(# legend.title = element_blank(),
    # axis.title.x = element_blank(),
    axis.text.x = element_text(face = "italic")) +
  labs(fill = "Method",
       x = "Channel",
       y = "Accuracy") +
  scale_fill_brewer(palette = "Pastel1") +
  coord_fixed(ratio = 5, ylim = c(0,1), expand = FALSE, clip = "on")

# combine stacked roc curves and bars
pdf("fig/Figure 1.pdf", height = 4, width = 12, onefile = FALSE)
egg::ggarrange(roc_stacked, task_bar, 
               nrow = 1, ncol = 2, 
               widths = c(1/3, 2/3),
               labels = c("a", "b"),
               label.args = list(gp = grid::gpar(fontface = "plain", cex = 1.3))) # get.gpar()
dev.off()

### noise dot chart 
fig <- read_delim("out/experiment 4/fig.csv", ";", 
                  trim_ws = TRUE, 
                  escape_double = FALSE, 
                  col_types = cols(X5 = col_skip(), 
                                   X6 = col_skip(), X7 = col_skip(), 
                                   X8 = col_skip(), X9 = col_skip(), 
                                   X10 = col_skip(), X11 = col_skip()))

# fig %>% 
#   ggplot(aes(x = log(noise), y = accuracy, 
#              fill = method, group = method)) + 
#   geom_line(aes(linetype = method, color = method)) +
#   geom_point(aes(shape = method, color = method)) +
#   geom_errorbar(aes(ymin = accuracy-sd, ymax = accuracy+sd, color = method), width=.1) +
#   theme_bw() +
#   scale_color_brewer(palette = "Dark2") +
#   coord_fixed(ratio = 8, xlim = c(-1.5, 1.5), ylim = c(0.69, 0.91), expand = FALSE, clip = "on") +
#   labs(x = "Log (Noise)",
#        y = "Accuracy")

noise_dot <- fig %>% ggplot(aes(x = factor(noise, levels = c(0.25, 0.5, 1, 2, 4)), 
                                y = accuracy,
                                ymin = accuracy-sd, 
                                ymax = accuracy+sd,
                                color = method)) + 
  geom_errorbar(position = position_dodge(width = 0.6)) + 
  geom_point(position = position_dodge(width = 0.6)) +
  scale_x_discrete(drop = FALSE) +
  theme_bw(plot.margin = margin(22, 5.5, 5.5, 5.5, "pt")) +
  coord_fixed(ratio = 6, ylim = c(0.5, 1), expand = FALSE, clip = "on") +
  scale_color_brewer("Method", palette = "Dark2") +
  labs(y = "Accuracy") +
  labs(x = expression("Sparsity" %<-% "Phenotypic information" %->% "Noise"))

# psd correction dot chart
fig <- read_delim("out/experiment 3/fig.csv", ";")

noise_corr <- fig %>% ggplot(aes(x = factor(corr, levels = c("clip", "flip", "shift")), 
                                 y = accuracy,
                                 ymin = accuracy-sd, 
                                 ymax = accuracy+sd,
                                 color = method)) + 
  geom_errorbar(position = position_dodge(width = 0.6)) + 
  geom_point(position = position_dodge(width = 0.6)) +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  theme(legend.position = "none",
        plot.margin = margin(22, 5.5, 5.5, 5.5, "pt")) +
  coord_fixed(ratio = 6, ylim = c(0.5, 1), expand = FALSE, clip = "on") +
  scale_color_brewer("Method", palette = "Dark2") +
  labs(y = "Accuracy",
       x = "Matrix correction")

# combine noise and psd correction plots
pdf("fig/Figure 2.pdf", height = 4, width = 12, onefile = FALSE)
egg::ggarrange(noise_corr, noise_dot, 
               nrow = 1, ncol = 2, 
               widths = c(1/3, 2/3),
               labels = c("a", "b"),
               label.args = list(gp = grid::gpar(fontface = "plain", cex = 1.3))) # get.gpar()
dev.off()

### visualize kernel matrices, expects objects from model script (Kb, Kt, Km, etc)
library(librarian)
librarian::shelf(tidyverse, 
                 ggcorrplot,
                 ggpubr,
                 cowplot)

# list kernel matrices of interest
# Kb: structural similarity
# Kt: task similarity
# Km: MTL kernel matrix
# Kp: phenotypic similarity
# Kmkl: MTMKL kernel matrix
mats <- list(Kb, Kt, Km, Kp, Kmkl)

# visualize
p <- list()
for(i in 1:length(mats)){
  mat <- mats[[i]]
  rownames(mat) <- NULL
  colnames(mat) <- NULL
  mat <- round(cor(mat), 3)
  p[[i]] <- ggcorrplot(mat, 
                       hc.order = TRUE,
                       hc.method = "average",
                       outline.color = NA,
                       legend.title = "  r",
                       colors = c("#6D9EC1", "white", "#E46726"))
  
  p[[i]] <- p[[i]] +
    theme_void() +
    theme(legend.position = "none") + # we will add the legend as a separate plot
    theme(plot.margin = margin(-0, -0, -0, -0, "cm")) # crop
}

l <- cowplot::get_legend(p[[i]] + theme(legend.position = "right"))
p[[i+1]] <- l

# arrange
pdf("fig/Figure 4.pdf", height = 4, width = 8, onefile = FALSE)
ggpubr::ggarrange(plotlist = p,
                  nrow = 2, ncol = 3, 
                  labels = c("a", "b", "c", "d", "e", ""),
                  font.label = list(face = "plain")
)
dev.off()

### dendrogram of weights for hierarchical decomposition MTMKL
# requires distance matrix G and list d from hierarchical MTMKL containing weights
# also computes refined pairwise task similarity γk,l (see Widmer et al. 2010)

# first note that distance matrix G from function dist.alignment contains square-roots of distances
G <- G^2 

df_2tasks <- setDT(CJ(names(G), names(G), unique = TRUE))
# df_2tasks <- df_2tasks[!duplicated(t(apply(df_2tasks, 1, sort))),]
mat <- matrix(data = NA, nrow = length(names(G)), ncol = length(names(G)))
rownames(mat) <- names(G)
colnames(mat) <- names(G)

# d <- mat_mkl[[4]]
for (i in 1:nrow(df_2tasks)){
  t <- as.character(df_2tasks[i,]) # get pairwise tasks
  tmp <- sapply(d$graph, function(x) t %in% x) # get indices of nodes containing task
  ind <- which(tmp[1,] & tmp[2,])
  wt <- sum(d$delta[ind]) # sum weight of task set
  mat[t[1], t[2]] <- wt
}

# TODO fix refined taskwise similarity from hierarchical decomposition
# image(cor(mat))
# image(cor(G))

# graph <- as.dist(G, diag = TRUE)
graph <- dist(as.matrix(G), diag = F)
graph <- as.matrix(graph)
colnames(graph) <- rownames(graph) <- names(G)
graph <- as.dist(graph)
graph <- hclust(graph, method = "average")
graph <- as.dendrogram(graph) 

lbl <- str_sort(names(G), numeric = TRUE) # more sensible label order

wt <- d$gamma %>%
  lapply(., function(x) x[[1]][[1]]) %>% # get phenotype weight for nodes
  lapply(., as.data.frame) %>%
  rbindlist(.) %>%
  unlist()

graph %>%
  set("nodes_pch", 1) %>%
  set("nodes_cex", 3*wt) %>%
  # set("labels", value = names(G)) %>%
  # assign_values_to_nodes_nodePar(value = 3*unlist(wt), nodePar = "cex") %>%
  # set("branches_lwd", 1+scale(d$delta[1:17])) %>%
  # dendextend::rotate(lbl) %>%
  # hang.dendrogram %>%
  plot()

### corrplot and distance matrix as network graph
# G <- as.matrix(G)
# rownames(G) <- colnames(G)
# G <- G[lbl,lbl] # fix order
# # G2 <- as.matrix(dist(G)) # probably not neccessary, as G is already a distance matrix
# 
# library(corrplot)
# corrplot(as.matrix(G), 
#          method = "color", 
#          cl.lim=c(0,1),
#          col = colorRampPalette(c("blue", "white", "red"))(200),
#          type = "full", 
#          diag = FALSE)
# 
# library(qgraph)
# qgraph(as.matrix(G), 
#        labels = colnames(G),
#        graph = "default",
#        layout = "spring", 
#        theme = "gray",
#        vsize = 10)