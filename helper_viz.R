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

# manual: labels
labels <- c("Uniform", "Global weights", "Hierarchical decomposition")
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
        axis.title.x = element_blank(),
        axis.text.x = element_text(face = "italic")) +
  labs(x = "Gene",
       y = "Accuracy") +
  scale_fill_brewer(palette = "Pastel1") +
  coord_fixed(ratio = 5, ylim = c(0,1), expand = FALSE, clip = "on")

### noise plot for figure 4, fig.csv is found in output folder
fig %>% 
  ggplot(aes(x = log(noise), y = accuracy, 
             fill = method, group = method)) + 
  geom_line(aes(linetype = method, color = method)) +
  geom_point(aes(shape = method, color = method)) +
  geom_errorbar(aes(ymin = accuracy-sd, ymax = accuracy+sd, color = method), width=.1) +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  coord_fixed(ratio = 8, xlim = c(-1.5, 1.5), ylim = c(0.69, 0.91), expand = FALSE, clip = "on") +
  labs(x = "Log (Noise)",
       y = "Accuracy")

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
  set("labels", names(G)) %>%
  assign_values_to_nodes_nodePar(value = 3*unlist(wt), nodePar = "cex") %>%
  set("branches_lwd", 1+scale(d$delta[1:17])) %>%
  dendextend::rotate(lbl) %>%
  # hang.dendrogram %>%
  plot()

### corrplot and distance matrix as network graph
G <- as.matrix(G)
rownames(G) <- colnames(G)
G <- G[lbl,lbl] # fix order
# G2 <- as.matrix(dist(G)) # probably not neccessary, as G is already a distance matrix

library(corrplot)
corrplot(as.matrix(G), 
         method = "color", 
         cl.lim=c(0,1),
         col = colorRampPalette(c("blue", "white", "red"))(200),
         type = "full", 
         diag = FALSE)

library(qgraph)
qgraph(as.matrix(G), 
       labels = colnames(G),
       graph = "default",
       layout = "spring", 
       theme = "gray",
       vsize = 10)


