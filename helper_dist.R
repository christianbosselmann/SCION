# SCN
# Functional variant prediction for voltage-gated sodium channels
# helper_dist.R
# script to generate a histogram of projections
# cf. Simple Method for Interpretation of High-Dimensional Nonlinear SVM Classification Models. January 2010. Source DBLP. Conference: Proceedings of The 2010 International Conference on Data Mining, DMIN 2010, July 12-15, 2010, Las Vegas, Nevada, USA

# pkg
library(librarian)
librarian::shelf(tidyverse,
                 ggplot2,
                 recipes,
                 gridExtra,
                 egg)

# helper functions
source("func.R")

# load raw predictions containg a column "dist" with distance to hyperplane for each observation
data_dist <- read_csv("out/dist/report_preds_2022-06-02 10:01:27.csv") %>%
  arrange(ind)

# load dataset with features
data_all <- read_csv("data/dat_prep.csv")

# plot histogram
plot_hist <- data_dist %>%
  ggplot(aes(x = dist, color = pred, fill = pred)) +
  geom_histogram(position = "identity", binwidth = 0.1, boundary = 0,
                 alpha = 0.3,
                 aes(y = ..density..)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  geom_density(alpha = 0) +
  geom_vline(xintercept=c(-1,-1), linetype = "dotted") +
  geom_vline(xintercept=c(1,1), linetype = "dotted") +
  geom_vline(xintercept=c(0,0)) +
  theme_bw() +
  scale_color_manual(name = "Predicted class", values = c("#44803F", "#FF5A33")) +
  scale_fill_manual(name = "Predicted class", values = c("#44803F", "#FF5A33")) +
  labs(x = "Distance from hyperplane", 
       y = "Density") 

# define set of observation for feature correlation
ind <- data_dist %>%
  filter(dist < 0.75, dist > -0.75) %>% # set bounds here
  select(ind) %>%
  unlist() %>%
  as.numeric()

# mark all observations that are in the set of interest (i.e. hard to predict)
data_all$is_hard <- NA
data_all$is_hard[ind] <- "Low"
data_all$is_hard[-ind] <- "High"

# simple feature correlation
pvalue <- vector(length = ncol(data_all))

for (i in 1:ncol(data_all)) {
  feature <- data_all[,i] %>%
    unlist()
  
  if(!is.numeric(feature)){next}
  
  pvalue[i] <- wilcox.test(feature ~ is_hard,
                           data = data_all)$p.value
}

# adjust p-values for multiple testing with Benjamini-Hochberg method
# then print significant vars
result <- cbind(colnames(data_all), pvalue) %>%
  as_tibble() %>%
  mutate(pvalue_adj = p.adjust(.$pvalue, method = "fdr")) %>%
  filter(pvalue_adj < 0.05, !pvalue_adj == 0) %>%
  print()

result <- result[c(1,3),] # manual filter due to redundant, i.e. correlated features

# for loop to iterate through all features
load("recipe.rds")

pretty_names <- c("Position on family alignment",
                  # "Relative accessible surface area",
                  "Accessible surface area" 
                  # "Relative sequence position"
                  )

plots <- list()
for (i in 1:length(result$V1)) {
  ft <- result$V1[[i]]
  
  # revert features of interest back into their original values
  data_all[[ft]] <- unNormalize(data_all[[ft]], dat_rec, ft)
  
  # plot confidence across the channel superfamily alignment
  plots[[i]] <- data_all %>%
    ggplot(aes_string(x = ft, fill = "is_hard")) +
    geom_histogram(position = "identity", 
                   bins = 30,
                   boundary = 0,
                   alpha = 0.5) +
    theme_bw() +
    scale_color_manual(name = "Confidence", values = c("#0571b0", "#ca0020")) +
    scale_fill_manual(name = "Confidence", values = c("#0571b0", "#ca0020")) +
    labs(x = pretty_names[[i]], # can also just be ft if the names don't need to be pretty
         y = "Number of variants") +
    theme(panel.grid.minor = element_blank()) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA))
  
  # need just one legend for the panel
  if(i != 2){plots[[i]] <- plots[[i]] + theme(legend.position = "none")}
}

# composit figure
bottom_row <- ggarrange(plots[[1]], plots[[2]], NULL, ncol = 3, labels = c("b", "c", ""), widths = c(4/10, 5/10, 1/10), font.label = list(size = 14, color = "black", face = "plain", family = NULL))
top_row <- ggarrange(plot_hist, NULL, ncol = 2, labels = c("a", ""), widths = c(9.2/10, 0.8/10), font.label = list(size = 14, color = "black", face = "plain", family = NULL))

pdf("fig/Figure 3.pdf", height = 8, width = 12, onefile = FALSE)
ggarrange(top_row, bottom_row, ncol = 1, heights = c(2/3, 1/3))
dev.off()