# packages
library(fmsb)
library(tidyverse)
library(hablar)

# helper functions
guide_axis_label_trans <- function(label_trans = identity, ...) {
  axis_guide <- guide_axis(...)
  axis_guide$label_trans <- rlang::as_function(label_trans)
  class(axis_guide) <- c("guide_axis_trans", class(axis_guide))
  axis_guide
}

guide_train.guide_axis_trans <- function(x, ...) {
  trained <- NextMethod()
  trained$key$.label <- x$label_trans(trained$key$.label)
  trained
}

# get data
data <- out %>%
  dplyr::select(prob.GOF, prob.LOF, prob.Neutral) %>%
  dplyr::rename(GOF = prob.GOF, LOF = prob.LOF, Neutral = prob.Neutral)

# add the max and min of each variable to show on the plot
data <- rbind(rep(1,5), rep(0,5), data)

# color vector
colors_border = c(rgb(12/255,45/255,72/255,0.9))
colors_in = c(rgb(20/255,93/255,160/255,0.4))

# set plot size
par(oma=c(0,1,1,1))
par(mar=c(0,0,0,0))

# plot with default options:
radarchart( data, axistype = 1, 
            # custom polygon
            pcol = colors_border, pfcol = colors_in, plwd = 4, plty = 1,
            # custom the grid
            cglcol = "darkgrey", cglty = 1, axislabcol = "darkgrey", caxislabels = c("0","","","","1"), cglwd = 0.8,
            # custom labels
            vlcex = 1
)

# add a legend
legend(x = 0.7, y = 1, legend = c("class probability"), bty = "n", pch = 15, col = colors_in, text.col = "darkgrey", cex = 1.2, pt.cex = 3)
legend(x = 0.9, y = 0.8, legend = paste("GOF:", round(out$prob.GOF, digits = 3), sep = " "), bty = "n", pch = NA, col = colors_in, text.col = "darkgrey", cex = 0.9, pt.cex = 3)
legend(x = 0.9, y = 0.7, legend = paste("LOF:", round(out$prob.LOF, digits = 3), sep = " "), bty = "n", pch = NA, col = colors_in, text.col = "darkgrey", cex = 0.9, pt.cex = 3)
legend(x = 0.9, y = 0.6, legend = paste("Neutral:", round(out$prob.Neutral, digits = 3), sep = " "), bty = "n", pch = NA, col = colors_in, text.col = "darkgrey", cex = 0.9, pt.cex = 3)

# record 
p <- recordPlot()

# another plot for confidence, i.e. distance from hyperplane
data2 <- data.frame(out$conf.LOF_Neutral, out$conf.LOF_GOF, out$conf.Neutral_GOF) %>%
  pivot_longer(cols = everything())

data2$fill <- NA
data2$fill[data2$name == "out.conf.LOF_Neutral"] <- ifelse(data2$value[1] > 0, "#FBB4AE", "#B3CDE3")
data2$fill[data2$name == "out.conf.LOF_GOF"] <- ifelse(data2$value[2] > 0, "#FBB4AE", "#CCEBC5")
data2$fill[data2$name == "out.conf.Neutral_GOF"] <- ifelse(data2$value[3] > 0, "#B3CDE3", "#CCEBC5")

data2$name <- c(0,1,2)
data2$name <- as.numeric(data2$name)

p2 <- ggplot(data2, aes(x = name, y = value, fill = fill, label = round(value, 2))) + 
  geom_bar(stat="identity", position = "identity", width = 0.4) +
  geom_text(size = 5, position = position_stack(vjust = 0.5)) +
  scale_x_continuous(breaks = seq(0, 2, by = 1), labels = c("Neutral", "GOF", "GOF"),
                     sec.axis = dup_axis(breaks = seq(0, 2, by = 1),
                                         labels = c("LOF", "LOF", "Neutral"))) +
  theme_bw() +
  geom_hline(yintercept = 0) +
  theme(axis.title.x = element_blank()) +
  labs(y="confidence (distance from hyperplane)") +
  scale_fill_identity() +
  coord_cartesian(ylim = c(-2,2)) +
  theme(text = element_text(size=15))

