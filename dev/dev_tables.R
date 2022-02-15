# SCN dev
# assemble result tables

# pkg
librarian::shelf(tidyverse,
                 grid, 
                 gridExtra,
                 cowplot)

# data
dirac <- read_csv("out/metrics dirac.csv") %>% round(3)
union <- read_csv("out/metrics union.csv") %>% round(3)
mtl <- read_csv("out/metrics mtl.csv") %>% round(3)
mkl <- read_csv("out/metrics jaccard semkl.csv") %>% round(3)

# set theme
cs_theme <- ttheme_default(core = list(fg_params = list(hjust = 0, x = 0.1, 
                                                       fontsize = 8)),
                          colhead = list(fg_params = list(fontsize = 9, 
                                                          fontface = "bold"))
)

# assemble
ls_metrics <- list(dirac, union, mtl, mkl)

d1 <- tableGrob(ls_metrics[[1]], rows = NULL, theme = cs_theme)
d2 <- tableGrob(ls_metrics[[2]], rows = NULL, theme = cs_theme)
d3 <- tableGrob(ls_metrics[[3]], rows = NULL, theme = cs_theme)
d4 <- tableGrob(ls_metrics[[4]], rows = NULL, theme = cs_theme)

# title
dd1 <- ggplot() + theme_bw() + theme(rect=element_rect(fill="transparent"),panel.border = element_blank()) + annotation_custom(tableGrob(ls_metrics[[1]], rows = NULL, theme = cs_theme)) + labs(title = 'Dirac')
dd2 <- ggplot() + theme_bw() + theme(rect=element_rect(fill="transparent"),panel.border = element_blank()) + annotation_custom(tableGrob(ls_metrics[[2]], rows = NULL, theme = cs_theme)) + labs(title = 'Union')
dd3 <- ggplot() + theme_bw() + theme(rect=element_rect(fill="transparent"),panel.border = element_blank()) + annotation_custom(tableGrob(ls_metrics[[3]], rows = NULL, theme = cs_theme)) + labs(title = 'MTL')
dd4 <- ggplot() + theme_bw() + theme(rect=element_rect(fill="transparent"),panel.border = element_blank()) + annotation_custom(tableGrob(ls_metrics[[4]], rows = NULL, theme = cs_theme)) + labs(title = 'MKL')

# display
p <- grid.arrange(dd1, dd2, dd3, dd4)

# save
save_plot("p.jpg", p)
