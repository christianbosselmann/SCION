# SCN
# Functional variant prediction for voltage-gated sodium channels
# helper_val.R
# this script helps prepare the validation data set for comparison with
# fuNCion by Heyne et al., DOI 10.1126/scitranslmed.aay6848 (Table S4)

# packages
library("librarian")
librarian::shelf(tidyverse,
                 data.table,
                 openxlsx,
                 quiet = TRUE)

# get data
dat_heyne <- read_tsv("data/aay6848_table_s4.txt") %>%
  select(obs, gene, refAA, pos, altAA) %>%
  rename(y = obs, gene = gene, aa1 = refAA, pos = pos, aa2 = altAA) %>%
  filter(y == "lof" | y == "gof") %>%
  filter(complete.cases(.)) # weird missing data in Heyne's Table S4

# rename label for consistency
dat_heyne$y <- toupper(dat_heyne$y)

# export
write_csv(dat_heyne, "data/dat_heyne.csv")

# get data
dat_heyne <- read_csv("data/dat_heyne.csv")
dat_scn <- read_csv("data/clean_tbl.csv") %>%
  select(-id, -pheno)

# remove Heyne's validation variants from our training data
dat_val <- anti_join(dat_scn, dat_heyne)
write_csv(dat_val, "data/dat_val.csv")

