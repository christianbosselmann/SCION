# this takes the supplement from DOI doi.org/10.1093/brain/awac006 and extracts
# SCN variants with their functional effect and main phenotype

# packages
library(librarian)
librarian::shelf(tidyverse,
                 pdftools,
                 rJava,
                 tabulizer,
                 data.table,
                 tesseract)

# scrape suppl pdf with tabula, select No. - Primary disease
raw_ext <- tabulizer::extract_areas("awac006_supplementary_data.pdf",
                                pages = 5:27,
                                output = "data.frame") 

# fix: first variant of a page ends up as colname
fix_names <- function(x) {rbind(colnames(x), x)}
raw_ext <- lapply(raw_ext, fix_names)

# fix page 12
split <- do.call(rbind, strsplit(raw_ext[[12]]$Zaharieva..2016.14.CMS, ' (?=[^ ]+$)', perl=TRUE))
raw_ext[[12]]$X <- split[,2]

# combine pages
raw_tbl <- rbindlist(raw_ext, use.names = FALSE)

# clean up
raw_tbl$X1 <- gsub("[^0-9+]", "", raw_tbl$X1) 

raw_tbl <- raw_tbl %>% # only select rows with variant ids
  filter(grepl("\\d", X))

raw_tbl <- raw_tbl %>%
  select(., -5) %>% # drop ref
  rename(id = X, gene = X.1, sub = X.2, y = Overall, pheno = Primary)

# filter only GOF and LOF, then fix label
raw_tbl <- raw_tbl %>%
  filter(y == "GoF" | y == "LoF")

raw_tbl$y[raw_tbl$y == "GoF"] <- "GOF"
raw_tbl$y[raw_tbl$y == "LoF"] <- "LOF"

# extract aa1, pos, aa2
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

raw_tbl$aa1 <- substr(raw_tbl$sub, 1, 1)
raw_tbl$pos <- str_extract(raw_tbl$sub, "[[:digit:]]+")
raw_tbl$aa2 <- substrRight(raw_tbl$sub, 1)
raw_tbl$sub <- NULL

# manual phenotype checks
# write_csv(raw_tbl, "raw_tbl.csv")
raw_tbl$pheno[raw_tbl$pheno == "without"] <- "NDD"
raw_tbl$pheno[raw_tbl$id == "139"] <- "fetal akinesia"
raw_tbl$pheno[raw_tbl$id == "154"] <- "HyperPP; PMC"
raw_tbl$pheno[raw_tbl$id == "226"] <- "HyperPP; PMC"
raw_tbl$pheno[raw_tbl$id == "378"] <- "LQT3; BrS"

# export
write_csv(raw_tbl, "clean_tbl.csv")
