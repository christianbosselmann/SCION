# SCN 
# functional variant prediction for voltage-gated sodium channels
# helper_merge.R
# this script reads a list of variants 
# {gene, variant position, aa substitution, functional effect}
# and merges in aa and str features

# packages
library("librarian")
librarian::shelf(tidyverse,
                 data.table,
                 openxlsx,
                 fastDummies,
                 hablar,
                 janitor,
                 tidymodels,
                 quiet = TRUE)

# read data
dat_raw <- read_csv("data/clean_tbl.csv") %>%
  select(-pheno) # remove clinical data for now

# remove redundant observations and find conflicts
dat_dupe <- dat_raw %>%
  unique() %>%
  get_dupes(aa1, pos, aa2, gene)

dat_raw <- anti_join(dat_raw, dat_dupe) 

# read feature tables
aa_feats <- read_tsv("features/aa_feats.tsv")
load("features/str_feats.rda")
cid_raw <- read_csv("features/cid.csv")

# merge aa
dat_raw <- dat_raw %>%
  merge(aa_feats, by = c("aa1", "aa2")) 

# unlist and merge str
str_feats <- do.call(rbind, Map(cbind, str_feats, id = names(str_feats))) %>%
  dplyr::rename(gene = "id")

dat_raw <- dat_raw %>%
  merge(str_feats, by = c("gene", "pos"))

# split, pivot and merge canonical alignment id (cid) and annotation features
cid_align <- cid_raw[,1:10] # split into domain annotation and actual alignment
cid_feats <- cid_raw[,10:15] # currently a magic number

cid_align <- cid_align %>%
  pivot_longer(!cid, names_to = "gene", values_to = "pos") %>%
  group_by(gene) %>%
  arrange(pos) %>%
  distinct(pos, .keep_all = TRUE)

dat_raw <- merge(dat_raw, cid_align)

dat_raw <- dat_raw %>%
  merge(cid_feats, by = c("cid"))

# drop unneccessary features
dat_prep <- dat_raw %>% 
  select(-one_of("pos", 
                 "id",
                 "aa1", 
                 "aa2"))

# fix data types
dat_prep <- hablar::retype(dat_prep)

# preprocessing
dat_rec <- recipe(y ~ ., dat_prep) %>%
  step_dummy(c("str_np_q3", 
               "str_np_q8", 
               "str_pbie", 
               "str_prhl", 
               "str_pito"), 
             one_hot = TRUE) %>%
  step_normalize(all_numeric())

dat_rec <- prep(dat_rec)

dat_prep <- bake(dat_rec, NULL)

# save preprocessing recipe
save(dat_rec, file = "recipe.rds") 

# export to data backup
write_csv(dat_prep, "data/dat_prep.csv")
