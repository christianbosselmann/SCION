# SCN
# Functional variant prediction for voltage-gated sodium channels
# helper_strf.R
# this script loads, slices and filters the raw feature output from 
# PredictProtein, NetSurfP, IUPred. The output is a list of named
# data.tables containing each channel's structural features.

# packages
library("librarian")
librarian::shelf(tidyverse,
                 data.table,
                 openxlsx,
                 quiet = TRUE)

# get vector of protein names
tmp_names <- list.dirs("features", full.names = FALSE, recursive = FALSE)
tmp_names_full <- list.dirs("features", full.names = TRUE, recursive = FALSE)

# empty list of data tables to store structural features
str_all <- list()

for(i in 1:length(tmp_names)){
  
  # get file path
  tmp_name <- tmp_names[[i]]
  tmp_path <- paste("features", tmp_name, "netsurfp.csv", sep = "/")
  
  # initialize empty data table and name accordingly
  str_all[[i]] <- data.table()
  names(str_all)[i] <- tmp_name
  
  # get seq length
  tmp_seqlength <- readr::read_delim(tmp_path, 
                                     delim = ",", 
                                     skip = 0, 
                                     col_names = TRUE)
  tmp_seqlength <- length(tmp_seqlength$seq)  
  
  # features/structure/predictprotein
  # documentation: https://rostlab.org/owiki/index.php/PredictProtein_-_Documentation
  # method: input canonical uniprot isoform, download all results
  str_all[[i]] <-
    readr::read_delim(paste(tmp_names_full[[i]], "query.consurf.grades", sep = "/"), 
                      delim = "\t", 
                      skip = 14, 
                      col_names = FALSE) %>%
    slice(1:tmp_seqlength) %>%
    select(X3) %>%
    rename(str_consurf = X3) %>%
    add_column(str_all[[i]])
  
  str_all[[i]] <-
    readr::read_delim(paste(tmp_names_full[[i]], "query.isis", sep = "/"), 
                      delim = " ", 
                      skip = 0, 
                      col_names = FALSE) %>%
    filter(grepl("^[0-9]*$", X1)) %>%
    select(X3) %>%
    rename(str_isis = X3) %>%
    add_column(str_all[[i]])
  
  str_all[[i]] <-
    readr::read_delim(paste(tmp_names_full[[i]], "query.mdisorder", sep = "/"),
                      delim = "\t", 
                      skip = 1, 
                      col_names = FALSE) %>%
    slice(1:tmp_seqlength) %>%
    select(X3, X5, X7) %>%
    rename(str_nors = X3, str_profbval = X5, str_ucon = X7) %>%
    add_column(str_all[[i]])
  
  # weird helper function
  # basically the PredictProtein output for the phdRdb query has a variable
  # header length, which is why we can't just read_delim and skip. First, we
  # read lines as character vectors, remove elements starting with '#' and then
  # let read_delim do the rest.
  tmp_tbl <- readLines(paste(tmp_names_full[[i]], "query.phdRdb", sep = "/")) %>%
    as.tibble(.) %>%
    filter(!grepl("(#)", .[[1]])) %>%
    pull(.)
  
  str_all[[i]] <- 
    readr::read_delim(tmp_tbl,
                      delim = "\t") %>% 
    slice(2:(tmp_seqlength+1)) %>%
    select(OtH, OtL, PRHL, PiTo) %>%
    rename(str_helix = OtH, str_loop = OtL, str_prhl = PRHL, str_pito = PiTo) %>%
    add_column(str_all[[i]])
  
  str_all[[i]] <-
    readr::read_delim(paste(tmp_names_full[[i]], "query.prof1Rdb", sep = "/"),
                      delim = "\t", 
                      skip = 78, 
                      col_names = TRUE) %>%
    select(PACC, PREL, Pbie) %>%
    rename(str_asa = PACC, str_rsa = PREL, str_pbie = Pbie) %>%
    add_column(str_all[[i]])
  
  # features/structure/iupred2a
  # documentation: DOI 10.1093/nar/gky384, DOI 10.1002/cpbi.99
  # method: input canonical uniprot isoform, select default long disorder, download html
  str_all[[i]] <-
    readr::read_delim(paste(tmp_names_full[[i]], "iupred.html", sep = "/"),
                      delim = "\t", 
                      skip = 7, 
                      col_names = TRUE) %>%
    slice(1:tmp_seqlength) %>%
    select("IUPRED SCORE", "ANCHOR SCORE") %>%
    rename(str_iupred = "IUPRED SCORE", str_anchor = "ANCHOR SCORE") %>%
    add_column(str_all[[i]])
  
  # features/structure/netsurfp
  # method: run netfsurp2.0 webservice (https://services.healthtech.dtu.dk/service.php?NetSurfP-2.0) and export all files as csv. rename accordingly.
  # documentation: DOI 10.1002/prot.25674
  str_all[[i]] <-
    readr::read_delim(paste(tmp_names_full[[i]], "netsurfp.csv", sep = "/"),
                      delim = ",", 
                      skip = 0, 
                      col_names = TRUE) %>%
    select("rsa", "asa", "q3", "q8") %>%
    rename(str_np_rsa = "rsa", str_np_asa = "asa", str_np_q3 = "q3", str_np_q8 = "q8") %>%
    add_column(str_all[[i]])
  
  # features/structure/perviewer
  # method: download gene family .txt file and rename accordingly
  # fix column shift error manually with excel
  # documentation: DOI 10.1101/gr.252601.119
  str_all[[i]] <- 
    readr::read_delim("features/per1.txt",
                      delim = "\t", 
                      skip = 0, 
                      col_names = TRUE) %>%
    filter(!grepl('-', !!as.symbol(tmp_name))) %>%
    select(Parazscore) %>%
    rename(str_paraz = Parazscore) %>%
    slice(1:tmp_seqlength) %>% # defensive coding to avoid some isoform weirdness
    add_column(str_all[[i]])
  
  # features/id/ecdf
  # method: absolute relative position of mutation, divided by total transcript length
  str_all[[i]]$pos <- seq_len(nrow(str_all[[i]]))
  str_all[[i]]$ecdf <- ecdf(str_all[[i]]$pos)(str_all[[i]]$pos)
}

# sanity check
# there should not be any missing data
if (sum(is.na(str_all)) > 0){
  stop("missing data, check pipeline.")
}

# fix name and export
str_feats <- str_all
save(str_feats, file = "features/str_feats.rda")
