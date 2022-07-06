# helper script
# takes SCN_Aligment.xlsx from scn-viewer.braodinstitute.org and converts it into a lookup table for the Shiny app

# pkg
library(zoo)
library(readxl)

# data
tbl_scn <- read_excel("data/SCN_Aligment.xlsx")

# save indices
ind <- which(is.na(tbl_scn$Index))

# fill NA indices ("indices of interest")
tbl_scn$Index <- na.locf(tbl_scn$Index, fromLast = FALSE)

# keep table of index and domain for NGLVieweR AlphaFold visualization
tbl_dom <- tbl_scn %>%
  filter(!is.na(Index)) %>%
  mutate(across(3:12, str_extract, "\\d+"))

# filter, remove redundant column and transform into a column vector with position of known analogous variants
tbl_scn <- tbl_scn[ind,] %>%
  select(-Functional_Effects) %>%
  transmute(cid = Index)

# export
write_csv(tbl_scn, "app/scnviewer_lookup.csv")
write_csv(tbl_dom, "app/nglviewer_lookup.csv")
