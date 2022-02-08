# SCN
# Functional variant prediction for voltage-gated sodium channels
# helper_msa.R
# this script generates a MSA of all .fasta files in the directory,
# displays and saves a phylogenetic tree und exports a distance matrix

# packages
library("librarian")
librarian::shelf(tidyverse,
                 data.table,
                 msa,
                 seqinr,
                 ape,
                 bios2mds,
                 binaryLogic,
                 fastDummies,
                 tools,
                 quiet = TRUE)

set.seed(777)

# get paths and fasta files
path <- list.files(path = "fasta", 
                   pattern = "\\.fasta$", 
                   full.names = TRUE)
fasta <- readAAStringSet(path)

# fix names
names <- basename(path)
names <- file_path_sans_ext(names)
names(fasta) <- names

# run alignment
aligned <- msa(fasta, method = "Muscle")

# get distance values and draw phylogenetic tree
aligned_conv <- msaConvert(aligned, type = "seqinr::alignment")
distance <- dist.alignment(aligned_conv, "identity", gap = FALSE) # tune here
tree <- nj(distance)

tiff("fig/tree.tiff", units="in", width=6, height=5, res=300)
plot(tree)
dev.off()

# export distance matrix for later use as similarity matrix
d <- data.frame(as.matrix(distance))
write_csv(d, "mat/distancematrix.csv")

# canonical id
aligned_cw <- msaConvert(aligned, "bios2mds::align")

cid <- data.frame()

for (i in seq_along(names)) {
  cid <- aligned_cw[names[i]] %>%
    unlist(use.names = FALSE) %>%
    rbind(cid)
}

colnames(cid) <- seq_along(cid)
cid <- t(cid)
colnames(cid) <- rev(names)
cid <- as.data.table(cid)

for (i in colnames(cid)) {
  cid[[i]] <- cumsum(cid[[i]] %in% LETTERS)
}

cid$cid <- rownames(cid)

# manual annotation of important structural and functional domains and motifs
# from DOI 10.1093/brain/awac006
# note: check this any time the MSA is re-run

cid$str_dom <- NA
cid$str_dom[1:128] <- "cytoplasmic"
cid$str_dom[129:147] <- "S1 of D1"
cid$str_dom[148:154] <- "extracellular"
cid$str_dom[155:175] <- "S2 of D1"
cid$str_dom[176:189] <- "cytoplasmic"
cid$str_dom[190:207] <- "S3 of D1"
cid$str_dom[208:213] <- "extracellular"
cid$str_dom[214:230] <- "S4 of D1"
cid$str_dom[231:249] <- "cytoplasmic"
cid$str_dom[250:269] <- "S5 of D1"
cid$str_dom[270:367] <- "extracellular"
cid$str_dom[368:392] <- "pore"
cid$str_dom[393:399] <- "extracellular"
cid$str_dom[400:420] <- "S6 of D1"
cid$str_dom[421:768] <- "cytoplasmic"
cid$str_dom[769:787] <- "S1 of D2"
cid$str_dom[788:798] <- "extracellular"
cid$str_dom[799:818] <- "S2 of D2"
cid$str_dom[819:832] <- "cytoplasmic"
cid$str_dom[833:852] <- "S3 of D2"
cid$str_dom[853:854] <- "extracellular"
cid$str_dom[855:872] <- "S4 of D2"
cid$str_dom[873:888] <- "cytoplasmic"
cid$str_dom[889:907] <- "S5 of D2"
cid$str_dom[908:936] <- "extracellular"
cid$str_dom[937:957] <- "pore"
cid$str_dom[958:970] <- "extracellular"
cid$str_dom[971:991] <- "S6 of D2"
cid$str_dom[992:1219] <- "cytoplasmic"
cid$str_dom[1220:1237] <- "S1 of D3"
cid$str_dom[1238:1250] <- "extracellular"
cid$str_dom[1251:1269] <- "S2 of D3"
cid$str_dom[1270:1283] <- "cytoplasmic"
cid$str_dom[1284:1302] <- "S3 of D3"
cid$str_dom[1303:1310] <- "extracellular"
cid$str_dom[1311:1329] <- "S4 of D3"
cid$str_dom[1330:1346] <- "cytoplasmic"
cid$str_dom[1347:1366] <- "S5 of D3"
cid$str_dom[1367:1418] <- "extracellular"
cid$str_dom[1419:1440] <- "pore"
cid$str_dom[1441:1457] <- "extracellular"
cid$str_dom[1458:1479] <- "S6 of D3"
cid$str_dom[1480:1542] <- "cytoplasmic"
cid$str_dom[1543:1560] <- "S1 of D4"
cid$str_dom[1561:1571] <- "extracellular"
cid$str_dom[1572:1590] <- "S2 of D4"
cid$str_dom[1591:1602] <- "cytoplasmic"
cid$str_dom[1603:1620] <- "S3 of D4"
cid$str_dom[1621:1633] <- "extracellular"
cid$str_dom[1634:1650] <- "S4 of D4"
cid$str_dom[1651:1669] <- "cytoplasmic"
cid$str_dom[1670:1687] <- "S5 of D4"
cid$str_dom[1688:1709] <- "extracellular"
cid$str_dom[1710:1732] <- "pore"
cid$str_dom[1733:1762] <- "extracellular"
cid$str_dom[1763:1785] <- "S6 of D4"
cid$str_dom[1786:2009] <- "cytoplasmic"

# binary encoding of ordered domains
encode_binary <- function(x, order = unique(x), name = "v_") {
  x <- as.numeric(factor(x, levels = order, exclude = NULL))
  x2 <- as.binary(x)
  maxlen <- max(sapply(x2, length))
  x2 <- lapply(x2, function(y) {
    l <- length(y)
    if (l < maxlen) {
      y <- c(rep(0, (maxlen - l)), y)
    }
    y
  })
  d <- as.data.frame(t(as.data.frame(x2)))
  rownames(d) <- NULL
  colnames(d) <- paste0(name, 1:maxlen)
  d
}

cid <- cbind(cid, encode_binary(cid[["str_dom"]], name = "str_dom"))
cid$str_dom <- NULL

# export
write_csv(cid, "features/cid.csv")
