#' description: HPO annotations for rare diseases [7377: OMIM; 47: DECIPHER; 3300 ORPHANET]
#' date: 2019-01-03
#' tracker: https://github.com/obophenotype/human-phenotype-ontology
#' HPO-version: http://purl.obolibrary.org/obo/hp/releases/2018-12-21/hp.owl
#' @params phenotype.hpoa
#' @returns phenotype.csv, a lookup table for use by the SCION app

df <- read_delim("app/phenotype.hpoa", skip = 4, delim = "\t") %>%
  select(`#DatabaseID`, DiseaseName, HPO_ID)

write_csv(df, "app/phenotype.csv")



