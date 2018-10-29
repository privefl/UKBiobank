ID <- bigsnpr::snp_attach("UKB-merged_sub1.rds")$fam$sample.ID

library(dplyr)

## Self-reported ancestry (https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=1001)
code_ancestry <- readr::read_tsv("coding1001.tsv")

## Cancer type (https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=3)
code_cancer <- readr::read_tsv("coding3.tsv")


df0 <- bigreadr::fread2(
  "ukb22544.csv",
  select = c("eid", "50-0.0", "34-0.0", "52-0.0", "22001-0.0", "21000-0.0",
             "2453-0.0", "20001-0.0", "21022-0.0", "189-0.0"),
  col.names = c("eid", "height", "year", "month", "sex", "pop",
                "has_cancer", "cancer_type", "age", "deprivation_index")
) %>%
  .[match(ID, .$eid), ] %>%
  mutate(
    sex  = factor(sex, levels = c(0, 1),  labels = c("Female", "Male")),
    pop  = factor(pop, levels = code_ancestry$coding,
                  labels = code_ancestry$meaning),
    has_cancer = as.logical(factor(has_cancer, levels = c(-3, -1, 0, 1),
                                   labels = c(NA, NA, FALSE, TRUE))),
    cancer_type = factor(cancer_type, levels = code_cancer$coding,
                         labels = code_cancer$meaning),
    date = (year - 1900) + (month - 0.5) / 12,
    year = NULL, month = NULL
  ) %>%
  as_tibble() %>%
  print()

saveRDS(df0, "pheno.rds")
