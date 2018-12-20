library(bigreadr)
library(dplyr, warn.conflicts = FALSE)

## Self-reported ancestry (https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=1001)
code_ancestry <- fread2("coding1001.tsv")

## Cancer type (https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=3)
code_cancer <- fread2("coding3.tsv")

csv <- "ukb22544.csv"
df0 <- fread2(
  csv,
  select = c("eid", "50-0.0", "34-0.0", "52-0.0", "22001-0.0", "21000-0.0",
             "22006-0.0", "2453-0.0", "20001-0.0", "21022-0.0", "189-0.0"),
  col.names = c("eid", "height", "year", "month", "sex", "pop", "is_caucasian",
                "has_cancer", "cancer_type", "age", "deprivation_index")
) %>%
  mutate(
    sex  = factor(sex, levels = c(0, 1),  labels = c("Female", "Male")),
    pop  = factor(pop, levels = code_ancestry$coding,
                  labels = code_ancestry$meaning),
    is_caucasian = as.logical(is_caucasian),
    has_cancer = as.logical(factor(has_cancer, levels = c(-3, -1, 0, 1),
                                   labels = c(NA, NA, FALSE, TRUE))),
    cancer_type = factor(cancer_type, levels = code_cancer$coding,
                         labels = code_cancer$meaning),
    date = (year - 1900) + (month - 0.5) / 12,
    year = NULL, month = NULL
  )
rownames(df0) <- df0$eid
head(df0)

# Principal Components
nPC <- 20
PC <- fread2(
  csv,
  select = paste0("22009-0.", 1:nPC),
  col.names = paste0("PC", 1:nPC)
)
rownames(PC) <- df0$eid
head(PC)
plot(PC2 ~ PC1, data = PC[sample(500e3, 20e3), ], pch = 20)

# Relatedness
rel <- fread2("ukb25589_rel_s488346.dat")
sum(ind.rm <- df0$eid %in% rel[rel$Kinship > 0.08, "ID2"])  ## 38,431

saveRDS(df0[!ind.rm, ], "pheno.rds")
saveRDS(PC[!ind.rm, ],  "PC.rds")
