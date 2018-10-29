ID <- bigsnpr::snp_attach("UKB-merged_sub1.rds")$fam$sample.ID

library(dplyr)
df.CAD <- bigreadr::fread2(
  "ukb22544.csv",
  select = c("eid", "6150-0.0", "6150-0.1", "6150-0.2", "6150-0.3", "6150-1.0", "6150-1.1",
             "6150-1.2", "6150-1.3", "6150-2.0", "6150-2.1", "6150-2.2", "6150-2.3")
) %>%
  .[match(ID, .$eid), ]

ind <- which(rowSums(is.na(df.CAD)) != ncol(df.CAD))

df.CAD[ind, ]

lvl.100605 <- c(-7,-3,1,2,3,4)
lbl.100605 <- c("None of the above","Prefer not to answer","Heart attack","Angina","Stroke","High blood pressure")
df.CAD <- bigreadr::fread2(
  "ukb22544.csv",
  select = c("eid", "6150-0.0", "20002-0.0"),
  col.names = c("eid", "CAD", "non_cancer")
) %>%
  .[match(ID, .$eid), ] %>%
  mutate(CAD = factor(CAD, lvl.100605, lbl.100605))
df.CAD

table(df.CAD$CAD)
table(df.CAD$`20002-0.0`)

ind.CAD <- which(df.CAD$CAD %in% c("None of the above", "Heart attack"))

library(bigsnpr)
ukb <- snp_attach("UKB-merged_sub1.rds")
G <- ukb$genotypes
CHR <- ukb$map$chromosome
POS <- ukb$map$physical.pos
rm(ukb)

PCs <- predict(readRDS("UKB_SVD_auto.rds"))
y.CAD <- (df.CAD$CAD == "Heart attack")
system.time(
  gwas.CAD <- big_univLogReg(G, y.CAD[ind.CAD], ind.train = ind.CAD,
                             covar.train = PCs[ind.CAD, ], ncores = nb_cores())
) # 13h

plot(gwas.CAD)
snp_qq(gwas.CAD) + xlim(2, NA)
snp_manhattan(snp_gc(gwas.CAD), CHR, POS, npoints = 20e3) +
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = 2)

ind.HBP <- which(df.CAD$CAD %in% c("None of the above", "High blood pressure"))
y.HBP <- (df.CAD$CAD == "High blood pressure")
system.time(
  gwas.HBP <- big_univLogReg(G, y.HBP[ind.HBP], ind.train = ind.HBP,
                         covar.train = PCs[ind.HBP, ],
                         ncores = nb_cores())
) # 41h

plot(gwas.HBP)
snp_qq(gwas.HBP) + xlim(2, NA)
snp_manhattan(snp_gc(gwas.HBP), CHR, POS, npoints = 20e3) +
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = 2)
