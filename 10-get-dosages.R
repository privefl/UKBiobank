sumstats_joined <- readRDS("breast-cancer/sumstats_all_chr.rds")
list_snp_id <- split(
  with(sumstats_joined, paste(chr, position_b37, a0, a1, sep = "_")),
  factor(sumstats_joined$chr, levels = 1:22, ordered = TRUE)
)
rm(sumstats_joined)

# system("./ukbgene imp -c1 -m ")
sample <- bigreadr::fread2("ukb25589_imp_chr1_v3_s487327.sample")
str(sample)
(N <- readBin("ukb_imp_chr1_v3.bgen", what = 1L, size = 4, n = 4)[4])
sample <- sample[-1, ]
sample$sex <- as.numeric(sample$sex)
stopifnot(nrow(sample) == N)

df0 <- readRDS("pheno.rds")
ind.indiv <- match(df0$eid[which(df0$is_caucasian)], sample$ID_2)
sum(is.na(ind.indiv))  # 766 -> withdrawals?
ind.indiv <- sort(na.omit(ind.indiv))

library(bigsnpr)
system.time(
  rds <- snp_readBGEN(
    bgenfiles = glue::glue("ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "UKB_imp_BC",
    ind_row = ind.indiv,
    bgi_dir = "ukb_imp_bgi",
    ncores = min(nb_cores(), 12)
  )
) # <10H

ukb_imp <- snp_attach(rds)
G <- ukb_imp$genotypes
dim(G)
G[, 1:5]

df_imp <- df0[match(sample$ID_2[ind.indiv], df0$eid), ]
table(df_imp$sex, sample$sex[ind.indiv])  ## verif -> OK
ukb_imp$fam <- df_imp
ukb_imp <- snp_save(ukb_imp)

# Verif with GWAS on height -> order seems OK
ind.nona <- which(!is.na(ukb_imp$fam$height))
system.time(
  gwas <- big_univLinReg(G, y.train = sample(ukb_imp$fam$height[ind.nona]),
                         ind.train = ind.nona, ind.col = seq_len(1e5),
                         ncores = 14)
)
plot(gwas)
library(ggplot2)
ind <- which(is.na(gwas$score))
snp_qq(gwas[-ind, ]) + xlim(1, NA)
