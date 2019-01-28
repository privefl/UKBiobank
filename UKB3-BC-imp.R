# ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/MichailidouK_29059683_GCST004988
sumstats <- bigreadr::fread2("breast-cancer/oncoarray_bcac_public_release_oct17.txt",
                             na.strings = "NULL")
sumstats_sub <- subset(sumstats, bcac_onco_icogs_gwas_P1df < 0.1, 3:6)
nrow(sumstats_sub) # 1,768,354
rm(sumstats)

library(doParallel)
registerDoParallel(cl <- makeCluster(12))
list_snp_id <- foreach(chr = 1:22) %dopar% {
  # cat("Processing chromosome", chr, "..\n")
  file <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
  infos_chr <- bigreadr::fread2(file, showProgress = FALSE, nThread = 1)
  ind <- which(sumstats_sub$chr == chr)
  sumstats_joined <- dplyr::inner_join(
    sumstats_sub[ind, ], infos_chr,
    by = c(position_b37 = "V3", a0 = "V4", a1 = "V5"),
    na_matches = "never"
  )
  with(sumstats_joined, paste(chr, position_b37, a0, a1, sep = "_"))
}
stopCluster(cl)

sum(lengths(list_snp_id))  # 1,698,015

sample <- bigreadr::fread2("ukb25589_imp_chr1_v3_s487327.sample")
str(sample)
(N <- readBin("data/ukb_imp_chr1_v3.bgen", what = 1L, size = 4, n = 4)[4])
sample <- sample[-1, ]
stopifnot(nrow(sample) == N)

df0 <- readRDS("pheno.rds")
ind.indiv <- match(df0$eid, sample$ID_2)
y <- ifelse(df0$has_breast_cancer, 1, ifelse(df0$has_cancer, NA, 0))

sub <- which(df0$sex == "Female" & df0$is_caucasian & !is.na(y) & !is.na(ind.indiv))
length(sub)  # 188,628
table(y[sub])
#      0      1
# 180094   8534


system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKB_imp_BC2",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = 10
  )
) # 5H


set.seed(1)
ind.train <- sort(sample(length(sub), 150e3))
ind.test <- setdiff(seq_along(sub), ind.train)
table(y[sub][ind.train])
#      0      1
# 143232   6768

library(bigsnpr)
ukb <- snp_attach("data/UKB_imp_BC2.rds")
G <- ukb$genotypes
CHR <- as.numeric(ukb$map$chromosome)
dim(G) # 188,628 x 1,698,015
file.size(G$backingfile) / 1024^3  # 298 GB

PC <- as.matrix(readRDS("PC.rds"))
plot(PC[sample(nrow(PC), 10e3), ], pch = 20)  # verif: all
plot(PC[sample(sub, 10e3), ], pch = 20)       # verif: caucasian only

system.time(
  mod <- big_spLogReg(G, y[sub][ind.train], ind.train,
                      covar.train = PC[sub[ind.train], ],
                      alphas = 10^(-(0:4)), return.all = TRUE,
                      ncores = nb_cores())
) # 12 min

preds <- sapply(mod, function(mods) {
  rowMeans(sapply(mods, function(mod) {
    predict(mod, G, ind.test, covar.row = PC[sub[ind.test], ])
  }))
})
apply(preds, 2, AUC, y[sub][ind.test])
# 0.6075806 0.5992255 0.5985956 0.5973514 0.5974256

sapply(mod, function(mods) {
  sum(rowMeans(sapply(mods, function(mod) mod$beta)) > 0)
})
# 14595  2435  1004   665   640

plot(mod[[1]][[1]]$loss.val, pch = 20)


sum(rowMeans(sapply(mod, function(x) x$beta)) > 0) # 397
pred <- predict(mod, G, ind.test, covar.row = PC[sub[ind.test], ])
AUC(pred, y[sub][ind.test])  # 58.7

