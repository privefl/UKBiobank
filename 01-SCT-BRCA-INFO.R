# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/MichailidouK_29059683_GCST004988/oncoarray_bcac_public_release_oct17%20(1).txt.gz",
#               destfile = "sumstats_BRCA.txt.gz")
# R.utils::gunzip("sumstats_BRCA.txt.gz")
library(bigreadr)
sumstats <- fread2("sumstats_BRCA.txt", na.strings = "NULL",
                   select = c("chr", "position_b37", "a0", "a1",
                              "bcac_onco_icogs_gwas_beta",
                              "bcac_onco_icogs_gwas_P1df",
                              "bcac_onco2_r2", "bcac_icogs2_r2"),
                   col.names = c("chr", "pos", "a0", "a1", "beta", "p", "INFO1", "INFO2"))
nrow(sumstats)  # 11,792,542

library(dplyr)
sumstats <- sumstats %>%
  mutate(chr = as.integer(chr),
         INFO = (INFO1 + INFO2) / 2, INFO1 = NULL, INFO2 = NULL,
         INFO = ifelse(INFO > 1 | INFO < 0.3, NA, INFO),
         beta = beta * INFO, p = p ** INFO, INFO = NULL) %>%
  filter(p < 0.1) %>%
  na.omit()
nrow(sumstats)  # 1,323,199
hist(sumstats$beta)
hist(sumstats$p)

info_snp_UKBB <- rbind_df(lapply(1:22, function(chr) {
  file <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
  df <- fread2(file, select = c(3:5, 8), col.names = c("pos", "a0", "a1", "info"))
  cbind.data.frame(chr = chr, df)
}))
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB)
# 1,323,199 variants in summary statistics.
# 183,778 ambiguous SNPs have been removed.
# 1,102,231 variants have been matched; 0 were flipped and 170 were reversed.
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB, strand_flip = FALSE)
# 1,281,721 variants have been matched; 0 were flipped and 170 were reversed.
info_snp <- subset(na.omit(info_snp), info > 0.3)

list_snp_id <- with(info_snp, split(paste(chr, pos, a0, a1, sep = "_"),
                                    factor(chr, levels = 1:22)))
beta <- info_snp$beta
lpval <- -log10(info_snp$p)
info <- info_snp$info


# subset samples
sample <- fread2("ukb25589_imp_chr1_v3_s487327.sample")
str(sample)
(N <- readBin("data/ukb_imp_chr1_v3.bgen", what = 1L, size = 4, n = 4)[4])
sample <- sample[-1, ]
stopifnot(nrow(sample) == N)

csv <- "ukb22544.csv"
df0 <- fread2(
  csv,
  select = c("eid", "22001-0.0", "22006-0.0"),
  col.names = c("eid", "sex", "is_caucasian")
)
ind.indiv <- match(df0$eid, sample$ID_2)
# related samples
df0$is_rel2 <- df0$eid %in% fread2("ukb25589_rel_s488346.dat")$ID2

df_cancer0 <- fread2(csv, select = paste0("2453-", 0:2, ".0"))
df_cancer1 <- fread2(csv, select = c(paste0("20001-0.", 0:5),
                                     paste0("20001-1.", 0:5),
                                     paste0("20001-2.", 0:5)))
df_cancer2 <- fread2(csv, select = c(paste0("40001-", 0:2, ".0"),
                                     paste0("40002-0.", 0:13),
                                     paste0("40002-1.", 0:13),
                                     paste0("40002-2.", 0:13),
                                     paste0("40006-", 0:31, ".0"),
                                     paste0("41202-0.", 0:379),
                                     paste0("41204-0.", 0:434)))
ind_BRCA <- sort(unique(unlist(c(
  lapply(df_cancer1,  function(x) which(x == 1002)),
  lapply(df_cancer2, function(x) which(substr(x, 1, 3) %in% c("C50", "D05")))
))))
table(df0$sex[ind_BRCA])
#     0     1
# 16681   116

y <- rep(NA, nrow(df0))
y[rowMeans(df_cancer0 == 0, na.rm = TRUE) == 1] <- 0
y[ind_BRCA] <- 1

sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian &
               df0$sex == 0 & !is.na(y))
length(sub)  # 169,969
table(y.sub <- y[sub])
#      0      1
# 158391  11578

NCORES <- 12
system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_BRCA_INFO",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 3.4H


library(bigsnpr)
ukbb <- snp_attach("data/UKBB_BRCA_INFO.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 203 GB
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos

set.seed(1); ind.train <- sort(sample(length(sub), 150e3))
# set.seed(2); ind.train <- c(sample(which(y.sub == 0), 2000), sample(which(y.sub == 1), 500))
ind.test <- setdiff(seq_along(sub), ind.train)

system.time(
  all_keep <- snp_grid_clumping(
    G, CHR, POS, lpS = lpval, ind.row = ind.train, infos.imp = info,
    grid.thr.imp = c(0.3, 0.6, 0.9, 0.95),
    grid.thr.r2 = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95),
    grid.base.size = c(50, 100, 200, 500),
    ncores = 6
  )
) # 2H

system.time(
  multi_PRS <- snp_grid_PRS(
    G, all_keep, betas = beta, lpS = lpval, ind.row = ind.train,
    n_thr_lpS = 50, backingfile = "data/UKBB_BRCA_INFO_scores", ncores = NCORES
  )
) # 1H

system.time(
  final_mod <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train], ncores = NCORES,
  )
) # 3H
mod <- final_mod$mod
plot(mod)
summary(mod)
# # A tibble: 3 x 6
#    alpha validation_loss intercept beta           nb_var message
#    <dbl>           <dbl>     <dbl> <list>          <int> <list>
# 1 0.0001           0.238     -1.91 <dbl [80,192]>  28320 <chr [10]>
# 2 0.01             0.239     -1.82 <dbl [80,192]>   6174 <chr [10]>
# 3 1                0.239     -1.79 <dbl [80,192]>   5391 <chr [10]>

new_beta <- final_mod$beta.G

length(ind <- which(new_beta != 0))  # 402,965
summary(new_beta)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -1.743e-01  0.000e+00  0.000e+00 -7.970e-06  0.000e+00  3.160e-01
summary(new_beta[which(sign(new_beta * beta) < 0)])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -5.629e-02 -5.340e-06  2.400e-07  1.451e-05  9.860e-06  3.087e-02

pred <- final_mod$intercept +
  big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind)

AUCBoot(pred, y.sub[ind.test])
# Before: 65.9 [64.4-67.4] / 62.9 [62.4-63.5]
# After:  66.0 [64.5-67.5]

library(tidyverse)

ggplot(data.frame(y = new_beta, x = beta)[ind, ]) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  geom_abline(slope = 0, intercept = 0, color = "blue") +
  geom_point(aes(x, y), size = 0.6) +
  theme_bigstatsr() +
  labs(x = "Effect sizes from GWAS", y = "Non-zero effect sizes from SCT")

grid2 <- attr(all_keep, "grid") %>%
  mutate(thr.lp = list(attr(multi_PRS, "grid.lpS.thr")), num = row_number()) %>%
  unnest()
s <- nrow(grid2)
grid2$auc <- big_apply(multi_PRS, a.FUN = function(X, ind, s, y.train) {
  # Sum over all chromosomes, for the same C+T parameters
  single_PRS <- rowSums(X[, ind + s * (0:21)])
  bigstatsr::AUC(single_PRS, y.train)
}, ind = 1:s, s = s, y.train = y.sub[ind.train],
a.combine = 'c', block.size = 1, ncores = NCORES)

std_prs <- grid2 %>%
  filter(thr.imp == 0.3, thr.r2 == 0.2, size == 500, thr.lp <= 8) %>%
  arrange(desc(auc)) %>%
  slice(1) %>%
  print()
#   size thr.r2 thr.imp num   thr.lp       auc
# 1  500    0.2     0.3  14 3.248484 0.6208515

ind.keep <- unlist(map(all_keep, std_prs$num))
# Verif on training set
AUC(
  snp_PRS(G, beta[ind.keep], ind.test = ind.train, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.train]
) # 0.6208515
# Eval on test set
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.test]
) # 62.1 [60.5-63.6] / 62.2 [61.6-62.7] -> 63.2 [61.7-64.8]
sum(lpval[ind.keep] > std_prs$thr.lp)  # 6256

max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
#    size thr.r2 thr.imp num   thr.lp       auc
# 1  2500    0.2    0.95 100 3.654696 0.6351843
# 2  5000    0.1    0.95  96 3.654696 0.6349136
# 3  1000    0.2    0.95  99 3.654696 0.6348574
# 4  2000    0.1    0.95  95 3.654696 0.6348562
# 5  2500    0.2    0.95 100 3.248484 0.6347538
# 6  2500    0.2    0.90  72 3.248484 0.6345844
# 7  2000    0.1    0.95  95 4.111703 0.6345656
# 8  5000    0.1    0.95  96 4.111703 0.6345364
# 9  5000    0.1    0.90  68 4.111703 0.6344550
# 10 2500    0.2    0.90  72 3.654696 0.6342989

ind.keep <- unlist(map(all_keep, max_prs$num))
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp),
  y.sub[ind.test]
) # 63.3 [61.7-64.8] / 63.4 [62.8-63.9] -> 63.7 [62.2-65.2]
sum(lpval[ind.keep] > max_prs$thr.lp)  # 2572

ggplot(grid2) +
  geom_point(aes(thr.lp, auc)) +
  facet_grid(thr.imp ~ thr.r2 + size) +
  scale_x_log10(limits = c(1, 8), breaks = c(1, 2, 5), minor_breaks = 1:10) +
  ylim(0.61, NA) +
  theme_bigstatsr(size.rel = 0.7) +
  labs(x = "-log10(p-value) threshold (log scale)", y = "AUC")


# lassosum
library(lassosum)
ukbb$map <- dplyr::mutate(ukbb$map, genetic.dist = 0, rsid = NULL,
                          chromosome = as.integer(chromosome))
ukbb$fam <- snp_fake(nrow(G), 1)$fam
bed <- snp_writeBed(ukbb, bedfile = "data/UKBB_BRCA_INFO_lassosum.bed",
                    ind.row = ind.train)
library(doParallel)
registerDoParallel(cl <- makeCluster(NCORES))
system.time(
  out <- lassosum.pipeline(
    cor = p2cor(10^-lpval, n = 137045 + 119078, sign = beta),
    snp = ukbb$map$marker.ID,
    A1 = ukbb$map$allele1,
    test.bfile = "data/UKBB_BRCA_INFO_lassosum",
    LDblocks = "EUR.hg19",
    cluster = cl,
    sample = 20e3,
    exclude.ambiguous = FALSE
  )
) # 79 min
stopCluster(cl)

v <- validate(out, pheno = y.sub[ind.train], validate.function = AUC)
length(ind <- which(v$best.beta != 0))  # 322,003 -> 89,877
pred_lassosum <- big_prodVec(G, v$best.beta[ind], ind.row = ind.test, ind.col = ind)
AUCBoot(pred_lassosum, y.sub[ind.test])
# Before: 57.9 [56.3-59.5] / 57.8 [57.3-58.4]
# After:  56.5 [54.9-58.1]
