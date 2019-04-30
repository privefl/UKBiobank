# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/WrayNR_29700475_GCST005839/MDD2018_ex23andMe.gz",
#               destfile = "sumstats_MDD.txt.gz")
# R.utils::gunzip("sumstats_MDD.txt.gz")
library(bigreadr)
sumstats <- fread2("sumstats_MDD.txt", fill = TRUE,
                   select = c("CHR", "BP", "A1", "A2", "OR", "P"),
                   col.names = c("chr", "pos", "a1", "a0", "or", "p"))
nrow(sumstats)  # 13,554,550
library(dplyr)
sumstats <- sumstats %>%
  filter(p < 0.1) %>%
  mutate(beta = log(or), or = NULL, chr = as.integer(chr))
nrow(sumstats)  # 1,610,995
hist(sumstats$beta)
hist(sumstats$p)


info_snp_UKBB <- rbind_df(lapply(1:22, function(chr) {
  file <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
  df <- fread2(file, select = c(3:5, 8), col.names = c("pos", "a0", "a1", "info"))
  cbind.data.frame(chr = chr, df)
}))
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB)
# 1,610,995 variants in summary statistics.
# 224,231 ambiguous SNPs have been removed.
# 1,201,894 variants have been matched; 102 were flipped and 521,969 were reversed.
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB, strand_flip = FALSE)
# 1,420,221 variants have been matched; 0 were flipped and 629,580 were reversed.
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
df0 <- fread2(csv, select = c("eid", "22006-0.0"),
              col.names = c("eid", "is_caucasian"))
ind.indiv <- match(df0$eid, sample$ID_2)
# related samples
df0$is_rel2 <- df0$eid %in% fread2("ukb25589_rel_s488346.dat")$ID2


df_illness <- fread2(csv, select = c(paste0("20002-0.", 0:28),
                                     paste0("20002-1.", 0:28),
                                     paste0("20002-2.", 0:28)))
df_ICD10 <- fread2(csv, select = c(paste0("40001-", 0:2, ".0"),
                                   paste0("40002-0.", 0:13),
                                   paste0("40002-1.", 0:13),
                                   paste0("40002-2.", 0:13),
                                   paste0("41202-0.", 0:379),
                                   paste0("41204-0.", 0:434)))
ind_MDD <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1286)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) %in% c("F32", "F33")))
))))
ind_psy <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x %in% 1286:1291)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 1) == "F"))
))))
y <- rep(0, nrow(df0)); y[ind_psy] <- NA; y[ind_MDD] <- 1

sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian & !is.na(y))
length(sub)  # 317,824
table(y.sub <- y[sub])
#      0      1
# 291690  26134


NCORES <- 12
system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_MDD",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 4H


library(bigsnpr)
ukbb <- snp_attach("data/UKBB_MDD.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 418 GB
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos

set.seed(1)
ind.train <- sort(sample(length(sub), 250e3))
ind.test <- setdiff(seq_along(sub), ind.train)

system.time(
  all_keep <- snp_grid_clumping(
    G, CHR, POS, lpS = lpval, ind.row = ind.train, infos.imp = info,
    grid.thr.imp = c(0.3, 0.6, 0.9, 0.95),
    grid.thr.r2 = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95),
    grid.base.size = c(50, 100, 200, 500),
    ncores = 6
  )
) # 13.4H

system.time(
  multi_PRS <- snp_grid_PRS(
    G, all_keep, betas = beta, lpS = lpval, ind.row = ind.train,
    n_thr_lpS = 50, backingfile = "data/UKBB_MDD_scores", ncores = NCORES
  )
) # 4.5H

system.time(
  final_mod <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train], ncores = NCORES, n.abort = 3
  )
) # 3.6H
mod <- final_mod$mod
plot(mod)
summary(mod)
# # A tibble: 5 x 6
#      alpha validation_loss intercept beta           nb_var message
#      <dbl>           <dbl>     <dbl> <list>          <int> <list>
#   1 0.0001           0.278      4.28 <dbl [58,380]>  26371 <chr [10]>
#   2 0.001            0.278      4.06 <dbl [58,380]>   9030 <chr [10]>
#   3 0.01             0.278      3.94 <dbl [58,380]>   4944 <chr [10]>
#   4 0.1              0.278      3.83 <dbl [58,380]>   3107 <chr [10]>
#   5 1                0.278      3.82 <dbl [58,380]>   2807 <chr [10]>

# Save
save(all_keep, multi_PRS, final_mod, file = "data/res_MDD.RData")

new_beta <- final_mod$beta.G

length(ind <- which(new_beta != 0))  # 504,886
summary(new_beta)
summary(new_beta[which(sign(new_beta * beta) < 0)])

pred <- final_mod$intercept +
  big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind)

AUCBoot(pred, y.sub[ind.test])  # 60.6 [59.8-61.3]


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
# 1  500    0.2     0.3  14 4.580883 0.7489942

ind.keep <- unlist(map(all_keep, std_prs$num))
# Verif on training set
AUC(
  snp_PRS(G, beta[ind.keep], ind.test = ind.train, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.train]
) # 0.7489942
# Eval on test set
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.test]
) # 75.4 [72.4-78.4]
sum(lpval[ind.keep] > std_prs$thr.lp)  # 1112

max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
#     size thr.r2 thr.imp num   thr.lp       auc
# 1  10000   0.01     0.9  58 4.580883 0.7730035
# 2  20000   0.01     0.9  59 4.580883 0.7730035
# 3  50000   0.01     0.9  60 4.580883 0.7730035
# 4   5000   0.01     0.9  57 4.580883 0.7722518
# 5   5000   0.01     0.9  57 5.149850 0.7716009
# 6  10000   0.01     0.9  58 5.149850 0.7713099
# 7  20000   0.01     0.9  59 5.149850 0.7713099
# 8  50000   0.01     0.9  60 5.149850 0.7713099
# 9   5000   0.01     0.9  57 7.316962 0.7691121
# 10  5000   0.01     0.9  57 5.789486 0.7686678

ind.keep <- unlist(map(all_keep, max_prs$num))
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp),
  y.sub[ind.test]
) # 76.9 [73.9-79.7]
sum(lpval[ind.keep] > max_prs$thr.lp)  # 267

ggplot(grid2) +
  geom_point(aes(thr.lp, auc)) +
  facet_grid(thr.imp ~ thr.r2 + size) +
  scale_x_log10(limits = c(1, 4), breaks = c(1, 2, 4), minor_breaks = 1:10) +
  ylim(0.52, NA) +
  theme_bigstatsr(size.rel = 0.7) +
  labs(x = "-log10(p-value) threshold (log scale)", y = "AUC")


ind3 <- intersect(which(abs(beta) > 8), ind)
gwas <- big_univLogReg(G, y.sub[ind.train], ind.train, ind3)
cbind(beta_GWAS = beta[ind3], pval_GWAS = 10^(-lpval[ind3]),
      af = big_scale()(G, ind.row = ind.train, ind.col = ind3)$center / 2,
      gwas, pval = predict(gwas, log10 = FALSE))
#    beta_GWAS pval_GWAS         af         estim      std.err niter       score      pval
# 1  -8.613200  0.034370 0.00091638 -2.924241e-01 2.262117e-01     5 -1.29270060 0.1961146
# 2  -8.111728  0.059710 0.00000916 -2.381605e+01 1.816082e+01     7 -1.31139708 0.1897236
# 3  -8.111728  0.014470 0.00000934  2.295782e+00 1.396240e+00     7  1.64426042 0.1001224
# 4  -8.383900  0.048800 0.00000148 -8.485884e+02 3.640731e+03    NA -0.23308185 0.8156979
# 5  -8.207900  0.014250 0.00055936  1.472254e-01 2.294435e-01     4  0.64166271 0.5210922
# 6   9.210340  0.084280 0.00438654  1.154473e-01 7.603308e-02     4  1.51838267 0.1289180
# 7   8.266100  0.055860 0.00008082 -2.605212e-01 9.480593e-01     5 -0.27479425 0.7834743
# 8  -9.210340  0.055620 0.00100462 -1.585354e-01 1.891760e-01     4 -0.83803121 0.4020132
# 9  -8.517193  0.080990 0.00004392 -1.320109e+00 2.056631e+00     5 -0.64187903 0.5209517
# 10 -8.628700  0.041030 0.00399966 -8.353803e-02 9.596195e-02     4 -0.87053289 0.3840093
# 11 -8.866900  0.021650 0.00295036 -8.937422e-02 1.025526e-01     4 -0.87149630 0.3834832
# 12 -9.210340  0.019190 0.00000814 -1.423166e+01 2.772481e+01     9 -0.51331861 0.6077285
# 13  8.262000  0.026560 0.00062238  1.921579e-01 2.894862e-01     4  0.66378948 0.5068251
# 14 -8.111728  0.032100 0.00000064 -7.188178e+02 3.685251e+03    NA -0.19505263 0.8453517
# 15  8.121100  0.046650 0.00306802 -7.772223e-03 9.983262e-02     4 -0.07785254 0.9379453
# 16  8.517193  0.050130 0.00260188  1.219611e-01 1.045409e-01     4  1.16663526 0.2433577
# 17 -8.517193  0.028230 0.00004034  9.229559e-01 7.870077e-01     6  1.17274068 0.2408998
# 18  8.257000  0.049410 0.00015974  5.142279e-01 3.565112e-01     5  1.44238928 0.1491926
# 19 -8.414900  0.018160 0.00204672  5.345157e-02 1.168808e-01     4  0.45731700 0.6474432
# 20 -9.210340  0.074780 0.00040210 -4.539010e-03 3.306237e-01     3 -0.01372863 0.9890465
# 21 -8.517193  0.006523 0.00436096  3.555955e-02 8.650870e-02     4  0.41105177 0.6810346
# 22 -8.456700  0.058450 0.00566180  8.875999e-02 6.855675e-02     4  1.29469376 0.1954259
# 23  8.111728  0.005863 0.99995660 -3.522933e-02 1.192194e+00     4 -0.02955000 0.9764259
# 24 -8.111728  0.099890 0.00017708 -5.316269e-01 5.344207e-01     5 -0.99477220 0.3198471
# 25 -8.111728  0.007289 0.00020222  2.313557e-01 3.593534e-01     5  0.64381092 0.5196980
# 26 -8.027300  0.096170 0.00371254 -1.489926e-01 1.030276e-01     4 -1.44614217 0.1481373
# 27 -8.036600  0.081220 0.00010902 -8.537821e-01 7.448401e-01     6 -1.14626220 0.2516867
# 28  8.654100  0.003695 0.00280562  1.649862e-03 1.040120e-01     3  0.01586224 0.9873443
# 29  8.517193  0.052420 0.00002000  1.868926e+00 1.337631e+00     6  1.39719052 0.1623563
# 30  9.054400  0.035150 0.00005730  9.110597e-02 7.857375e-01     4  0.11594963 0.9076925
# 31 -9.172600  0.099850 0.00104842 -3.276111e-01 2.033252e-01     5 -1.61126649 0.1071217
# 32 -9.210340  0.035940 0.00000494 -1.918516e+01 2.793304e+01     8 -0.68682673 0.4921919
# 33 -9.210340  0.035610 0.00000602 -1.531764e+01 1.851051e+01     7 -0.82751077 0.4079476
