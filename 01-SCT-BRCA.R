# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/MichailidouK_29059683_GCST004988/oncoarray_bcac_public_release_oct17%20(1).txt.gz",
#               destfile = "sumstats_BRCA.txt.gz")
# R.utils::gunzip("sumstats_BRCA.txt.gz")
library(bigreadr)
sumstats <- fread2("sumstats_BRCA.txt",
                   select = c("chr", "position_b37", "a0", "a1",
                              "bcac_onco_icogs_gwas_beta",
                              "bcac_onco_icogs_gwas_P1df"),
                   col.names = c("chr", "pos", "a0", "a1", "beta", "p"))
nrow(sumstats)  # 11,792,542

library(tidyverse)
sumstats <- sumstats %>%
  filter(p < 0.1, !paste(a0, a1) %in% c("A T", "T A", "C G", "G C"))
nrow(sumstats)  # 1,523,756
hist(sumstats$beta)
hist(sumstats$p)

# augment dataset to match reverse alleles
ACTG <- setNames(c("A", "C", "T", "G"), c("T", "G", "A", "C"))
sumstats2 <- sumstats
sumstats2$a0 <- ACTG[sumstats$a0]
sumstats2$a1 <- ACTG[sumstats$a1]
sumstats2 <- rbind(sumstats, na.omit(sumstats2))
sumstats3 <- sumstats2
sumstats3$a0 <- sumstats2$a1
sumstats3$a1 <- sumstats2$a0
sumstats3$beta <- -sumstats2$beta
sumstats3 <- rbind(sumstats2, sumstats3)

# match variants with UKBB
library(doParallel)
NCORES <- 12
registerDoParallel(cl <- makeCluster(NCORES))
info_snp <- foreach(chr = 1:22, .packages = "dplyr") %dopar% {

  paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt") %>%
    bigreadr::fread2() %>%
    inner_join(
      sumstats3[sumstats3$chr == chr, ], ., na_matches = "never",
      by = c(pos = "V3", a0 = "V4", a1 = "V5")
    ) %>%
    arrange(pos) %>%
    transmute(
      id    = paste(chr, pos, a0, a1, sep = "_"),
      beta  = beta,
      lpval = -log10(p),
      info  = V8
    ) %>%
    na.omit() %>%
    filter(info > 0.3)
}
stopCluster(cl)

list_snp_id <-  lapply(info_snp, function(df) df$id)
beta <-  unlist(lapply(info_snp, function(df) df$beta))
lpval <- unlist(lapply(info_snp, function(df) df$lpval))
info <-  unlist(lapply(info_snp, function(df) df$info))
length(beta)  # 1,460,446


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

system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_BRCA",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 3.8H


library(bigsnpr)
ukbb <- snp_attach("data/UKBB_BRCA.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 231 GB
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos

set.seed(1)
ind.train <- sort(sample(length(sub), 140e3))
ind.test <- setdiff(seq_along(sub), ind.train)

system.time(
  all_keep <- snp_grid_clumping(
    G, CHR, POS, lpS = lpval, ind.row = ind.train, infos.imp = info,
    grid.thr.imp = c(0.3, 0.6, 0.9, 0.95),
    grid.thr.r2 = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95),
    grid.base.size = c(50, 100, 200, 500),
    ncores = NCORES
  )
) # 2.8H
plot(lengths(all_keep[[1]]), pch = 20)

system.time(
  multi_PRS <- snp_grid_PRS(
    G, all_keep, betas = beta, lpS = lpval, ind.row = ind.train,
    n_thr_lpS = 50, backingfile = "data/UKBB_BRCA_scores", ncores = NCORES
  )
) # 50 min

system.time(
  final_mod <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train], ncores = NCORES, n.abort = 5
  )
) # 100 min
mod <- final_mod$mod
plot(mod)
summary(mod)
# # A tibble: 3 x 6
#      alpha validation_loss intercept beta           nb_var message
#      <dbl>           <dbl>     <dbl> <list>          <int> <list>
#   1 0.0001           0.239     -2.15 <dbl [80,416]>  25341 <chr [10]>
#   2 0.01             0.239     -2.07 <dbl [80,416]>   5237 <chr [10]>
#   3 1                0.239     -2.06 <dbl [80,416]>   3799 <chr [10]>

new_beta <- final_mod$beta.G

length(ind <- which(new_beta != 0))  # 582,878
summary(new_beta)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -1.615e-01  0.000e+00  0.000e+00 -9.490e-06  0.000e+00  2.485e-01
summary(new_beta[which(sign(new_beta * beta) < 0)])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -2.177e-02 -1.857e-05  2.800e-07  2.849e-05  3.421e-05  2.477e-02

pred <- final_mod$intercept +
  big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind)

AUCBoot(pred, y.sub[ind.test])  # 65.9 [64.7-67.1]


# Save
save(all_keep, multi_PRS, final_mod, file = "data/res_BRCA.RData")


library(tidyverse)

ind2 <- sort(sample(ind, 10e3))
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
# 1  500    0.2     0.3  14 3.248484 0.6214921

ind.keep <- unlist(map(all_keep, std_prs$num))
# Verif on training set
AUC(
  snp_PRS(G, beta[ind.keep], ind.test = ind.train, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.train]
)
# Eval on test set
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.test]
) # 62.4 [61.2-63.7]
sum(lpval[ind.keep] > std_prs$thr.lp)  # 5903

max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
#    size thr.r2 thr.imp num   thr.lp       auc
# 1  5000    0.1    0.95  96 3.654696 0.6349646
# 2  2000    0.1    0.95  95 3.654696 0.6348287
# 3  2500    0.2    0.95 100 3.654696 0.6346257
# 4  2000    0.1    0.95  95 4.111703 0.6346209
# 5  5000    0.1    0.95  96 4.111703 0.6346123
# 6  5000    0.1    0.90  68 4.111703 0.6343516
# 7  1000    0.2    0.95  99 3.654696 0.6341778
# 8  2500    0.2    0.90  72 3.248484 0.6340949
# 9  2000    0.1    0.90  67 4.111703 0.6339956
# 10 2500    0.2    0.95 100 3.248484 0.6339942

ind.keep <- unlist(map(all_keep, max_prs$num))
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp),
  y.sub[ind.test]
) # 63.5 [62.3-64.8]
sum(lpval[ind.keep] > max_prs$thr.lp)  # 1985

ggplot(grid2) +
  geom_point(aes(thr.lp, auc)) +
  facet_grid(thr.imp ~ thr.r2 + size) +
  scale_x_log10(limits = c(1, 8), breaks = c(1, 2, 5), minor_breaks = 1:10) +
  ylim(0.61, NA) +
  theme_bigstatsr(size.rel = 0.7) +
  labs(x = "-log10(p-value) threshold (log scale)", y = "AUC")

# Check large effects
ind3 <- intersect(which(abs(beta) > 2), ind)
gwas <- big_univLogReg(G, y.sub[ind.train], ind.train, ind3)
cbind(beta_GWAS = beta[ind3], pval_GWAS = 10^(-lpval[ind3]),
      af = big_scale()(G, ind.row = ind.train, ind.col = ind3)$center / 2,
      gwas, pval = predict(gwas, log10 = FALSE))
#   beta_GWAS pval_GWAS           af       estim    std.err niter      score      pval
# 1   -2.6335  0.004775 3.189286e-04  -0.9875785  0.7397290     5 -1.3350545 0.1818585
# 2   -2.9353  0.007172 5.000000e-05  -3.8826748  6.4837872     9 -0.5988282 0.5492874
# 3   -2.2807  0.011160 8.300000e-05  -2.7064779  2.9339564     8 -0.9224670 0.3562850
# 4   -2.1569  0.004768 8.153571e-05  -0.2359545  1.2740710     4 -0.1851973 0.8530743
# 5   -2.0555  0.004252 1.500000e-05 -32.1185585 31.5562989     9 -1.0178177 0.3087646
# 6   -2.5219  0.003705 1.543571e-04  -1.9503401  1.6900563     6 -1.1540089 0.2484965
# 7   -2.0048  0.004386 6.222143e-04  -0.2669400  0.3894919     5 -0.6853544 0.4931203
