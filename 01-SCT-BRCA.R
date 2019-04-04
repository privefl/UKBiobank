download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/MichailidouK_29059683_GCST004988/oncoarray_bcac_public_release_oct17%20(1).txt.gz",
              destfile = "sumstats_BRCA.txt.gz")
R.utils::gunzip("sumstats_BRCA.txt.gz")
library(bigreadr)
sumstats <- fread2("sumstats_BRCA.txt", na.strings = "NULL",
                   select = c("chr", "position_b37", "a0", "a1",
                              "bcac_onco_icogs_gwas_beta",
                              "bcac_onco_icogs_gwas_P1df"))
nrow(sumstats)  # 11,792,542

# augment dataset to match reverse alleles
sumstats2 <- sumstats
sumstats2$a1 <- sumstats$a2
sumstats2$a2 <- sumstats$a1
sumstats2$beta <- -sumstats$beta
sumstats2 <- rbind(sumstats, sumstats2)

# match variants with UKBB
library(doParallel)
NCORES <- 12
registerDoParallel(cl <- makeCluster(NCORES))
info_snp <- foreach(chr = 1:22, .packages = "dplyr") %dopar% {

  paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt") %>%
    bigreadr::fread2() %>%
    inner_join(
      sumstats2[sumstats2$chr == chr, ], ., na_matches = "never",
      by = c(position_b37 = "V3", a0 = "V4", a1 = "V5")
    ) %>%
    arrange(position_b37) %>%
    transmute(
      id    = paste(chr, position_b37, a0, a1, sep = "_"),
      beta  = bcac_onco_icogs_gwas_beta,
      lpval = -log10(bcac_onco_icogs_gwas_P1df),
      info  = V8
    ) %>%
    na.omit() %>%
    filter(lpval > 1 & info > 0.3)
}
stopCluster(cl)

list_snp_id <-  lapply(info_snp, function(df) df$id)
beta <-  unlist(lapply(info_snp, function(df) df$beta))
lpval <- unlist(lapply(info_snp, function(df) df$lpval))
info <-  unlist(lapply(info_snp, function(df) df$info))
length(beta)  # 1,693,570

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
df_cancer2 <- fread2(csv, select = c(paste0("40006-", 0:31, ".0"),
                                     paste0("40001-", 0:2, ".0"),
                                     paste0("40002-0.", 0:13),
                                     paste0("40002-1.", 0:13),
                                     paste0("40002-2.", 0:13)))
ind_BRCA <- sort(unique(unlist(c(
  lapply(df_cancer1,  function(x) which(x == 1002)),
  lapply(df_cancer2, function(x) which(substr(x, 1, 3) %in% c("C50", "D05")))
))))
table(df0$sex[ind_BRCA])
#     0     1
# 15411    87

y <- rep(NA, nrow(df0))
y[rowMeans(df_cancer0 == 0, na.rm = TRUE) == 1] <- 0
y[ind_BRCA] <- 1

sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian &
               df0$sex == 0 & !is.na(y))
length(sub)  # 169,853
table(y.sub <- y[sub])
#      0      1
# 159163  10690

system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_BRCA",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 5H


library(bigsnpr)
ukbb <- snp_attach("data/UKBB_BRCA.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 268 GB
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos

set.seed(1)
ind.train <- sort(sample(length(sub), 150e3))
ind.test <- setdiff(seq_along(sub), ind.train)

system.time(
  all_keep <- snp_grid_clumping(
    G, CHR, POS, lpS = lpval, ind.row = ind.train, infos.imp = info,
    grid.thr.imp = c(0.3, 0.6, 0.9, 0.95),
    grid.thr.r2 = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95),
    grid.base.size = c(50, 100, 200, 500),
    ncores = NCORES
  )
) # 5H
plot(lengths(all_keep[[1]]), pch = 20)

system.time(
  multi_PRS <- snp_grid_PRS(
    G, all_keep, betas = beta, lpS = lpval, ind.row = ind.train,
    grid.lpS.thr = seq_log(1, 0.999 * max(lpval), 50),
    backingfile = "data/UKBB_BRCA_scores", ncores = NCORES
  )
) # 2.4H

nPC <- 20
PC <- fread2(csv, select = paste0("22009-0.", 1:nPC),
             col.names = paste0("PC", 1:nPC))
PC_sub <- as.matrix(PC[sub, ])

system.time(
  final_mod <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train], covar.train = PC_sub[ind.train, ], ncores = NCORES
  )
) # 5H
mod <- final_mod$mod
plot(mod)
summary(mod)
# # A tibble: 5 x 6
#      alpha validation_loss intercept beta           nb_var message
#      <dbl>           <dbl>     <dbl> <list>          <int> <list>
#   1 0.0001           0.224     -2.14 <dbl [80,380]>  23407 <chr [10]>
#   2 0.001            0.224     -1.95 <dbl [80,380]>   8922 <chr [10]>
#   3 0.01             0.225     -1.84 <dbl [80,380]>   4271 <chr [10]>
#   4 0.1              0.225     -1.83 <dbl [80,380]>   3858 <chr [10]>
#   5 1                0.225     -1.82 <dbl [80,380]>   3608 <chr [10]>

new_beta <- final_mod$beta.G

length(ind <- which(new_beta != 0))  # 658,065
summary(new_beta)
summary(new_beta[which(sign(new_beta * beta) < 0)])

pred <- final_mod$intercept +
  big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind) +
  PC_sub[ind.test, ] %*% final_mod$beta.covar

AUCBoot(pred, y.sub[ind.test])  # 65.6 [64.0-67.1]


max(aucs <- apply(multi_PRS, 2, AUC, y.sub[ind.train]))  # 64.0

library(tidyverse)

ind2 <- sort(sample(ind, 10e3))
ggplot(data.frame(y = new_beta, x = beta)[ind, ]) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  geom_abline(slope = 0, intercept = 0, color = "blue") +
  geom_point(aes(x, y), size = 0.6) +
  theme_bigstatsr() +
  labs(x = "Effect sizes from GWAS", y = "Non-zero effect sizes from SCT")

grid2 <- attr(all_keep, "grid") %>%
  mutate(thr.lp = list(attr(multi_PRS, "grid.lpS.thr"))) %>%
  unnest() %>%
  mutate(auc = aucs)

grid2 %>%
  filter(thr.imp == 0.3, thr.r2 == 0.2, size == 500, thr.lp <= 8) %>%
  arrange(desc(auc)) %>%
  slice(1)
#   size thr.r2 thr.imp   thr.lp       auc
# 1  500    0.2     0.3 3.248033 0.6245214

grid2 %>% arrange(desc(auc)) %>% slice(1:10)
#     size thr.r2 thr.imp   thr.lp       auc
# 1   5000   0.10    0.90 3.248033 0.6397758
# 2   2000   0.10    0.90 3.248033 0.6395280
# 3  10000   0.05    0.90 3.248033 0.6395090
# 4  10000   0.05    0.95 3.248033 0.6394382
# 5   5000   0.10    0.90 4.110972 0.6393490
# 6   4000   0.05    0.90 3.248033 0.6393441
# 7   2000   0.10    0.90 4.110972 0.6393047
# 8   5000   0.10    0.90 3.654117 0.6392969
# 9  10000   0.05    0.90 4.110972 0.6392300
# 10  4000   0.05    0.90 4.110972 0.6391322

ggplot(grid2) +
  geom_point(aes(thr.lp, auc)) +
  facet_grid(thr.imp ~ thr.r2 + size) +
  scale_x_log10(limits = c(1, 10), breaks = c(1, 2, 5), minor_breaks = 1:10) +
  ylim(0.61, NA) +
  geom_hline(yintercept = max(aucs), color = "blue", linetype = 2) +
  geom_hline(yintercept = AUC(pred, y.sub[ind.test]), color = "red", linetype = 2) +
  theme_bigstatsr(size.rel = 0.7) +
  labs(x = "-log10(p-value) threshold (log scale)", y = "AUC")


ind3 <- which(abs(beta) > 2)
gwas <- big_univLogReg(G, y.sub[ind.train], ind.train, ind3)
cbind(beta = beta[ind3], maf = snp_MAF(G, ind.row = ind.train, ind.col = ind3),
      gwas, pval = predict(gwas, log10 = FALSE))
#       beta          maf        estim    std.err niter       score       pval
# 1  -2.6335 3.065333e-04  -1.62947420  0.9379711     6 -1.73723282 0.08234608
# 2  -2.9353 4.130000e-05  -3.49990545  5.8237872     8 -0.60096726 0.54786179
# 3  -2.2807 9.290000e-05  -2.74233804  2.8699579     8 -0.95553249 0.33930848
# 4  -2.1291 3.073667e-04   0.01915371  0.4714915     4  0.04062366 0.96759592
# 5  -2.1569 7.990000e-05  -0.16540450  1.2701422     4 -0.13022518 0.89638827
# 6  -2.0555 1.663333e-05 -30.22087424 31.2814228     9 -0.96609654 0.33399588
# 7  -2.5219 1.401333e-04  -1.62227389  1.6787426     6 -0.96636248 0.33386283
# 8  -2.0048 6.308333e-04  -0.17502723  0.3706654     4 -0.47219733 0.63678595
# 9  -2.3273 5.773333e-05  -0.45063707  1.5811693     5 -0.28500242 0.77564229
# 10 -2.9494 9.463333e-05  -0.33016997  1.1110281     5 -0.29717517 0.76633278
# 11 -2.5205 1.568667e-04  -0.26062116  0.7867266     5 -0.33127286 0.74043839
