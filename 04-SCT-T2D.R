# Download from http://diagram-consortium.org/downloads.html
# DIAGRAM 1000G GWAS meta-analysis Stage 1 Summary statistics
# Published in Scott et al (2017)
# unzip("METAANALYSIS_DIAGRAM_SE1.zip")
library(bigreadr)
sumstats <- fread2("METAANALYSIS_DIAGRAM_SE1.txt", select = c(1:4, 6))
sumstats <- tidyr::separate(sumstats, "Chr:Position", c("chr", "pos"), convert = TRUE)
names(sumstats) <- c("chr", "pos", "a1", "a0", "beta", "p")
nrow(sumstats)  # 12,056,346
hist(sumstats$beta)
hist(sumstats$p)

info_snp_UKBB <- rbind_df(lapply(1:22, function(chr) {

  file <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")

  df <- fread2(file, select = c(3:5, 8), col.names = c("pos", "a0", "a1", "info"))

  cbind.data.frame(chr = chr, df)

}))
info_snp <- bigsnpr::snp_match(subset(sumstats, p < 0.1), info_snp_UKBB)
# 1,408,672 variants in summary statistics.
# 215,821 ambiguous SNPs have been removed.
# 1,145,260 variants have been matched; 38 were flipped and 499,125 were reversed.
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB, strand_flip = FALSE)
# 1,350,844 variants have been matched; 0 were flipped and 602,001 were reversed.
info_snp <- subset(na.omit(info_snp), info > 0.3)

list_snp_id <- with(info_snp, split(paste(chr, pos, a0, a1, sep = "_"),
                                    factor(chr, levels = 1:22)))
beta <- info_snp$beta
lpval <- -log10(info_snp$p)
info <- info_snp$info

# subset samples
library(bigreadr)
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
ind_diabetes <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x %in% 1220:1223)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) %in% paste0("E", 10:14)))
))))
ind_TD1 <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1222)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "E10"))
))))
ind_TD2 <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1223)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "E11"))
))))
y <- rep(0, nrow(df_illness))
y[ind_diabetes] <- NA
y[ind_TD2] <- 1
y[ind_TD1] <- NA

sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian & !is.na(y))
length(sub)  # 328,723
table(y.sub <- y[sub])
#      0      1
# 314547  14176

NCORES <- 12
system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_T2D",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 3.6H

library(bigsnpr)
ukbb <- snp_attach("data/UKBB_T2D.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 348 GB
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos

set.seed(1)
ind.train <- sort(sample(length(sub), 250e3))
ind.test <- setdiff(seq_along(sub), ind.train)

system.time(
  all_keep <- snp_grid_clumping(
    G, CHR, POS, lpS = lpval, ind.row = sort(sample(ind.train, 20e3)),
    grid.thr.imp = c(0.3, 0.6, 0.9, 0.95), infos.imp = info,
    grid.thr.r2 = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95),
    grid.base.size = c(50, 100, 200, 500),
    ncores = 6  # use less cores because of swapping if not enough memory
  )
) # 4H

system.time(
  multi_PRS <- snp_grid_PRS(
    G, all_keep, betas = beta, lpS = lpval, ind.row = ind.train,
    n_thr_lpS = 50, backingfile = "data/UKBB_T2D_scores", ncores = NCORES
  )
) # 2.2H

system.time(
  final_mod <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train], ncores = NCORES, n.abort = 5
  )
) # 5H
mod <- final_mod$mod
plot(mod)
summary(mod)
# # A tibble: 3 x 6
#    alpha validation_loss intercept beta           nb_var message
#    <dbl>           <dbl>     <dbl> <list>          <int> <list>
# 1 0.0001           0.173     -2.04 <dbl [65,520]>  23813 <chr [10]>
# 2 0.01             0.172     -2.01 <dbl [65,520]>   7280 <chr [10]>
# 3 1                0.172     -1.99 <dbl [65,520]>   5317 <chr [10]>

new_beta <- final_mod$beta.G

length(ind <- which(new_beta != 0))  # 477,267
summary(new_beta)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -0.3631526  0.0000000  0.0000000  0.0000209  0.0000000  0.2320248
summary(new_beta[which(sign(new_beta * beta) < 0)])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -4.759e-02 -2.303e-04 -1.078e-05 -8.076e-05  1.185e-04  3.865e-02

pred <- final_mod$intercept +
  big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind)

AUCBoot(pred, y.sub[ind.test])  # 63.7 [62.8-64.7]


library(tidyverse)

ind2 <- sample(ind, 10e3)
ggplot(data.frame(y = new_beta, x = beta)[ind, ]) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  geom_abline(slope = 0, intercept = 0, color = "blue") +
  geom_point(aes(x, y), size = 0.8) +
  theme_bigstatsr() +
  labs(x = "Effect sizes from GWAS", y = "Non-zero effect sizes from SCT")

grid2 <- attr(all_keep, "grid") %>%
  mutate(thr.lp = list(attr(multi_PRS, "grid.lpS.thr")), num = row_number()) %>%
  unnest()
s <- nrow(grid2)
grid2$auc <- big_apply(multi_PRS, a.FUN = function(X, ind, s, y.train) {
  # Sum over all chromosomes, for the same C+T parameters
  single_PRS <- rowSums(X[, ind + s * (0:21)])
  bigstatsr::AUC(single_PRS, 1 - y.train)
}, ind = 1:s, s = s, y.train = y.sub[ind.train],
a.combine = 'c', block.size = 1, ncores = NCORES)

std_prs <- grid2 %>%
  filter(thr.imp == 0.3, thr.r2 == 0.2, size == 500, thr.lp <= 8) %>%
  arrange(desc(auc)) %>%
  slice(1) %>%
  print()
#   size thr.r2 thr.imp num   thr.lp      auc
# 1  500    0.2     0.3  14 2.263405 0.566519

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
) # 56.9 [56.2-57.6]
sum(lpval[ind.keep] > std_prs$thr.lp)  # 3288

max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
#    size thr.r2 thr.imp num   thr.lp       auc
# 1  2500    0.2    0.90  72 3.645079 0.5703098
# 2  2500    0.2    0.30  16 3.645079 0.5702892
# 3  2500    0.2    0.60  44 3.645079 0.5702892
# 4  2500    0.2    0.95 100 3.645079 0.5702163
# 5  5000    0.1    0.95  96 4.176706 0.5698633
# 6  2500    0.2    0.95 100 3.405206 0.5698004
# 7  5000    0.1    0.95  96 4.470926 0.5697556
# 8  2500    0.2    0.90  72 3.405206 0.5697458
# 9  2500    0.2    0.95 100 3.901849 0.5697330
# 10 5000    0.1    0.95  96 3.901849 0.5697202

ind.keep <- unlist(map(all_keep, max_prs$num))
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp),
  y.sub[ind.test]
) # 57.3 [56.7-58.0]
sum(lpval[ind.keep] > max_prs$thr.lp)  # 368

ggplot(grid2) +
  geom_point(aes(thr.lp, auc)) +
  facet_grid(thr.imp ~ thr.r2 + size) +
  scale_x_log10(limits = c(2, 10), breaks = c(1, 2, 5), minor_breaks = 1:10) +
  ylim(0.57, NA) +
  geom_hline(yintercept = max(1 - aucs), color = "blue", linetype = 2) +
  geom_hline(yintercept = AUC(pred, y.sub[ind.test]), color = "red", linetype = 2) +
  theme_bigstatsr(size.rel = 0.7) +
  labs(x = "-log10(p-value) threshold (log scale)", y = "AUC")


# Redo GWAS quickly to check
bmi_age <- as.matrix(fread2(csv, select = c("21001-0.0", "21003-0.0")))[sub, ]
ind.train2 <- ind.train[-which(is.na(bmi_age[ind.train, ]), arr.ind = TRUE)[, "row"]]
system.time(
  gwas <- big_univLinReg(G, y.sub[ind.train2], ind.train = ind.train2, ind.col = 1:10000,
                         covar.train = cbind(COVAR, bmi_age)[ind.train2, ], ncores = NCORES)
)

qplot(beta[1:10000], gwas$estim) +
  theme_bigstatsr() +
  ylim(-2, 2) +
  geom_abline(slope = -1, intercept = 0, color = "red")

qplot(lpval[1:10000], -predict(gwas), alpha = I(0.2)) +
  theme_bigstatsr() +
  ylim(NA, 5) +
  geom_abline(slope = -1, intercept = 0, color = "red") +
  labs(x = "-log10(p-values) of initial GWAS", y = "-log10(p-values) of new GWAS")
