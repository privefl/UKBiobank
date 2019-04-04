# Download from http://diagram-consortium.org/downloads.html
# DIAGRAM 1000G GWAS meta-analysis Stage 1 Summary statistics
# Published in in Scott et al (2017)
unzip("METAANALYSIS_DIAGRAM_SE1.zip")
library(bigreadr)
sumstats <- fread2("METAANALYSIS_DIAGRAM_SE1.txt", select = c(1:4, 6))
sumstats <- tidyr::separate(sumstats, "Chr:Position", c("chr", "pos"), convert = TRUE)
names(sumstats) <- c("chr", "pos", "a1", "a2", "beta", "pval")
nrow(sumstats)  # 12,056,346
hist(sumstats$beta)
hist(sumstats$pval)

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
    bigreadr::fread2(showProgress = FALSE, nThread = 1) %>%
    inner_join(
      sumstats2[sumstats2$chr == chr, ], ., na_matches = "never",
      by = c(pos = "V3", a1 = "V4", a2 = "V5")
    ) %>%
    arrange(pos) %>%
    transmute(
      id    = paste(chr, pos, a1, a2, sep = "_"),
      beta  = -beta,
      lpval = -log10(pval),
      info  = V8
    ) %>%
    na.omit() %>%
    filter(lpval > 1, info > 0.3)
}
stopCluster(cl)

list_snp_id <-  lapply(info_snp, function(df) df$id)
beta <-  unlist(lapply(info_snp, function(df) df$beta))
lpval <- unlist(lapply(info_snp, function(df) df$lpval))
info <-  unlist(lapply(info_snp, function(df) df$info))
length(beta)  # 1,340,712

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


df_illness <- fread2(csv, select = c(paste0("20002-0.", 0:28),
                                     paste0("20002-1.", 0:28),
                                     paste0("20002-2.", 0:28)))
df_ICD10 <- fread2(csv, select = c(paste0("41202-0.", 0:379),
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

sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian & !is.na(y) &
               !is.na(df0$sex))
length(sub)  # 328,726
table(y.sub <- y[sub])
#      0      1
# 314570  14156

system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_T2D",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 4H

library(bigsnpr)
ukbb <- snp_attach("data/UKBB_T2D.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 410 GB
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos

set.seed(1)
ind.train <- sort(sample(length(sub), 300e3))
ind.test <- setdiff(seq_along(sub), ind.train)

system.time(
  all_keep <- snp_grid_clumping(
    G, CHR, POS, lpS = lpval, ind.row = sort(sample(ind.train, 20e3)),
    grid.thr.imp = c(0.3, 0.6, 0.9, 0.95), infos.imp = info,
    grid.thr.r2 = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95),
    grid.base.size = c(50, 100, 200, 500),
    ncores = 6  # use less cores because of swapping if not enough memory
  )
) # 6H
plot(lengths(all_keep[[1]]), pch = 20)

system.time(
  multi_PRS <- snp_grid_PRS(
    G, all_keep, betas = beta, lpS = lpval, ind.row = ind.train,
    grid.lpS.thr = seq_log(1, 0.999 * max(lpval), 50),
    backingfile = "data/UKBB_T2D_scores", ncores = NCORES
  )
) # 4.5H

nPC <- 20
PC <- fread2(csv, select = paste0("22009-0.", 1:nPC),
             col.names = paste0("PC", 1:nPC))
COVAR <- as.matrix(cbind(df0$sex, PC)[sub, ])

system.time(
  final_mod <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train], covar.train = COVAR[ind.train, ], ncores = NCORES
  )
) # 13H
mod <- final_mod$mod
plot(mod)
summary(mod)
# # A tibble: 5 x 6
#      alpha validation_loss intercept beta           nb_var message
#      <dbl>           <dbl>     <dbl> <list>          <int> <list>
#   1 0.0001           0.171     -3.71 <dbl [62,349]>  34485 <chr [10]>
#   2 0.001            0.171     -3.67 <dbl [62,349]>  11379 <chr [10]>
#   3 0.01             0.171     -3.65 <dbl [62,349]>   6876 <chr [10]>
#   4 0.1              0.171     -3.67 <dbl [62,349]>   5795 <chr [10]>
#   5 1                0.171     -3.65 <dbl [62,349]>   5952 <chr [10]>

new_beta <- final_mod$beta.G

length(ind <- which(new_beta != 0))  # 520,027
summary(new_beta)
summary(new_beta[which(sign(new_beta * beta) < 0)])

pred <- final_mod$intercept +
  big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind) +
  COVAR[ind.test, ] %*% final_mod$beta.covar

AUCBoot(pred, y.sub[ind.test])  # 65.6 [64.0-67.1]


max(aucs <- apply(multi_PRS, 2, AUC, y.sub[ind.train]))  # 47.1

library(tidyverse)

ind2 <- sample(ind, 10e3)
ggplot(data.frame(y = new_beta, x = beta)[ind, ]) +
  geom_abline(slope = -1, intercept = 0, color = "red") +
  geom_abline(slope = 0, intercept = 0, color = "blue") +
  geom_point(aes(x, y), size = 0.8) +
  theme_bigstatsr() +
  labs(x = "Effect sizes from GWAS", y = "Non-zero effect sizes from SCT")

grid2 <- attr(all_keep, "grid") %>%
  mutate(thr.lp = list(attr(multi_PRS, "grid.lpS.thr"))) %>%
  unnest() %>%
  mutate(auc = 1 - aucs)

grid2 %>%
  filter(thr.imp == 0.3, thr.r2 == 0.2, size == 500, thr.lp <= 8) %>%
  arrange(desc(auc)) %>%
  slice(1)
#   size thr.r2 thr.imp   thr.lp       auc
# 1  500    0.2     0.3 5.818992 0.5901172

grid2 %>% arrange(desc(auc)) %>% slice(1:10)
#     size thr.r2 thr.imp   thr.lp       auc
# 1   1000   0.05    0.95 6.384147 0.5955242
# 2   2000   0.05    0.95 6.384147 0.5955242
# 3   4000   0.05    0.95 6.384147 0.5955242
# 4  10000   0.05    0.95 6.384147 0.5955242
# 5   5000   0.01    0.95 5.818992 0.5953669
# 6  10000   0.01    0.95 5.818992 0.5953669
# 7  20000   0.01    0.95 5.818992 0.5953669
# 8  50000   0.01    0.95 5.818992 0.5953669
# 9   5000   0.01    0.95 6.384147 0.5951054
# 10 10000   0.01    0.95 6.384147 0.5951054

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
