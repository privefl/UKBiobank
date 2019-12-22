# download.file("https://portals.broadinstitute.org/collaboration/giant/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz",
#               destfile = "sumstats_BMI.txt.gz")
# R.utils::gunzip("sumstats_BMI.txt.gz")
library(bigreadr)
sumstats <- fread2("sumstats_BMI.txt", select = c(1:3, 5, 7),
                   col.names = c("snp", "a1", "a0", "beta", "p"))
nrow(sumstats)  # 2,554,637

sumstats <- subset(sumstats, p < 0.5)
nrow(sumstats)  # 1,310,033
hist(sumstats$beta)
hist(sumstats$p)

info_snp_UKBB <- rbind_df(lapply(1:22, function(chr) {
  file <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
  df <- fread2(file, select = c(2:5, 8), col.names = c("snp", "pos", "a0", "a1", "info"))
  cbind.data.frame(chr = chr, df)
}))
sumstats2 <- merge(sumstats, info_snp_UKBB, by = "snp", suffixes = c("", ".2"))
info_snp <- bigsnpr::snp_match(sumstats2[2:7], info_snp_UKBB)
# 1,312,161 variants in summary statistics.
# 203,136 ambiguous SNPs have been removed.
# 1,108,965 variants have been matched; 394 were flipped and 533,570 were reversed.
info_snp <- bigsnpr::snp_match(sumstats2[2:7], info_snp_UKBB, strand_flip = FALSE)
# 1,311,678 variants have been matched; 0 were flipped and 632,162 were reversed.
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
  select = c("eid", "21001-0.0", "22006-0.0"),
  col.names = c("eid", "BMI", "is_caucasian")
)
ind.indiv <- match(df0$eid, sample$ID_2)
# related samples
df0$is_rel2 <- df0$eid %in% fread2("ukb25589_rel_s488346.dat")$ID2



sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian & !is.na(df0$BMI))
length(y.sub <- df0$BMI[sub])  # 334,540
hist(y.sub)
hist(log(y.sub))

NCORES <- 15
system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_BMI",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 5.5H

library(bigsnpr)
ukbb <- snp_attach("data/UKBB_BMI.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 269 GB
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos

set.seed(1)
ind.train <- sort(sample(length(sub), 300e3))
ind.test <- setdiff(seq_along(sub), ind.train)

#### SCT ####

system.time(
  all_keep <- snp_grid_clumping(
    G, CHR, POS, lpS = lpval, ind.row = ind.train, infos.imp = info,
    grid.thr.imp = c(0.3, 0.6, 0.9, 0.95),
    grid.thr.r2 = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95),
    grid.base.size = c(50, 100, 200, 500),
    ncores = 8
  )
) # 3H

system.time(
  multi_PRS <- snp_grid_PRS(
    G, all_keep, betas = beta, lpS = lpval, ind.row = ind.train,
    n_thr_lpS = 50, backingfile = "data/UKBB_BMI_scores", ncores = NCORES
  )
) # 27 min

system.time(
  final_mod <- snp_grid_stacking(
    multi_PRS, log(y.sub[ind.train]), ncores = NCORES, n.abort = 5
  )
) # 87.5H
mod <- final_mod$mod
plot(mod)
summary(mod)
# # A tibble: 3 x 6
#    alpha validation_loss intercept beta           nb_var message
#    <dbl>           <dbl>     <dbl> <list>          <int> <list>
# 1 0.0001          0.0249      3.37 <dbl [80,640]>  73398 <chr [10]>
# 2 0.01            0.0249      3.37 <dbl [80,640]>  34252 <chr [10]>
# 3 1               0.0249      3.38 <dbl [80,640]>  29705 <chr [10]>

# Save
save(all_keep, multi_PRS, final_mod, file = "data/res_BMI.RData")

new_beta <- final_mod$beta.G

length(ind <- which(new_beta != 0))  # 545,179
summary(new_beta)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -0.0105007  0.0000000  0.0000000  0.0000004  0.0000000  0.0117095
summary(new_beta[which(sign(new_beta * beta) < 0)])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -2.946e-03 -1.389e-05 -9.400e-08 -1.526e-06  1.272e-05  5.308e-03

pred <- final_mod$intercept +
  big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind)

cor(exp(pred), y.sub[ind.test])  # 30.8 -> 9.5% variance


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
grid2$cor <- big_apply(multi_PRS, a.FUN = function(X, ind, s, y.train) {
  # Sum over all chromosomes, for the same C+T parameters
  single_PRS <- rowSums(X[, ind + s * (0:21)])
  cor(exp(single_PRS), y.train)
}, ind = 1:s, s = s, y.train = y.sub[ind.train],
a.combine = 'c', block.size = 1, ncores = NCORES)

std_prs <- grid2 %>%
  filter(thr.imp == 0.3, thr.r2 == 0.2, size == 500, thr.lp <= 8) %>%
  arrange(desc(cor)) %>%
  slice(1) %>%
  print()
#   size thr.r2 thr.imp num   thr.lp       cor
# 1  500    0.2     0.3  14 2.298557 0.2084777

ind.keep <- unlist(map(all_keep, std_prs$num))
# Verif on training set
cor(
  snp_PRS(G, beta[ind.keep], ind.test = ind.train, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.train]
) # 0.2285089 -> not the same?!
# Eval on test set
cor(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.test]
) # 0.2323813
sum(lpval[ind.keep] > std_prs$thr.lp)  # 3956

max_prs <- grid2 %>% arrange(desc(cor)) %>% slice(1:10) %>% print() %>% slice(1)
#    size thr.r2 thr.imp num   thr.lp       cor
# 1  5000    0.1    0.95  96 2.298557 0.2150117
# 2  2000    0.1    0.95  95 2.298557 0.2149480
# 3  1000    0.1    0.95  94 2.298557 0.2147743
# 4   500    0.1    0.95  93 2.298557 0.2142485
# 5  5000    0.1    0.90  68 2.298557 0.2137942
# 6  2000    0.1    0.90  67 2.298557 0.2137314
# 7  1000    0.1    0.90  66 2.298557 0.2136516
# 8  5000    0.1    0.60  40 2.298557 0.2132900
# 9  5000    0.1    0.30  12 2.298557 0.2132774
# 10 2000    0.1    0.60  39 2.298557 0.2132288

ind.keep <- unlist(map(all_keep, max_prs$num))
cor(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp),
  y.sub[ind.test]
) # 0.231033
sum(lpval[ind.keep] > max_prs$thr.lp)  # 3166

ggplot(grid2) +
  geom_point(aes(thr.lp, cor)) +
  facet_grid(thr.imp ~ thr.r2 + size) +
  scale_x_log10(limits = c(1, 8), breaks = c(1, 2, 5), minor_breaks = 1:10) +
  ylim(0.1, NA) +
  theme_bigstatsr(size.rel = 0.7) +
  labs(x = "-log10(p-value) threshold (log scale)", y = "Correlation")


#### prioPLR ####
conf <- pmax(lpval, 1) * info
system.time(
  mod <- big_spLinReg(G, y.sub[ind.train], ind.train, #pf.X = 1 / conf,
                      K = 10, alphas = 0.01, # lambda.min = 1e-5,
                      dfmax = 300e3, n.abort = 3, ncores = NCORES)
) # 20H
plot(mod)
summary(mod)
# # A tibble: 1 x 6
#   alpha validation_loss intercept beta              nb_var message
#   <dbl>           <dbl>     <dbl> <list>             <int> <list>
# 1     1            20.2      34.5 <dbl [1,306,039]>  94015 <chr [10]>
summary(mod)$message  ## No more improvement

pred <- predict(mod, G, ind.test)
cor(pred, y.sub[ind.test])         # 0.3600648 / 0.3601948 -> 0.3494391 without pf.X
cor(predict(mod[3], G, ind.test), y.sub[ind.test])  #
plot(pred, y.sub[ind.test], pch = 20); abline(0, 1, col = "red", lwd = 2)
plot(pred_log, log(y.sub[ind.test]), pch = 20); abline(0, 1, col = "red", lwd = 2)

new_beta <- summary(mod, best.only = TRUE)$beta[[1]]
ind <- which(new_beta != 0)
plot(beta[attr(mod, "ind.col")[ind]], new_beta[ind], pch = 20); abline(0, 1, col = "red")
plot(beta[attr(mod, "ind.col")[ind]], new_beta[ind], pch = 20, ylim = c(-1, 1)); abline(0, 1, col = "red")

snp_MAF(G, ind.col = attr(mod, "ind.col")[new_beta > 80]) # 2e-07
plot(G[,  attr(mod, "ind.col")[new_beta > 80]], y.sub)
big_univLinReg(G, y.sub, ind.col = attr(mod, "ind.col")[new_beta > 80])
#      estim  std.err    score
# 1 244.7094 67.28775 3.636759

# Verif
# pred2 <- summary(mod, best.only = TRUE)$intercept +
#   big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = attr(mod, "ind.col")[ind])
# all.equal(pred, pred2, check.attributes = FALSE)

pred_log <- predict(mod_log, G, ind.test)
plot(pred, exp(pred_log), pch = 20); abline(0, 1, col = "red")


save(mod, mod_log, file = "BMI.RData")


pred_PRS <- snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
                    lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp)
pred_PRS_train <- snp_PRS(G, beta[ind.keep], ind.test = ind.train, ind.keep = ind.keep,
                          lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp)
plot(pred_PRS, pred, pch = 20)

summary(lm(y.sub[ind.train] ~ pred_PRS_train))$coef
pred_PRS2 <- 1.909265 * pred_PRS + 28.531365
plot(pred_PRS2, pred, pch = 20); abline(0, 1, col = "red", lwd = 2)
sapply(c(0.5, 0.3, 0.2, 0.15, 0.1, 0.05, 0), function(a) {
  cor((a * pred_PRS2 + pred) / (1 + a), y.sub[ind.test])
})
a <- 0.1; cor((a * pred_PRS2 + pred) / (1 + a), y.sub[ind.test])
cor(pred, y.sub[ind.test])
