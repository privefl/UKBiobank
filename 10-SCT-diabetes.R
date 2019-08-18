#### T1D SUMSTATS ####

# urls <- gsubfn::strapply(
#   readLines("https://datadryad.org//resource/doi:10.5061/dryad.ns8q3"),
#   "<a href=\"(/bitstream/handle/10255/dryad\\.[0-9]+/meta_chr_[0-9]+\\?sequence=1)\">",
#   simplify = 'c')
#
# sumstats <- purrr::map_dfr(urls, ~ {
#   download.file(paste0("https://datadryad.org", .x),
#                 destfile = (tmp <- tempfile(fileext = ".txt")))
#   sumstats <- bigreadr::fread2(
#     tmp, select = c("chromosome", "position", "a0", "a1", "beta.meta", "p.meta"),
#     col.names = c("chr", "pos", "a0", "a1", "beta", "p"))
#   na.omit(sumstats)
# })
#
# saveRDS(sumstats, "sumstats_T1D.rds")

sumstats_T1D <- readRDS("sumstats_T1D.rds")
nrow(sumstats_T1D)  # 8,996,866
hist(sumstats_T1D$p)
hist(sumstats_T1D$beta)

sumstats_T1D <- subset(sumstats_T1D, p < 0.1)

library(bigreadr)
info_snp_UKBB <- rbind_df(lapply(1:22, function(chr) {
  file <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
  df <- fread2(file, select = c(3:5, 8), col.names = c("pos", "a0", "a1", "info"))
  cbind.data.frame(chr = chr, df)
}))
info_snp_T1D <- bigsnpr::snp_match(sumstats_T1D, info_snp_UKBB)
# 1,241,172 variants in summary statistics.
# 172,374 ambiguous SNPs have been removed.
# 1,051,668 variants have been matched; 76 were flipped and 522 were reversed.
info_snp_T1D <- bigsnpr::snp_match(sumstats_T1D, info_snp_UKBB, strand_flip = FALSE)
# 1,222,801 variants have been matched; 0 were flipped and 471 were reversed.
info_snp_T1D <- subset(na.omit(info_snp_T1D), info > 0.6)


#### T2D SUMSTATS ####

# Download from http://diagram-consortium.org/downloads.html
# DIAGRAM 1000G GWAS meta-analysis Stage 1 Summary statistics
# Published in Scott et al (2017)
# unzip("METAANALYSIS_DIAGRAM_SE1.zip")
library(bigreadr)
sumstats_T2D <- fread2("METAANALYSIS_DIAGRAM_SE1.txt", select = c(1:4, 6))
sumstats_T2D <- tidyr::separate(sumstats_T2D, "Chr:Position", c("chr", "pos"), convert = TRUE)
names(sumstats_T2D) <- c("chr", "pos", "a1", "a0", "beta", "p")
nrow(sumstats_T2D)  # 12,056,346
hist(sumstats_T2D$beta)
hist(sumstats_T2D$p)
sumstats_T2D <- subset(sumstats_T2D, p < 0.1)

info_snp_T2D <- bigsnpr::snp_match(sumstats_T2D, info_snp_UKBB)
# 1,408,672 variants in summary statistics.
# 215,821 ambiguous SNPs have been removed.
# 1,145,260 variants have been matched; 38 were flipped and 499,125 were reversed.
info_snp_T2D <- bigsnpr::snp_match(sumstats_T2D, info_snp_UKBB, strand_flip = FALSE)
# 1,350,844 variants have been matched; 0 were flipped and 602,001 were reversed.
info_snp_T2D <- subset(na.omit(info_snp_T2D), info > 0.6)

info_snp0 <- merge(info_snp_T1D, info_snp_T2D,
                   by = c("chr", "pos", "a0", "a1", "info"),
                   suffixes = c("_T1D","_T2D"))
str(info_snp0)
plot(info_snp0[c("beta_T1D", "beta_T2D")], pch = 20); abline(0, 1, col = "red")

library(dplyr)
info_snp <- rbind(cbind(info_snp_T1D, id = 1), cbind(info_snp_T2D, id = 2)) %>%
  arrange(chr, pos)
str(info_snp)

list_snp_id <- with(info_snp, split(paste(chr, pos, a0, a1, sep = "_"),
                                    factor(chr, levels = 1:22)))
beta <- info_snp$beta
lpval <- -log10(info_snp$p)
info <- info_snp$info
# Infinite values because of large effects on chromosome 6
lpval <- pmin(lpval, -log10(.Machine$double.xmin) + abs(beta))

#### PHENOTYPES ####
library(bigreadr)
sample <- fread2("ukb25589_imp_chr1_v3_s487327.sample")
str(sample)
(N <- readBin("data/ukb_imp_chr1_v3.bgen", what = 1L, size = 4, n = 4)[4])
sample <- sample[-1, ]
stopifnot(nrow(sample) == N)

csv <- "ukb22544.csv"
df0 <- fread2(csv, select = c("eid", "22006-0.0", "21001-0.0", "31-0.0",
                              paste0("2976-", 0:2, ".0")),
              col.names = c("eid", "is_caucasian", "BMI", "sex",
                            paste0("diag_age", 0:2)))
diag_age <- df0[5:7] %>%
  mutate_all(~ ifelse(. > 0, ., NA)) %>%
  rowMeans(na.rm = TRUE)

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
ind_TD1 <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1222)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "E10"))
))))
ind_TD2 <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1223)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "E11"))
))))
y <- rep(NA, nrow(df0))
y[ind_TD2] <- 1
y[ind_TD1] <- 0
y[intersect(ind_TD1, ind_TD2)] <- NA

sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian & !is.na(y) &
               diag_age > 0)
length(sub)  # 10,288
table(y.sub <- y[sub])
#      0    1
#    642 9646

library(bigsnpr)
NCORES <- nb_cores()
system.time(
  rds <- snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_diabetes",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 7H

ukbb <- snp_attach("data/UKBB_diabetes.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 24 GB
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos

set.seed(1); ind.train <- sort(sample(length(sub), 8e3))
ind.test <- setdiff(seq_along(sub), ind.train)

system.time(
  all_keep <- snp_grid_clumping(
    G, CHR, POS, lpS = lpval, ind.row = ind.train, infos.imp = info,
    grid.thr.imp = c(0.6, 0.9, 0.95),
    grid.thr.r2 = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95),
    grid.base.size = c(50, 200, 500),
    groups = split(rows_along(info_snp), info_snp$id),
    ncores = NCORES
  )
) # 2.4H

system.time(
  multi_PRS <- snp_grid_PRS(
    G, all_keep, betas = beta, lpS = lpval, ind.row = ind.train,
    n_thr_lpS = 50, backingfile = "data/UKBB_diabetes_scores", ncores = NCORES
  )
) # 2 min

system.time(
  final_mod <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train], ncores = NCORES, K = 4
  )
) # 4 min
mod <- final_mod$mod
plot(mod)
summary(mod)
# A tibble: 3 x 6
# alpha validation_loss intercept beta           nb_var message
#    <dbl>           <dbl>     <dbl> <list>          <int> <list>
# 1 0.0001           0.189      3.46 <dbl [61,404]>  22919 <chr [4]>
# 2 0.01             0.190      3.63 <dbl [61,404]>   3763 <chr [4]>
# 3 1                0.190      3.82 <dbl [61,404]>   1278 <chr [4]>

new_beta <- final_mod$beta.G

length(ind <- which(new_beta != 0))  # 824,130

pred <- final_mod$intercept +
  big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind)

AUCBoot(pred, y.sub[ind.test])  # 80.3 [76.4-83.9]

library(tidyverse)

ggplot(data.frame(pheno = factor(y.sub[ind.test], levels = 0:1, labels = c("T1D", "T2D")),
                  pred = 1 / (1 + exp(-pred)))) +
  theme_bigstatsr() +
  geom_density(aes(pred, fill = pheno), alpha = 0.3) +
  theme(legend.position = c(0.2, 0.7)) +
  labs(fill = NULL, x = "Probability of type 2 (PRS)") +
  ggtitle("AUC: 80.3 [76.4-83.9]")

ggplot(data.frame(pheno = factor(y.sub[ind.test], levels = 0:1, labels = c("T1D", "T2D")),
                  pred = 1 / (1 + exp(-pred)),
                  age = diag_age[sub][ind.test])) +
  theme_bigstatsr() +
  geom_point(aes(pred, age, color = pheno), alpha = 0.7) +
  theme(legend.position = c(0.1, 0.55)) +
  labs(color = NULL, x = "Probability of type 2 (PRS)", y = "Age at diagnosis") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

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


max_prs <- grid2 %>% arrange(desc(pmax(auc, 1 - auc))) %>%
  slice(1:10) %>% print() %>% slice(1)
#    size thr.r2 grp.num thr.imp num   thr.lp       auc
# 1  5000    0.1       1    0.90  68 3.624586 0.2090946
# 2  2000    0.1       1    0.90  67 3.624586 0.2091041
# 3  2000    0.1       1    0.60  11 3.624586 0.2091687
# 4  5000    0.1       1    0.60  12 3.624586 0.2096240
# 5  2000    0.1       1    0.60  11 4.074776 0.2097023
# 6  5000    0.1       1    0.60  12 4.074776 0.2098944
# 7  2000    0.1       1    0.95 123 3.624586 0.2100633
# 8  5000    0.1       1    0.95 124 3.624586 0.2104783
# 9  2000    0.1       1    0.60  11 4.580883 0.2106477
# 10 1000    0.1       1    0.95 122 3.624586 0.2107607

ind.keep <- unlist(map(all_keep, max_prs$num))
AUCBoot(
  snp_PRS(G, -beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp),
  y.sub[ind.test]
) # 76.8 [72.8-80.5]
sum(lpval[ind.keep] > max_prs$thr.lp)  # 1220

ggplot(grid2) +
  geom_point(aes(thr.lp, pmax(auc, 1 - auc), color = as.factor(grp.num))) +
  facet_grid(thr.imp ~ thr.r2 + size) +
  scale_x_log10(limits = c(2, 10), breaks = c(1, 2, 5), minor_breaks = 1:10) +
  ylim(0.50, NA) +
  theme_bigstatsr(size.rel = 0.7) +
  labs(x = "-log10(p-value) threshold (log scale)", y = "AUC", color = "Type")


#### TRY CLUSTERING ####

ms <- big_scale()(multi_PRS)

svd <- big_randomSVD(multi_PRS, fun.scaling = big_scale(),
                     ind.col = which(ms$scale > 0),
                     ncores = nb_cores())
apply(svd$u, 2, AUC, y.sub[ind.train])
plot(svd)
plot(svd, type = "scores", scores = 1:2) +
  aes(color = as.factor(y.sub[ind.train] + 1), alpha = I(0.6)) +
  labs(color = "Type")
plot(svd, type = "scores", scores = 3:4) +
  aes(color = as.factor(y.sub[ind.train] + 1), alpha = I(0.6)) +
  labs(color = "Type")

#### USING AGE AT DIAGNOSIS ####

age.train <- diag_age[sub][ind.train]
system.time(
  final_mod2 <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train], ncores = NCORES, K = 4,
    covar.train = cbind(age.train, log1p(age.train)), pf.covar = c(0, 0)
  )
) # 4 min
mod2 <- final_mod2$mod
plot(mod2)
summary(mod2)
# ## A tibble: 3 x 6
#    alpha validation_loss intercept beta           nb_var message
#    <dbl>           <dbl>     <dbl> <list>          <int> <list>
# 1 0.0001           0.137    -0.117 <dbl [80,949]>  22854 <chr [5]>
# 2 0.01             0.138     0.108 <dbl [80,949]>   2819 <chr [5]>
# 3 1                0.138     0.135 <dbl [80,949]>   1188 <chr [5]>
## Adding log(age):
# A tibble: 3 x 6
#    alpha validation_loss intercept beta           nb_var message
#    <dbl>           <dbl>     <dbl> <list>          <int> <list>
# 1 0.0001           0.128      7.17 <dbl [80,950]>  20818 <chr [5]>
# 2 0.01             0.129      7.49 <dbl [80,950]>   2915 <chr [5]>
# 3 1                0.130      7.53 <dbl [80,950]>   1216 <chr [5]>

new_beta2 <- final_mod2$beta.G

length(ind2 <- which(new_beta2 != 0))  # 779,802 / 726,426

pred2 <- final_mod2$intercept +
  big_prodVec(G, new_beta2[ind2], ind.row = ind.test, ind.col = ind2)

AUCBoot(pred2, y.sub[ind.test])  # 79.4 [75.5-83.2] / 79.5 [75.5-83.2]


age.test <- diag_age[sub][ind.test]
pred.covar <- cbind(age.test, log1p(age.test)) %*% final_mod2$beta.covar
pred3 <- pred2 + pred.covar

curve(crossprod(rbind(x, log1p(x)), final_mod2$beta.covar), from = 0.5, to = 80,
      xlab = "Age of diagnosis", ylab = "Effect of age of diagnosis")

AUCBoot(pred3, y.sub[ind.test])  # 90.5 [87.3-93.3] / 90.9 [87.7-93.7]
AUCBoot(pred.covar, y.sub[ind.test])

ggplot(data.frame(pheno = factor(y.sub[ind.test], levels = 0:1, labels = c("T1D", "T2D")),
                  pred = 1 / (1 + exp(-pred3)))) +
  theme_bigstatsr() +
  geom_density(aes(pred, fill = pheno), alpha = 0.3) +
  theme(legend.position = c(0.2, 0.7)) +
  labs(fill = NULL, x = "Probability of type 2 (PRS + diag_age + log1p(diag_age))") +
  ggtitle("AUC: 91.0 [87.8-93.9]")

ggplot(data.frame(pheno = factor(y.sub[ind.test], levels = 0:1, labels = c("T1D", "T2D")),
                  pred = 1 / (1 + exp(-pred3)),
                  age = diag_age[sub][ind.test])) +
  theme_bigstatsr() +
  geom_point(aes(pred, age, color = pheno), alpha = 0.7) +
  theme(legend.position = c(0.15, 0.8)) +
  labs(x = "Probability of type 2 (PRS + diag_age + log1p(diag_age))",
       y = "Age at diagnosis", color = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggplot(data.frame(pheno = factor(y.sub[ind.test], levels = 0:1, labels = c("T1D", "T2D")),
                  pred = 1 / (1 + exp(-pred3)),
                  BMI = df0$BMI[sub][ind.test])) +
  theme_bigstatsr() +
  geom_point(aes(pred, BMI, color = pheno), alpha = 0.7) +
  theme(legend.position = c(0.1, 0.8)) +
  labs(color = "Type", x = "Probability of type 2 (PRS + diag_age)", y = "BMI")

ggplot(data.frame(pheno = factor(y.sub[ind.test], levels = 0:1, labels = c("T1D", "T2D")),
                  pred = 1 / (1 + exp(-pred3)),
                  age = diag_age[sub][ind.test],
                  BMI = df0$BMI[sub][ind.test],
                  sex = factor(df0$sex[sub][ind.test], levels = 0:1,
                               labels = c("Female", "Male")))) +
  theme_bigstatsr() +
  geom_point(aes(age, BMI, color = pheno), alpha = 0.7) +
  theme(legend.position = "top") +
  labs(color = NULL, x = "Age at diagnosis", y = "BMI") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  facet_wrap(~ sex)

# df <- data.frame(age = diag_age[sub][ind.test], y = y.sub[ind.test], PRS = pred,
#                  PRS_prob = 1 / (1 + exp(-pred)))
# summary(glm(y ~ poly(age, 2) + log(age) + PRS + PRS_prob, data = df, family = "binomial"))
# summary(glm(y ~ (log(age) + age) * PRS, data = df, family = "binomial"))

ggplot(data.frame(pheno = factor(y.sub[ind.test], levels = 0:1, labels = c("T1D", "T2D")),
                  pred = 1 / (1 + exp(-pred3)),
                  age = diag_age[sub][ind.test],
                  BMI = df0$BMI[sub][ind.test])) +
  theme_bigstatsr() +
  geom_density(aes(age, fill = pheno), alpha = 0.6) +
  theme(legend.position = c(0.1, 0.8)) +
  labs(color = "Type", x = "Age at diagnosis")

AUC(diag_age[sub][ind.test], y.sub[ind.test])  # 0.88

sub2 <- which(diag_age[sub] < 30)
ind.train2 <- intersect(ind.train, sub2)
table(y.sub[ind.train2])
# 0 -> 252 ;  1 -> 189
system.time(
  final_mod3 <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train2], ind.train = match(ind.train2, ind.train),
    ncores = NCORES, K = 3
  )
) # 4 min
mod3 <- final_mod3$mod
plot(mod3)
summary(mod3)

new_beta3 <- final_mod3$beta.G

length(ind3 <- which(new_beta3 != 0))  # 998,445

pred3 <- final_mod3$intercept +
  big_prodVec(G, new_beta3[ind3], ind.row = ind.test, ind.col = ind3)

AUCBoot(pred3, y.sub[ind.test])  # 78.2 [74.4-82.0]

bool <- diag_age[sub][ind.test] < 30
AUCBoot(pred3[bool], y.sub[ind.test][bool])  # 87.2 [80.3-93.1]

library(ggplot2)
ggplot(data.frame(pheno = factor(y.sub[ind.test], levels = 0:1, labels = c("T1D", "T2D")),
                  pred = 1 / (1 + exp(-pred3)),
                  age = diag_age[sub][ind.test])[bool, ]) +
  theme_bigstatsr() +
  geom_point(aes(pred, age, color = pheno), alpha = 0.9) +
  labs(color = "Type", x = "Probability of type 2", y = "Age at diagnosis")

#### XGBoost ####

library(xgboost)
pred_all <- final_mod$intercept + big_prodVec(G, new_beta[ind], ind.col = ind)
mat <- cbind(PRS = pred_all, PRS_prob = 1 / (1 + exp(-pred_all)),
             BMI = df0$BMI[sub], age = diag_age[sub])
plot(as.data.frame(mat))
dtrain <- xgb.DMatrix(mat[ind.train, ], label = y.sub[ind.train])
dtest <- xgb.DMatrix(mat[ind.test, ], label = y.sub[ind.test])
watchlist <- list(train = dtrain, eval = dtest)
param <- list(max_depth = 2, eta = 0.3, silent = 1, nthread = 2,
              objective = "binary:logistic", eval_metric = "auc")
bst <- xgb.train(param, dtrain, nrounds = 60, watchlist)
xgboost::xgb.importance(model = bst)
# xgboost::xgb.ggplot.importance(.Last.value)

pred.train <- predict(bst, dtrain)
qplot(pred.train, fill = as.factor(y.sub[ind.train]), geom = "density", alpha = I(0.5))


#### PLR ####
age.train <- diag_age[sub][ind.train]
system.time(
  mod_PLR <- big_spLogReg(G, y.sub[ind.train], ind.train = ind.train, pf.X = 1 / beta^2,
                          alphas = 10^(-(0:4)), K = 5, ncores = nb_cores(),
                          covar.train = cbind(age.train, log1p(age.train)), pf.covar = c(0, 0),
                          nlam.min = 25)
)
plot(mod_PLR)
summary(mod_PLR)
pred_PLR <- predict(mod_PLR, G, ind.test, covar.row = matrix(0, length(ind.test), 2))
AUCBoot(pred_PLR, y.sub[ind.test])
# 77.9 [73.9-81.7] / 78.4 [74.4-82.1] / 78.7 [74.9-82.4]


age.test <- diag_age[sub][ind.test]
pred_PLR2 <- predict(mod_PLR, G, ind.test, covar.row = cbind(age.test, log1p(age.test)))
AUCBoot(pred_PLR2, y.sub[ind.test]) # 91.1 [88.0-93.9]
library(ggplot2)
ggplot(data.frame(pheno = factor(y.sub[ind.test], levels = 0:1, labels = c("T1D", "T2D")),
                  pred = 1 / (1 + exp(-pred_PLR)),
                  age = diag_age[sub][ind.test])) +
  theme_bigstatsr() +
  geom_point(aes(pred, age, color = pheno), alpha = 0.7) +
  theme(legend.position = c(0.2, 0.8)) +
  labs(color = "Type", x = "Probability of type 2", y = "Age at diagnosis")

table(pred = pred_PLR > 0.7, pheno = y.sub[ind.test])
beta_PLR <- summary(mod_PLR, best.only = TRUE)$beta[[1]]
sum(beta_PLR != 0) # 5736
