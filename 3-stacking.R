sample <- bigreadr::fread2("ukb25589_imp_chr1_v3_s487327.sample")
str(sample)
(N <- readBin("data/ukb_imp_chr1_v3.bgen", what = 1L, size = 4, n = 4)[4])
sample <- sample[-1, ]
stopifnot(nrow(sample) == N)

df0 <- readRDS("pheno.rds")
ind.indiv <- match(df0$eid, sample$ID_2)
y <- df0$height
sub <- which(df0$is_caucasian & !is.na(y) & !is.na(ind.indiv))
length(sub)  # 374,131

set.seed(1)
ind.train <- sort(sample(length(sub), 350e3))
ind.test <- setdiff(seq_along(sub), ind.train)

library(bigsnpr)
ukb <- snp_attach("data/UKB_imp_height.rds")
G <- ukb$genotypes
CHR <- as.numeric(ukb$map$chromosome)
POS <- ukb$map$physical.pos
A1 <- ukb$map$allele1
A2 <- ukb$map$allele2
dim(G) # 374,131 x 656,060
file.size(G$backingfile) / 1024^3  # 229 GB

PC <- as.matrix(readRDS("PC.rds"))
COVAR <- cbind(PC, df0$date, df0$sex)[sub, ]

# system.time(
#   gwas <- big_univLinReg(G, y[sub][ind.train], ind.train,
#                          covar.train = COVAR[ind.train, ],
#                          ncores = nb_cores())
# ) # 94 min
# saveRDS(gwas, "gwas_height_train.rds")
gwas <- readRDS("gwas_height_train.rds")
library(ggplot2)
plot(gwas) + scale_x_continuous(breaks = 0:10 / 10)
snp_qq(gwas) + xlim(1, NA)
snp_manhattan(gwas, CHR, POS) +
  geom_hline(yintercept = -log10(5e-8), linetype = 2, color = "red")
lpval <- -predict(gwas)
betas <- gwas$estim


grid <- expand.grid(
  thr.clmp = c(0.05, 0.2, 0.5, 0.8, 0.95),
  chr = unique(CHR)
)

n_thr_pval <- 100
scores <- FBM(length(ind.test), nrow(grid) * n_thr_pval,
              backingfile = "height_stacking")$save()
dim(scores)

set.seed(1)
ind.val <- sort(sample(ind.test, 10e3))
ind.test2 <- setdiff(ind.test, ind.val)


library(doParallel)
registerDoParallel(cl <- makeCluster(11, outfile = ""))
foreach(ic = rows_along(grid)) %dopar% {

  chr      <- grid[ic, "chr"]
  thr.clmp <- grid[ic, "thr.clmp"]

  if (!identical(chr, grid[ic - 1, "chr"]))
    cat("Starting with chromosome", chr, "\n")

  ind.keep <- bigsnpr::snp_clumping(
    G, CHR, S = lpval, thr.r2 = thr.clmp,
    ind.row = ind.val,
    exclude = which(CHR != chr),
    infos.pos = POS, is.size.in.bp = TRUE
  )

  lp.thr <- exp(seq(log(0.1), log(0.99 * max(lpval[ind.keep])),
                    length.out = n_thr_pval))

  prs <- bigsnpr::snp_PRS(G, ind.test = ind.test,
                          betas.keep = betas[ind.keep], ind.keep = ind.keep,
                          lpS.keep = lpval[ind.keep], thr.list = lp.thr)

  scores[, (ic - 1) * n_thr_pval + 1:n_thr_pval] <- prs

  NULL
}
stopCluster(cl)

scores[, 1]
scores[, ncol(scores)]
corr <- big_apply(scores, function(X, ind, y.test) {
  apply(X[, ind], 2, cor, y.test)
}, a.combine = 'c', ncores = 10, y.test = y[sub][ind.test])

plot(corr, col = rep(1:22, 500), pch = 20)


summary(mod.base <- lm(height ~ date + sex, data = df0[sub[ind.val], ]))
pred.base <- predict(mod.base, df0[sub, ])
system.time(
  prs_final <- big_spLinReg(scores, y[sub][ind.val], match(ind.val, ind.test),
                            base.train = pred.base[ind.val],
                            ncores = nb_cores(), alphas = 10^(-(0:8 / 2)))
) # 17 min

str(prs_final)

preds_final <- pred.base[ind.test2] +
  predict(prs_final, scores, match(ind.test2, ind.test))
cor(preds_final, y[sub][ind.test2])
library(dplyr)
cbind.data.frame(pred = preds_final, df0[sub[ind.test2], c("sex", "height")]) %>%
  group_by(sex) %>%
  summarize(cor(pred, height), RMSE = sqrt(mean((pred - height)^2)))
# A tibble: 2 x 3
#   sex    `cor(pred, height)`  RMSE
#   <fct>                <dbl> <dbl>
# 1 Female               0.584  5.05
# 2 Male                 0.580  5.47

cbind.data.frame(pred = preds_final, df0[sub[ind.test2], c("sex", "height")]) %>%
  ggplot(aes(pred, height, color = sex)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  coord_equal() +
  theme_bigstatsr()

# Try another strategy
size <- 500
scores2 <- scores[, 1:500]
for (i in 1:21) inplace::`%+<-%`(scores2, scores[, 1:500 + i * 500])
corr2 <- apply(scores2, 2, cor, y[sub][ind.test])
plot(corr2, pch = 20)
X <- cbind(scores2[match(ind.val, ind.test), ],
           df0$date[sub[ind.val]],
           as.integer(df0$sex[sub[ind.val]]))
mod <- glmnet::glmnet(X, y[sub][ind.val], alpha = 0.01, lower.limits = 0,
                      penalty.factor = c(rep(1, 500), rep(0, 2)))

X2 <- cbind(scores2[match(ind.test2, ind.test), ],
           df0$date[sub[ind.test2]],
           as.integer(df0$sex[sub[ind.test2]]))
preds_glmnet <- predict(mod, X2)
corr3 <- apply(preds_glmnet, 2, cor, y[sub[ind.test2]])
plot(corr3, pch = 20)

cbind.data.frame(pred = preds_glmnet[, 80], df0[sub[ind.test2], c("sex", "height")]) %>%
  group_by(sex) %>%
  summarize(cor(pred, height), RMSE = sqrt(mean((pred - height)^2)))
# # A tibble: 2 x 3
#   sex    `cor(pred, height)`  RMSE
#   <fct>                <dbl> <dbl>
# 1 Female               0.570  5.10
# 2 Male                 0.569  5.50

library(ggplot2)
cbind.data.frame(pred = preds_glmnet[, 80], df0[sub[ind.test2], c("sex", "height")]) %>%
  ggplot(aes(pred, height, color = sex)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  coord_equal() +
  theme_bigstatsr()
