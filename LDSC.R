# use a df for (ld, w_ld, N, M)?
# directy use (N / M * ld)?
WEIGHTS <- function(hsq, ld, w_ld, N, M, int) {

  h <- max(0, min(hsq, 1))
  ld <- pmax(ld, 1)
  w_ld <- pmax(w_ld, 1)

  het_w <- 1 / (2 * (int + h * N / M * ld)^2)   # heteroscedasticity weights
  oc_w <- 1 / w_ld                              # overcounting weights
  het_w * oc_w                                  # merged weights
}


# ld <- bigreadr::fread2("tables/LDSCORE.1000G_EUR.tab.gz")

# w <- runif(10)
# x <- rnorm(10)
# y <- rnorm(10)
#
WEIGHT <- function(w, x) {
  w1 <- sqrt(w)
  w_norm <- w1 / sum(w1)
  x * w_norm
}
# xw <- WEIGHT(w, x)
# yw <- WEIGHT(w, y)
# summary(lm(yw ~ xw + 0))
# summary(lm(y ~ x + 0, weights = w))

#### step 1 ####

ld <- myld
chi2 <- gwas0$score^2                                 # 1.058188 / 0.4503596
# chi2 <- gwas$score^2                                  # 1.160764 / 0.8286497
# chi2 <- gwas_bolt$CHISQ_BOLT_LMM_INF                  # 1.080658 / 0.5951337
# chi2 <- gwas_bolt$CHISQ_BOLT_LMM                      # 0.981509 / 0.6645476
# chi2 <- gwas0_g$score^2; ld <- ld[CHR != CHR[ind_g]]  # 1.055765 / -0.002249073
# lambdaGC(gwas0_g) = 1.04112
M <- ncol(G)
N <- nrow(G)
(hsq1 <- M / N * mean(chi2 - 1) / mean(ld))  # 0.33 -> 0.45 with `median(ld)`
(lam_rel <- bigsnpr:::getLambdaGC(gwas0_g)^2)

# MASS::rlm(X2 ~ L2 + offset(off) + 0, data = data.frame(X2 = chi2, L2 = ld, off = 1))
# lmfit <- lm(X2 ~ L2 + offset(off) + 0, data = data.frame(X2 = chi2, L2 = ld, off = 1))
# summary(lmfit)  # 0.111887   0.001157   96.71   <2e-16 ***
# library(sandwich)
# coeftest(lmfit, vcov = vcovHC(lmfit))  # 0.111887   0.002146  52.138 < 2.2e-16 ***
#
# lmfit2 <- lm(X2 ~ L2, data = data.frame(X2 = chi2, L2 = ld))
# summary(lmfit2)
# #             Estimate Std. Error t value Pr(>|t|)
# # (Intercept) 1.880556   0.039380   47.75   <2e-16 ***
# # L2          0.084989   0.001668   50.96   <2e-16 ***
# library(lmtest)
# bptest(lmfit2)
# coeftest(lmfit2, vcov = vcovHC(lmfit2))
# #              Estimate Std. Error t value  Pr(>|t|)
# # (Intercept) 1.8805565  0.0539749  34.841 < 2.2e-16 ***
# # L2          0.0849888  0.0034653  24.526 < 2.2e-16 ***
# lmrobfit <- robustbase::lmrob(X2 ~ L2, data = data.frame(X2 = chi2, L2 = ld))
# summary(lmrobfit)$coef
# #                 Estimate   Std. Error    t value     Pr(>|t|)
# # (Intercept) 0.7299110467 0.0055786019 130.841215 0.000000e+00
# # L2          0.0007829041 0.0001771272   4.420011 9.876282e-06

THR <- 30
x1 <- ld[chi2 < THR]
x1_int <- cbind(x1, 1)
y1 <- chi2[chi2 < THR]

(hsq_up <- hsq1)
reg_int <- 1

for (i in 1:10) {
  initial_w1 <- WEIGHTS(hsq_up, x1, x1, N, M, reg_int)
  xw1 <- WEIGHT(initial_w1, x1_int)
  yw1 <- WEIGHT(initial_w1, y1)
  coefw <- summary(lm(yw1 ~ 0 + xw1))$coef
  reg_coef <- coefw[1, 1]
  reg_int <- coefw[2, 1]
  # coefw <- summary(lm(y1 ~ x1, weights = initial_w1))$coef
  # reg_int <- coefw[1, 1]
  # reg_coef <- coefw[2, 1]
  print(hsq_up <- M / N * reg_coef)
}

# initial_w1 <- WEIGHTS(hsq_up, x1, x1, N, M, reg_int)
# coefw <- summary(lm(y1 ~ x1, weights = initial_w1))$coef
# reg_int <- coefw[1, 1]
# reg_coef <- coefw[2, 1]
# (hsq_up <- M / N * reg_coef)


blocks <- bigstatsr:::CutBySize(length(x1), nb = 200)
n_blocks <- nrow(blocks)
xtx_block_values <- xty_block_values <- list()
for (i in 1:n_blocks) {
  s <- bigstatsr:::seq2(blocks[i, ])
  X <- xw1[s, ]
  xtx_block_values[[i]] <- crossprod(X)
  xty_block_values[[i]] <- crossprod(X, yw1[s])
}
xtx <- Reduce('+', xtx_block_values)
xty <- Reduce('+', xty_block_values)
est <- solve(xtx, xty)

delete_values <- matrix(NA_real_, nrow(blocks), 2)
for (i in rows_along(blocks)) {
  delete_values[i, ] <-
    solve(xtx - xtx_block_values[[i]], xty - xty_block_values[[i]])
}
delete_values
pseudovalues <- sweep((1 - n_blocks) * delete_values, 2, n_blocks * est, '+')

jknife_cov <- cov(pseudovalues) / n_blocks
jknife_var <- diag(jknife_cov)
jknife_se <- sqrt(jknife_var)
jknife_est <- colMeans(pseudovalues)

(step1_int <- jknife_est[2])
((lam_rel - 1) * h2 + 1)
(step1_int_se <- jknife_se[2])
library(lmtest)
lmfit0 <- lm(yw1 ~ 0 + xw1)
coeftest(lmfit0, vcov = vcovHC(lmfit0))
lmfit0 <- lm(y1 ~ 1, weights = initial_w1)
coeftest(lmfit0, vcov = vcovHC(lmfit0))

#### step 2 ####

x <- ld
y <- chi2
yp <- chi2 - step1_int

(hsq_up <- hsq1)
reg_int <- 1

for (i in 1:5) {
  initial_w2 <- WEIGHTS(hsq_up, x, x, N, M, reg_int)
  xw2 <- WEIGHT(initial_w2, x)
  yw2 <- WEIGHT(initial_w2, yp)
  reg_coef2 <- lm(yw2 ~ 0 + xw2)$coef[[1]]
  reg_int <- step1_int
  print(hsq_up <- M / N * reg_coef2)
}


xtx_block_values <- xty_block_values <- list()
for (i in 1:n_blocks) {
  s <- bigstatsr:::seq2(blocks[i, ])
  X <- xw2[s]
  xtx_block_values[[i]] <- crossprod(X)
  xty_block_values[[i]] <- crossprod(X, yw2[s])
}
xtx <- Reduce('+', xtx_block_values)
xty <- Reduce('+', xty_block_values)
est <- drop(solve(xtx, xty))

delete_values <- rep(NA_real_, nrow(blocks))
for (i in rows_along(blocks)) {
  delete_values[i] <-
    solve(xtx - xtx_block_values[[i]], xty - xty_block_values[[i]])
}
delete_values
pseudovalues <- n_blocks * est - (n_blocks - 1) * delete_values

jknife_cov <- cov(as.matrix(pseudovalues)) / n_blocks
jknife_var <- diag(jknife_cov)
jknife_se <- sqrt(jknife_var)
jknife_est <- mean(pseudovalues)


(int <- step1_int)              # 1.5247 /  1.0436
(int_se <- step1_int_se)        # 0.0453 /  0.0097
(h2_est <- M / N * jknife_est)  # 0.3072 / -0.0021
(h2_se <- M / N * jknife_se)    # 0.0297 /  0.0047

(int - 1) / (mean(chi2) - 1)    # 0.2254

gwas0_gc3 <- structure(gwas0, transfo = function(x) {
  attr(gwas0, "transfo")(x) / sqrt(step1_int)
  # should be h2g for binary traits
})

cowplot::plot_grid(
  snp_qq(gwas0_gc[ is_odd_chr, ]) + xlim(1, NA),
  snp_qq(gwas0_gc[!is_odd_chr, ]) + xlim(1, NA),
  snp_qq(gwas0_gc3[ is_odd_chr, ]) + xlim(1, NA),
  snp_qq(gwas0_gc3[!is_odd_chr, ]) + xlim(1, NA),
  scale = 0.95, labels = c("A", "B"), label_size = 20
)
