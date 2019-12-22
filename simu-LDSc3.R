library(bigsnpr)
# celiac <- snp_attach("data/celiacQC.rds")
# subset(celiac, ind.row = which(celiac$fam$affection == 1))
celiac <- snp_attach("data/celiacQC_sub1.rds")
snp_writeBed(celiac, "data/celiacQC_sub1.bed")

G <- celiac$genotypes
CHR <- celiac$map$chromosome
POS <- celiac$map$physical.pos
# SD <- big_scale()(G)$scale
NCORES <- nb_cores()

cor.scaling <- function(X, ind.row, ind.col) {
  ms <- big_scale(center = TRUE, scale = TRUE)(X, ind.row, ind.col)
  ms$scale <- ms$scale * sqrt(length(ind.row) - 1)
  ms
}
svd2 <- snp_autoSVD(G, CHR, POS, ncores = NCORES, fun.scaling = cor.scaling)
plot(svd2, type = "scores")
ind.keep <- attr(svd2, "subset")


library(doParallel)
registerDoParallel(cl <- makeCluster(22, outfile = ""))
myld <- foreach(ind = split(cols_along(G), CHR), .combine = 'c') %dopar% {
  print(min(ind))
  cor <- bigsnpr::snp_cor(G, size = 1000, ind.col = ind)
  Matrix::colSums(cor^2)
}
stopCluster(cl)

plot(svd2$v[, 2], myld[ind.keep])


GtU <- big_apply(G, function(X, ind) {
  print(min(ind))
  sc <- big_scale()(X, ind.col = ind)$scale
  sweep(crossprod(X[, ind], svd2$u), 1, sc, '/')
}, a.combine = "rbind")
N <- nrow(G)
eta <- (GtU^2 / N) - (1 - (GtU^2 / N)) / (N - 2)
plot(eta[, c(1, 10)])
p <- eta %*% (svd2$d / sum(svd2$d))
plot(p)

pcadapt <- snp_pcadapt(G, svd2$u, ncores = NCORES)
plot(pcadapt)
snp_qq(pcadapt) + ggplot2::xlim(1, NA)
R2 <- qchisq(predict(pcadapt, log10 = FALSE), df = 10, lower.tail = FALSE) / N
plot(p, R2, pch = 20, xlab = "PC score", ylab = "pcadapt score")
abline(0, 1, col = "red")
chi2 <- qchisq(predict(pcadapt, log10 = FALSE), df = 1, lower.tail = FALSE)
# plot(p, chi2, pch = 20, xlab = "PC score", ylab = "pcadapt X2(1)")
plot(myld, chi2, pch = 20, xlab = "LD score", ylab = "pcadapt X2(1)")
summary(lm(myld ~ R2 + p))

ms <- big_scale()(G)
# set.seed(2)

h2 <- 0.5; m <- 50e3
gamma <- 0.1
is_odd_chr <- (CHR %% 2 == 1)
set <- sort(sample(which(is_odd_chr), size = m))
effects <- rnorm(m, sd = sqrt(h2 / m))
y1 <- big_prodVec(G, effects, ind.col = set,
                  center = ms$center[set],
                  scale = ms$scale[set])
y2 <- y1 + sqrt(gamma) * scale(svd2$u[, 1])
hist(y <- y2 + rnorm(nrow(G), sd = sqrt(1 - h2 - gamma)))
c(var(y1), var(y2), var(y))

system.time(
  gwas <- big_univLinReg(G, y, ncores = NCORES)
) # 6 sec with no covar
X2 <- gwas$score^2
summary(lm(X2 ~ myld))
summary(lm(X2 ~ myld + p))
summary(lm(X2 ~ myld + R2))
summary(lm(X2 ~ myld + R2 + p))

y1 <- X2
x1 <- myld
x1_int <- cbind(x1, 1)
M <- ncol(G); N <- nrow(G)
(hsq_up <- M / N * mean(x1 - 1) / mean(y1))  # 0.33 -> 0.45 with `median(ld)`
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
(step1_int <- reg_int)

x <- myld
y <- X2
yp <- X2 - step1_int

reg_int <- 1

for (i in 1:5) {
  initial_w2 <- WEIGHTS(hsq_up, x, x, N, M, reg_int)
  xw2 <- WEIGHT(initial_w2, x)
  yw2 <- WEIGHT(initial_w2, yp)
  reg_coef2 <- lm(yw2 ~ 0 + xw2)$coef[[1]]
  reg_int <- step1_int
  print(hsq_up <- M / N * reg_coef2)
}

c(reg_int, hsq_up)
