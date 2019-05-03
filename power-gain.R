
beta <- summary(mod, best.only = TRUE)$beta[[1]]
beta[CHR == 1] <- 0
ind <- which(head(beta, -20) != 0)
pred <- pred.base[ind.train] +
  big_prodVec(G, beta[ind], ind.train, ind) +
  drop(PC[sub[ind.train], ] %*% tail(beta, 20))


system.time(
  gwas2 <- big_univLinReg(G, y[sub][ind.train], ind.train, which(CHR == 1),
                          covar.train = cbind(COVAR[sub[ind.train], ], pred),
                          ncores = nb_cores())
) # 6 min

plot(gwas2)
library(ggplot2)
snp_qq(gwas2)
snp_qq(gwas) + xlim(1, NA)
plot(gwas)
plot(-predict(gwas)[CHR == 1], -predict(gwas2), pch = 20); abline(0, 1, col = "red")
plot(gwas$estim[CHR == 1], gwas2$estim, pch = 20); abline(0, 1, col = "red")

pred2 <- pred.base[ind.train] +
  predict(mod, G, ind.train, covar.row = PC[sub[ind.train], ])
system.time(
  gwas3 <- big_univLinReg(G, y[sub][ind.train], ind.train, which(CHR == 1),
                          covar.train = cbind(COVAR[sub[ind.train], ], pred2),
                          ncores = nb_cores())
) # 6 min

plot(gwas3)
library(ggplot2)
snp_qq(gwas3)
snp_qq(gwas) + xlim(1, NA)
plot(gwas)
plot(-predict(gwas)[CHR == 1], -predict(gwas3), pch = 20); abline(0, 1, col = "red")
plot(gwas$estim[CHR == 1], gwas3$estim, pch = 20); abline(0, 1, col = "red")

which(-predict(gwas3) > 15)
plot(as.factor(G[, 27858]), y[sub], pch = 20)
plot(as.factor(G[ind.train, 27858]), y[sub][ind.train] - pred2, pch = 20)

beta <- summary(mod, best.only = TRUE)$beta[[1]]
