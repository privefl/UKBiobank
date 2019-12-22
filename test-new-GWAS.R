library(bigsnpr)
simu <- snp_attach("data/ukbb4simu.rds")
G <- simu$genotypes
CHR <- as.integer(simu$map$chromosome)
POS <- simu$map$physical.pos
INFO <- readRDS("data/ukbb4simu_info.rds")$info
SD <- readRDS("data/ukbb4simu_stats.rds")$scale
PC <- readRDS("PC_sub.rds")
# load("data/ukbb4simu_ind.RData")
NCORES <- nb_cores()

# G.train <- snp_attach("data/ukbb4simu_train.rds")$genotypes
# G.test  <- snp_attach("data/ukbb4simu_test.rds")$genotypes

#### SIMU PHENO ####

set.seed(1)
h2 <- 0.5; M <- 10e3
set <- sort(sample(ncol(G), size = M))
effects <- rnorm(M, sd = sqrt(h2 / M))
y <- big_prodVec(G, effects / SD[set], ind.col = set)
y <- (y - mean(y)) / sd(y) * sqrt(h2)         ## make sure that var(y) = h2
y <- y + rnorm(nrow(G), sd = sqrt(1 - h2))

#### STANDARD GWAS ####

# system.time(
#   gwas <- big_univLinReg(G, y, covar.train = PC, ncores = NCORES)
# ) # 8230
# saveRDS(gwas, "gwas1.rds")
gwas <- readRDS("gwas1.rds")
plot(gwas) + ggplot2::scale_y_log10()
hist(beta_gwas <- gwas$estim)
lpval <- -predict(gwas)
snp_manhattan(gwas, CHR, POS, npoints = 20e3)

#### PLR ####
system.time(
  mod <- big_spLinReg(G, y, covar.train = PC, ncores = nb_cores(),
                      alphas = c(1, 0.01, 0.0001), n.abort = 2)
) # 31H
summary(mod)
summary(mod)$message
plot(mod)
beta <- head(summary(mod, best.only = TRUE)$beta[[1]], -10)
poly <- lapply(split(cols_along(G), CHR), function(ind) {
  ind2 <- ind[beta[ind] != 0]
  big_prodVec(G, beta[ind2], ind.col = ind2)
})

#### GWAS WITH POLY ####
system.time(
  gwas2 <- big_univLinReg(G, y, covar.train = cbind(PC, Reduce('+', poly[-1])),
                          ind.col = which(CHR == 1), ncores = NCORES)
) # 65 sec

plot(gwas2) + ggplot2::scale_y_log10()
snp_qq(gwas[rows_along(gwas2), ])
snp_qq(gwas2)

library(ggplot2)
ggplot(data.frame(pval1 = -predict(gwas[rows_along(gwas2), ]), pval2 = -predict(gwas2)),
       aes(pval1, pval2)) +
  theme_bigstatsr() +
  geom_point(aes(color = which(CHR == 1) %in% set),
             alpha = 0.4 + 0.6 * (which(CHR == 1) %in% set)) +
  geom_abline(linetype = 3) +
  geom_smooth(method = "lm", linetype = 2) +
  labs(x = "P-value (standard GWAS)", y = "P-value (GWAS with polygenic effect)",
       color = "Causal effect?") +
  theme(legend.position = c(0.25, 0.75))

snp_manhattan(gwas2, CHR[CHR == 1], POS[CHR == 1],
              ind.highlight = intersect(which(CHR == 1), set)) +
  aes(alpha = I(0.2 + 0.8 * (which(CHR == 1) %in% set)))
