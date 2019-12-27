# download.file("https://github.com/privefl/bigsnpr/raw/master/data-raw/public-data.zip",
#               "public-data.zip")
# unzip("public-data.zip")
sumstats <- bigreadr::fread2("tmp-data/public-data-sumstats.txt")
hist(sumstats$p)

test.bed <- "tmp-data/public-data.bed"
library(bigsnpr)
snp_readBed(test.bed)
obj.bigsnp <- snp_attach("tmp-data/public-data.rds")
G <- obj.bigsnp$genotypes
dim(G) # 559 131276
table(CHR <- obj.bigsnp$map$chromosome)
POS <- obj.bigsnp$map$physical.pos
lpval <- -log10(sumstats$p)
NCORES <- nb_cores()

POS2 <- POS + 1e9 * CHR
K <- snp_cor(G, infos.pos = POS2)

library(Matrix)
apply2_sp <- function(X, FUN) {
  res <- numeric(ncol(X))
  X2 <- as(X, "dgTMatrix")
  tmp <- tapply(X2@x, X2@j, FUN)
  res[as.integer(names(tmp)) + 1] <- tmp
  res
}
ld <- apply2_sp(K, crossprod)
chi2 <- qchisq(-lpval * log(10), df = 1, lower.tail = FALSE, log.p = TRUE)
beta <- sumstats$beta

library(dplyr)
data.frame(ld, chi2) %>%
  group_by(cut(ld, c(-Inf, quantile(ld, 1:19 / 20), Inf))) %>%
  summarise_all(mean) %>%
  plot(chi2 ~ ld, data = .)
source("https://raw.githubusercontent.com/privefl/UKBiobank/ed1339b7b5971d63fc1d0f56d68bab4db64b6b26/LDSC-FUN.R")
M <- ncol(G); N <- 10e3
(ldsc <- LDSC(chi2 = chi2, ld = ld, M = M, N = N))
slp <- ldsc[[3]] * N / M

diag(K) <- 1 + 1 / slp
# diag(K) <- 50
new_beta <- as.vector(Matrix::solve(K, beta))
plot(beta, new_beta, pch = 20); abline(0, 1, col = "red")

pred <- big_prodVec(G, new_beta)
y <- obj.bigsnp$fam$affection
cor(pred, y)^2  # 0.007585999 / 0.03039965
AUC(pred, y - 1) #  / 61.9%


# remotes::install_github("baolinwu/MTAR")
library(MASS)
MTAR::SHvr(Z = sqrt(chi2), ld, N)




corr <- sapply(10^(-6:0), function(alpha) {

  K <- snp_cor(G, infos.pos = POS[ind.chr], ind.col = ind.chr, alpha = alpha)

  ld <- apply2_sp(K, crossprod)
  chi2 <- qchisq(-lpval * log(10), df = 1, lower.tail = FALSE, log.p = TRUE)[ind.chr]
  beta <- sumstats$beta[ind.chr]

  M <- length(ind.chr); N <- 10e3
  (ldsc <- LDSC(chi2 = chi2, ld = ld, M = M, N = N))
  slp <- ldsc[[3]] * N / M

  diag(K) <- 1 + 1 / slp
  new_beta <- as.vector(Matrix::solve(K, beta))

  pred <- big_prodVec(G, new_beta, ind.col = ind.chr)
  y <- obj.bigsnp$fam$affection - 1
  AUC(pred, y)
})
corr
# 0.6194811 0.6193553 0.6192138 0.6190723 0.6189937 0.6189780 0.6189465

# Clumping
all_keep <- snp_grid_clumping(G, CHR, POS, lpS = lpval, ncores = NCORES,
                              grid.base.size = 200)
attr(all_keep, "grid")
# Thresholding
beta <- sumstats$beta
tmp <- tempfile(tmpdir = "tmp-data")
multi_PRS <- snp_grid_PRS(G, all_keep, beta, lpval,
                          backingfile = tmp,
                          n_thr_lpS = 50, ncores = NCORES)
dim(multi_PRS)  ## 559 x 1050

library(tidyverse)
grid2 <- attr(all_keep, "grid") %>%
  mutate(thr.lp = list(attr(multi_PRS, "grid.lpS.thr")), num = row_number()) %>%
  unnest()

## Warning: `cols` is now required.
## Please use `cols = c(thr.lp)`

s <- nrow(grid2)
grid2$auc <- big_apply(multi_PRS, a.FUN = function(X, ind, s, y.train) {
  # Sum over all chromosomes, for the same C+T parameters
  print(single_PRS <- rowSums(X[, ind, drop = FALSE]))  ## replace by 0:21 in real data
  bigstatsr::AUC(single_PRS, y.train)
}, ind = 1:s, s = s, y.train = y,
a.combine = 'c', block.size = 1, ncores = NCORES)

max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
#     size thr.r2 grp.num thr.imp num    thr.lp       auc
# 1  20000   0.01       1       1   1 0.2076388 0.5946541
# 2  20000   0.01       1       1   1 0.2252008 0.5944025
# 3   4000   0.05       1       1   2 1.7143646 0.5942767
# 4  20000   0.01       1       1   1 0.0999900 0.5942138
# 5  20000   0.01       1       1   1 0.1084471 0.5940881
# 6  20000   0.01       1       1   1 0.5072046 0.5940881
# 7  20000   0.01       1       1   1 0.1627513 0.5940252
# 8  20000   0.01       1       1   1 0.8953906 0.5940252
# 9  20000   0.01       1       1   1 0.1176195 0.5940094
# 10 20000   0.01       1       1   1 0.1383572 0.5940094

file.remove(paste0(tmp, c(".bk", ".rds")))
