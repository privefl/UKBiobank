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

library(dplyr)
readRDS("tmp-data/public-data-pheno.rds")[-1] %>%
  as_tibble() %>%
  dplyr::arrange(desc(effects)) %>%
  mutate(chr = CHR[set])

ind.chr <- which(CHR == 2)
K <- snp_cor(G, infos.pos = POS[ind.chr], ind.col = ind.chr)

library(Matrix)
apply2_sp <- function(X, FUN) {
  res <- numeric(ncol(X))
  X2 <- as(X, "dgTMatrix")
  tmp <- tapply(X2@x, X2@j, FUN)
  res[as.integer(names(tmp)) + 1] <- tmp
  res
}
ld <- apply2_sp(K, crossprod)
chi2 <- qchisq(-lpval * log(10), df = 1, lower.tail = FALSE, log.p = TRUE)[ind.chr]
beta <- sumstats$beta[ind.chr]

library(dplyr)
data.frame(ld, chi2) %>%
  group_by(cut(ld, c(-Inf, quantile(ld, 1:19 / 20), Inf))) %>%
  summarise_all(mean) %>%
  plot(chi2 ~ ld, data = .)
source("https://raw.githubusercontent.com/privefl/UKBiobank/ed1339b7b5971d63fc1d0f56d68bab4db64b6b26/LDSC-FUN.R")
M <- length(ind.chr); N <- 10e3
(ldsc <- LDSC(chi2 = chi2, ld = ld, M = M, N = N))
slp <- ldsc[[3]] * N / M

diag(K) <- 1 + 1 / slp
new_beta <- as.vector(Matrix::solve(K, beta))
plot(beta, new_beta, pch = 20); abline(0, 1, col = "red")

pred <- big_prodVec(G, new_beta, ind.col = ind.chr)
y <- obj.bigsnp$fam$affection - 1
cor(pred, y)^2  # 0.007585999 / 0.03039965 /
AUC(pred, y) # 55.9% / 61.9% / 59.2%

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
all_keep <- snp_grid_clumping(G, CHR, POS, lpS = lpval, ncores = NCORES)
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
  single_PRS <- rowSums(X[, ind + s * (0:2)])  ## replace by 0:21 in real data
  bigstatsr::AUC(single_PRS, y.train)
}, ind = 1:s, s = s, y.train = y,
a.combine = 'c', block.size = 1, ncores = NCORES)

max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
#     size thr.r2 grp.num thr.imp num   thr.lp       auc
# 1     52   0.95       1       1  25 4.187722 0.6708255
# 2    105   0.95       1       1  26 4.187722 0.6708255
# 3    210   0.95       1       1  27 4.187722 0.6708255
# 4    526   0.95       1       1  28 4.187722 0.6708255
# 5     62   0.80       1       1  21 4.187722 0.6594182
# 6    125   0.80       1       1  22 4.187722 0.6594182
# 7    250   0.80       1       1  23 4.187722 0.6594182
# 8    625   0.80       1       1  24 4.187722 0.6594182
# 9   5000   0.01       1       1   1 4.187722 0.6528066
# 10 10000   0.01       1       1   2 4.187722 0.6528066

file.remove(paste0(tmp, c(".bk", ".rds")))
