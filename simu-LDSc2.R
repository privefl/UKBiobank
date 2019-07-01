library(bigsnpr)
celiac <- snp_attach("data/celiacQC_sub1.rds")
G <- celiac$genotypes
CHR <- celiac$map$chromosome
POS <- celiac$map$physical.pos

dim(G)
N <- nrow(G)
NCORES <- nb_cores()

svd <- snp_autoSVD(G, CHR, POS, fun.scaling = big_scale(), ncores = NCORES)
ind.keep <- attr(svd, "subset")

system.time(K <- big_tcrossprodSelf(G, big_scale(), ind.col = ind.keep)) # 6 min
system.time(eigs <- eigen(K[]))# 6 min
all.equal(tcrossprod(sweep(eigs$vectors, 2, sqrt(eigs$values), '*')), K[])
rbind(svd$d^2, head(eigs$values, 10))

ms <- big_scale()(G)
GtU <- big_cprodMat(G, eigs$vectors / N, center = ms$center, scale = ms$scale,
                    ncores = NCORES)
dim(GtU)
ld <- drop(GtU^2 %*% eigs$values)

library(doParallel)
registerDoParallel(cl <- makeCluster(NCORES, outfile = ""))
myld <- foreach(ind = split(cols_along(G), CHR), .combine = 'c') %dopar% {
  print(CHR[min(ind)])
  cor <- bigsnpr::snp_cor(G, size = 1000, ind.col = ind, infos.pos = POS[ind],
                          alpha = 0.1) # change alpha? use 1?
  Matrix::colSums(cor^2)
}
stopCluster(cl)
# cor <- bigsnpr::snp_cor(G, size = 1000, alpha = 1) # change alpha? use 1?
# myld <- Matrix::colSums(cor^2)
plot(ld, myld, pch = 20); abline(0, 1, col = "red")

nPC <- 6
GtU_part <- big_cprodMat(G, eigs$vectors[, 1:nPC] / N, ncores = NCORES,
                         center = ms$center, scale = ms$scale)
dim(GtU_part)
ld_part <- drop(GtU_part^2 %*% eigs$values[1:nPC])
plot(ld, ld_part)
lm(ld ~ ld_part)
ld2 <- drop((svd$v[, 1:nPC] / N)^2 %*% svd$d[1:nPC]^4)
all.equal(ld2, ld_part[ind.keep])
cor(ld_part, ld)   # 0.98
cor(myld, ld)      # 0.08
cor(myld, ld_part) # 0.05
cor(ld - ld_part, myld) # 0.19
cor(ld - ld_part, ld)   # 0.04
cor(ld - ld_part, ld_part) # -0.15

plot(ld - ld_part, myld, pch = 20); abline(0, 1, col = "red")
plot(ld - ld_part, ld, pch = 20); abline(0, 1, col = "red")


# Simulate some phenotypes
M <- ncol(G)
h2 <- 0.5; m <- 10e3
gamma <- 0.2
set <- sample(M, size = m)
effects <- rnorm(m, sd = sqrt(h2 / m))
y1 <- big_prodVec(G, effects, ind.col = set,
                  center = ms$center[set],
                  scale = ms$scale[set])
y2 <- y1 + sqrt(gamma) * scale(svd$u[, 1])
hist(y <- y2 + rnorm(nrow(G), sd = sqrt(1 - h2 - gamma)))
c(var(y1), var(y2), var(y))

X2 <-  big_univLinReg(G, y)$score^2
# X2 <-  big_univLinReg(G, y, covar.train = svd$u, ncores = 4)$score^2

# Implementation of LDSC reg in R -> should verify, but based on
# https://helda.helsinki.fi/bitstream/handle/10138/273501/MT_ldsc_hautakangas.pdf
source('LDSC-FUN.R')

LDSC(chi2 = X2, ld = myld,         M = M, N = N)
LDSC(chi2 = X2, ld = myld,         M = M, N = N, THR = 100)
LDSC(chi2 = X2, ld = ld,           M = M, N = N)
# LDSC(chi2 = X2, ld = ld - ld_part, M = M, N = N)
# LDSC(chi2 = X2, ld = myld - ld_part, M = M, N = N)
# LDSC(chi2 = X2, ld = ld_part, M = M, N = N)

LDSC2(chi2 = X2, ld = myld,         ld_part = ld_part, M = M, N = N)
LDSC2(chi2 = X2, ld = myld,         ld_part = ld_part, M = M, N = N, THR = 100)
LDSC2(chi2 = X2, ld = ld,           ld_part = ld_part, M = M, N = N)
# LDSC2(chi2 = X2, ld = ld - ld_part, ld_part = ld_part, M = M, N = N)
# LDSC2(chi2 = X2, ld = myld - ld_part, ld_part = ld_part, M = M, N = N)
