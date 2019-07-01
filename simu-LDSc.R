# devtools::install_github("privefl/bigsnpr")
library(bigsnpr)

# Read data in my format (matrix format)
snp_readBed("data/celiacQC_sub1.bed")
celiac <- snp_attach("data/celiacQC_sub1.rds")
G <- celiac$genotypes
CHR <- celiac$map$chromosome
POS <- celiac$map$physical.pos
NCORES <- nb_cores()
N <- nrow(G); M <- ncol(G)

# Compute PCA with same scaling as correlation
cor.scaling <- function(X, ind.row, ind.col) {
  ms <- big_scale(center = TRUE, scale = TRUE)(X, ind.row, ind.col)
  ms$scale <- ms$scale * sqrt(length(ind.row) - 1)
  ms
}
svd <- snp_autoSVD(G, CHR, POS, ncores = NCORES, fun.scaling = cor.scaling)
plot(svd)
plot(svd, type = "scores")
ind.keep <- attr(svd, "subset")

# Compute LD scores
library(doParallel)
registerDoParallel(cl <- makeCluster(NCORES, outfile = ""))
myld <- foreach(ind = split(cols_along(G), CHR), .combine = 'c') %dopar% {
  print(CHR[min(ind)])
  cor <- bigsnpr::snp_cor(G, size = 1000, ind.col = ind, infos.pos = POS[ind],
                          alpha = 0.1) # change alpha? use 1?
  Matrix::colSums(cor^2)
}
stopCluster(cl)


ms <- big_scale()(G)

# Simulate some phenotypes
h2 <- 0.5; m <- 50e3
gamma <- 0.0
set <- sample(M, size = m)
effects <- rnorm(m, sd = sqrt(h2 / m))
y1 <- big_prodVec(G, effects, ind.col = set,
                  center = ms$center[set],
                  scale = ms$scale[set])
y2 <- y1 + sqrt(gamma) * scale(svd$u[, 1])
hist(y <- y2 + rnorm(nrow(G), sd = sqrt(1 - h2 - gamma)))
c(var(y1), var(y2), var(y))

system.time(
  gwas <- big_univLinReg(G, y, ncores = NCORES)
) # 6 sec with no covar


# Implementation of LDSC reg in R -> should verify, but based on
# https://helda.helsinki.fi/bitstream/handle/10138/273501/MT_ldsc_hautakangas.pdf
source('LDSC-FUN.R')

LDSC(chi2 = gwas$score^2, ld = myld, M = M, N = N)
