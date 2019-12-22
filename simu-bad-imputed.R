library(bigsnpr)
simu <- snp_attach("data/ukbb4simu.rds")
G <- simu$genotypes
CHR <- as.integer(simu$map$chromosome)
POS <- simu$map$physical.pos
info_final <- readRDS("data/ukbb4simu_info.rds")
INFO <- info_final$info
SD <- readRDS("data/ukbb4simu_stats.rds")$scale
PC <- readRDS("PC_sub.rds")
load("data/ukbb4simu_ind.RData")
NCORES <- nb_cores()

set <- snp_clumping(G, CHR, infos.pos = POS, exclude = which(INFO < 1),
                    ncores = NCORES)
hist(INFO[set])

set.seed(1)
h2 <- 0.5
M <- length(set)
effects <- rnorm(M, sd = sqrt(h2 / M))
y <- big_prodVec(G, effects / SD[set], ind.col = set)
y <- (y - mean(y)) / sd(y) * sqrt(h2)         ## make sure that var(y) = h2
y <- y + rnorm(nrow(G), sd = sqrt(1 - h2))



G.train <- snp_attach("data/ukbb4simu_train.rds")$genotypes
system.time(
  gwas <- big_univLinReg(G.train, y[ind.train],
                         covar.train = PC[ind.train, ],
                         ncores = NCORES)
) # 209 / 29
plot(gwas)
snp_manhattan(gwas, CHR, POS, npoints = 20e3, ind.highlight = set)
plot(gwas[set, ])
plot(INFO[set], -predict(gwas[set, ]), pch = 20)
summary(lm(-predict(gwas[set, ]) ~ INFO[set]))
plot(abs(effects), -predict(gwas[set, ]), pch = 20)
summary(lm(-predict(gwas[set, ]) ~ abs(effects)))



library(bigreadr)
sample <- fread2("ukb25589_imp_chr1_v3_s487327.sample")
(N <- readBin("data/ukb_imp_chr1_v3.bgen", what = 1L, size = 4, n = 4)[4])
sample <- sample[-1, ]
stopifnot(nrow(sample) == N)

csv <- "ukb22544.csv"
df0 <- fread2(csv, select = c("eid", "22006-0.0"),
              col.names = c("eid", "is_caucasian"))
# sample still in data
ind.indiv <- match(df0$eid, sample$ID_2)
# related samples
rel <- fread2("ukb25589_rel_s488346.dat")
df0$is_rel2 <- df0$eid %in% rel$ID2
# + keep caucasian only
sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian)
length(sub)  # 335,609

list_snp <- split(with(info_final, paste(chr, pos, a1, a2, sep = "_"))[set], CHR[set])
# Read all data as random hard calls
system.time(
  snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = lapply(list_snp, rep, each = 20),
    backingfile = "data/ukbb4simu_bad20",
    ind_row = ind.indiv[sub][ind.train],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES,
    read_as = "random"
  )
) # 4H

set2 <- snp_clumping(G, CHR, infos.pos = POS, exclude = which(INFO > 0.6),
                     ncores = NCORES)
dim(G.train2 <- snp_attach("data/ukbb4simu_bad20.rds")$genotypes)
G.train2[, 1:3]
system.time(
  gwas2 <- big_univLinReg(G.train2, y[ind.train],
                          covar.train = PC[ind.train, ],
                          ncores = NCORES)
) # 209 / 29
plot(gwas2)
snp_manhattan(gwas2, rep(CHR[set2], each = 20), rep(POS[set2], each = 20), npoints = 20e3)

gwas2[1:20, ]
gwas[set2[1], ]
beta <- matrix(gwas2$estim, nrow = 20)
apply(beta[, 1:5], 2, dplyr::cummean)
beta_MI <- colMeans(beta)
std <- matrix(gwas2$std.err, nrow = 20)
var_W <- colMeans(std^2)
var_B <- colSums(sweep(beta, 2, beta_MI, '-')^2) / 19
var_MI <- var_W + 21 / 20 * var_B
df1 <- 19 * (1 + 20 / 21 * var_W / var_B)^2
df2 <- 19 * (1 + 1 / 20 * var_W / var_B)^2

lpval <- pt(abs(beta_MI / sqrt(var_MI)), df = df1, lower.tail = FALSE, log.p = TRUE)
lpval <- (log(2) + lpval) / log(10)
hist(10^lpval)

# ?mitools::MIcombine
# test <- sapply(1:20, function(k) {
#   tmp <- mitools::MIcombine(as.list(beta[, k]), as.list(std[, k]^2))
#   stopifnot(all.equal(c(beta_MI[k], var_MI[k]),
#                       c(tmp$coefficients, tmp$variance)))
#   c(df1[k], tmp$df, tmp$missinfo)
# })
# plot(test[1, ], log(test[2, ]), pch = 20)
# abline(lm(log(test[2, ]) ~ test[1, ]), col = "red")
# summary(lm(test[2, ] ~ poly(test[1, ], 4)))
# hist(pchisq((beta_MI / sqrt(var_MI))[]^2, df = 1, lower.tail = FALSE))
# lpval <- log10(pchisq((beta_MI / sqrt(var_MI))[]^2, df = 1, lower.tail = FALSE))
# df1 <- 19 * (1 + 20 / 21 * var_W / var_B)^2
# df2 <- 19 * (1 + 1 / 20 * var_W / var_B)^2
# df3 <- 19 * (1 + 1 / 20 * var_W / var_MI)^2
# all.equal(test[2, ], df1[1:20])
# plot(test[3, ], INFO[set2[1:20]])

plot(gwas$estim[set2], beta_MI, pch = 20); abline(0, 1, col = "red")
plot(var_B, var_W, pch = 20)
plot(gwas$std.err[set2]^2, var_MI, pch = 20); abline(0, 1, col = "red")
plot(-predict(gwas[set2, ]), -lpval, pch = 20); abline(0, 1, col = "red")
plot(INFO[set2] * gwas$estim[set2], beta_MI, pch = 20); abline(0, 1, col = "red")
plot(-predict(gwas[set2, ]) * INFO[set2], -lpval, pch = 20); abline(0, 1, col = "red")

library(ggplot2)
cowplot::plot_grid(
  qplot(gwas$estim[set2], beta_MI, color = INFO[set2]) +
    theme_bigstatsr() +
    scale_colour_viridis_c() +
    geom_abline(color = "red") +
    geom_abline(slope = 0.6, linetype = 2, color = "red") +
    geom_abline(slope = 0.3, linetype = 2, color = "red") +
    labs(x = "Effect size GWAS", y = "Effect size MI", color = "INFO") +
    theme(legend.position = c(0.8, 0.2)),


  qplot(-predict(gwas[set2, ]), -lpval, color = INFO[set2]) +
    theme_bigstatsr() +
    scale_colour_viridis_c() +
    geom_abline(color = "red") +
    geom_abline(slope = 0.6, linetype = 2, color = "red") +
    geom_abline(slope = 0.3, linetype = 2, color = "red") +
    labs(x = "-log10(pval) GWAS", y = "-log10(pval) MI", color = "INFO") +
    theme(legend.position = "none"),

  scale = 0.95

)

