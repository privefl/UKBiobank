library(doParallel)
registerDoParallel(cl <- makeCluster(12))
system.time({
  list_snp_id <- foreach(chr = 1:22) %dopar% {
    # cat("Processing chromosome", chr, "..\n")
    mfi <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
    infos_chr <- bigreadr::fread2(mfi, showProgress = FALSE)
    infos_chr_sub <- subset(infos_chr, V6 > 0.01, V8 > 0.99)
    bim <- paste0("data/ukb_snp_bim/ukb_snp_chr", chr, "_v2.bim")
    map_chr <- bigreadr::fread2(bim)
    joined <- dplyr::inner_join(
      infos_chr_sub, map_chr,
      by = c("V3" = "V4", "V4" = "V5", "V5" = "V6")
    )
    with(joined, paste(chr, V3, V4, V5, sep = "_"))
  }
}) # 80 sec
stopCluster(cl)

sum(lengths(list_snp_id))  # 656,060


# subset samples
library(bigreadr)
sample <- fread2("ukb25589_imp_chr1_v3_s487327.sample")
str(sample)
(N <- readBin("data/ukb_imp_chr1_v3.bgen", what = 1L, size = 4, n = 4)[4])
sample <- sample[-1, ]
stopifnot(nrow(sample) == N)

csv <- "ukb22544.csv"
df0 <- fread2(csv, select = "eid")
ind.indiv <- match(df0$eid, sample$ID_2)
sub <- which(!is.na(ind.indiv))

system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_full",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = 22
  )
) # 2.2H



library(bigsnpr)
NCORES <- nb_cores()
simu <- snp_attach("data/UKBB_full.rds")
G <- simu$genotypes
CHR <- as.integer(simu$map$chromosome)
POS <- simu$map$physical.pos
SD <- big_parallelize(G, function(X, ind) {
  big_scale()(X, ind.col = ind)$scale
}, p.combine = 'c', ncores = NCORES)
nPC <- 10
PC_all <- fread2(csv, select = paste0("22009-0.", 1:nPC),
                 col.names = paste0("PC", 1:nPC))
PC <- as.matrix(PC_all)[sub, ]

#### SIMU PHENO ####
set.seed(1)
h2 <- 0.5; M <- 10e3
is_odd_chr <- (CHR %% 2 == 1)
set <- sort(sample(which(is_odd_chr), size = M))
effects <- rnorm(M, sd = sqrt(h2 / M))
y <- big_prodVec(G, effects / SD[set], ind.col = set)
y1 <- (y - mean(y)) / sd(y) * sqrt(h2)         ## make sure that var(y) = h2
hist(y <- y1 + rnorm(nrow(G), sd = sqrt(1 - h2)))
# y <- (y > 1) + 0L

#### STANDARD GWAS ####
system.time(
  gwas0 <- big_univLinReg(G, y, covar.train = PC, ncores = NCORES)
  # gwas0 <- big_univLogReg(G, y, covar.train = PC, ncores = NCORES)
) # 1 min -> 17 min
library(ggplot2)
plot(gwas0) + scale_y_log10()
snp_manhattan(gwas0, CHR, POS, npoints = 20e3)
cowplot::plot_grid(
  snp_qq(gwas0[ is_odd_chr, ]) + xlim(1, NA),
  snp_qq(gwas0[!is_odd_chr, ]) + xlim(1, NA),
  scale = 0.95, labels = c("A", "B"), label_size = 20
)
# 1.03 if thr = 1

#### PLR ####
system.time(
  mod <- big_spLinReg(G, y, covar.train = PC, ncores = nb_cores(),
                      alphas = c(1, 0.1), n.abort = 3, nlam.min = 20)
) # 99 min
# saveRDS(mod, "mod-simu-GWAS.rds")
# mod <- readRDS("mod-simu-GWAS.rds")
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
  gwas_all_chr <- lapply(1:22, function(chr) {
    big_univLinReg(G, y, ind.col = which(CHR == chr),
                   covar.train = cbind(PC, Reduce('+', poly[-chr])),
                   ncores = NCORES)
  })
) # 6 min
gwas <- do.call("rbind", gwas_all_chr)

library(ggplot2)
snp_manhattan(gwas0, CHR, POS, npoints = 20e3)
snp_manhattan(gwas,  CHR, POS, npoints = 20e3)
cowplot::plot_grid(
  snp_qq(gwas[ is_odd_chr, ]) + xlim(1, NA),
  snp_qq(gwas[!is_odd_chr, ]) + xlim(1, NA),
  scale = 0.95, labels = c("A", "B"), label_size = 20
)

#### LDSC ####
# download.file("https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2",
#               destfile = "eur_w_ld_chr.tar.bz2")
# untar("eur_w_ld_chr.tar.bz2")
ldsc <- fread2(glue::glue("eur_w_ld_chr/{chr}.l2.ldscore.gz", chr = 1:22))

library(dplyr)
map_ldsc <- simu$map %>%
  select(CHR = chromosome, BP = physical.pos) %>%
  mutate(X20 = gwas0$score^2, X2 = gwas$score^2, myld = myld) %>%
  mutate(CHR = as.integer(CHR)) %>%
  left_join(ldsc)
map_ldsc
summary(mylm <- lm(X2 ~ L2, data = map_ldsc))
mylm$coefficients[["L2"]] * ncol(G) / nrow(G)
# plot(map_ldsc[c("L2", "X2")], pch = 20, col = scales::alpha("black", 0.2))
# abline(mylm, col = "red")

summary(mylm0 <- lm(X20 ~ L2, data = map_ldsc))
mylm0$coefficients[["L2"]] * ncol(G) / nrow(G)
# plot(map_ldsc[c("L2", "X20")], pch = 20, col = scales::alpha("black", 0.2),
#      xlab = "LD score", ylab = "X2 stat from standard GWAS")
# abline(mylm0, col = "red")

gwas_gc <- structure(gwas0, transfo = function(x) {
  attr(gwas, "transfo")(x) / sqrt(lm(X2 ~ L2, data = map_ldsc)$coeff[[1]])
})
gwas0_gc <- structure(gwas0, transfo = function(x) {
  attr(gwas0, "transfo")(x) / sqrt(lm(X20 ~ L2, data = map_ldsc)$coeff[[1]])
})
cowplot::plot_grid(
  snp_qq(gwas_gc[ is_odd_chr, ]) + xlim(1, NA),
  snp_qq(gwas_gc[!is_odd_chr, ]) + xlim(1, NA),
  snp_qq(gwas0_gc[ is_odd_chr, ]) + xlim(1, NA),
  snp_qq(gwas0_gc[!is_odd_chr, ]) + xlim(1, NA),
  scale = 0.95
)

#### Estim pred ####
set.seed(1)
ind.train <- sample(nrow(G), 50e3)
ind.test <- setdiff(rows_along(G), ind.train)
system.time(
  mod2 <- big_spLinReg(G, y[ind.train], ind.train = ind.train,
                       covar.train = PC[ind.train, ],
                       alphas = c(1, 0.1, 0.01, 0.001),
                       ncores = nb_cores(),
                       dfmax = 100e3, n.abort = 3)
) # 7H
plot(mod2)
summary(mod2)
pred <- predict(mod2, G, ind.test, covar.row = PC[ind.test, ])
cor(pred, y[ind.test])
# M = 50K -> r=0.343 -> v=11.8% -> 13.4% power increase
# M = 10K -> r=0.427 -> v=18.2% -> 22.3% power increase
