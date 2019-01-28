sumstats_joined <- readRDS("breast-cancer/sumstats_all_chr.rds")
list_snp_id <- split(
  with(sumstats_joined, paste(chr, position_b37, a0, a1, sep = "_")),
  factor(sumstats_joined$chr, levels = 1:22, ordered = TRUE)
)
rm(sumstats_joined)

# system("./ukbgene imp -c1 -m")
sample <- bigreadr::fread2("ukb25589_imp_chr1_v3_s487327.sample")
str(sample)
(N <- readBin("data/ukb_imp_chr1_v3.bgen", what = 1L, size = 4, n = 4)[4])
sample <- sample[-1, ]
stopifnot(nrow(sample) == N)

df0 <- readRDS("pheno.rds")
ind.indiv <- match(df0$eid, sample$ID_2)
y <- df0$has_cancer
y[y] <- NA  ## set to missing all types of cancer
y[df0$cancer_type == "breast cancer"] <- 1  ## keep only BC
sub <- which(df0$sex == "Female" & df0$is_caucasian & !is.na(y) & !is.na(ind.indiv))
table(y[sub])
sub2 <- sort(c(sub[y[sub] == 1], sample(sub[y[sub] == 0], 10e3)))

library(bigsnpr)
system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKB_imp_BC3",
    ind_row = ind.indiv[sub2],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = 10
  )
)

ukb_imp <- snp_attach("data/UKB_imp_BC3.rds")
G <- ukb_imp$genotypes
dim(G)
G[, 1:5]

set.seed(1)
ind.train <- sort(sample(nrow(G), 15e3))
ind.test <- setdiff(rows_along(G), ind.train)

y2 <- y[sub2]

system.time(
  mod <- big_spLogReg(G, y2[ind.train], ind.train, ncores = nb_cores(),
                      alphas = 10^(-(0:4)), return.all = TRUE)
)
attr(mod, "alpha")
str(mod)

pred <- predict(mod, G, ind.test)
AUC(pred, y2[ind.test]) # 60.8
str(mod)

str(mod[[1]])
attr(mod[[1]], "alpha")
plot(mod[[1]][[1]]$loss.val, pch = 20)
plot(mod[[1]][[2]]$loss.val, pch = 20)
plot(mod[[1]][[3]]$loss.val, pch = 20)
plot(mod[[1]][[2]]$loss.val, pch = 20)
plot(mod[[1]][[2]]$loss.val, pch = 20)



system.time(
  mod2 <- big_spLogReg(G, y2[ind.train], ind.train, ncores = 10,
                       alphas = 10^(-4), return.all = TRUE, dfmax = 200e3)
)
plot(mod2[[1]][[1]]$loss.val, pch = 20)
lapply(2:10, function(i) {
  points(mod2[[1]][[i]]$loss.val, pch = 20, col = i)
})

pred2 <- rowMeans(sapply(mod2[[1]], function(mod) {
  predict(mod, G, ind.test)
}))
bigstatsr::AUC(pred2, y2[ind.test]) # 62.1

library(biglasso)
X <- big.matrix(length(ind.train), ncol(G),
                backingfile = "BM_BC.bk", backingpath = "data")
big_apply(G, function(G, ind) {
  X[, ind] <- G[ind.train, ind]
  NULL
})
mod3 <- biglasso(X, y2[ind.train], penalty = "enet", family = "binomial",
                 ncores = 10, alpha = 1e-4, dfmax = 20e3, verbose = TRUE)
str(mod3)
plot(mod3$loss, pch = 20)
library(Matrix)
colSums(mod3$beta != 0)
preds <- big_prodMat(G, mod3$beta[-1, ], ind.test)
aucs <- apply(preds, 2, bigstatsr::AUC, y2[ind.test])
max(aucs)  # 0.6194623
plot(aucs, pch = 20)
plot(mod3$loss, pch = 20)

str(mod[[1]][[1]])
sum(mod[[1]][[1]]$beta != 0)

coef <- -1 / log10(sumstats_joined$bcac_icogs2_P1df_Wald)
hist(coef); summary(coef)

mod4 <- biglasso(X, y2[ind.train], penalty = "enet", family = "binomial",
                 ncores = 10, alpha = 1e-4, dfmax = 20e3, verbose = TRUE,
                 penalty.factor = coef)
