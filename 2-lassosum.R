sample <- bigreadr::fread2("ukb25589_imp_chr1_v3_s487327.sample")
str(sample)
(N <- readBin("data/ukb_imp_chr1_v3.bgen", what = 1L, size = 4, n = 4)[4])
sample <- sample[-1, ]
stopifnot(nrow(sample) == N)

df0 <- readRDS("pheno.rds")
ind.indiv <- match(df0$eid, sample$ID_2)
y <- df0$height
sub <- which(df0$is_caucasian & !is.na(y) & !is.na(ind.indiv))
length(sub)  # 374,131

set.seed(1)
ind.train <- sort(sample(length(sub), 350e3))
ind.test <- setdiff(seq_along(sub), ind.train)

library(bigsnpr)
ukb <- snp_attach("data/UKB_imp_height.rds")
G <- ukb$genotypes
CHR <- as.numeric(ukb$map$chromosome)
POS <- ukb$map$physical.pos
A1 <- ukb$map$allele1
A2 <- ukb$map$allele2
dim(G) # 374,131 x 656,060
file.size(G$backingfile) / 1024^3  # 229 GB

PC <- as.matrix(readRDS("PC.rds"))
COVAR <- cbind(PC, df0$date, df0$sex)[sub, ]

system.time(
  gwas <- big_univLinReg(G, y[sub][ind.train], ind.train,
                         covar.train = COVAR[ind.train, ],
                         ncores = nb_cores())
) # 94 min
# saveRDS(gwas, "gwas_height_train.rds")
plot(gwas)
snp_manhattan(gwas, CHR, POS, npoints = 20e3)

summary(mod.base <- lm(height ~ date + sex, data = df0[sub[ind.train], ]))
# 154.6 cm  +  1 cm every 6.4 year since 1900  +  13.3 cm for males
pred.base <- predict(mod.base, df0[sub, ])
y2 <- y[sub] - pred.base

# devtools::install_github("tshmak/lassosum")
library(lassosum)
cor <- p2cor(p = predict(gwas, log10 = FALSE),
             n = length(ind.train),
             sign = gwas$estim)

ukb$map$genetic.dist <- 0
ukb$fam <- data.frame(affection = y2, "family.ID" = NA, "sample.ID" = NA,
                      "paternal.ID" = NA, "maternal.ID" = NA, "sex" = NA)
unlink("lassosum/*")
set.seed(1)
size_val <- c(1, 2, 5, 10) * 1000
for (i in seq_along(size_val)) {
  ind.val <- sort(sample(ind.test, size_val[i]))
  snp_writeBed(ukb, paste0("lassosum/val", i, ".bed"), ind.row = ind.val)
}

snp_writeBed(ukb, "lassosum/test1.bed", ind.row = ind.test.split[[1]])
# ind.test.split <- split(ind.test, sample(1:2, length(ind.test), TRUE))
# snp_writeBed(ukb, "lassosum/test1.bed", ind.row = ind.test.split[[1]])
# ukb$fam$affection[] <- 9999      ## just to be sure that this is not used
# snp_writeBed(ukb, "lassosum/test2.bed", ind.row = ind.test.split[[2]])

library(doParallel)
registerDoParallel(cl <- makeCluster(nb_cores()))
system.time(
  out <- lassosum.pipeline(cor = cor, chr = CHR, pos = POS, A1 = A2, A2 = A1,
                           test.bfile = "lassosum/val1",
                           cluster = cl,
                           s = c(0.2, 0.1, 0.05),
                           exclude.ambiguous = FALSE,
                           LDblocks = "EUR.hg19")
) # 98 min
stopCluster(cl)

v <- validate(out, covar = COVAR[ind.val1, ])
v$best.validation.result  # 0.588

ind <- which(v$best.beta != 0)
pred <- big_prodVec(G, -v$best.beta[ind],
                    ind.row = ind.test,
                    ind.col = ind)
plot(pred, y2[ind.test], pch = 20)
cor(pred, y2[ind.test])
