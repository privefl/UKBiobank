#### MAF > 0.01 & imputation score > 95% ####

library(doParallel)
registerDoParallel(cl <- makeCluster(12))
system.time({
  list_snp_id <- foreach(chr = 1:22) %dopar% {
    # cat("Processing chromosome", chr, "..\n")
    file <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
    infos_chr <- bigreadr::fread2(file, showProgress = FALSE)
    infos_chr_sub <- subset(infos_chr, V6 > 0.01 & V8 > 0.95)
    with(infos_chr_sub, paste(chr, V3, V4, V5, sep = "_"))
  }
})
stopCluster(cl)
sum(lengths(list_snp_id)) # 8.35M

sample <- bigreadr::fread2("ukb25589_imp_chr1_v3_s487327.sample")
str(sample)
(N <- readBin("data/ukb_imp_chr1_v3.bgen", what = 1L, size = 4, n = 4)[4])
sample <- sample[-1, ]
stopifnot(nrow(sample) == N)

df0 <- readRDS("pheno.rds")
ind.indiv <- match(df0$eid, sample$ID_2)
y <- df0$height

sub <- which(df0$is_caucasian & !is.na(y) & !is.na(ind.indiv))

#### Quick GWAS (on training set) to subset variants ####

set.seed(1)
ind.test <- sort(sample(length(sub), 20e3))
ind.train <- setdiff(seq_along(sub), ind.test)

summary(mod.base <- lm(height ~ date + sex, data = df0[sub[ind.train], ]))

system.time(
  lpval <- bigsnpr::snp_assocBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    y_row = mod.base$residuals,
    ind_row = ind.indiv[sub][ind.train],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = 12
  )
) # 10H
sum(lengths(lpval))

hist(-unlist(lpval))

#### Read variants as bigSNP ####

list_snp_id2 <- lapply(seq_along(list_snp_id), function(i) {
  list_snp_id[[i]][-lpval[[i]] > 1]
})
sum(lengths(list_snp_id2))

system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id2,
    backingfile = "data/UKB_imp_height2",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = 12
  )
) # 8.3H

#### Predictive models ####

library(bigsnpr)
ukb <- snp_attach("data/UKB_imp_height2.rds")
G <- ukb$genotypes
dim(G) # 374,131 x 2,929,097
file.size(G$backingfile) / 1024^3  # 1 TB
CHR <- as.numeric(ukb$map$chromosome)

PC <- as.matrix(readRDS("PC.rds"))
plot(PC[sample(nrow(PC), 10e3), ], pch = 20)  # verif: all
plot(PC[sample(sub, 10e3), ], pch = 20)       # verif: caucasian only


summary(mod.base <- lm(height ~ date + sex, data = df0[sub[ind.train], ]))
# 154.6 cm  +  1 cm every 6.4 year since 1900  +  13.3 cm for males
pred.base <- predict(mod.base, df0[sub, ])

chrs <- list(1:2, 3:4, 5:6, 7:9, 10:12, 13:16, 17:22)
sum(lengths(chrs))  # 22
sapply(chrs, function(chrs) sum(CHR %in% chrs))
# [1] 486965 402204 425944 448007 440921 364630 360426

mods <- vector("list", length(chrs))
times <- rep(NA, length(chrs))

library(doParallel)
registerDoParallel(cl <- makeCluster(2))
foreach(chr = 1:22) %dopar% {
  time <- system.time(
    mod <- bigstatsr::big_spLinReg(G, y[sub][ind.train], ind.train,
                                   ind.col = which(CHR == chr),
                                   covar.train = PC[sub[ind.train], ],
                                   base.train = pred.base[ind.train],
                                   alphas = c(1, 0.5, 0.1),
                                   ncores = 10)
  )[3]
  saveRDS(list(time, mod), paste0("mod_height_chr", chr, ".rds"))
}
stopCluster(cl)


times <- sapply(1:22, function(chr) {
  readRDS(paste0("mod_height_chr", chr, ".rds"))[[1]]
})
# 38h in total

preds <- sapply(1:22, function(chr) {
  mod <- readRDS(paste0("mod_height_chr", chr, ".rds"))[[2]]
  pred.base + predict(mod, G, covar.row = PC[sub, ])
})
apply(preds[ind.test, ], 2, function(pred) cor(pred, y[sub[ind.test]]))

df.sub <- df0[sub, ]
df.sub$PC <- I(PC[sub, ])
df.sub$pred <- I(preds)
summary(mylm <- lm(height ~ date + sex + PC + pred, data = df.sub[ind.train, ]))
df.sub$pred_final <- predict(mylm, df.sub)
cor(df.sub$height[ind.test], df.sub$pred_final[ind.test])

library(dplyr)
df.sub[ind.test, ] %>%
  group_by(sex) %>%
  summarise(cor = cor(pred_final, height))

# another strategy
lpval2 <- lapply(seq_along(list_snp_id), function(i) {
  lpval[[i]][-lpval[[i]] > 1]
}) %>%
  unlist()

ind.keep <- snp_clumping(G, CHR,
                         ind.row = sort(sample(ind.train, 20e3)),
                         S = -lpval2,
                         thr.r2 = 0.95,
                         ncores = 12)


which(-lpval2 > 3)
big_cor(G, ind.col = 1:)[]
