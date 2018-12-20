#### MAF > 0.01 & imputation score > 99% ####

library(foreach)
system.time({
  list_snp_id <- foreach(chr = 1:22) %do% {
    cat("Processing chromosome", chr, "..\n")
    file <- paste0("ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
    infos_chr <- bigreadr::fread2(file, showProgress = FALSE)
    infos_chr_sub <- subset(infos_chr, V6 > 0.01 & V8 > 0.99)
    with(infos_chr_sub, paste(chr, V3, V4, V5, sep = "_"))
  }
})
sum(lengths(list_snp_id)) # 5.4M

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

system.time(
  lpval <- bigsnpr::snp_assocBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    y_row = y[sub][ind.train],
    ind_row = ind.indiv[sub][ind.train],
    bgi_dir = "ukb_imp_bgi",
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
    backingfile = "data/UKB_imp_height",
    ind_row = ind.indiv[sub],
    bgi_dir = "ukb_imp_bgi",
    ncores = 12
  )
) # 5H

#### Predictive models ####

library(bigsnpr)
ukb <- snp_attach("data/UKB_imp_height.rds")
G <- ukb$genotypes
dim(G) # 374,131 x 1,501,593
file.size(G$backingfile) / 1024^3  # 523 GB

PC <- as.matrix(readRDS("PC.rds"))
plot(PC[sample(nrow(PC), 10e3), ], pch = 20)  # verif: all
plot(PC[sample(sub, 10e3), ], pch = 20)       # verif: caucasian only


summary(mod.base <- lm(height ~ date + sex, data = df0[sub[ind.train], ]))
# 154.6 cm  +  1 cm every 6.4 year since 1900  +  13.3 cm for males
pred.base <- predict(mod.base, df0[sub, ])

ind.train2 <- sample(ind.train, 5e3)
system.time(
  logit <- big_spLinReg(G, y[sub][ind.train2], ind.train2,
                        covar.train = PC[sub[ind.train2], ],
                        base.train = pred.base[ind.train2],
                        dfmax = 10e3,
                        ncores = 3)
)
