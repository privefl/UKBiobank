library(doParallel)
registerDoParallel(cl <- makeCluster(12))
system.time({
  list_snp_id <- foreach(chr = 1:22) %dopar% {
    # cat("Processing chromosome", chr, "..\n")
    mfi <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
    infos_chr <- bigreadr::fread2(mfi, showProgress = FALSE)
    infos_chr_sub <- subset(infos_chr, V6 > 0.01)  ## MAF > 1%
    bim <- paste0("data/ukb_snp_bim/ukb_snp_chr", chr, "_v2.bim")
    map_chr <- bigreadr::fread2(bim)
    joined <- dplyr::inner_join(
      infos_chr_sub, map_chr,
      by = c("V3" = "V4", "V4" = "V5", "V5" = "V6")
    )
    with(joined, paste(chr, V3, V4, V5, sep = "_"))
  }
}) # < 2 min
stopCluster(cl)

sum(lengths(list_snp_id))  # 656,060

sample <- bigreadr::fread2("ukb25589_imp_chr1_v3_s487327.sample")
str(sample)
(N <- readBin("data/ukb_imp_chr1_v3.bgen", what = 1L, size = 4, n = 4)[4])
sample <- sample[-1, ]
stopifnot(nrow(sample) == N)

df0 <- readRDS("pheno.rds")
ind.indiv <- match(df0$eid, sample$ID_2)
y <- ifelse(df0$has_breast_cancer, 1, ifelse(df0$has_cancer, NA, 0))

sub <- which(df0$sex == "Female" & df0$is_caucasian & !is.na(y) & !is.na(ind.indiv))
length(sub)  # 188,628
table(y[sub])
#      0      1
# 180094   8534


system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKB_imp_BC",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = 10
  )
) # 2H


set.seed(1)
ind.train <- sort(sample(length(sub), 150e3))
ind.test <- setdiff(seq_along(sub), ind.train)
table(y[sub][ind.train])
#      0      1
# 143232   6768

library(bigsnpr)
ukb <- snp_attach("data/UKB_imp_BC.rds")
G <- ukb$genotypes
CHR <- as.numeric(ukb$map$chromosome)
dim(G) # 188,628 x 656,060
file.size(G$backingfile) / 1024^3  # 115 GB

PC <- as.matrix(readRDS("PC.rds"))
plot(PC[sample(nrow(PC), 10e3), ], pch = 20)  # verif: all
plot(PC[sample(sub, 10e3), ], pch = 20)       # verif: caucasian only

# Fit L1-penalized logistic regression
system.time(
  mod <- big_spLogReg(G, y[sub][ind.train], ind.train,
                      covar.train = PC[sub[ind.train], ],
                      ncores = 10)
) # 13 min
summary(mod, best.only = TRUE)$nb_var  # 2653
summary(mod, best.only = TRUE)$validation_loss  # 0.1819475
pred <- predict(mod, G, ind.test, covar.row = PC[sub[ind.test], ])
AUC(pred, y[sub][ind.test])  # 59.8

# Fit logistic regression with elastic net penalization
# (using the same sets as before)
system.time(
  mod2 <- big_spLogReg(G, y[sub][ind.train], ind.train,
                       covar.train = PC[sub[ind.train], ],
                       alphas = c(0.5, 0.1, 0.01, 0.001),
                       ind.sets = attr(mod, "ind.sets"),
                       ncores = nb_cores())
) # 45 min
plot(mod2)
summary(mod2)$validation_loss  # 0.1821156 0.1820076 0.1819554 0.1819486
pred2 <- predict(mod2, G, ind.test, covar.row = PC[sub[ind.test], ])
AUC(pred2, y[sub][ind.test])  # 59.8


#### C+T ####
res <- pkg.paper.PRS::PRS(G, CHR, ukb$map$physical.pos, y[sub], PC[sub, ],
                          ind.train, ind.test)
# 15h for the GWAS part
res$AUC <- sapply(res$pred, AUC, y[sub[ind.test]])
saveRDS(res, "BC-paper2-gwas.rds")
res
# # A tibble: 9 x 6
#   method        pred           thr.r2 set             timing   AUC
#   <chr>         <list>          <dbl> <list>           <dbl> <dbl>
# 1 PRS-all       <dbl [38,628]>   0.05 <int [159,453]> 51236. 0.543
# 2 PRS-stringent <dbl [38,628]>   0.05 <int [10]>      51236. 0.578
# 3 PRS-max       <dbl [38,628]>   0.05 <int [21]>      51236. 0.589
# 4 PRS-all       <dbl [38,628]>   0.2  <int [293,226]> 51236. 0.546
# 5 PRS-stringent <dbl [38,628]>   0.2  <int [12]>      51236. 0.580
# 6 PRS-max       <dbl [38,628]>   0.2  <int [30]>      51236. 0.589
# 7 PRS-all       <dbl [38,628]>   0.8  <int [575,360]> 51236. 0.552
# 8 PRS-stringent <dbl [38,628]>   0.8  <int [22]>      51236. 0.574
# 9 PRS-max       <dbl [38,628]>   0.8  <int [3,151]>   51236. 0.582
