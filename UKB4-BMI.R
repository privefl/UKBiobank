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
}) # 80 sec
stopCluster(cl)

sum(lengths(list_snp_id))  # 656,060

sample <- bigreadr::fread2("ukb25589_imp_chr1_v3_s487327.sample")
str(sample)
(N <- readBin("data/ukb_imp_chr1_v3.bgen", what = 1L, size = 4, n = 4)[4])
sample <- sample[-1, ]
stopifnot(nrow(sample) == N)

csv <- "ukb22544.csv"
library(bigreadr)
df0 <- fread2(csv, select = c("eid", "21001-0.0", "22006-0.0"),
              col.names = c("eid", "BMI", "is_caucasian"))
ind.indiv <- match(df0$eid, sample$ID_2)
# related samples
df0$is_rel2 <- df0$eid %in% fread2("ukb25589_rel_s488346.dat")$ID2

y <- df0$BMI
sub <- which(df0$is_caucasian & !is.na(y) & !is.na(ind.indiv) & !df0$is_rel2)
length(sub)  # 334,540


system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKB_imp_BMI",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = 12
  )
) # 2.5H


set.seed(1)
ind.train <- sort(sample(length(sub), 300e3))
ind.test <- setdiff(seq_along(sub), ind.train)

library(bigsnpr)
ukb <- snp_attach("data/UKB_imp_BMI.rds")
G <- ukb$genotypes
dim(G) # 334,540 x 656,060
file.size(G$backingfile) / 1024^3  # 204 GB

system.time(
  mod <- big_spLinReg(G, log(y[sub][ind.train]), ind.train,
                      dfmax = 200e3, n.abort = 2, ncores = 12)
) # 4h
plot(mod)
summary(mod)
summary(mod)$message

pred <- exp(predict(mod, G, ind.test))
library(ggplot2)
qplot(pred, y[sub[ind.test]], alpha = I(0.2)) +
  theme_bigstatsr() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  labs(x = "Predicted BMI with lasso", y = "True BMI")
cor(pred, y[sub[ind.test]])
