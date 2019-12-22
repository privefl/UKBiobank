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

df0 <- readRDS("pheno.rds")
ind.indiv <- match(df0$eid, sample$ID_2)
y <- df0$height

set.seed(1)
sub <- sample(which(df0$is_caucasian & !is.na(y) & !is.na(ind.indiv)), 20e3)
length(sub)  # 20,000

system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKB_imp_height_20K",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = 12
  )
) # 2H


library(bigsnpr)
ukb <- snp_attach("data/UKB_imp_height_20K.rds")
G <- ukb$genotypes
CHR <- as.numeric(ukb$map$chromosome)
POS <- ukb$map$physical.pos
dim(G) # 374,131 x 656,060
file.size(G$backingfile) / 1024^3  # 12 GB

PC <- as.matrix(readRDS("PC.rds"))
plot(PC[sample(nrow(PC), 10e3), ], pch = 20)  # verif: all
plot(PC[sample(sub, 10e3), ], pch = 20)       # verif: caucasian only

COVAR <- cbind(PC, df0$date, df0$sex)
system.time(
  gwas <- big_univLinReg(G, y[sub],  covar.train = COVAR[sub, ],
                         ncores = nb_cores())
) # 36 sec

library(ggplot2)
snp_manhattan(gwas, CHR, POS) + ggtitle(NULL)

lpval <- -predict(gwas)
system.time(
  ind.keep <- snp_clumping(G, CHR, infos.pos = POS, S = lpval, ncores = nb_cores())
) # 88 sec

snp_manhattan(gwas, CHR, POS) + ggtitle(NULL) +
  aes(alpha = I(ifelse(cols_along(G) %in% ind.keep, 1, 0)))

ind.keep2 <- ind.keep[lpval[ind.keep] > 3]

snp_manhattan(gwas, CHR, POS) + ggtitle(NULL) +
  aes(alpha = I(ifelse(cols_along(G) %in% ind.keep2, 1, 0)))

