infos_chr1 <- bigreadr::fread2("ukb_imp_mfi/ukb_mfi_chr1_v3.txt")

sumstats_chr1 <- subset(sumstats, chr == 1)

str(infos_chr1)
head(infos_chr1$V3)
hist(infos_chr1$V8)  ## r2_impute?
hist(infos_chr1$V6)  ## MAF?
head(sumstats_chr1$position_b37)


find_chr1 <- match(sumstats_chr1$position_b37, infos_chr1$V3)
mean(!is.na(find_chr1)) # 98.8%

library(magrittr)
library(hexbin)
maf_chr1 <- sumstats_chr1$bcac_icogs2_eaf_controls %>%
  pmin(., 1 - .)
plot(hexbin(maf_chr1, infos_chr1$V6[find_chr1], xbins = 100))



plot(hexbin(sumstats_chr1$bcac_icogs2_r2, infos_chr1$V8[find_chr1], xbins = 100))
cor(sumstats_chr1$bcac_icogs2_r2, infos_chr1$V8[find_chr1],
    use = "pairwise.complete.obs", method = "spearman")

A1 <- sumstats_chr1[c("a0", "a1")]
A2 <- infos_chr1[find_chr1, 4:5]
cbind(A1, A2)

library(bigsnpr)
same_ref(A1$a0, A1$a1, A2$V4, A2$V5)

alleles <- cbind.data.frame(A1$a0, A1$a1, A2$V4, A2$V5)

alleles[1:60, ]
mean(A1$a0 == A2$V4 & A1$a1 == A2$V5, na.rm = TRUE) # >99%
ind <- which(!(A1$a0 == A2$V4 & A1$a1 == A2$V5))
alleles[ind, ][1:140, ]
