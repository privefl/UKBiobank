sumstats <- bigreadr::fread2("breast-cancer/oncoarray_bcac_public_release_oct17.txt",
                             na.strings = "NULL")
dim(sumstats)  ## 11,792,542 x 66
format(object.size(sumstats), unit = "GB")  # 7.3 Gb

pval_all <- sumstats$bcac_onco_icogs_gwas_P1df
hist(pval_all)

sub <- which(pval_all < 0.0001)
plot(-log10(pval_all[sub]), pch = 20, col = sumstats$chr[sub])

pval_onco2 <- sumstats$bcac_onco2_P1df_Wald
plot(-log10(pval_onco2[sub]), -log10(pval_all[sub]), pch = 20)   ## 10^-150
r2_onco2 <- sumstats$bcac_onco2_r2
hist(r2_onco2)  ## des < 0 OMG

pval_icogs2 <- sumstats$bcac_icogs2_P1df_Wald
plot(-log10(pval_icogs2[sub]), -log10(pval_all[sub]), pch = 20)  ## 10^-120
r2_icogs2 <- sumstats$bcac_icogs2_r2
hist(r2_icogs2)
af_icogs2 <- sumstats$bcac_icogs2_eaf_controls
maf_icogs2 <- pmin(af_icogs2, 1 - af_icogs2)
hist(log10(maf_icogs2))
beta_icogs2 <- sumstats$bcac_icogs2_beta
summary(abs(beta_icogs2))

library(hexbin)
plot(hexbin(r2_icogs2, maf_icogs2, xbins = 20))

length(which(pval_icogs2 < 0.1)) # 1.5M
length(which(pval_icogs2 < 0.1 & r2_icogs2 > 0.3)) # 1.5M
length(which(pval_icogs2 < 0.1 & r2_icogs2 > 0.5)) # 1.4M
length(which(pval_icogs2 < 0.1 & r2_icogs2 > 0.8)) # 840K
length(which(pval_icogs2 < 0.1 & r2_icogs2 > 0.9)) # 590K

sub2 <- which(maf_icogs2 < 0.002)
plot(maf_icogs2[sub2], beta_icogs2[sub2], pch = 20)
sub3 <- which(abs(beta_icogs2) > 1)
maf_icogs2[sub3]
pval_icogs2[sub3]
