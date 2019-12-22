# download.file("https://data.broadinstitute.org/alkesgroup/BOLT-LMM/downloads/BOLT-LMM_v2.3.2.tar.gz",
#               destfile = "BOLT-LMM_v2.3.2.tar.gz")
# untar("BOLT-LMM_v2.3.2.tar.gz")

# Source L23-76 of simu-rel.R
# simu$map <- dplyr::mutate(simu$map, genetic.dist = 0, rsid = NULL,
#                           chromosome = as.integer(chromosome))
# simu$fam <- snp_fake(nrow(G), 1)$fam
# bed <- snp_writeBed(simu, bedfile = "data/UKBB_full.bed")

pcfile <- bigreadr::fwrite2(
  cbind.data.frame(FID = 0, IID = paste0("ind_", seq_along(y)), y, PC),
  sep = " ", file = "pheno-covar.txt")
readLines(pcfile, n = 5)

# system("BOLT-LMM_v2.3.2/bolt -h")
system.time(
  system(paste(
    "BOLT-LMM_v2.3.2/bolt",
    "--bfile data/UKBB_full",
    "--phenoFile=pheno-covar.txt --phenoCol=y",
    "--covarFile=pheno-covar.txt --qCovarCol=PC{1:10}",
    "--lmm --verboseStats",
    "--LDscoresFile=BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz",
    "--statsFile bolt-res3.txt",
    "--numThreads=20"
  ))
)

# Time for estimating mixture parameters = 4259.84 sec
#
# Calibration stats: mean and lambdaGC (over SNPs used in GRM)
# (note that both should be >1 because of polygenicity)
# Mean LINREG: 1.35209 (656060 good SNPs)   lambdaGC: 1.2963
# Mean BOLT_LMM_INF: 1.23721 (656060 good SNPs)   lambdaGC: 1.18338
# Note that LINREG may be confounded by a factor of 1.09739
#
# === Streaming genotypes to compute and write assoc stats at all SNPs ===
#
# Time for streaming genotypes and writing output = 537.981 sec
#
# Total elapsed time for analysis = 7166.07 sec

# M= 10K -> 10.8% variance explained in CV

str(gwas_bolt <- bigreadr::fread2("bolt-res3.txt"))
library(ggplot2)
qplot(-log10(gwas_bolt$P_LINREG), -predict(gwas0)) + xlim(1, NA) +
  theme_bigstatsr() +
  geom_abline(color = "red")  # small differences
qplot(-predict(gwas0_gc), -log10(gwas_bolt$P_BOLT_LMM_INF)) + xlim(1, NA) +
  theme_bigstatsr() +
  labs(x = "-log10(p) (Standard GWAS + correction for relatedness)",
       y = "-log10(p) (BOLT-LMM)") +
  geom_abline(color = "blue", slope = 1 / (1 - 0.027)) +
  geom_abline(color = "red")  # no power gain?

gwas_bolt2 <- structure(
  data.frame(score = gwas_bolt$CHISQ_BOLT_LMM_INF),
  class = c("mhtest", "data.frame"),
  transfo = identity,
  predict = function(xtr) {
    pchisq(xtr, df = 1, lower.tail = FALSE, log.p = TRUE) / log(10)
  }
)
plot(gwas_bolt2)
snp_manhattan(gwas_bolt2, CHR, POS, npoints = 20e3)
snp_qq(gwas_bolt2) + xlim(1, NA)

cowplot::plot_grid(
  snp_qq(gwas_bolt2[ is_odd_chr, , drop = FALSE]) + xlim(1, NA),
  snp_qq(gwas_bolt2[!is_odd_chr, , drop = FALSE]) + xlim(1, NA),
  snp_qq(gwas[ is_odd_chr, ]) + xlim(1, NA),
  snp_qq(gwas[!is_odd_chr, ]) + xlim(1, NA),
  snp_qq(gwas0[ is_odd_chr, ]) + xlim(1, NA),
  snp_qq(gwas0[!is_odd_chr, ]) + xlim(1, NA),
  scale = 0.95, labels = LETTERS[1:6], label_size = 20, ncol = 2
)

plot(bolt ~ lm, data = data.frame(bolt = gwas_bolt2$score, lm = gwas0$score^2))
# 29% increase
