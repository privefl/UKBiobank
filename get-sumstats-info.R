setwd("tmp-data")

library(bigreadr)

not_duplicated <- function(sumstats) {
  sumstats[!vctrs::vec_duplicate_detect(sumstats[c("chr", "pos")]), ]
}

#### BRCA ####
# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/MichailidouK_29059683_GCST004988/oncoarray_bcac_public_release_oct17%20(1).txt.gz",
#               destfile = "sumstats_BRCA.txt.gz")
# R.utils::gunzip("sumstats_BRCA.txt.gz")
sumstats_brca <- fread2(
  "sumstats_BRCA.txt", na.strings = "NULL",
  select = c("chr", "position_b37", "bcac_onco2_r2", "bcac_icogs2_r2"),
  col.names = c("chr", "pos", "info_onco_brca", "info_icogs_brca"))
str(sumstats_brca)

merged <- not_duplicated(sumstats_brca)

#### PRCA ####
# download.file("http://practical.icr.ac.uk/blog/wp-content/uploads/uploadedfiles/oncoarray/MetaSummaryData/meta_v3_onco_euro_overall_ChrAll_1_release.zip",
#               destfile = "sumstats_PRCA.zip")
# unzip("sumstats_PRCA.zip")
sumstats_prca <- fread2(
  "meta_v3_onco_euro_overall_ChrAll_1_release.txt",
  select = c("Chr", "position", "OncoArray_imputation_r2"),
  col.names = c("chr", "pos", "info_onco_prca")
)
str(sumstats_prca)  # 20,370,946

merged <- dplyr::full_join(merged, not_duplicated(sumstats_prca))

#### MDD ####
# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/WrayNR_29700475_GCST005839/MDD2018_ex23andMe.gz",
#               destfile = "sumstats_MDD.txt.gz")
# R.utils::gunzip("sumstats_MDD.txt.gz")
# sumstats_mdd <- readr::read_tsv("sumstats_MDD.txt") # Warning: 122 parsing failures.
# sumstats_mdd <- dplyr::select(sumstats_mdd,
#                               c(chr = "CHR", pos = "BP", info_mdd = "INFO"))
# sumstats_mdd  # 13,554,550
#
# merged <- dplyr::full_join(merged, not_duplicated(sumstats_mdd))

#### CAD ####
# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/NikpayM_26343387_GCST003116/cad.add.160614.website.txt",
#               destfile = "sumstats_CAD.txt")
sumstats_cad <- fread2("sumstats_CAD.txt",
                   select = c("chr", "bp_hg19", "median_info"),
                   col.names = c("chr", "pos", "info_cad"))

merged <- dplyr::full_join(merged, not_duplicated(sumstats_cad))

#### Intelligence ####
# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/SavageJE_29942086_GCST006250/harmonised/29942086-GCST006250-EFO_0004337-Build37.f.tsv.gz",
#               destfile = "sumstats_intel.tsv.gz")
# R.utils::gunzip("sumstats_intel.tsv.gz")
# sumstats_intel <- fread2("sumstats_intel.tsv",
#                        select = c("chromosome", "base_pair_location", "mininfo"),
#                        col.names = c("chr", "pos", "info_intel"))
# sumstats_intel$chr <- as.integer(sumstats_intel$chr)
#
# merged <- dplyr::full_join(merged, not_duplicated(sumstats_intel))

#### T1D ####
# https://datadryad.org//resource/doi:10.5061/dryad.ns8q3
# R.utils::gunzip("65009084.tar.gz"); untar("65009084.tar")
sumstats_t1d <- fread2(paste0("meta_chr_", 1:22),
                       select = c("chromosome", "position", "info_score.I", "info_score.A"),
                       col.names = c("chr", "pos", "info_illu_t1d", "info_affy_t1d"))

merged <- dplyr::full_join(merged, not_duplicated(sumstats_t1d))

#### ADHD ####
# https://ipsych.dk/en/research/downloads/data-download-agreement-adhd-european-ancestry-gwas-june-2017/
# R.utils::gunzip("adhd_eur_jun2017.gz")
sumstats_adhd <- fread2("adhd_eur_jun2017",
                         select = c("CHR", "BP", "INFO"),
                         col.names = c("chr", "pos", "info_adhd"))

merged <- dplyr::full_join(merged, not_duplicated(sumstats_adhd))

#### Physical activity UKBB ####
# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/KlimentidisYC_29899525_GCST006097/Klimentidis_29899525_MVPA_Model1_BOLTLMM_500K.txt.gz",
#               destfile = "sumstats_UKBB.txt.gz")
# R.utils::gunzip("sumstats_UKBB.txt.gz")
sumstats_ukbb <- fread2("sumstats_UKBB.txt",
                        select = c("CHR", "BP", "INFO", "SNP"),
                        col.names = c("chr", "pos", "info_ukbb", "snp_id"))

merged <- dplyr::full_join(merged, not_duplicated(sumstats_ukbb))
mean(substr(merged$snp_id, 1, 2) == "rs", na.rm = TRUE)  # 1

#### ####
# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/FadistaJ_29379196_GCST005289/Hirschsprung_GWAS_summaryStats.txt.gz",
#               destfile = "sumstats_dk.txt.gz")
sumstats_dk <- fread2("sumstats_dk.txt.gz",
                      select = c("CHR", "BP"),
                      col.names = c("chr", "pos"))
sumstats_dk$info_dk <- 0.9

merged <- dplyr::full_join(merged, not_duplicated(sumstats_dk))

round(100 * cor(merged[3:10], use = "pairwise.complete.obs"), 1)
#                 info_onco_brca info_icogs_brca info_onco_prca info_cad
# info_onco_brca           100.0            62.3           99.0     75.3
# info_icogs_brca           62.3           100.0           62.3     46.1
# info_onco_prca            99.0            62.3          100.0     74.4
# info_cad                  75.3            46.1           74.4    100.0
# info_illu_t1d             82.7            52.4           82.2     83.7
# info_affy_t1d             74.2            52.0           73.9     74.3
# info_adhd                 74.8            33.4           74.4     70.5
# info_ukbb                 65.1            42.1           72.6     59.3
#                 info_illu_t1d info_affy_t1d info_adhd info_ukbb
# info_onco_brca           82.7          74.2      74.8      65.1
# info_icogs_brca          52.4          52.0      33.4      42.1
# info_onco_prca           82.2          73.9      74.4      72.6
# info_cad                 83.7          74.3      70.5      59.3
# info_illu_t1d           100.0          81.7      74.9      55.7
# info_affy_t1d            81.7         100.0      65.0      47.3
# info_adhd                74.9          65.0     100.0      40.7
# info_ukbb                55.7          47.3      40.7     100.0
colMeans(merged[3:10], na.rm = TRUE)
# info_onco_brca info_icogs_brca  info_onco_prca        info_cad   info_illu_t1d
#      0.9054941       0.7671034       0.7427423       0.9136934       0.8942405
#  info_affy_t1d       info_adhd       info_ukbb
#      0.8252165       0.9684464       0.9407504

merged_nona <- na.omit(merged)
merged_nona$snp_id
str(merged_nona)  # 4.3M
hist(merged_nona$info_illu_t1d)
merged_nona$info_min_all <- apply(merged_nona[3:10], 1, min)
merged_nona$info_mean_all <- rowMeans(merged_nona[3:10])
hist(merged_nona$info_min_all)

sapply(setNames(nm = seq(0.5, 1, by = 0.01)), function(thr) {
  sum(merged_nona$info_min_all >= thr)
})
#     0.5    0.51    0.52    0.53    0.54    0.55    0.56    0.57    0.58    0.59     0.6
# 4130614 4105505 4079272 4051637 4021407 3989708 3955934 3918975 3881674 3842299 3799976
#    0.61    0.62    0.63    0.64    0.65    0.66    0.67    0.68    0.69     0.7    0.71
# 3755214 3708109 3660611 3610927 3557704 3504224 3448578 3390500 3329844 3267896 3203855
#    0.72    0.73    0.74    0.75    0.76    0.77    0.78    0.79     0.8    0.81    0.82
# 3138485 3071489 3002852 2932242 2860682 2786816 2709615 2631038 2550091 2465408 2376319
#    0.83    0.84    0.85    0.86    0.87    0.88    0.89     0.9    0.91    0.92    0.93
# 2287060 2198145 2105843 2013789 1913469 1808736 1697296 1579028 1454460 1320057 1174549
#    0.94    0.95    0.96    0.97    0.98    0.99       1
# 1016258  840767  652274  448905  234998   44129       0


map_hapmap3 <- bigreadr::fread2("ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3/plink_format/draft_2/hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2")
merged_nona_hapmap3 <-
  dplyr::left_join(map_hapmap3[2], merged_nona, by = c(V2 = "snp_id"))
na.omit(merged_nona_hapmap3)
merged_hapmap3 <- dplyr::left_join(map_hapmap3[2], merged, by = c(V2 = "snp_id"))
str(na.omit(merged_hapmap3))  # 739,701
str(merged_hapmap3)

hist(merged_nona_hapmap3$info_min_all, freq = FALSE,
     breaks = seq(0.2, 1, by = 0.02),
     col = scales::alpha("blue", 0.4))
hist(pmax(0, merged_nona$info_min_all), freq = FALSE, add = TRUE,
     breaks = seq(0, 1, by = 0.02),
     col = scales::alpha("red", 0.4))

logit <- function(x) {
  x <- pmin(pmax(1e-5, x), 1 - 1e-5)
  log(x / (1 - x))
}

hist(logit(merged$info_ukbb))
hist(logit(merged$info_illu_t1d))

cor(do.call("cbind", lapply(merged[3:10], logit)), use = "pairwise.complete.obs")
