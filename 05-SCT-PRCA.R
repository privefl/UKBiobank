download.file("http://practical.icr.ac.uk/blog/wp-content/uploads/uploadedfiles/oncoarray/MetaSummaryData/meta_v3_onco_euro_overall_ChrAll_1_release.zip",
              destfile = "sumstats_PRCA.zip")
unzip("sumstats_PRCA.zip")
library(bigreadr)
sumstats <- fread2("meta_v3_onco_euro_overall_ChrAll_1_release.txt",
                   select = c("Chr", "position", "Allele1", "Allele2", "Freq1",
                              "Effect", "Pvalue", "OncoArray_imputation_r2"))
nrow(sumstats)  # 20,370,946
library(dplyr)
sumstats <- sumstats %>%
  filter(pmin(Freq1, 1 - Freq1) > 1e-3) %>%
  mutate(Freq1 = NULL, Allele1 = toupper(Allele1), Allele2 = toupper(Allele2))
nrow(sumstats)  # 15,975,602

sumstats2 <- sumstats
sumstats2$Allele1 <- sumstats$Allele2
sumstats2$Allele2 <- sumstats$Allele1
sumstats2$Effect <- -sumstats$Effect
sumstats2 <- rbind(sumstats, sumstats2)

# match variants with UKBB
library(doParallel)
NCORES <- 12
registerDoParallel(cl <- makeCluster(NCORES))
info_snp <- foreach(chr = 1:22, .packages = "dplyr") %dopar% {
  paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt") %>%
    bigreadr::fread2(showProgress = FALSE, nThread = 1) %>%
    dplyr::inner_join(
      sumstats2[sumstats2$Chr == chr, ], ., na_matches = "never",
      by = c(position = "V3", Allele1 = "V4", Allele2 = "V5")
    ) %>%
    arrange(position) %>%
    transmute(
      id    = paste(chr, position, Allele1, Allele2, sep = "_"),
      beta  = -Effect,
      lpval = -log10(Pvalue),
      info  = pmin(OncoArray_imputation_r2, V8)
    ) %>%
    na.omit() %>%
    subset(lpval > 1 & info > 0.3)
}
stopCluster(cl)

list_snp_id <-  lapply(info_snp, function(df) df$id)
beta <-  unlist(lapply(info_snp, function(df) df$beta))
lpval <- unlist(lapply(info_snp, function(df) df$lpval))
info <-  unlist(lapply(info_snp, function(df) df$info))
length(beta)  # 1,496,403

# subset samples
sample <- fread2("ukb25589_imp_chr1_v3_s487327.sample")
str(sample)
(N <- readBin("data/ukb_imp_chr1_v3.bgen", what = 1L, size = 4, n = 4)[4])
sample <- sample[-1, ]
stopifnot(nrow(sample) == N)

csv <- "ukb22544.csv"
df0 <- fread2(
  csv,
  select = c("eid", "22001-0.0", "22006-0.0"),
  col.names = c("eid", "sex", "is_caucasian")
)
ind.indiv <- match(df0$eid, sample$ID_2)
# related samples
df0$is_rel2 <- df0$eid %in% fread2("ukb25589_rel_s488346.dat")$ID2

df_cancer0 <- fread2(csv, select = paste0("2453-", 0:2, ".0"))
df_cancer1 <- fread2(csv, select = c(paste0("20001-0.", 0:5),
                                     paste0("20001-1.", 0:5),
                                     paste0("20001-2.", 0:5)))
df_cancer2 <- fread2(csv, select = c(paste0("40006-", 0:31, ".0"),
                                     paste0("40001-", 0:2, ".0"),
                                     paste0("40002-0.", 0:13),
                                     paste0("40002-1.", 0:13),
                                     paste0("40002-2.", 0:13)))
ind_PRCA <- sort(unique(unlist(c(
  lapply(df_cancer1,  function(x) which(x == 1044)),
  lapply(df_cancer2, function(x) which(x %in% c("C61", "D075")))
))))
table(df0$sex[ind_PRCA])
# 0    1
# 2 8847

y <- rep(NA, nrow(df0))
y[rowMeans(df_cancer0 == 0, na.rm = TRUE) == 1] <- 0
y[ind_PRCA] <- 1

sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian &
               df0$sex == 1 & !is.na(y))
length(sub)  # 147,921
table(y.sub <- y[sub])
#      0      1
# 141591   6330

system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_PRCA",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 3.6H


library(bigsnpr)
ukbb <- snp_attach("data/UKBB_PRCA.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 206 GB
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos

set.seed(1)
ind.train <- sort(sample(length(sub), 120e3))
ind.test <- setdiff(seq_along(sub), ind.train)

system.time(
  all_keep <- snp_grid_clumping(
    G, CHR, POS, lpS = lpval, ind.row = ind.train, infos.imp = info,
    grid.thr.imp = c(0.3, 0.6, 0.9, 0.95),
    grid.thr.r2 = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95),
    grid.base.size = c(50, 100, 200, 500),
    ncores = NCORES
  )
) # 7.3H -> swapping -> should have use less cores
plot(lengths(all_keep[[1]]), pch = 20)

system.time(
  multi_PRS <- snp_grid_PRS(
    G, all_keep, betas = beta, lpS = lpval, ind.row = ind.train,
    grid.lpS.thr = seq_log(1, 0.999 * max(lpval), 50),
    backingfile = "data/UKBB_PRCA_scores", ncores = NCORES
  )
) # 1.9H

nPC <- 20
PC <- fread2(csv, select = paste0("22009-0.", 1:nPC),
             col.names = paste0("PC", 1:nPC))
PC_sub <- as.matrix(PC[sub, ])

system.time(
  final_mod <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train], covar.train = PC_sub[ind.train, ], ncores = NCORES
  )
) # 2.4H
mod <- final_mod$mod
plot(mod)
summary(mod)
# # A tibble: 5 x 6
#      alpha validation_loss intercept beta           nb_var message
#      <dbl>           <dbl>     <dbl> <list>          <int> <list>
#   1 0.0001           0.166    -0.627 <dbl [73,436]>  21181 <chr [10]>
#   2 0.001            0.166    -0.484 <dbl [73,436]>   8058 <chr [10]>
#   3 0.01             0.166    -0.397 <dbl [73,436]>   4515 <chr [10]>
#   4 0.1              0.166    -0.352 <dbl [73,436]>   3553 <chr [10]>
#   5 1                0.166    -0.349 <dbl [73,436]>   3326 <chr [10]>

new_beta <- final_mod$beta.G

length(ind <- which(new_beta != 0))  # 461,733
summary(new_beta)
summary(new_beta[which(sign(new_beta * beta) < 0)])

pred <- final_mod$intercept +
  big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind) +
  PC_sub[ind.test, ] %*% final_mod$beta.covar

AUCBoot(pred, y.sub[ind.test])  # 67.2 [65.7-68.7]

aucs <- apply(multi_PRS, 2, AUC, y.sub[ind.train])

library(tidyverse)

ind2 <- sort(sample(ind, 10e3))
ggplot(data.frame(y = new_beta, x = beta)[ind, ]) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  geom_abline(slope = 0, intercept = 0, color = "blue") +
  geom_point(aes(x, y), size = 0.6) +
  theme_bigstatsr() +
  labs(x = "Effect sizes from GWAS", y = "Non-zero effect sizes from SCT")

grid2 <- attr(all_keep, "grid") %>%
  mutate(thr.lp = list(attr(multi_PRS, "grid.lpS.thr")), num = row_number()) %>%
  unnest() %>%
  mutate(auc = aucs)

std_prs <- grid2 %>%
  filter(thr.imp == 0.3, thr.r2 == 0.2, size == 500, thr.lp <= 8) %>%
  arrange(desc(auc)) %>%
  slice(1) %>%
  print()
#   size thr.r2 thr.imp num   thr.lp       auc
# 1  500    0.2     0.3  14 4.456008 0.6639146

ind.keep <- unlist(map(all_keep, std_prs$num))
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.test]
) # 64.9 [63.2-66.5]
sum(lpval[ind.keep] > std_prs$thr.lp)  # 1683

max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
#     size thr.r2 thr.imp num   thr.lp       auc
# 1   2000   0.05     0.9  62 4.456008 0.6786133
# 2   1000   0.05     0.9  61 4.456008 0.6784085
# 3   4000   0.05     0.9  63 4.456008 0.6783587
# 4  10000   0.05     0.9  64 4.456008 0.6783587
# 5   5000   0.01     0.9  57 5.516359 0.6775051
# 6  10000   0.01     0.9  58 5.516359 0.6775051
# 7  20000   0.01     0.9  59 5.516359 0.6775051
# 8  50000   0.01     0.9  60 5.516359 0.6775051
# 9   5000   0.01     0.9  57 4.456008 0.6774659
# 10 10000   0.01     0.9  58 4.456008 0.6774067

ind.keep <- unlist(map(all_keep, max_prs$num))
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp),
  y.sub[ind.test]
) # 65.2 [63.6-66.8]
sum(lpval[ind.keep] > max_prs$thr.lp)  # 610

ggplot(grid2) +
  geom_point(aes(thr.lp, auc)) +
  facet_grid(thr.imp ~ thr.r2 + size) +
  scale_x_log10(limits = c(1, 10), breaks = c(1, 2, 5), minor_breaks = 1:10) +
  ylim(0.62, NA) +
  theme_bigstatsr(size.rel = 0.7) +
  labs(x = "-log10(p-value) threshold (log scale)", y = "AUC")


ind3 <- which(abs(beta) > 10)
gwas <- big_univLogReg(G, y.sub[ind.train], ind.train, ind3)
cbind(beta_GWAS = beta[ind3], pval_GWAS = 10^(-lpval[ind3]),
      af = big_scale()(G, ind.row = ind.train, ind.col = ind3)$center / 2,
      gwas, pval = predict(gwas, log10 = FALSE))
#    beta_GWAS pval_GWAS           af         estim      std.err niter         score       pval
# 1   -18.3006 0.0409000 7.916667e-06 -1.114748e+01 2.632595e+01     7 -4.234405e-01 0.67197388
# 2   -16.6150 0.0259600 4.708333e-06  1.780807e+01 8.954184e+00     6  1.988799e+00 0.04672339
# 3   -16.6150 0.0259600 4.708333e-06  1.780807e+01 8.954184e+00     6  1.988799e+00 0.04672339
# 4   -18.4922 0.0610000 6.275000e-05 -3.584283e-01 2.581716e+00     4 -1.388333e-01 0.88958185
# 5   -20.0049 0.0976600 2.916667e-07 -7.663363e+02 6.604582e+03    NA -1.160310e-01 0.90762797
# 6   -13.3498 0.0035520 5.916667e-06 -1.146015e+01 3.175944e+01     7 -3.608424e-01 0.71821730
# 7   -20.1500 0.0611500 1.683333e-05 -3.052031e+00 9.004343e+00     7 -3.389510e-01 0.73464666
# 8   -20.1371 0.0513400 3.208333e-06  2.799164e+01 2.220995e+01     6  1.260320e+00 0.20755386
# 9   -20.2094 0.0149800 2.295833e-05 -4.502343e+01 7.346198e+01    10 -6.128806e-01 0.53995529
# 10  -11.7501 0.0030230 8.750000e-07 -8.251973e+02 6.265142e+03    NA -1.317125e-01 0.89521172
# 11  -20.0370 0.0945400 2.720833e-05 -3.569095e-01 2.628259e+00     5 -1.357969e-01 0.89198182
# 12  -20.0798 0.0694400 2.750000e-06 -9.641282e+02 8.145099e+03    NA -1.183691e-01 0.90577519
# 13  -20.0855 0.0483300 5.250000e-06  6.842568e+00 1.938198e+01     5  3.530377e-01 0.72406018
# 14  -16.4113 0.0652200 4.170833e-05  8.911422e-02 2.083875e+00     4  4.276370e-02 0.96588990
# 15  -20.0437 0.0001334 8.333333e-07 -8.811054e+02 6.767293e+03    NA -1.302006e-01 0.89640775
# 16  -20.2844 0.0007550 8.070833e-04 -4.827150e-02 4.718294e-01     4 -1.023071e-01 0.91851292
# 17  -20.0834 0.0807000 2.533333e-05 -9.495788e+00 1.726726e+01     8 -5.499303e-01 0.58236718
# 18  -20.1468 0.0457300 1.190000e-04 -3.474032e+00 3.488454e+00     6 -9.958656e-01 0.31931543
# 19  -20.0996 0.0146000 1.208333e-05  3.208018e+03 6.523114e+07     6  4.917925e-05 0.99996076
# 20  -20.0677 0.0180800 2.500000e-07 -6.903006e+02 6.386592e+03    NA -1.080859e-01 0.91392754
# 21  -20.0265 0.0979400 5.970833e-05  9.267366e-01 9.213898e-01     6  1.005803e+00 0.31451033
# 22  -20.1598 0.0310700 1.583333e-06 -5.827571e+00 3.645173e+01     6 -1.598709e-01 0.87298277
# 23  -12.3406 0.0307600 7.495833e-05  2.286299e+00 9.302098e-01    NA  2.457831e+00 0.01397789
# 24  -17.5909 0.0295900 3.208333e-06 -3.261151e+02 3.567589e+03    NA -9.141050e-02 0.92716642
# 25  -20.2353 0.0051480 1.833333e-05 -2.320849e+01 2.127291e+01     6 -1.090988e+00 0.27527824
# 26  -16.6329 0.0451600 1.733333e-05  2.347977e+00 2.092458e+00     6  1.122114e+00 0.26181386
# 27  -20.0541 0.0908900 2.750000e-06 -8.223249e+02 6.603018e+03    NA -1.245377e-01 0.90088952
# 28  -12.3663 0.0686900 5.750000e-04  1.519177e-01 4.759680e-01     4  3.191763e-01 0.74959284
# 29  -20.0854 0.0818800 1.422083e-04 -2.331458e+00 2.619950e+00     6 -8.898865e-01 0.37352683
# 30  -12.8060 0.0749900 2.008333e-05 -1.634974e+01 1.842828e+01     6 -8.872092e-01 0.37496630
# 31  -10.2548 0.0494400 1.758333e-05 -2.101394e+00 9.308812e+00     6 -2.257424e-01 0.82140175
# 32  -20.0600 0.0852700 2.791667e-06 -1.676997e+00 1.631476e+01     5 -1.027901e-01 0.91812955
# 33  -19.2376 0.0138300 9.762500e-05 -1.480131e+01 1.214701e+01     8 -1.218515e+00 0.22302829
# 34  -20.1045 0.0159400 1.625000e-06  2.527350e+01 3.764214e+01     6  6.714150e-01 0.50195620
# 35  -20.0780 0.0821600 3.541667e-06 -3.422472e+01 6.884415e+01     7 -4.971334e-01 0.61909501
# 36  -20.1866 0.0544900 4.125000e-05 -1.084944e+01 1.538298e+01     8 -7.052887e-01 0.48063057
# 37  -10.8455 0.0507300 1.000833e-04 -2.177248e+00 3.139957e+00     6 -6.934006e-01 0.48805819
# 38  -20.1452 0.0717200 2.333333e-06  4.658841e+01 2.475281e+01     5  1.882146e+00 0.05981616
# 39  -19.0779 0.0976600 1.416667e-06 -9.707272e+02 8.418890e+03    NA -1.153035e-01 0.90820459
# 40  -20.1832 0.0048290 3.506250e-04 -8.044312e-02 5.835241e-01     4 -1.378574e-01 0.89035312
# 41  -19.4530 0.0539700 9.767083e-04 -6.729510e-01 4.956811e-01     5 -1.357629e+00 0.17458141
# 42  -20.0847 0.0838200 1.250000e-07 -1.671252e+03 2.009160e+06    19 -8.318162e-04 0.99933631
# 43  -20.0598 0.0731300 1.379167e-05  5.206318e-01 7.991851e+00     4  6.514533e-02 0.94805829
# 44  -20.1889 0.0660600 3.112500e-05 -2.776870e-01 3.138487e+00     5 -8.847798e-02 0.92949678
# 45  -20.1271 0.0595600 9.933333e-05 -1.787639e+01 3.258004e+01    11 -5.486913e-01 0.58321731
# 46  -20.1356 0.0236000 1.434125e-03  5.756712e-01 3.035800e-01     5  1.896275e+00 0.05792368
# 47  -13.7656 0.0814700 3.041667e-06  2.231468e+01 1.940135e+01     6  1.150161e+00 0.25007754
# 48  -14.2224 0.0495300 5.500000e-06 -6.131212e+01 7.376851e+01     8 -8.311421e-01 0.40589335
# 49  -19.0200 0.0487600 1.875000e-06 -1.224031e+01 4.759521e+01     6 -2.571752e-01 0.79704353
# 50  -19.1616 0.0496600 1.875000e-06 -1.224031e+01 4.759521e+01     6 -2.571752e-01 0.79704353
# 51  -17.6729 0.0253300 5.625000e-06 -7.680014e+02 5.944457e+03    NA -1.291962e-01 0.89720240
# 52  -20.1866 0.0595900 1.870833e-05 -2.838930e+01 6.847626e+01    10 -4.145861e-01 0.67844493
# 53  -20.1836 0.0082580 2.666667e-06 -9.336175e+00 4.588885e+01     7 -2.034519e-01 0.83878182
# 54  -20.1303 0.0592500 1.208333e-06 -8.475698e+02 7.493811e+03    NA -1.131026e-01 0.90994918
# 55  -20.2112 0.0309700 1.208333e-06  7.075493e+01 5.191219e+01     6  1.362973e+00 0.17289092
# 56  -20.0926 0.0231400 1.625000e-05 -6.636759e-01 4.750666e+00     5 -1.397017e-01 0.88889571
# 57  -20.0938 0.0227500 2.217917e-04 -2.598800e-02 7.505464e-01     4 -3.462544e-02 0.97237842
# 58  -20.1367 0.0219800 3.804167e-05  7.373114e-01 1.506695e+00     6  4.893569e-01 0.62458904
# 59  -20.0351 0.0703300 5.000000e-05  8.918601e-01 1.263982e+00     6  7.055958e-01 0.48043954
# 60  -20.1770 0.0216200 2.916667e-06  7.537592e+00 2.524572e+01     5  2.985691e-01 0.76526885
# 61  -20.1251 0.0320200 3.458333e-06  7.283398e+00 1.441364e+01     6  5.053129e-01 0.61333910
# 62  -20.0800 0.0414500 4.791667e-06  5.090344e+00 1.379978e+01     6  3.688714e-01 0.71222359
# 63  -20.1029 0.0765300 8.879167e-05 -4.192614e+00 4.570170e+00     7 -9.173868e-01 0.35893999
# 64  -10.8130 0.0525800 5.583333e-06 -9.186973e+02 7.560468e+03    NA -1.215133e-01 0.90328449
# 65  -10.8237 0.0443700 2.916667e-07 -6.238713e+02 6.380383e+03    NA -9.777960e-02 0.92210731
# 66  -20.0628 0.0312600 4.166667e-08 -1.606971e+03 1.456908e+06    17 -1.103001e-03 0.99911993
# 67  -12.3087 0.0338400 9.000000e-06 -6.166627e+00 2.061072e+01     8 -2.991951e-01 0.76479118
# 68  -11.6432 0.0081540 4.333333e-06 -2.587701e+01 4.766300e+01     5 -5.429160e-01 0.58718763
# 69  -15.5527 0.0225400 1.837500e-05 -9.581848e+02 5.824406e+03    NA -1.645120e-01 0.86932809
# 70  -10.5940 0.0129700 1.208333e-06 -1.683480e+02 2.074980e+03    NA -8.113231e-02 0.93533673
# 71  -16.4361 0.0438000 1.025000e-05  1.871634e+00 6.796387e+00     5  2.753867e-01 0.78301919
# 72  -20.0240 0.0859900 2.166667e-06 -3.051091e+01 7.242026e+01     6 -4.213035e-01 0.67353351
# 73  -20.0540 0.0809800 3.375000e-06 -1.283294e+01 3.719543e+01     6 -3.450139e-01 0.73008394
# 74  -20.1307 0.0706300 2.708333e-06 -8.568958e+02 7.342301e+03    NA -1.167067e-01 0.90709247
# 75  -10.4862 0.0069630 2.858333e-05  1.686414e+00 1.115964e+00     6  1.511172e+00 0.13074451
# 76  -16.8825 0.0466500 7.248333e-04  9.592938e-02 4.956639e-01     4  1.935371e-01 0.84653832
# 77  -20.1123 0.0008016 2.875000e-06 -9.443819e+02 7.379998e+03    NA -1.279651e-01 0.89817662
# 78  -20.0921 0.0540900 7.833333e-06 -2.405113e+01 3.723064e+01     5 -6.460038e-01 0.51827692
# 79  -20.1045 0.0100700 3.916667e-05  2.289969e+00 1.026168e+00     8  2.231574e+00 0.02564314
# 80  -20.1897 0.0773000 4.166667e-07 -7.249988e+02 7.589044e+03    NA -9.553229e-02 0.92389204
# 81  -20.0606 0.0733600 8.137500e-05  1.529175e+00 1.067063e+00     6  1.433068e+00 0.15183830
# 82  -10.4140 0.0978800 2.029167e-05 -1.711239e+00 5.452795e+00     6 -3.138279e-01 0.75365173
# 83  -16.3895 0.0291500 2.166667e-05 -1.020024e+03 7.824641e+03    NA -1.303604e-01 0.89628128
# 84  -20.1105 0.0734400 4.750000e-06 -8.587098e+02 8.550134e+03    NA -1.004323e-01 0.92000111
# 85  -20.1096 0.0353600 3.041667e-06 -8.680055e+02 7.314546e+03    NA -1.186684e-01 0.90553806
# 86  -20.0622 0.0705900 6.625000e-06 -9.134078e+02 6.887117e+03    NA -1.326256e-01 0.89448952
