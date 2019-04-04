# urls <- gsubfn::strapply(
#   readLines("https://datadryad.org//resource/doi:10.5061/dryad.ns8q3"),
#   "<a href=\"(/bitstream/handle/10255/dryad\\.[0-9]+/meta_chr_[0-9]+\\?sequence=1)\">",
#   simplify = 'c')
#
# sumstats <- purrr::map_dfr(urls, ~ {
#   download.file(paste0("https://datadryad.org", .x),
#                 destfile = (tmp <- tempfile(fileext = ".txt")))
#   sumstats <- bigreadr::fread2(
#     tmp, select = c("chromosome", "position", "a0", "a1", "beta.meta", "p.meta"))
#   na.omit(sumstats)
# })
#
# saveRDS(sumstats, "sumstats_T1D.rds")

sumstats <- readRDS("sumstats_T1D.rds")
nrow(sumstats)  # 8,996,866
hist(sumstats$p.meta)
hist(sumstats$beta.meta)
rle(sumstats$chromosome)

# augment dataset to match reverse alleles
sumstats2 <- sumstats
sumstats2$a0 <- sumstats$a1
sumstats2$a1 <- sumstats$a0
sumstats2$beta.meta <- -sumstats$beta.meta
sumstats2 <- rbind(sumstats, sumstats2)

# match variants with UKBB
library(doParallel)
NCORES <- 12
registerDoParallel(cl <- makeCluster(NCORES))
info_snp <- foreach(chr = 1:22, .packages = "dplyr") %dopar% {

  paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt") %>%
    bigreadr::fread2(showProgress = FALSE, nThread = 1) %>%
    dplyr::inner_join(
      sumstats2[sumstats2$chromosome == chr, ], ., na_matches = "never",
      by = c(position = "V3", a0 = "V4", a1 = "V5")
    ) %>%
    arrange(position) %>%
    transmute(
      id = paste(chromosome, position, a0, a1, sep = "_"),
      beta = beta.meta,
      lpval = pmin(-log10(p.meta), -log10(.Machine$double.xmin) + abs(beta.meta)),
      info = V8
    ) %>%
    na.omit() %>%
    filter(lpval > 1 & info > 0.3)
}
stopCluster(cl)

list_snp_id <-  lapply(info_snp, function(df) df$id)
beta <-  unlist(lapply(info_snp, function(df) df$beta))
lpval <- unlist(lapply(info_snp, function(df) df$lpval))
info <-  unlist(lapply(info_snp, function(df) df$info))
length(beta)  # 1,222,701

# subset samples
library(bigreadr)
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

# Non-cancer illness codes (https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=6)
df_illness <- fread2(csv, select = c(paste0("20002-0.", 0:28),
                                     paste0("20002-1.", 0:28),
                                     paste0("20002-2.", 0:28)))
df_ICD10 <- fread2(csv, select = c(paste0("41202-0.", 0:379),
                                   paste0("41204-0.", 0:434)))
ind_diabetes <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x %in% 1220:1223)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) %in% paste0("E", 10:14)))
))))
ind_TD1 <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1222)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "E10"))
))))
ind_TD2 <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1223)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "E11"))
))))
y <- rep(0, nrow(df_illness))
y[ind_diabetes] <- NA
y[ind_TD1] <- 1
y[ind_TD2] <- NA

sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian & !is.na(y) &
               !is.na(df0$sex))
length(sub)  # 315,338
table(y.sub <- y[sub])
#      0      1
# 314570    768

system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_T1D",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 3.5H

library(bigsnpr)
ukbb <- snp_attach("data/UKBB_T1D.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 359 GB
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos

set.seed(1)
ind.train <- sort(sample(length(sub), 250e3))
ind.test <- setdiff(seq_along(sub), ind.train)

system.time(
  all_keep <- snp_grid_clumping(
    G, CHR, POS, lpS = lpval, ind.row = ind.train, infos.imp = info,
    grid.thr.imp = c(0.3, 0.6, 0.9, 0.95),
    grid.thr.r2 = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95),
    grid.base.size = c(50, 100, 200, 500),
    ncores = NCORES
  )
) # 9H -> swapping -> need to use less cores
plot(lengths(all_keep[[1]]), pch = 20)

system.time(
  multi_PRS <- snp_grid_PRS(
    G, all_keep, betas = beta, lpS = lpval, ind.row = ind.train,
    grid.lpS.thr = seq_log(1, 0.999 * max(lpval), 50),
    backingfile = "data/UKBB_T1D_scores", ncores = NCORES
  )
) # 3H

nPC <- 20
PC <- fread2(csv, select = paste0("22009-0.", 1:nPC),
             col.names = paste0("PC", 1:nPC))
COVAR <- as.matrix(cbind(df0$sex, PC)[sub, ])

system.time(
  final_mod <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train], covar.train = COVAR[ind.train, ], ncores = NCORES
  )
) # 2.3H
mod <- final_mod$mod
plot(mod)
summary(mod)
# # A tibble: 5 x 6
#      alpha validation_loss intercept beta           nb_var message
#      <dbl>           <dbl>     <dbl> <list>          <int> <list>
#   1 0.0001          0.0153     -6.69 <dbl [55,349]>   5882 <chr [10]>
#   2 0.001           0.0153     -6.73 <dbl [55,349]>   2359 <chr [10]>
#   3 0.01            0.0153     -6.79 <dbl [55,349]>   2118 <chr [10]>
#   4 0.1             0.0153     -6.76 <dbl [55,349]>   1603 <chr [10]>
#   5 1               0.0153     -6.75 <dbl [55,349]>   1646 <chr [10]>

new_beta <- final_mod$beta.G

length(ind <- which(new_beta != 0))  # 160,756
summary(new_beta)
summary(new_beta[which(sign(new_beta * beta) < 0)])

pred <- final_mod$intercept +
  big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind) +
  COVAR[ind.test, ] %*% final_mod$beta.covar

AUCBoot(pred, y.sub[ind.test])  # 79.4 [75.5-83.2]


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
# 1  500    0.2     0.3  14 2.269113 0.7542697

ind.keep <- unlist(map(all_keep, std_prs$num))
# Verif on training set
AUC(
  snp_PRS(G, beta[ind.keep], ind.test = ind.train, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.train]
)
# Eval on test set
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.test]
) # 73.7 [69.4-77.8]
sum(lpval[ind.keep] > std_prs$thr.lp)  # 13,169

max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
#     size thr.r2 thr.imp num   thr.lp       auc
# 1   5000   0.01     0.9  57 5.148873 0.7740869
# 2  10000   0.01     0.9  58 5.148873 0.7733324
# 3  20000   0.01     0.9  59 5.148873 0.7733324
# 4  50000   0.01     0.9  60 5.148873 0.7733324
# 5   5000   0.01     0.9  57 6.507068 0.7730739
# 6   5000   0.01     0.9  57 5.788269 0.7730179
# 7   5000   0.01     0.9  57 7.315128 0.7728356
# 8   5000   0.01     0.9  57 8.223534 0.7725093
# 9  10000   0.01     0.9  58 7.315128 0.7719940
# 10 20000   0.01     0.9  59 7.315128 0.7719940

ind.keep <- unlist(map(all_keep, max_prs$num))
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp),
  y.sub[ind.test]
) # 76.6 [72.6-80.5]
sum(lpval[ind.keep] > max_prs$thr.lp)  # 156

ggplot(grid2) +
  geom_point(aes(thr.lp, auc)) +
  facet_grid(thr.imp ~ thr.r2 + size) +
  scale_x_log10(limits = c(1, 420)) +
  ylim(0.72, NA) +
  theme_bigstatsr(size.rel = 0.7) +
  labs(x = "-log10(p-value) threshold (log scale)", y = "AUC")


ind3 <- intersect(which(abs(beta) > 2), ind)
gwas <- big_univLogReg(G, y.sub[ind.train], ind.train, ind3)
cbind(beta_GWAS = beta[ind3], pval_GWAS = 10^(-lpval[ind3]),
      af = big_scale()(G, ind.row = ind.train, ind.col = ind3)$center / 2,
      gwas, pval = predict(gwas, log10 = FALSE))
#    beta_GWAS     pval_GWAS         af      estim    std.err niter      score          pval
# 1     -2.001  2.251049e-02 0.99975364  2.3287190 5.84184181     7  0.3986275  6.901677e-01
# 2     -2.169  2.392526e-79 0.06927200 -0.8535355 0.16711651     6 -5.1074274  3.265746e-07
# 3     -3.677  1.387642e-11 0.00306040 -4.5634947 2.55744684     7 -1.7843947  7.435954e-02
# 4      2.100 1.767439e-310 0.69515640  1.2119391 0.09025996     6 13.4272072  4.189145e-41
# 5      2.095 1.787905e-310 0.75285852  1.0006257 0.09333209     6 10.7211325  8.100531e-27
# 6     -2.071  3.110612e-88 0.08276406 -0.7658117 0.14766578     6 -5.1861152  2.147260e-07
# 7     -2.123  5.370921e-83 0.06626514 -0.8967162 0.17864911     6 -5.0194271  5.182581e-07
# 8     -2.669  3.235930e-27 0.00978682 -1.0801250 0.52659284     6 -2.0511578  4.025158e-02
# 9     -2.281  2.244712e-11 0.01285704 -1.9022608 0.68049480     6 -2.7954083  5.183420e-03
# 10    -2.006  5.695059e-99 0.08703216 -0.9611371 0.15661173     6 -6.1370697  8.405749e-10
# 11     2.683  5.610216e-68 0.00072592  2.0885938 0.85729681     8  2.4362552  1.484021e-02
# 12     2.394 8.981408e-311 0.28604648  1.3098143 0.06164685     6 21.2470605 3.509616e-100
# 13     2.316 1.074842e-310 0.10616746  1.0776314 0.06724236     7 16.0260785  8.402040e-58
# 14    -2.158  4.619196e-44 0.03250884 -1.2091644 0.29091046     6 -4.1564833  3.231838e-05
# 15    -2.412  2.190292e-10 0.00817600 -0.6826325 0.44880196     6 -1.5210106  1.282572e-01
# 16    -2.089  3.527391e-03 0.99937586  2.2394513 3.48901633     7  0.6418575  5.209657e-01
