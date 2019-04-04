download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/NikpayM_26343387_GCST003116/cad.add.160614.website.txt",
              destfile = "sumstats_CAD.txt")
library(bigreadr)
sumstats <- fread2("sumstats_CAD.txt",
                   select = c("chr", "bp_hg19", "noneffect_allele",
                              "effect_allele", "beta", "p_dgc"))
nrow(sumstats)  # 9,455,778
hist(sumstats$beta)
hist(sumstats$p_dgc)

# augment dataset to match reverse alleles
sumstats2 <- sumstats
sumstats2$noneffect_allele <- sumstats$effect_allele
sumstats2$effect_allele <- sumstats$noneffect_allele
sumstats2$beta <- -sumstats$beta
sumstats2 <- rbind(sumstats, sumstats2)

# match variants with UKBB
library(doParallel)
NCORES <- 12
registerDoParallel(cl <- makeCluster(NCORES))
info_snp <- foreach(chr = 1:22, .packages = "dplyr") %dopar% {

  paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt") %>%
    bigreadr::fread2() %>%
    inner_join(
      sumstats2[sumstats2$chr == chr, ], ., na_matches = "never",
      by = c(bp_hg19 = "V3", noneffect_allele = "V4", effect_allele = "V5")
    ) %>%
    arrange(bp_hg19) %>%
    transmute(
      id    = paste(chr, bp_hg19, noneffect_allele, effect_allele, sep = "_"),
      beta  = beta,
      lpval = -log10(p_dgc),
      info  = V8
    ) %>%
    na.omit() %>%
    filter(lpval > 1 & info > 0.3)
}
stopCluster(cl)

list_snp_id <-  lapply(info_snp, function(df) df$id)
beta <-  unlist(lapply(info_snp, function(df) df$beta))
lpval <- unlist(lapply(info_snp, function(df) df$lpval))
info <-  unlist(lapply(info_snp, function(df) df$info))
length(beta)  # 897,949

# subset samples
sample <- fread2("ukb25589_imp_chr1_v3_s487327.sample")
str(sample)
(N <- readBin("data/ukb_imp_chr1_v3.bgen", what = 1L, size = 4, n = 4)[4])
sample <- sample[-1, ]
stopifnot(nrow(sample) == N)

csv <- "ukb22544.csv"
df0 <- fread2(csv, select = c("eid", "22006-0.0"),
              col.names = c("eid", "is_caucasian"))
ind.indiv <- match(df0$eid, sample$ID_2)
# related samples
df0$is_rel2 <- df0$eid %in% fread2("ukb25589_rel_s488346.dat")$ID2


df_heart <- fread2(csv, select = c(paste0("6150-0.", 0:3),
                                   paste0("6150-1.", 0:3),
                                   paste0("6150-2.", 0:3)))
df_illness <- fread2(csv, select = c(paste0("20002-0.", 0:28),
                                     paste0("20002-1.", 0:28),
                                     paste0("20002-2.", 0:28)))
df_ICD10 <- fread2(csv, select = c(paste0("40001-", 0:2, ".0"),
                                   paste0("40002-0.", 0:13),
                                   paste0("40002-1.", 0:13),
                                   paste0("40002-2.", 0:13),
                                   paste0("41202-0.", 0:379),
                                   paste0("41204-0.", 0:434)))
ind_CAD <- sort(unique(unlist(c(
  lapply(df_heart,   function(x) which(x == 1)),
  lapply(df_illness, function(x) which(x == 1075)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 3) %in% paste0("I", 21:23))),
  lapply(df_ICD10,   function(x) which(x == "I252"))
))))
ind_heart <- sort(unique(unlist(c(
  lapply(df_heart,   function(x) which(x %in% 1:3)),
  lapply(df_illness, function(x) which(x %in% 1074:1080)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 1) == "I"))
))))
y <- rep(0, nrow(df0)); y[ind_heart] <- NA; y[ind_CAD] <- 1

sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian & !is.na(y))
length(sub)  # 238,190
table(y.sub <- y[sub])
#      0      1
# 225927  12263

system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_CAD",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 2.7H


library(bigsnpr)
ukbb <- snp_attach("data/UKBB_CAD.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 199 GB
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos

set.seed(1)
ind.train <- sort(sample(length(sub), 200e3))
ind.test <- setdiff(seq_along(sub), ind.train)

system.time(
  all_keep <- snp_grid_clumping(
    G, CHR, POS, lpS = lpval, ind.row = ind.train, infos.imp = info,
    grid.thr.imp = c(0.3, 0.6, 0.9, 0.95),
    grid.thr.r2 = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95),
    grid.base.size = c(50, 100, 200, 500),
    ncores = 6
  )
) # 1.8H
plot(lengths(all_keep[[1]]), pch = 20)

system.time(
  multi_PRS <- snp_grid_PRS(
    G, all_keep, betas = beta, lpS = lpval, ind.row = ind.train,
    grid.lpS.thr = seq_log(1, 0.999 * max(lpval), 50),
    backingfile = "data/UKBB_CAD_scores", ncores = NCORES
  )
) # 2.6H

system.time(
  final_mod <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train], ncores = NCORES, n.abort = 3
  )
) # 4.6H
mod <- final_mod$mod
plot(mod)
summary(mod)
# # A tibble: 5 x 6
#      alpha validation_loss intercept beta           nb_var message
#      <dbl>           <dbl>     <dbl> <list>          <int> <list>
#   1 0.0001           0.197     -3.46 <dbl [65,940]>  28921 <chr [10]>
#   2 0.001            0.197     -3.42 <dbl [65,940]>  11020 <chr [10]>
#   3 0.01             0.197     -3.44 <dbl [65,940]>   6675 <chr [10]>
#   4 0.1              0.197     -3.42 <dbl [65,940]>   6164 <chr [10]>
#   5 1                0.197     -3.42 <dbl [65,940]>   5886 <chr [10]>

new_beta <- final_mod$beta.G

length(ind <- which(new_beta != 0))  # 315,373
summary(new_beta)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# -0.138676  0.000000  0.000000  0.000016  0.000000  0.326784
summary(new_beta[which(sign(new_beta * beta) < 0)])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -4.674e-02 -1.193e-04 -1.930e-06 -1.685e-05  9.601e-05  3.452e-02

pred <- final_mod$intercept +
  big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind)

AUCBoot(pred, y.sub[ind.test])  # 64.0 [62.7-65.2]


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
#   size thr.r2 thr.imp num   thr.lp      auc
# 1  500    0.2     0.3  14 3.370896 0.603575

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
) # 59.9 [58.6-61.2]
sum(lpval[ind.keep] > std_prs$thr.lp)  # 1184

max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
#    size thr.r2 thr.imp num   thr.lp       auc
# 1   526   0.95    0.95 112 1.453393 0.6190787
# 2   210   0.95    0.95 111 1.453393 0.6190581
# 3   105   0.95    0.95 110 1.453393 0.6188843
# 4   526   0.95    0.95 112 1.323693 0.6188750
# 5   210   0.95    0.95 111 1.323693 0.6188560
# 6   105   0.95    0.95 110 1.323693 0.6187532
# 7   250   0.80    0.95 107 1.453393 0.6185265
# 8   625   0.80    0.95 108 1.453393 0.6184453
# 9   625   0.80    0.95 108 1.595802 0.6184363
# 10  526   0.95    0.95 112 1.205568 0.618408

ind.keep <- unlist(map(all_keep, max_prs$num))
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp),
  y.sub[ind.test]
) # 61.1 [59.9-62.4]
sum(lpval[ind.keep] > max_prs$thr.lp)  # 87,572

ggplot(grid2) +
  geom_point(aes(thr.lp, auc)) +
  facet_grid(thr.imp ~ thr.r2 + size) +
  scale_x_log10(limits = c(1, 7), breaks = c(1, 2, 5), minor_breaks = 1:10) +
  ylim(0.58, NA) +
  theme_bigstatsr(size.rel = 0.7) +
  labs(x = "-log10(p-value) threshold (log scale)", y = "AUC")
