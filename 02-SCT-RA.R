download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/OkadaY_24390342_GCST002318/RA_GWASmeta_TransEthnic_v2.txt.gz",
              destfile = "sumstats_RA.txt.gz")
R.utils::gunzip("sumstats_RA.txt.gz")
library(bigreadr)
sumstats <- fread2("sumstats_RA.txt", select = 2:7,
                   col.names = c("chr", "pos", "a1", "a2", "or", "p"))
nrow(sumstats)  # 9,739,303

# augment dataset to match reverse alleles
sumstats2 <- sumstats
sumstats2$a1 <- sumstats$a2
sumstats2$a2 <- sumstats$a1
sumstats2$or <- 1 / sumstats$or
sumstats2 <- rbind(sumstats, sumstats2)

# match variants with UKBB
library(doParallel)
NCORES <- 12
registerDoParallel(cl <- makeCluster(NCORES))
info_snp <- foreach(chr = 1:22) %dopar% {
  # cat("Processing chromosome", chr, "..\n")
  file <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
  infos_chr <- bigreadr::fread2(file, showProgress = FALSE, nThread = 1)
  sumstats_joined <- dplyr::inner_join(
    sumstats2[sumstats2$chr == chr, ], infos_chr, na_matches = "never",
    by = c(pos = "V3", a1 = "V4", a2 = "V5")
  )
  df <- na.omit(dplyr::transmute(
    sumstats_joined,
    id = paste(chr, pos, a1, a2, sep = "_"),
    beta = log(or),
    lpval = -log10(p),
    info = V8
  ))
  subset(df, lpval > 1 & info > 0.3)
}
stopCluster(cl)

list_snp_id <-  lapply(info_snp, function(df) df$id)
beta <-  unlist(lapply(info_snp, function(df) -df$beta))
lpval <- unlist(lapply(info_snp, function(df) df$lpval))
info <-  unlist(lapply(info_snp, function(df) df$info))
length(beta)  # 1,010,055

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

# Non-cancer illness codes (https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=6)
df_illness <- fread2(csv, select = c(paste0("20002-0.", 0:32),
                                     paste0("20002-1.", 0:32),
                                     paste0("20002-2.", 0:32)))
ind_RA <- sort(unique(unlist(
  lapply(df_illness, function(x) which(x == 1464))
)))
ind_muscu <- sort(unique(unlist(
  lapply(df_illness, function(x) which(x %in% c(1295, 1464:1467, 1477, 1538)))
)))
y <- rep(0, nrow(df0)); y[ind_muscu] <- NA; y[ind_RA] <- 1

sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian & !is.na(y) &
               !is.na(df0$sex))
length(sub)  # 295,009
table(y.sub <- y[sub])
#      0      1
# 291155   3854

system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_RA",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 2.9H


library(bigsnpr)
ukbb <- snp_attach("data/UKBB_RA.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 278 GB
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
) # 9.3H -> swapping? too many cores?
plot(lengths(all_keep[[1]]), pch = 20)

system.time(
  multi_PRS <- snp_grid_PRS(
    G, all_keep, betas = beta, lpS = lpval, ind.row = ind.train,
    grid.lpS.thr = exp(seq(log(1), log(0.999 * max(lpval)), length.out = 50)),
    backingfile = "data/UKBB_RA_scores", ncores = NCORES
  )
) # 2.7H

nPC <- 20
PC <- fread2(csv, select = paste0("22009-0.", 1:nPC),
             col.names = paste0("PC", 1:nPC))
COVAR <- as.matrix(cbind(df0$sex, PC)[sub, ])

system.time(
  final_mod <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train], covar.train = COVAR[ind.train, ], ncores = NCORES
  )
) # 3.4H
mod <- final_mod$mod
plot(mod)
summary(mod)
# # A tibble: 5 x 6
#      alpha validation_loss intercept beta           nb_var message
#      <dbl>           <dbl>     <dbl> <list>          <int> <list>
#   1 0.0001          0.0676     -4.21 <dbl [63,525]>   9712 <chr [10]>
#   2 0.001           0.0676     -4.23 <dbl [63,525]>   3829 <chr [10]>
#   3 0.01            0.0676     -4.25 <dbl [63,525]>   2489 <chr [10]>
#   4 0.1             0.0676     -4.26 <dbl [63,525]>   2000 <chr [10]>
#   5 1               0.0676     -4.26 <dbl [63,525]>   2040 <chr [10]>

new_beta <- final_mod$beta.G

length(ind <- which(new_beta != 0))  # 276,688
summary(new_beta)
summary(new_beta[which(sign(new_beta * beta) < 0)])

pred <- final_mod$intercept +
  big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind) +
  COVAR[ind.test, ] %*% final_mod$beta.covar

AUCBoot(pred, y.sub[ind.test])  # 65.5 [63.1-67.8]


max(aucs <- apply(multi_PRS, 2, AUC, y.sub[ind.train]))  # 61.6

library(tidyverse)

ind2 <- sort(sample(ind, 10e3))
ggplot(data.frame(y = new_beta, x = beta)[ind, ]) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  geom_abline(slope = 0, intercept = 0, color = "blue") +
  geom_point(aes(x, y), size = 0.6) +
  theme_bigstatsr() +
  labs(x = "Effect sizes from GWAS", y = "Non-zero effect sizes from SCT")

grid2 <- attr(all_keep, "grid") %>%
  filter(chr == 1) %>%
  mutate(thr.lp = list(attr(multi_PRS, "grid.lpS.thr")), chr = NULL) %>%
  unnest() %>%
  mutate(auc = aucs)

grid2 %>%
  filter(thr.imp == 0.3, thr.r2 == 0.2, size == 500, thr.lp <= 8) %>%
  arrange(desc(auc)) %>%
  slice(1)
#   size thr.r2 thr.imp   thr.lp       auc
# 1  500    0.2     0.3 2.200396 0.6124643

grid2 %>% arrange(desc(auc)) %>% slice(1:10)

ggplot(grid2) +
  geom_point(aes(thr.lp, auc)) +
  facet_grid(thr.imp ~ thr.r2 + size) +
  scale_x_log10(limits = c(1, 10), breaks = c(1, 2, 5), minor_breaks = 1:10) +
  ylim(0.58, NA) +
  geom_hline(yintercept = max(aucs), color = "blue", linetype = 2) +
  geom_hline(yintercept = AUC(pred, y.sub[ind.test]), color = "red", linetype = 2) +
  theme_bigstatsr(size.rel = 0.7) +
  labs(x = "-log10(p-value) threshold (log scale)", y = "AUC")
