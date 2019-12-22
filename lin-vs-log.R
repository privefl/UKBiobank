# urls <- gsubfn::strapply(
#   readLines("https://datadryad.org//resource/doi:10.5061/dryad.ns8q3"),
#   "<a href=\"(/bitstream/handle/10255/dryad\\.[0-9]+/meta_chr_[0-9]+\\?sequence=1)\">",
#   simplify = 'c')
#
# sumstats <- purrr::map_dfr(urls, ~ {
#   download.file(paste0("https://datadryad.org", .x),
#                 destfile = (tmp <- tempfile(fileext = ".txt")))
#   sumstats <- bigreadr::fread2(
#     tmp, select = c("chromosome", "position", "a0", "a1", "beta.meta", "p.meta"),
#     col.names = c("chr", "pos", "a0", "a1", "beta", "p"))
#   na.omit(sumstats)
# })
#
# saveRDS(sumstats, "sumstats_T1D.rds")

sumstats <- readRDS("sumstats_T1D.rds")
nrow(sumstats)  # 8,996,866
hist(sumstats$p)
hist(sumstats$beta)

sumstats <- subset(sumstats, p < 0.1)

library(bigreadr)
info_snp_UKBB <- rbind_df(lapply(1:22, function(chr) {
  file <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
  df <- fread2(file, select = c(3:5, 8), col.names = c("pos", "a0", "a1", "info"))
  cbind.data.frame(chr = chr, df)
}))
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB)
# 1,241,172 variants in summary statistics.
# 172,374 ambiguous SNPs have been removed.
# 1,051,668 variants have been matched; 76 were flipped and 522 were reversed.
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB, strand_flip = FALSE)
# 1,222,801 variants have been matched; 0 were flipped and 471 were reversed.
info_snp <- subset(na.omit(info_snp), info > 0.3)

list_snp_id <- with(info_snp, split(paste(chr, pos, a0, a1, sep = "_"),
                                    factor(chr, levels = 1:22)))
beta <- info_snp$beta
lpval <- -log10(info_snp$p)
info <- info_snp$info

# Infinite values because of large effects on chromosome 6
lpval <- pmin(lpval, -log10(.Machine$double.xmin) + abs(beta))

# subset samples
library(bigreadr)
sample <- fread2("ukb25589_imp_chr1_v3_s487327.sample")
str(sample)
(N <- readBin("data/ukb_imp_chr1_v3.bgen", what = 1L, size = 4, n = 4)[4])
sample <- sample[-1, ]
stopifnot(nrow(sample) == N)

csv <- "data/ukb22544.csv"
df0 <- fread2(csv, select = c("eid", "22006-0.0"),
              col.names = c("eid", "is_caucasian"))
ind.indiv <- match(df0$eid, sample$ID_2)
# related samples
df0$is_rel2 <- df0$eid %in% fread2("ukb25589_rel_s488346.dat")$ID2

df_illness <- fread2(csv, select = c(paste0("20002-0.", 0:28),
                                     paste0("20002-1.", 0:28),
                                     paste0("20002-2.", 0:28)))
df_ICD10 <- fread2(csv, select = c(paste0("40001-", 0:2, ".0"),
                                   paste0("40002-0.", 0:13),
                                   paste0("40002-1.", 0:13),
                                   paste0("40002-2.", 0:13),
                                   paste0("41202-0.", 0:379),
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

sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian & !is.na(y))
length(sub)  # 315,318
table(y.sub <- y[sub])
#      0      1
# 314547    771

NCORES <- 12
system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_T1D",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 4H

library(bigsnpr)
ukbb <- snp_attach("data/UKBB_T1D.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 359 GB
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos

ind.train <- sort(sample(length(sub), 200e3))
ind.test <- setdiff(seq_along(sub), ind.train)

system.time(
  lin <- big_spLinReg(G, y.sub[ind.train], ind.train = ind.train, K = 4,
                      ncores = nb_cores())
)
lin_pred <- predict(lin, G, ind.test)
AUCBoot(lin_pred, y.sub[ind.test]) # 75.7 [72.7-78.5]

system.time(
  log <- big_spLogReg(G, y.sub[ind.train], ind.train = ind.train, K = 4,
                      ncores = nb_cores())
)
log_pred <- predict(log, G, ind.test)
AUCBoot(log_pred, y.sub[ind.test]) # 76.7 [73.8-79.6]

# remotes::install_github("privefl/paper2-PRS/pkg.paper.PRS")
y.simu <- pkg.paper.PRS::get_pheno(G, h2 = 0.4, M = 5000,
                                   effects.dist = "laplace", K = 0.02)$pheno

ind.train <- sort(sample(nrow(G), 250e3))
ind.test <- setdiff(rows_along(G), ind.train)

system.time(
  lin <- big_spLinReg(G, y.simu[ind.train], ind.train = ind.train, K = 5,
                      ncores = nb_cores())
)
lin_pred <- predict(lin, G, ind.test)
AUCBoot(lin_pred, y.simu[ind.test]) # 82.7 [82.5-82.9]

system.time(
  log <- big_spLogReg(G, y.simu[ind.train], ind.train = ind.train, K = 5,
                      ncores = nb_cores())
)
log_pred <- predict(log, G, ind.test, proba = FALSE)
AUCBoot(log_pred, y.simu[ind.test]) # 82.7 [82.5-83.0]

plot(lin_pred, log_pred, pch = 20); abline(MASS::rlm(log_pred ~ lin_pred), col = "red")
# Useful in the tails?
