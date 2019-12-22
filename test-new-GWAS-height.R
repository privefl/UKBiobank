library(doParallel)
registerDoParallel(cl <- makeCluster(12))
system.time({
  list_snp_id <- foreach(chr = 1:22) %dopar% {
    # cat("Processing chromosome", chr, "..\n")
    mfi <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
    infos_chr <- bigreadr::fread2(mfi, showProgress = FALSE)
    infos_chr_sub <- subset(infos_chr, V6 > 0.01)  ## MAF > 1%
    bim <- paste0("data/ukb_snp_chr", chr, "_v2.bim")
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

library(dplyr)
csv <- "ukb22544.csv"
df0 <- bigreadr::fread2(
  csv,
  select = c("eid", "50-0.0", "34-0.0", "52-0.0", "22001-0.0", "22006-0.0"),
  col.names = c("eid", "height", "year", "month", "sex", "is_caucasian")
) %>%
  mutate(
    sex  = factor(sex, levels = c(0, 1),  labels = c("Female", "Male")),
    is_caucasian = as.logical(is_caucasian),
    date = (year - 1900) + (month - 0.5) / 12,
    year = NULL, month = NULL
  )
ind.indiv <- match(df0$eid, sample$ID_2)
y <- df0$height

sub <- which(df0$is_caucasian & !is.na(y) & !is.na(ind.indiv) &
               !is.na(df0$sex) & !is.na(df0$date))
length(sub)  # 408,036

library(bigsnpr)
NCORES <- nb_cores()
system.time(
  rds <- snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKB_imp_height",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 2H

library(bigsnpr)
ukb <- snp_attach("data/UKB_imp_height.rds")
G <- ukb$genotypes
CHR <- as.numeric(ukb$map$chromosome)
POS <- ukb$map$physical.pos
dim(G) # 408,036 x 656,060
file.size(G$backingfile) / 1024^3  # 249 GB

nPC <- 10
PC <- bigreadr::fread2(
  csv,
  select = paste0("22009-0.", 1:nPC),
  col.names = paste0("PC", 1:nPC)
)
COVAR <- cbind(as.matrix(PC), df0$sex, df0$date)[sub, ]

# system.time(
#   mod <- big_spLinReg(G, y[sub], alphas = c(1, 0.1),
#                       covar.train = COVAR, pf.covar = rep(0, nPC + 2),
#                       dfmax = Inf, n.abort = 2, ncores = NCORES)
# ) # 63H
# saveRDS(mod, "elnet-mod-height-full.rds")
mod <- readRDS("elnet-mod-height-full.rds")
plot(mod)
summary(mod)
# # A tibble: 2 x 6
# alpha validation_loss intercept beta            nb_var message
#   <dbl>           <dbl>     <dbl> <list>           <int> <list>
# 1   0.1            23.6      142. <dbl [656,072]> 193173 <chr [10]>
# 2   1              23.5      141. <dbl [656,072]> 163303 <chr [10]>
summary(mod)$message  # all "No more improvement"

beta <- head(summary(mod, best.only = TRUE)$beta[[1]], ncol(G))
poly <- lapply(split(cols_along(G), CHR), function(ind) {
  ind2 <- ind[beta[ind] != 0]
  big_prodVec(G, beta[ind2], ind.col = ind2)
})

system.time(
  gwas0 <- big_univLinReg(G, y[sub], covar.train = COVAR, ncores = NCORES)
) # 499 sec

system.time(
  gwas_all_chr <- lapply(1:22, function(chr) {
    big_univLinReg(G, y[sub], ind.col = which(CHR == chr),
                   covar.train = cbind(COVAR, Reduce('+', poly[-chr])),
                   ncores = NCORES)
  })
) # 1108 sec
gwas <- do.call("rbind", gwas_all_chr)
attributes(gwas) <- attributes(gwas_all_chr[[1]])

library(ggplot2)
plot(gwas0) + scale_y_log10()
plot(gwas) + scale_y_log10()
snp_manhattan(gwas0, CHR, POS, npoints = 20e3)
snp_manhattan(gwas,  CHR, POS, npoints = 20e3)
snp_qq(gwas0) + xlim(1, NA)
snp_qq(gwas) + xlim(1, NA)
qplot(gwas0$score^2, gwas$score^2) +
  xlim(4, NA) +
  theme_bigstatsr() +
  geom_abline(color = "blue") +
  geom_abline(slope = 2, color = "red") +
  labs(x = "X2 with standard GWAS", y = "X2 when including polygenic effect")
qplot(gwas0$score^2, gwas$score^2, alpha = I(0.05)) +
  xlim(NA, 10) + ylim(NA, 20) +
  theme_bigstatsr() +
  geom_abline(color = "blue") +
  geom_abline(slope = 2, color = "red") +
  geom_smooth(color = "green") +
  labs(x = "X2 with standard GWAS", y = "X2 when including polygenic effect")
