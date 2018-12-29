library(doParallel)
registerDoParallel(cl <- makeCluster(12))
system.time({
  list_snp_id <- foreach(chr = 1:22) %dopar% {
    # cat("Processing chromosome", chr, "..\n")
    mfi <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
    infos_chr <- bigreadr::fread2(mfi, showProgress = FALSE)
    infos_chr_sub <- subset(infos_chr, V6 > 0.01)  ## MAF > 1%
    bim <- paste0("data/ukb_snp_bim/ukb_snp_chr", chr, "_v2.bim")
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

df0 <- readRDS("pheno.rds")
ind.indiv <- match(df0$eid, sample$ID_2)
y <- df0$height

sub <- which(df0$is_caucasian & !is.na(y) & !is.na(ind.indiv))
length(sub)  # 374,131


system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKB_imp_height",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = 10
  )
) # 2H


set.seed(1)
ind.train <- sort(sample(length(sub), 350e3))
ind.test <- setdiff(seq_along(sub), ind.train)

library(bigsnpr)
ukb <- snp_attach("data/UKB_imp_height.rds")
G <- ukb$genotypes
CHR <- as.numeric(ukb$map$chromosome)
dim(G) # 374,131 x 656,060
file.size(G$backingfile) / 1024^3  # 229 GB

PC <- as.matrix(readRDS("PC.rds"))
plot(PC[sample(nrow(PC), 10e3), ], pch = 20)  # verif: all
plot(PC[sample(sub, 10e3), ], pch = 20)       # verif: caucasian only


summary(mod.base <- lm(height ~ date + sex, data = df0[sub[ind.train], ]))
# 154.6 cm  +  1 cm every 6.4 year since 1900  +  13.3 cm for males
pred.base <- predict(mod.base, df0[sub, ])

chrs <- list(1, 1:2, 1:6, 1:22)
times <- rep(NA, 4)
preds <- matrix(NA, length(ind.test), 4)
for (i in 1:4) {
  times[i] <- system.time(
    mod <- big_spLinReg(G, y[sub][ind.train], ind.train,
                        ind.col = which(CHR %in% chrs[[i]]),
                        covar.train = PC[sub[ind.train], ],
                        base.train = pred.base[ind.train],
                        ncores = 10)
  )[3]
  preds[, i] <- pred.base[ind.test] +
    predict(mod, G, ind.test, covar.row = PC[sub[ind.test], ])
}

library(dplyr)
res <- foreach(i = 1:4, .combine = "rbind") %do% {
  df0 %>%
    slice(sub[ind.test]) %>%
    mutate(pred = preds[, i]) %>%
    group_by(sex) %>%
    summarise(cor = cor(pred, height)) %>%
    cbind.data.frame(nvar = sum(CHR %in% chrs[[i]]),
                     time = times[i])
}
saveRDS(res, "height-paper2.rds")
res
#      sex       cor   nvar      time
# 1 Female 0.2334127  51572   386.673
# 2   Male 0.2271432  51572   386.673
# 3 Female 0.3204401 104235  7709.831
# 4   Male 0.3106693 104235  7709.831
# 5 Female 0.4488728 276683 40238.400
# 6   Male 0.4523325 276683 40238.400
# 7 Female 0.6545584 656060 75378.401
# 8   Male 0.6551584 656060 75378.401

library(ggplot2)
ggplot(res, aes(nvar, time / 3600)) +
  theme_bigstatsr() +
  geom_point(size = 3) +
  geom_line(size = 1.5) +
    labs(x = "Number of variables (350K individuals)",
         y = "Fitting time (in hours)")

ggsave("height-paper2.pdf", width = 887, height = 755, scale = 1/100)


#### C+T ####
COVAR <- cbind(PC, df0$date, df0$sex)
system.time(
  res2 <- pkg.paper.PRS::PRS(G, CHR, ukb$map$physical.pos, y[sub], COVAR[sub, ],
                             ind.train, ind.test, family = "gaussian")
) # 69 min for the GWAS part

cor_by_sex <- lapply(res2$pred, function(pred) {
  # should have been learned on training set, but let's be optimistic for C+T
  mylm <- lm(y ~ pred + COVAR, data.frame(pred, y = y[sub][ind.test],
                                          COVAR = I(COVAR[sub[ind.test], ])))
  pred2 <- predict(mylm)
  tapply(seq_along(pred2), df0$sex[sub[ind.test]], function(ind) {
    cor(pred2[ind], y[sub][ind.test][ind])
  })
})

res2$cor_female <- purrr::map_dbl(cor_by_sex, "Female")
res2$cor_male   <- purrr::map_dbl(cor_by_sex, "Male")
saveRDS(res2, "height-paper2-gwas.rds")
res2
# # A tibble: 9 x 7
#   method        pred          thr.r2 set           timing cor_female cor_male
#   <chr>         <list>         <dbl> <list>         <dbl>      <dbl>    <dbl>
# 1 PRS-all       <dbl [24,131…   0.05 <int [158,38…  4148.      0.528    0.537
# 2 PRS-stringent <dbl [24,131…   0.05 <int [692]>    4148.      0.454    0.457
# 3 PRS-max       <dbl [24,131…   0.05 <int [13,715…  4148.      0.551    0.557
# 4 PRS-all       <dbl [24,131…   0.2  <int [292,93…  4148.      0.537    0.551
# 5 PRS-stringent <dbl [24,131…   0.2  <int [1,055]>  4148.      0.426    0.429
# 6 PRS-max       <dbl [24,131…   0.2  <int [45,570…  4148.      0.549    0.561
# 7 PRS-all       <dbl [24,131…   0.8  <int [575,15…  4148.      0.484    0.494
# 8 PRS-stringent <dbl [24,131…   0.8  <int [2,635]>  4148.      0.327    0.327
# 9 PRS-max       <dbl [24,131…   0.8  <int [575,15…  4148.      0.484    0.494
