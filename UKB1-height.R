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
POS <- ukb$map$physical.pos
dim(G) # 374,131 x 656,060
file.size(G$backingfile) / 1024^3  # 229 GB

PC <- as.matrix(readRDS("PC.rds"))
plot(PC[sample(nrow(PC), 10e3), ], pch = 20)  # verif: all
plot(PC[sample(sub, 10e3), ], pch = 20)       # verif: caucasian only


summary(mod.base <- lm(height ~ date + sex, data = df0[sub[ind.train], ]))
# 154.6 cm  +  1 cm every 6.4 year since 1900  +  13.3 cm for males
pred.base <- predict(mod.base, df0[sub, ])

system.time(
  mod <- big_spLinReg(G, y[sub][ind.train], ind.train,
                      covar.train = PC[sub[ind.train], ],
                      base.train = pred.base[ind.train],
                      dfmax = Inf,
                      ncores = 10)
) # 21h with dfmax = 50e3  //  61h with dfmax = Inf
# saveRDS(mod, "lasso-mod-height-full.rds")
plot(mod)
## For dfmax = 50e3:
summary(mod)
# # A tibble: 1 x 6
#   alpha validation_loss intercept beta            nb_var message
#   <dbl>           <dbl>     <dbl> <list>           <int> <list>
# 1     1            24.4    -0.285 <dbl [656,080]> 115997 <chr [10]>
summary(mod)$message  # all "Too many variables"
## For dfmax = Inf:
summary(mod)
# # A tibble: 1 x 6
#   alpha validation_loss intercept beta            nb_var message
#   <dbl>           <dbl>     <dbl> <list>           <int> <list>
# 1     1            24.3   -0.0899 <dbl [656,080]> 148341 <chr [10]>
summary(mod)$message  # all "No more improvement"

pred <- pred.base[ind.test] +
  predict(mod, G, ind.test, covar.row = PC[sub[ind.test], ])

library(dplyr)
df0 %>%
  slice(sub[ind.test]) %>%
  mutate(pred = pred) %>%
  group_by(sex) %>%
  summarise(cor(pred, height))
# # A tibble: 2 x 2
#   sex    `cor(pred, height)`
#   <fct>                <dbl>
# 1 Female               0.656  ->  0.657
# 2 Male                 0.657  ->  0.658

COVAR <- cbind(PC, df0$date, df0$sex)
system.time(
  gwas <- big_univLinReg(G, y[sub][ind.train], ind.train,
               covar.train = COVAR[sub[ind.train], ],
               ncores = nb_cores())
) # 68 min

snp_manhattan(gwas, CHR, POS, npoints = 20e3,
              ind.highlight = which(beta[cols_along(G)] != 0))

lpval <- -predict(gwas)
hist(log(lpval))
beta <- summary(mod, best.only = TRUE)$beta[[1]]
hist(log(lpval[which(beta[cols_along(G)] != 0)]))

ind.max <- order(predict(gwas))[1:100e3]
snp_manhattan(gwas, CHR, POS, npoints = 20e3, ind.highlight = ind.max)
system.time(
  mod2 <- big_spLinReg(G, y[sub][ind.train], ind.train,
                       ind.col = ind.max,
                       covar.train = PC[sub[ind.train], ],
                       base.train = pred.base[ind.train],
                       dfmax = Inf,
                       ncores = 10)
) # 48.5h

summary(mod2)
# # A tibble: 1 x 6
#   alpha validation_loss intercept beta            nb_var message
#   <dbl>           <dbl>     <dbl> <list>           <int> <list>
# 1     1            21.7     -5.59 <dbl [100,020]>  70162 <chr [10]>
summary(mod2)$message  # all "No more improvement"

pred2 <- pred.base[ind.test] +
  predict(mod2, G, ind.test, covar.row = PC[sub[ind.test], ])

library(dplyr)
df0 %>%
  slice(sub[ind.test]) %>%
  mutate(pred = pred2) %>%
  group_by(sex) %>%
  summarise(cor(pred, height))
# # A tibble: 2 x 2
#   sex    `cor(pred, height)`
#   <fct>                <dbl>
# 1 Female               0.634
# 2 Male                 0.643

#### C+T ####
res_CT <- lapply(c(0.05, 0.2, 0.8), function(thr.r2) {
  ind.keep <- snp_clumping(G, infos.chr = CHR, ind.row = ind.test,
                           thr.r2 = thr.r2, S = abs(gwas$score), size = 500,
                           is.size.in.bp = TRUE, infos.pos = POS, ncores = nb_cores())
  thrs <- c(0, -log10(5e-08), exp(seq(log(0.1), log(100), length.out = 100)))
  lpS <- -stats::predict(gwas)
  prs <- snp_PRS(G, betas.keep = gwas$estim[ind.keep],
                 ind.test = ind.test, ind.keep = ind.keep, lpS.keep = lpS[ind.keep],
                 thr.list = thrs)
  ind.best <- which.max(apply(prs, 2, cor, y[sub][ind.test]))
  methods <- c("PRS-all", "PRS-stringent", "PRS-max")
  indices <- c(1:2, ind.best)
  lapply(1:3, function(i) {
    k <- indices[i]
    tibble(
      method = methods[i],
      pred = list(prs[, k]),
      thr.r2 = thr.r2,
      set = list(intersect(ind.keep, which(lpS > thrs[k])))
    )
  }) %>% bind_rows()
}) %>% bind_rows()

cor_by_sex <- lapply(res_CT$pred, function(pred) {
  # should have been learned on training set, but let's be optimistic for C+T
  mylm <- lm(y ~ pred + COVAR, data.frame(pred, y = y[sub][ind.test],
                                          COVAR = I(COVAR[sub[ind.test], ])))
  pred2 <- predict(mylm)
  tapply(seq_along(pred2), df0$sex[sub[ind.test]], function(ind) {
    cor(pred2[ind], y[sub][ind.test][ind])
  })
})

res_CT$cor_female <- purrr::map_dbl(cor_by_sex, "Female")
res_CT$cor_male   <- purrr::map_dbl(cor_by_sex, "Male")
saveRDS(res_CT, "height-paper2-gwas.rds")
res_CT
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
