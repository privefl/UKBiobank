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
y <- df0$has_HBP

sub <- which(df0$is_caucasian & !is.na(y) & !is.na(ind.indiv))
length(sub)  # 349,586
table(y[sub])
#  FALSE   TRUE
# 260330  89256



system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKB_imp_HBP",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = 22
  )
) # 2.5H


set.seed(1)
ind.train <- sort(sample(length(sub), 320e3))
ind.test <- setdiff(seq_along(sub), ind.train)

library(bigsnpr)
ukb <- snp_attach("data/UKB_imp_HBP.rds")
G <- ukb$genotypes
CHR <- as.numeric(ukb$map$chromosome)
POS <- ukb$map$physical.pos
dim(G) # 349,586 x 656,060
file.size(G$backingfile) / 1024^3  # 214 GB

PC <- as.matrix(readRDS("PC.rds"))
plot(PC[sample(nrow(PC), 10e3), ], pch = 20)  # verif: all
plot(PC[sample(sub, 10e3), ], pch = 20)       # verif: caucasian only


summary(mod.base <- glm(has_HBP ~ sex, data = df0[sub[ind.train], ], family = "binomial"))
pred.base <- predict(mod.base, df0[sub, ])

system.time(
  mod <- big_spLogReg(G, y[sub][ind.train] + 0, ind.train,
                      covar.train = PC[sub[ind.train], ],
                      base.train = pred.base[ind.train],
                      dfmax = Inf,
                      ncores = 10)
) # 7h with dfmax = Inf
# saveRDS(mod, "lasso-mod-height-full.rds")
plot(mod)
## For dfmax = 50e3:
summary(mod)
# # A tibble: 1 x 6
#   alpha validation_loss intercept beta            nb_var message
#   <dbl>           <dbl>     <dbl> <list>           <int> <list>
# 1     1           0.545    -0.954 <dbl [656,080]>  47762 <chr [10]>
summary(mod)$message  # all "No more improvement"

# WARNING: do not sum things on different scales
pred <- pred.base[ind.test] +
  predict(mod, G, ind.test, covar.row = PC[sub[ind.test], ], proba = FALSE)

library(dplyr)
df0 %>%
  slice(sub[ind.test]) %>%
  mutate(pred = pred) %>%
  group_by(sex) %>%
  summarise(AUC(pred, has_HBP + 0))
# # A tibble: 2 x 2
#   sex    `AUC(pred, has_HBP + 0)`
#   <fct>                <dbl>
# 1 Female               0.658
# 2 Male                 0.620


system.time(
  mod2 <- big_spLogReg(G, y[sub][ind.train] + 0, ind.train,
                       covar.train = PC[sub[ind.train], ],
                       base.train = pred.base[ind.train],
                       alpha = c(0.5, 0.1, 0.01),
                       ind.sets = attr(mod, "ind.sets"),
                       dfmax = Inf, n.abort = 3,
                       ncores = nb_cores())
) # 7h
plot(mod2)
## For dfmax = 50e3:
summary(mod2)
# A tibble: 3 x 6
#   alpha validation_loss intercept beta            nb_var message
#   <dbl>           <dbl>     <dbl> <list>           <int> <list>
# 1  0.01           0.545    -0.977 <dbl [656,080]>  63712 <chr [10]>
# 2  0.1            0.545    -0.992 <dbl [656,080]>  50370 <chr [10]>
# 3  0.5            0.545    -0.984 <dbl [656,080]>  49195 <chr [10]>
summary(mod2)$message  # all "No more improvement"

# WARNING: do not sum things on different scales
pred2 <- pred.base[ind.test] +
  predict(mod2, G, ind.test, covar.row = PC[sub[ind.test], ], proba = FALSE)

library(dplyr)
df0 %>%
  slice(sub[ind.test]) %>%
  mutate(pred = pred2) %>%
  group_by(sex) %>%
  summarise(AUC(pred, has_HBP + 0))
# # A tibble: 2 x 2
#   sex    `AUC(pred, has_HBP + 0)`
#   <fct>                     <dbl>
# 1 Female                    0.659
# 2 Male                      0.621

COVAR <- cbind(PC, df0$sex, df0$age)
system.time(
  gwas <- big_univLogReg(G, y[sub][ind.train] + 0, ind.train,
                         covar.train = COVAR[sub[ind.train], ],
                         ncores = nb_cores())
) # 68 min

snp_manhattan(gwas, CHR, POS, npoints = 20e3,
              ind.highlight = which(beta[cols_along(G)] != 0))



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
  ind.best <- which.max(apply(prs, 2, AUC, y[sub][ind.test]))
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
