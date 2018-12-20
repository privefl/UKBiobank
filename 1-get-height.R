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
    logit <- big_spLinReg(G, y[sub][ind.train], ind.train,
                          ind.col = which(CHR %in% chrs[[i]]),
                          covar.train = PC[sub[ind.train], ],
                          base.train = pred.base[ind.train],
                          ncores = 10)
  )[3]
  preds[, i] <- pred.base[ind.test] +
    predict(logit, G, ind.test, covar.row = PC[sub[ind.test], ])
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
ggplot(res) +
  theme_bigstatsr() +
  geom_point(aes(nvar, cor, color = sex), size = 3) +
  geom_line(aes(nvar, cor, color = sex), size = 1.5) +
  labs(x = "Number of variables (350K individuals)",
       y = "Correlation between predicted and true heights",
       color = "Sex") +
  theme(legend.position = c(0.2, 0.85)) +
  geom_point(aes(nvar, time / 3600 / 40), size = 3) +
  geom_line(aes(nvar, time / 3600 / 40), size = 1.5) +
  scale_y_continuous(sec.axis = sec_axis(~.*40, name = "Fitting time (in hours)"))

ggsave("height-paper2.pdf", width = 887, height = 755, scale = 1/100)
