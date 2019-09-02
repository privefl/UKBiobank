readLines("ukb23144.txt", n = 4)
library(bigreadr)
HLA <- fread2("ukb23144.txt", skip = 1, sep = "\t", sep2 = ",")

csv <- "ukb22544.csv"
df0 <- fread2(csv, select = c("eid", "22006-0.0"),
              col.names = c("eid", "is_caucasian"))
ind.indiv <- match(df0$eid, HLA$V2)
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
length(sub)  # 316,038
table(y.sub <- y[sub])
#      0      1
# 315267    771


writeLines(HLA$V3[ind.indiv[sub]], con = (tmp <- tempfile()))
alleles <- big_read(tmp, select = seq_along(scan(text = HLA$V3[[1]], sep = ",")))

gwas <- big_univLogReg(alleles, y.sub, ncores = nb_cores())
plot(gwas)
plot(gwas, type = "Manhattan")
plot(gwas, type = "Volcano")
gwas[which(predict(gwas, log10 = FALSE) < 0.01), ]

set.seed(1)
ind.train <- sort(sample(length(sub), 200e3))
ind.test <- setdiff(seq_along(sub), ind.train)

mod <- big_spLogReg(alleles, y.sub[ind.train], ind.train = ind.train,
                    alphas = 10^(-(0:4)), K = 5, ncores = nb_cores())
plot(mod)
summary(mod)

pred <- predict(mod, alleles, ind.test)
AUCBoot(pred, y.sub[ind.test])  # 74.6 [71.2-78.0]

library(ggplot2)
ggplot(data.frame(pred = pred, pheno = y.sub[ind.test])) +
  theme_bigstatsr() +
  geom_density(aes(pred, fill = as.factor(pheno)), alpha = 0.3) +
  scale_x_continuous(limits = c(0, 0.03)) +
  theme(legend.position = c(0.7, 0.8)) +
  labs(x = "Prediction", fill = "Phenotype")
