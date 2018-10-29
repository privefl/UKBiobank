
df0 <- readRDS("pheno.rds")

head(sort(table(df0$cancer_type), decreasing = TRUE))


filter(df0, cancer_type == "breast cancer", sex == "Male") %>%
  print(n = Inf)


y <- df0$has_cancer
y[y] <- NA
y[df0$cancer_type == "breast cancer"] <- 1
ind <- which(df0$sex == "Female" & !is.na(y))
system.time(
  gwas.breast <- big_univLogReg(G, y[ind], ind.train = ind,
                                covar.train = PCs[ind, ],
                                ncores = nb_cores())
) # 9h

plot(gwas.breast)
snp_qq(gwas.breast) + xlim(2, NA)
snp_manhattan(snp_gc(gwas.breast), CHR, POS, npoints = 20e3) +
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = 2)

set.seed(1)
ind.train <- sample(ind, 200e3)
ind.test <- setdiff(ind, ind.train)
system.time(
  prs.breast <- big_spLogReg(G, y[ind.train], ind.train = ind.train,
                             covar.train = PCs[ind.train, ],
                             alphas = c(1, 0.5, 0.1, 0.01, 0.001),
                             ncores = nb_cores())
) # 3h
str(prs.breast)
sum(rowMeans(sapply(prs.breast, function(x) x$beta.X)) > 0) # 252 predictors

preds <- predict(prs.breast, G)
table(y[ind.train])  # 8328 cases
AUC(preds[ind.test], y[ind.test])  # 60.7

df2 <- data.frame(y, pred = preds, age = df0$age, PC = I(PCs))
df3 <- df2[ind.train, ]
summary(myglm <- glm(y ~ pred + poly(age, 2) + PC, data = df3, family = "binomial"))
preds2 <- predict(myglm, df2[ind.test, ], type = "response")
AUCBoot(preds2, y[ind.test]) # 67.9 [66.5 - 69.4]
qplot(preds2, fill = as.factor(y[ind.test]), geom = "density", alpha = I(0.4))
# all < 20%  ->  because of unbalanced?

SNPs_breast_cancer <- readxl::read_excel("breast-cancer/SNPs-breast-cancer.xlsx")
plot(-log10(as.numeric(SNPs_breast_cancer$onco_icogs_gwas_P1df)), pch = 20)
map <- bigsnpr::snp_attach("UKB-merged_sub1.rds")$map
df3 <- tidyr::separate(SNPs_breast_cancer, phase_3_1kg_id, c("rsID", "pos", "A0", "A1"))
sum(!is.na(ind <- match(map$marker.ID, df3$rsID))) # 1505
plot(-predict(gwas.breast), -log10(as.numeric(SNPs_breast_cancer$onco_icogs_gwas_P1df[ind])), pch = 20)
beta_gwas <- (-1)^(map$allele1 != SNPs_breast_cancer$a0[ind]) * as.numeric(SNPs_breast_cancer$onco_icogs_gwas_beta[ind])
plot(gwas.breast$estim, beta_gwas, pch = 20); abline(0, 1, col = "red")


tmp <- bigreadr::fread2("breast-cancer/nature24284-s5.csv", dec = ",")
as.numeric(tmp$onco_icogs_gwas_P1df)


## BALANCED
# One run
ind.train.control <- which(y[ind.train] == 0)
ind.train2 <- sort(c(ind.train[y[ind.train] == 1],
                     ind.train[sample(ind.train.control, 20e3)]))
table(y[ind.train2])
system.time(
  prs.breast2 <- big_spLogReg(G, y[ind.train2], ind.train = ind.train2,
                              covar.train = PCs[ind.train2, ],
                              alphas = c(1, 0.5, 0.1),
                              ncores = nb_cores())
) # 1h
str(prs.breast2)
sum(rowMeans(sapply(prs.breast2, function(x) x$beta.X)) > 0) # 252 -> 401

preds <- predict(prs.breast2, G, covar.row = PCs)
AUC(preds[ind.test], y[ind.test])  # 59.4

df2 <- data.frame(y, pred = preds, age = df0$age, PC = I(PCs))
df3 <- df2[ind.train2, ]
summary(myglm <- glm(y ~ pred + poly(age, 2) + PC, data = df3, family = "binomial"))
preds2 <- predict(myglm, df2[ind.test, ], type = "response")
AUCBoot(preds2, y[ind.test]) # 67.0 [65.6 - 68.5]
qplot(preds2, fill = as.factor(y[ind.test]), geom = "density", alpha = I(0.4))
# all < 20%  ->  because of unbalanced?

pROC::ggroc(pROC::roc(y[ind.test], preds2)) +
  theme_bigstatsr()

# Many runs
set.seed(1)
control_sets <- split(ind.train.control, sample(rep_len(1:8, length(ind.train.control))))

summary(myglm <- glm(y ~ age + PC, data = df2[ind.train, ], family = "binomial"))
base <- predict(myglm, df2)
AUC(base[ind.test], y[ind.test]) # 65.0
qplot(1 / (1 + exp(-base[ind.test])), fill = as.factor(y[ind.test]),
      geom = "density", alpha = I(0.4)) +
  theme_bigstatsr() +
  labs(fill = "Case?")

library(doParallel)
test <- foreach(set = control_sets) %do% {

  ind.train2 <- sort(c(ind.train[y[ind.train] == 1], ind.train[set]))
  bigstatsr::big_spLogReg(G, y[ind.train2], ind.train = ind.train2,
                          covar.train = PCs[ind.train2, ],
                          base.train = base[ind.train2],
                          alphas = c(1, 0.5, 0.1),
                          ncores = nb_cores())
}

K_scores <- foreach(mod = test, .combine = "cbind") %do% {
  ind.col <- attr(mod, "ind.col")
  foreach(obj = mod, .combine = "cbind") %do% {
    beta.X <- obj$beta.X
    ind.nozero <- which(beta.X != 0)
    big_prodVec(G, beta.X[ind.nozero], ind.col = ind.col[ind.nozero]) +
      obj$intercept + drop(PCs %*% obj$beta.covar)
  }
}


aucs <- apply(K_scores, 2, function(pred) AUC(pred[ind.test], y[ind.test]))
summary(aucs)
hist(aucs)
glm <- glmnet::glmnet(K_scores[ind.train, ], y[ind.train], alpha = 0.0001,
                      lower.limits = 0, family = "binomial")
glm$beta
predf <- predict(glm, K_scores[ind.test, ])
max(aucsf <- apply(predf, 2, AUC, target = y[ind.test]))
plot(aucsf[-1], pch = 20)
# alpha = 0.01   -> 61.6 // 61.5
# alpha = 0.001  -> 61.8 // 61.5
# alpha = 0.0001 -> 61.8 // 61.5

qplot(1 / (1 + exp(-predf[, 40])), fill = as.factor(y[ind.test]),
      geom = "density", alpha = I(0.4)) +
  theme_bigstatsr() +
  labs(fill = "Case?")

glm <- glmnet::glmnet(cbind(K_scores, df0$age, PCs)[ind.train, ], y[ind.train],
                      alpha = 0.001, family = "binomial",
                      lower.limits = c(rep(0, 80), rep(-Inf, 11)),
                      penalty.factor = c(rep(1, 80), rep(0, 11)))

predf <- predict(glm, cbind(K_scores, df0$age, PCs)[ind.test, ])
max(aucsf <- apply(predf, 2, AUC, target = y[ind.test]))
plot(aucsf[-1], pch = 20)
# alpha = 0.01   -> 68.3
# alpha = 0.001  -> 68.3
# alpha = 0.0001 -> 68.2

qplot(1 / (1 + exp(-predf[, 40])), fill = as.factor(y[ind.test]),
      geom = "density", alpha = I(0.4)) +
  theme_bigstatsr() +
  labs(fill = "Case?")

AUCBoot(predf[, 40], y[ind.test])
