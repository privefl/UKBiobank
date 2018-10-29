library(bigsnpr)
ukb <- snp_attach("UKB-merged_sub1.rds")
G <- ukb$genotypes
CHR <- ukb$map$chromosome
POS <- ukb$map$physical.pos
rm(ukb)

COV <- predict(readRDS("UKB_SVD_auto.rds"))

df0 <- readRDS("pheno.rds")
library(hexbin)
plot(hexbin(df0$date, df0$height))

ind.noNA <- with(df0, which(!is.na(height) & !is.na(sex) & !is.na(date) & pop == "British" & height > 140 & height < 200))

mod.date <- lm(height ~ date * sex, data = df0[ind.noNA, ])
summary(mod.date)
mod.date <- lm(height ~ date + sex, data = df0[ind.noNA, ])
summary(mod.date)
# 154.6 cm  +  1 cm every 6.4 year  +  13.3 cm for males

df.test <- df0[ind.test, ]
df.test$pred <- predict(mod.date, df.test)
ggplot(df.test, aes(pred, height, color = sex)) +
  geom_point(alpha = 0.4) +
  # geom_smooth(linetype = 2, method = "lm") +
  theme_bigstatsr() +
  labs(color = "Sex") +
  scale_color_viridis_d() +
  geom_abline(color = "red")

# PCs <- setNames(as.data.frame(COV), paste0("PC", 1:10))
mod.PC <- lm(height ~ date + sex + COV, data = df0,
             subset = height > 140 & height < 200)
summary(mod.PC)  ## all PCs (except #4) are significant
base <- predict(mod.PC, df0)
stopifnot(length(base) == nrow(G))
hist(base - df0$height)

# system.time(
#   gwas_sex <- big_univLogReg(G, as.integer(sex[ind.noNA]) - 1, ind.noNA,
#                              ncores = nb_cores())
# ) # 87 min
# library(ggplot2)
# plot(gwas_sex)
# snp_qq(gwas_sex) + xlim(2, NA)
# snp_manhattan(snp_gc(gwas_sex), CHR, POS, npoints = 20e3) +
#   geom_hline(yintercept = -log10(5e-8), color = "red", linetype = 2)

# ind.test <- sample(ind.noNA, 20e3)
# ind.train <- setdiff(ind.noNA, ind.test)

COV2 <- cbind(COV, df0$date, as.integer(df0$sex))
system.time(
  prs <- big_spLinReg(G, df0$height[ind.train],
                      ind.train = ind.train,
                      # covar.train = COV2[ind.train, ],
                      base.train = base[ind.train],
                      alphas = c(1, 0.1, 0.01),
                      ncores = nb_cores())
) # 21h -> 17.5h without covars
str(prs)  ## alpha = 1

preds.test <- base[ind.test] + predict(prs, G, ind.test)

library(ggplot2)

ggplot(data.frame(pred = preds.test, true = df0$height[ind.test], sex = df0$sex[ind.test]),
       aes(pred, true, color = sex)) +
  geom_point(alpha = 0.4) +
  geom_smooth(linetype = 2, method = "lm") +
  theme_bigstatsr() +
  labs(color = "Sex") +
  scale_color_viridis_d() +
  geom_abline(color = "red")


preds.train <- base[ind.train] +
  predict(prs, G, ind.train)

mylm <- lm(height ~ pred * sex + date + COV,
           data = cbind(df0[ind.train, ], pred = preds.train, COV = I(COV[ind.train, ])),
           subset = height > 140 & height < 200)
summary(mylm)

preds.test.adj <- predict(mylm, cbind(df0[ind.test, ], pred = preds.test,
                                      COV = I(COV[ind.test, ])))

ggplot(data.frame(pred = preds.test.adj, height = df0$height[ind.test], sex = df0$sex[ind.test]),
       aes(pred, height, color = sex)) +
  geom_point(alpha = 0.4) +
  geom_smooth(linetype = 2, method = "lm") +
  theme_bigstatsr() +
  labs(color = "Sex") +
  scale_color_viridis_d() +
  geom_abline(color = "red")

library(dplyr)
data.frame(pred = preds.test.adj, true = df0$height[ind.test], sex = df0$sex[ind.test]) %>%
  group_by(sex) %>%
  summarize(cor = cor(pred, true))

hist(preds.test.adj - df0$height[ind.test], breaks = 20)
# 0.844 (Female: 0.655 // Male: 0.65)

resid <- rbind(
  data.frame(eps = df0$height[ind.test] - preds.test.adj, model = "final"),
  data.frame(eps = df.test$height - df.test$pred, model = "baseline")
)
ggplot(resid) +
  geom_density(aes(x = eps, fill = model), alpha = 0.4) +
  bigstatsr::theme_bigstatsr() +
  labs(x = "residuals")

resid %>%
  group_by(model) %>%
  summarize(MAE = mean(abs(eps)))
  summarise(resids = list(quantile(eps, c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99)))) %>%
  pull(resids)

beta.X <- rowMeans(sapply(prs, function(x) x$beta.X))
beta.X[beta.X == 0] <- NA
plot(beta.X, pch = 20)
hist(beta.X, breaks = 50) # Laplace
mean(!is.na(beta.X)) # 18% (-> 104K SNPs)
rowMeans(sapply(prs, function(x) x$beta.covar))


ind.other <- with(df0, which(!is.na(height) & !is.na(sex) & !is.na(date) & pop != "British" & height > 140 & height < 200))
preds.other <- base[ind.other] + predict(prs, G, ind.other, covar.row = COV2[ind.other, ])
preds.other.adj <- predict(mylm, cbind(df0[ind.other, ], pred = preds.other,
                                      COV = I(COV[ind.other, ])))

ggplot(data.frame(pred = preds.other.adj, true = df0$height[ind.other], sex = df0$sex[ind.other]),
       aes(pred, true, color = sex)) +
  geom_point(alpha = 0.4) +
  geom_smooth(linetype = 2, method = "lm") +
  theme_bigstatsr() +
  labs(color = "Sex") +
  scale_color_viridis_d() +
  geom_abline(color = "red")



RMSE <- function(pred, true) {
  sqrt(mean((pred - true)^2))
}
MAE <- function(pred, true) {
  mean(abs(pred - true))
}

cbind(pred = preds.other.adj, df0[ind.other, ]) %>%
  group_by(sex, pop) %>%
  summarise(RMSE = RMSE(pred, height), MAE = MAE(pred, height)) %>%
  arrange(sex, desc(RMSE)) %>%
  print(n = Inf)

cbind(pred = preds.test.adj, df0[ind.test, ]) %>%
  group_by(sex, pop) %>%
  summarise(RMSE = RMSE(pred, height), MAE = MAE(pred, height)) %>%
  arrange(sex, desc(RMSE)) %>%
  print(n = Inf)
