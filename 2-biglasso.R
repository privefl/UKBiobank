library(bigsnpr)
ukb <- snp_attach("UKB-merged_sub1.rds")
G <- ukb$genotypes
CHR <- ukb$map$chromosome
POS <- ukb$map$physical.pos
NCORES <- nb_cores()

COV <- readRDS("UKB_SVD_auto.rds")$u

df0 <- bigreadr::fread2("ukb22544.csv")
sex <- factor(df0[match(ukb$fam$sample.ID, df0$eid), "22001-0.0"],
              levels = c(0, 1), labels = c("Female", "Male"))

lvl.1001 <- c(-3,-1,1,2,3,4,5,6,1001,1002,1003,2001,2002,2003,2004,3001,3002,3003,3004,4001,4002,4003)
lbl.1001 <- c("Prefer not to answer","Do not know","White","Mixed","Asian or Asian British","Black or Black British","Chinese","Other ethnic group","British","Irish","Any other white background","White and Black Caribbean","White and Black African","White and Asian","Any other mixed background","Indian","Pakistani","Bangladeshi","Any other Asian background","Caribbean","African","Any other Black background")
pop0 <- ordered(df0[match(ukb$fam$sample.ID, df0$eid), "21000-0.0"],
                levels = lvl.1001, labels = lbl.1001)
height <- df0[match(ukb$fam$sample.ID, df0$eid), "50-0.0"]
ind.noNA <- which(!is.na(height) & !is.na(sex) & pop0 == "British" & height > 130)
rm(df0, ukb)  ## make session load faster

# New dataset
set.seed(1)
ind <- sort(sample(ind.noNA, size = 25e3))
# G2 <- big_copy(G, ind.row = ind, type = "double", backingfile = "ukb-biglasso")$save()
# G2$add_columns(11)
# G2[, tail(cols_along(G2), 11)] <- cbind(COV, sex)[ind, ]
G2 <- big_attach("ukb-biglasso.rds")
stopifnot(G2$is_saved)
stopifnot(all.equal(dim(G2), c(length(ind), ncol(G) + 11)))

set.seed(1)
ind.train <- sort(sample(25e3, 20e3))
ind.test <- setdiff(seq_len(25e3), ind.train)


system.time(
  test <- biglasso::biglasso(G2$bm(), height[ind], row.idx = ind.train,
                             ncores = nb_cores())
)
#     user   system  elapsed
# 3042.715   57.977  677.233

plot(AIC(test), pch = 20)
points(BIC(test), pch = 20, col = "red")
which.min(BIC(test))


str(test)
plot(test$loss, pch = 20)
preds <- predict(test, G2$bm(), row.idx = ind.test)


library(ggplot2)
qplot(preds[, 100], height[ind[ind.test]], color = sex[ind[ind.test]]) +
  theme_bigstatsr() +
  geom_abline(slope = 1, color = "red") +
  # coord_equal() +
  labs(y = "True height", x = "Predicted height", color = "Sex") +
  scale_color_viridis_d()

library(dplyr)
data_frame(pred = preds[, 100], true = height[ind[ind.test]], sex = sex[ind[ind.test]]) %>%
  group_by(sex) %>%
  summarize(cor = cor(pred, true))
#   sex      cor
#   <fct>  <dbl>
# 1 Female 0.124
# 2 Male   0.143

library(Matrix)
nb.pred <- colSums(test$beta != 0)
plot(nb.pred, pch = 20)
beta <- test$beta[, 100]
summary(beta[beta != 0])
head(sort(abs(beta), decreasing = TRUE), 10)



set.seed(1)
ind.train <- sort(sample(25e3, 20e3))
ind.test <- setdiff(seq_len(25e3), ind.train)

mod <- lm(height ~ sex, data = data.frame(sex, height), subset = ind[ind.train])
base <- predict(mod, data.frame(sex = sex[ind]))
COV.ind <- cbind(COV, sex)[ind, ]
# debugonce(big_spLinReg)
system.time(
  prs <- big_spLinReg(G, height[ind[ind.train]],
                      ind.train = ind[ind.train],
                      covar.train = COV.ind[ind.train, ],
                      base.train = base[ind.train],
                      alphas = 1,
                      ncores = 10)
)

pred <- base + predict(prs, G, ind.row = ind, covar.row = COV.ind)

data_frame(pred = pred[ind.test], true = height[ind[ind.test]], sex = sex[ind[ind.test]]) %>%
  group_by(sex) %>%
  summarize(cor = cor(pred, true))

qplot(pred[ind.test], height[ind[ind.test]], color = sex[ind[ind.test]]) +
  theme_bigstatsr() +
  geom_abline(slope = 1, color = "red") +
  # coord_equal() +
  labs(y = "True height", x = "Predicted height", color = "Sex") +
  scale_color_viridis_d()

adj <- lm(y ~ pred * sex, data = data.frame(y = height[ind[ind.train]], pred = pred[ind.train], sex = sex[ind[ind.train]]))
summary(adj)
pred.adj <- predict(adj, data.frame(pred = pred[ind.test], sex = sex[ind[ind.test]]))


qplot(pred.adj, height[ind[ind.test]], color = sex[ind[ind.test]], alpha = I(0.5)) +
  theme_bigstatsr() +
  geom_abline(slope = 1, color = "red") +
  # coord_equal() +
  labs(y = "True height", x = "Predicted height", color = "Sex") +
  scale_color_viridis_d()


## Prediction on African population
levels(pop0)
ind.Af <- which(pop0 == "African" & !is.na(sex) & !is.na(height))

# G.Af <- big_copy(G, ind.row = ind.Af, backingfile = "ukb-af")$save()
G.Af <- big_attach("ukb-af.rds")
set.seed(1)
ind.Af.train <- sort(sample(length(ind.Af), 2000))
ind.Af.test <- setdiff(seq_along(ind.Af), ind.Af.train)

COV.Af <- cbind(COV, sex)[ind.Af, ]
pred.Af <- predict(mod, data.frame(sex = sex[ind.Af])) +
  predict(prs, G.Af, covar.row = COV.Af)

data_frame(pred = pred.Af, true = height[ind.Af], sex = sex[ind.Af]) %>%
  group_by(sex) %>%
  summarize(cor = cor(pred, true))
# 0.13 and 0.10

mod.Af1 <- big_spLinReg(G.Af, height[ind.Af[ind.Af.train]],
                        ind.train = ind.Af.train,
                        covar.train = COV.Af[ind.Af.train, ],
                        base.train = pred.Af[ind.Af.train],
                        alphas = c(1, 0.1, 0.001),
                        ncores = nb_cores())
str(mod.Af1)

pred.Af1 <- pred.Af[ind.Af.test] +
  predict(mod.Af1, G.Af, ind.row = ind.Af.test, covar.row = COV.Af[ind.Af.test, ])
plot(pred.Af1, height[ind.Af[ind.Af.test]], pch = 20); abline(0, 1, col = "red")

data_frame(pred = pred.Af1, true = height[ind.Af[ind.Af.test]],
           sex = sex[ind.Af[ind.Af.test]]) %>%
  group_by(sex) %>%
  summarize(cor = cor(pred, true))
# 0.06 & 0.15



pred.Af.adj <- predict(adj, data.frame(pred = pred.Af, sex = sex[ind.Af]))
qplot(pred.Af.adj, height[ind.Af], color = sex[ind.Af]) +
  theme_bigstatsr() +
  geom_abline(slope = 1, color = "red") +
  # coord_equal() +
  labs(y = "True height", x = "Predicted height", color = "Sex") +
  scale_color_viridis_d()

data_frame(pred = pred.Af.adj, true = height[ind.Af], sex = sex[ind.Af]) %>%
  group_by(sex) %>%
  summarize(cor = cor(pred, true))

preds2 <- bind_rows(
  data_frame(pred = pred.Af.adj, true = height[ind.Af],
             sex = sex[ind.Af], pop = "African"),
  data_frame(pred = pred.adj, true = height[ind[ind.test]],
             sex = sex[ind[ind.test]], pop = "British")
)

ggplot(preds2) +
  geom_density(aes(pred), alpha = 0.4, fill = "blue") +
  geom_density(aes(true), alpha = 0.4, fill = "red") +
  facet_grid(pop ~ sex) +
  theme_bigstatsr()

gwas <- big_univLinReg(as_FBM(COV[c(ind.Af, ind[ind.test]), ]),
               y.train = preds2$pred - preds2$true)
plot(gwas, type = "Manhattan")
