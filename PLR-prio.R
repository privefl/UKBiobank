#### BRCA ####
conf <- pmax(lpval, 1) * info
system.time(
  mod <- big_spLogReg(G, y.sub[ind.train], ind.train, pf.X = 1 / conf,
                      K = 10, alphas = 10^(-(0:4)), dfmax = 200e3, n.abort = 3,
                      ncores = NCORES)
) # 7.3H
plot(mod)
summary(mod)
# # A tibble: 5 x 6
#    alpha validation_loss intercept beta              nb_var message
#    <dbl>           <dbl>     <dbl> <list>             <int> <list>
# 1 0.0001           0.240     -3.26 <dbl [1,697,409]>  60764 <chr [10]>
# 2 0.001            0.241     -3.74 <dbl [1,697,409]>  15862 <chr [10]>
# 3 0.01             0.241     -3.92 <dbl [1,697,409]>   4780 <chr [10]>
# 4 0.1              0.241     -3.97 <dbl [1,697,409]>   3497 <chr [10]>
# 5 1                0.241     -3.98 <dbl [1,697,409]>   3300 <chr [10]>

pred <- predict(mod, G, ind.test)
AUCBoot(pred, y.sub[ind.test])                          # 65.3 [63.8-66.8]
AUCBoot(predict(mod[5], G, ind.test), y.sub[ind.test])  # 64.3 [62.8-65.8]

new_beta <- summary(mod, best.only = TRUE)$beta[[1]]
ind <- which(new_beta != 0)
plot(beta[ind], new_beta[ind], pch = 20); abline(0, 1, col = "red")

plot(1 / (1 + exp(-pred)), pred0, pch = 20, col = scales::alpha("black", 0.3),
     xlab = "Prediction with SCT", ylab = "Prediction with prioPLR")
abline(0, 1, col = "red")


library(ggplot2)
ggplot(data.frame(pred = pred, pheno = y.sub[ind.test])) +
  theme_bigstatsr() +
  geom_density(aes(pred, fill = as.factor(pheno)), alpha = 0.3)

#### T1D ####
conf <- pmax(lpval, 1) * info
system.time(
  mod <- big_spLogReg(G, y.sub[ind.train], ind.train, pf.X = 1 / conf,
                      K = 5, alphas = 10^(-(0:4)), ncores = NCORES)
) # 55 min
plot(mod)
summary(mod)
# # A tibble: 5 x 6
#    alpha validation_loss intercept beta              nb_var message
#    <dbl>           <dbl>     <dbl> <list>             <int> <list>
# 1 0.0001          0.0154     -7.99 <dbl [1,222,687]>   2019 <chr [5]>
# 2 0.001           0.0154     -7.93 <dbl [1,222,687]>    506 <chr [5]>
# 3 0.01            0.0154     -7.90 <dbl [1,222,687]>    291 <chr [5]>
# 4 0.1             0.0154     -7.79 <dbl [1,222,687]>    280 <chr [5]>
# 5 1               0.0154     -7.82 <dbl [1,222,687]>    280 <chr [5]>

pred <- predict(mod, G, ind.test)
AUCBoot(pred, y.sub[ind.test]) # 78.1 [75.1-81.1]
AUCBoot(predict(mod[5], G, ind.test), y.sub[ind.test]) # 77.7 [74.7-80.6]


#### T2D ####
conf <- pmax(lpval, 1) * info
system.time(
  mod <- big_spLogReg(G, y.sub[ind.train], ind.train, pf.X = 1 / conf,
                      K = 10, alphas = 10^(-(0:4)), ncores = NCORES)
) # 7.3H
plot(mod)
summary(mod)
# # A tibble: 5 x 6
#    alpha validation_loss intercept beta              nb_var message
#    <dbl>           <dbl>     <dbl> <list>             <int> <list>
# 1 0.0001           0.173      23.8 <dbl [1,339,724]> 103247 <chr [10]>
# 2 0.001            0.173      24.0 <dbl [1,339,724]>  14775 <chr [10]>
# 3 0.01             0.174      23.3 <dbl [1,339,724]>   4531 <chr [10]>
# 4 0.1              0.174      17.9 <dbl [1,339,724]>   3053 <chr [10]>
# 5 1                0.174      18.1 <dbl [1,339,724]>   2982 <chr [10]>

pred <- predict(mod, G, ind.test)
AUCBoot(pred, y.sub[ind.test])                          # 65.3 [64.4-66.2]
AUCBoot(predict(mod[5], G, ind.test), y.sub[ind.test])  # 63.4 [62.4-64.3]


#### PRCA ####
conf <- pmax(lpval, 1) * info
system.time(
  mod <- big_spLogReg(G, y.sub[ind.train], ind.train, pf.X = 1 / conf,
                      K = 10, alphas = 10^(-(0:4)), ncores = NCORES)
) # 2.9H
plot(mod)
summary(mod)
# # A tibble: 5 x 6
#    alpha validation_loss intercept beta              nb_var message
#    <dbl>           <dbl>     <dbl> <list>             <int> <list>
# 1 0.0001           0.172     -1.13 <dbl [1,309,368]>  52382 <chr [10]>
# 2 0.001            0.173     -1.69 <dbl [1,309,368]>  11483 <chr [10]>
# 3 0.01             0.173     -2.00 <dbl [1,309,368]>   4699 <chr [10]>
# 4 0.1              0.173     -2.14 <dbl [1,309,368]>   3648 <chr [10]>
# 5 1                0.173     -2.16 <dbl [1,309,368]>   3581 <chr [10]>

pred <- predict(mod, G, ind.test)
AUCBoot(pred, y.sub[ind.test])                          # 70.2 [68.7-71.6]
AUCBoot(predict(mod[5], G, ind.test), y.sub[ind.test])  # 68.7 [68.2-71.1]


#### Depression ####
conf <- pmax(lpval, 1) * info
system.time(
  mod <- big_spLogReg(G, y.sub[ind.train], ind.train, pf.X = 1 / conf,
                      K = 10, alphas = 10^(-(0:4)), dfmax = 200e3, n.abort = 3,
                      ncores = NCORES)
) # 7.3H
plot(mod)
summary(mod)
# # A tibble: 5 x 6
#    alpha validation_loss intercept beta              nb_var message
#    <dbl>           <dbl>     <dbl> <list>             <int> <list>
# 1 0.0001           0.278      5.31 <dbl [1,410,615]> 107416 <chr [10]>
# 2 0.001            0.279      5.29 <dbl [1,410,615]>  17547 <chr [10]>
# 3 0.01             0.279      1.82 <dbl [1,410,615]>   5266 <chr [10]>
# 4 0.1              0.279      1.30 <dbl [1,410,615]>   4015 <chr [10]>
# 5 1                0.279      1.31 <dbl [1,410,615]>   3930 <chr [10]>

pred <- predict(mod, G, ind.test)
AUCBoot(pred, y.sub[ind.test])                          # 56.2 [55.0-57.5]
AUCBoot(predict(mod[5], G, ind.test), y.sub[ind.test])  # 54.4 [53.1-55.6]

new_beta <- summary(mod, best.only = TRUE)$beta[[1]]
ind <- which(new_beta != 0)
plot(beta[ind], new_beta[ind], pch = 20); abline(0, 1, col = "red")
# some laaaaarge effects


#### CAD ####
conf <- pmax(lpval, 1) * info
system.time(
  mod <- big_spLogReg(G, y.sub[ind.train], ind.train, pf.X = 1 / conf,
                      K = 10, alphas = 10^(-(0:4)), dfmax = 200e3, n.abort = 3,
                      ncores = NCORES)
) # 7.3H
plot(mod)
summary(mod)
# # A tibble: 5 x 6
#    alpha validation_loss intercept beta              nb_var message
#    <dbl>           <dbl>     <dbl> <list>             <int> <list>
# 1 0.0001           0.398     -2.51 <dbl [1,013,109]> 307973 <chr [10]>
# 2 0.001            0.399     -2.24 <dbl [1,013,109]>  70832 <chr [10]>
# 3 0.01             0.400     -2.24 <dbl [1,013,109]>  21589 <chr [10]>
# 4 0.1              0.400     -2.33 <dbl [1,013,109]>  15609 <chr [10]>
# 5 1                0.400     -2.33 <dbl [1,013,109]>  15413 <chr [10]>

pred <- predict(mod, G, ind.test)
AUCBoot(pred, y.sub[ind.test])                          # 63.5 [62.8-64.1]
AUCBoot(predict(mod[5], G, ind.test), y.sub[ind.test])  # 62.6 [62.0-63.3]

new_beta <- summary(mod, best.only = TRUE)$beta[[1]]
ind <- which(new_beta != 0)
plot(beta[ind], new_beta[ind], pch = 20); abline(0, 1, col = "red")


#### Asthma ####
conf <- pmax(lpval, 1) * info
system.time(
  mod <- big_spLogReg(G, y.sub[ind.train], ind.train, pf.X = 1 / conf,
                      K = 10, alphas = 10^(-(0:4)), dfmax = 200e3, n.abort = 3,
                      ncores = NCORES)
) # 7.3H
plot(mod)
summary(mod)
# # A tibble: 5 x 6
#    alpha validation_loss intercept beta              nb_var message
#    <dbl>           <dbl>     <dbl> <list>             <int> <list>
# 1 0.0001           0.398     -2.51 <dbl [1,013,109]> 307973 <chr [10]>
# 2 0.001            0.399     -2.24 <dbl [1,013,109]>  70832 <chr [10]>
# 3 0.01             0.400     -2.24 <dbl [1,013,109]>  21589 <chr [10]>
# 4 0.1              0.400     -2.33 <dbl [1,013,109]>  15609 <chr [10]>
# 5 1                0.400     -2.33 <dbl [1,013,109]>  15413 <chr [10]>

pred <- predict(mod, G, ind.test)
AUCBoot(pred, y.sub[ind.test])                          # 63.5 [62.8-64.1]
AUCBoot(predict(mod[5], G, ind.test), y.sub[ind.test])  # 62.6 [62.0-63.3]

new_beta <- summary(mod, best.only = TRUE)$beta[[1]]
ind <- which(new_beta != 0)
plot(beta[ind], new_beta[ind], pch = 20); abline(0, 1, col = "red")
