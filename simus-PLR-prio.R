library(bigsnpr)
# snp_readBed("simus/data_train.bed")
train <- snp_attach("simus/data_train.rds")
G <- train$genotypes
y <- train$fam$affection

sumstats <- bigreadr::fread2("simus/sumstats.txt")
CHR <- sumstats$chromosome
POS <- sumstats$physical.pos
LPVAL <- -log10(sumstats$pval)
pf <- 1 / pmax(LPVAL, 1)
hist(pf)

system.time(
  mod <- big_spLogReg(G, y, alphas = 10^(-(0:4)), ncores = nb_cores(),
                      dfmax = Inf, nlam.min = 30, n.abort = 3)
) # 10 min
plot(mod)
summary(mod)
# # A tibble: 5 x 6
#    alpha validation_loss intercept beta            nb_var message
#    <dbl>           <dbl>     <dbl> <list>           <int> <list>
# 1 0.0001           0.457     -1.97 <dbl [686,066]> 227815 <chr [10]>
# 2 0.001            0.456     -1.40 <dbl [686,066]>  71018 <chr [10]>
# 3 0.01             0.463     -1.41 <dbl [686,066]>  19455 <chr [10]>
# 4 0.1              0.466     -1.65 <dbl [686,066]>   6925 <chr [10]>
# 5 1                0.466     -1.87 <dbl [686,066]>   4893 <chr [10]>

# snp_readBed("simus/data_test.bed")
test <- snp_attach("simus/data_test.rds")
pred.test <- predict(mod, test$genotypes)
AUCBoot(pred.test, test$fam$affection)  # 71.3 [68.5-74.3]


system.time(
  mod2 <- big_spLogReg(G, y, alphas = 10^(-(0:4)), ncores = nb_cores(),
                       pf.X = pf, dfmax = Inf, nlam.min = 30, n.abort = 3)
) # 10 min
plot(mod2)
summary(mod2)
# # A tibble: 5 x 6
#    alpha validation_loss intercept beta            nb_var message
#    <dbl>           <dbl>     <dbl> <list>           <int> <list>
# 1 0.0001           0.419   -0.0354 <dbl [686,066]> 112615 <chr [10]>
# 2 0.001            0.417    0.210  <dbl [686,066]>  31689 <chr [10]>
# 3 0.01             0.426    0.141  <dbl [686,066]>   9420 <chr [10]>
# 4 0.1              0.428   -0.172  <dbl [686,066]>   3616 <chr [10]>
# 5 1                0.429   -0.341  <dbl [686,066]>   2264 <chr [10]>

# snp_readBed("simus/data_test.bed")
test <- snp_attach("simus/data_test.rds")
pred2.test <- predict(mod2, test$genotypes)
AUCBoot(pred2.test, test$fam$affection)  # 77.1 [74.5-79.7]
