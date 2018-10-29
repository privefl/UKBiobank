library(bigsnpr)
ukb <- snp_attach("UKB-merged_sub1.rds")
G <- ukb$genotypes
CHR <- ukb$map$chromosome
POS <- ukb$map$physical.pos
dim(G) # 449926 x 564148
# big_counts(G)[4, ] # OK
(NCORES <- nb_cores())

ind.rm <- snp_indLRLDR(CHR, POS)
system.time(
  ind.keep <- snp_clumping(G, CHR, thr.r2 = 0.1, exclude = ind.rm,
                           ncores = NCORES)
) # 195K
# utilisateur     système      écoulé
#   27488.856     220.668   29665.410
#    user   system  elapsed
# 402.985  203.360 6352.371
saveRDS(ind.keep, "keep_pruning.rds")

system.time(
  obj.svd <- big_randomSVD(G, snp_scaleBinom(), ind.col = ind.keep,
                           verbose = TRUE, ncores = NCORES)
)
#  user   system  elapsed
# 5.275   10.437 1700.229

obj.svd2 <- snp_autoSVD(G, CHR, POS, k = 20, thr.r2 = 0.1, ncores = NCORES)
# saveRDS(obj.svd2, "UKB_SVD_auto.rds")
p <- plot(obj.svd2, type = "scores", scores = 1:2)
# system.time({png("test.png"); print(p); dev.off()})  # 36 sec
# system.time(ggplot2::ggsave("test.png", p))  # 70 sec
# system.time(ggplot2::ggsave("test.svg", p))  # 19 sec
## no overlap :O (expect for chr 8)
attr(obj.svd2, "lrldr")
subset(LD.wiki34, Chr %in% c(1, 3, 6, 8))


# saveRDS(obj.svd, "UKB_SVD.rds")
obj.svd <- readRDS("UKB_SVD.rds")
plot(obj.svd)  #  K = 5
# Compare with https://biobank.ctsu.ox.ac.uk/crystal/docs/genotyping_qc.pdf (p11)
plot(obj.svd, type = "scores")
plot(obj.svd, type = "scores", scores = 3:4)
plot(obj.svd, type = "scores", scores = 5:6)
plot(obj.svd, type = "loadings", loadings = 1:10, coeff = 0.5)
ggplot2::ggsave("tmp.png")

# loadings <- bigreadr::fread2("snp_pca_map.txt")
# match(loadings$V2, ukb$map$marker.ID)  ## not all in common

obj.pcadapt <- snp_pcadapt(G, obj.svd$u[, 1:2], ncores = NCORES)
saveRDS(obj.pcadapt, "UKB_pcadapt.rds")
plot(obj.pcadapt)
ggplot2::ggsave("tmp.png")
snp_qq(obj.pcadapt) + xlim(2, NA)
ggplot2::ggsave("tmp.png")
obj.pcadapt.gc <- snp_gc(obj.pcadapt)
plot(obj.pcadapt.gc)
snp_manhattan(obj.pcadapt.gc, CHR, POS, npoints = 20e3) +
  ggplot2::geom_hline(yintercept = -log10(5e-8), color = "red")
ggplot2::ggsave("tmp.png")


## Need ethnicity (Data-Field 21000)
# First dataset?
system("./ukbmd5 ukb22544.enc")
system("./ukbunpack ukb22544.enc key-22544.txt")
system("./ukbconv ukb22544.enc_ukb csv")
# system("./ukbconv ukb22544.enc_ukb R") ## to get categories
# Second dataset?
system("./ukbmd5 ukb23144.enc")
system("./ukbunpack ukb23144.enc key-23144.txt")
system("./ukbconv ukb23144.enc_ukb csv")
df0 <- bigreadr::fread2("ukb22544.csv")
df <- bigreadr::fread2("ukb23144.csv")

dim(df0)
dim(df)
c(2977063, 3158119, 4956959) %in% ukb$fam$sample.ID
which(ukb$fam$sample.ID <= 0)

names(df) %in% names(df0)
tmp_21000 <- df0[grep("21000", names(df0))]
head(tmp_21000)
table(tmp_21000$`21000-0.0`)

lvl.1001 <- c(-3,-1,1,2,3,4,5,6,1001,1002,1003,2001,2002,2003,2004,3001,3002,3003,3004,4001,4002,4003)
lbl.1001 <- c("Prefer not to answer","Do not know","White","Mixed","Asian or Asian British","Black or Black British","Chinese","Other ethnic group","British","Irish","Any other white background","White and Black Caribbean","White and Black African","White and Asian","Any other mixed background","Indian","Pakistani","Bangladeshi","Any other Asian background","Caribbean","African","Any other Black background")
pop0 <- ordered(df0[match(ukb$fam$sample.ID, df0$eid), "21000-0.0"],
                levels = lvl.1001, labels = lbl.1001)

lbl.1001.2 <- c(NA, NA,"White",NA,"Asian","Black","Asian",NA,"White","White","White","White and Black","White and Black","White and Asian",NA,"Asian","Asian","Asian","Asian","Black","Black","Black")
pop <- ordered(df0[match(ukb$fam$sample.ID, df0$eid), "21000-0.0"],
               levels = lvl.1001, labels = lbl.1001.2)

library(ggplot2)
plot(obj.svd, type = "scores", scores = 3:4) +
  aes(color = pop) +
  labs(color = "Ancestry") +
  theme(legend.position = c(0.8, 0.8))

ggsave("PC-3-4.png", width = 8, height = 8)

ind.BW <- which(pop %in% c("Black", "White", "White and Black"))
system.time(
  obj.svd3 <- snp_autoSVD(G, CHR, POS, thr.r2 = 0.1, ncores = NCORES,
                          ind.row = ind.BW)
) # 5h

plot(obj.svd3)
plot(obj.svd3, type = "loadings", loadings = 1:10, coeff = 0.5)
p <- plot(obj.svd3, type = "scores") +
  aes(color = pop[ind.BW]) +
  labs(color = "Ancestry") +
  theme(legend.position = c(0.8, 0.8))

ggsave("PC-BW-1-2.png", p, width = 8, height = 8)

p <- plot(obj.svd3, type = "scores", scores = 3:4) +
  aes(color = pop0[ind.BW]) +
  labs(color = "Ancestry") +
  theme(legend.position = c(0.8, 0.2))

ggsave("PC-BW-3-4.png", p, width = 8, height = 8)


obj.pcadapt3 <- snp_pcadapt(G, obj.svd3$u[, 1:3], ind.row = ind.BW, ncores = NCORES)
plot(obj.pcadapt3)
snp_qq(obj.pcadapt3) + xlim(2, NA)
obj.pcadapt3.gc <- snp_gc(obj.pcadapt3)
plot(obj.pcadapt3.gc)
snp_manhattan(obj.pcadapt3.gc, CHR, POS, npoints = 20e3) +
  ggplot2::geom_hline(yintercept = -log10(5e-8), color = "red")



ind.AW <- which(pop %in% c("Asian", "White", "White and Asian"))
system.time(
  obj.svd4 <- snp_autoSVD(G, CHR, POS, thr.r2 = 0.1, ncores = NCORES,
                          ind.row = ind.AW)
) # 7h

lbl.1001.3 <- c(NA, NA,"White",NA,"Asian","Black","Asian",NA,"White","White","White","White and Black","White and Black","White and Asian",NA,"Indian","Indian","Indian","Asian (other)","Black","Black","Black")
pop3 <- ordered(df0[match(ukb$fam$sample.ID, df0$eid), "21000-0.0"],
                levels = lvl.1001, labels = lbl.1001.3)

table(pop3)

plot(obj.svd4)
plot(obj.svd4, type = "loadings", loadings = 1:10, coeff = 0.5)
p <- plot(obj.svd4, type = "scores") +
  aes(color = pop3[ind.AW]) +
  labs(color = "Ancestry") +
  theme(legend.position = c(0.2, 0.2))

ggsave("PC-AW-1-2.png", p, width = 8, height = 8)

p <- plot(obj.svd4, type = "scores", scores = 3:4) +
  aes(color = pop0[ind.AW]) +
  labs(color = "Ancestry") +
  theme(legend.position = c(0.8, 0.8))

ggsave("PC-AW-3-4.png", p, width = 8, height = 8)


obj.pcadapt4 <- snp_pcadapt(G, obj.svd4$u[, 1:3], ind.row = ind.AW, ncores = NCORES)
plot(obj.pcadapt4)
snp_qq(obj.pcadapt4) + xlim(2, NA)
obj.pcadapt4.gc <- snp_gc(obj.pcadapt4)
plot(obj.pcadapt4.gc)
snp_manhattan(obj.pcadapt4.gc, CHR, POS, npoints = 20e3) +
  ggplot2::geom_hline(yintercept = -log10(5e-8), color = "red")


ind.A <- which(pop %in% c("Asian", "Any other Asian background", "White and Asian", "Indian", "Pakistani", "Bangladeshi"))
system.time(
  obj.svd5 <- snp_autoSVD(G, CHR, POS, thr.r2 = 0.1, ncores = NCORES,
                          ind.row = ind.A)
) # 7h

plot(obj.svd5)
plot(obj.svd5, type = "loadings", loadings = 1:10, coeff = 0.5)
p <- plot(obj.svd5, type = "scores") +
  aes(color = pop0[ind.A]) +
  labs(color = "Ancestry") +
  theme(legend.position = c(0.75, 0.2))

ggsave("PC-A-1-2.png", p, width = 8, height = 8)

p <- plot(obj.svd5, type = "scores", scores = 3:4) +
  aes(color = pop0[ind.A]) +
  labs(color = "Ancestry") +
  theme(legend.position = c(0.25, 0.8))

ggsave("PC-A-3-4.png", p, width = 8, height = 8)


ind.B <- which(pop0 %in% c("Black or Black British","White and Black Caribbean","White and Black African","Caribbean","African","Any other Black background"))
system.time(
  obj.svd6 <- snp_autoSVD(G, CHR, POS, thr.r2 = 0.1, ncores = NCORES,
                          ind.row = ind.B)
) # 2h

plot(obj.svd6)
plot(obj.svd6, type = "loadings", loadings = 1:10, coeff = 0.5)
p <- plot(obj.svd6, type = "scores") +
  aes(color = pop0[ind.B]) +
  labs(color = "Ancestry") +
  theme(legend.position = c(0.22, 0.85))

ggsave("PC-B-1-2.png", p, width = 8, height = 8)

p <- plot(obj.svd6, type = "scores", scores = 3:4) +
  aes(color = pop0[ind.B]) +
  labs(color = "Ancestry") +
  theme(legend.position = c(0.22, 0.15))

ggsave("PC-B-3-4.png", p, width = 8, height = 8)

height <- df0[match(ukb$fam$sample.ID, df0$eid), "50-0.0"]
ind.noNA <- which(!is.na(height))
system.time(
  gwas <- big_univLinReg(G, height[ind.noNA], ind.train = ind.noNA,
                         covar.train = obj.svd$u[ind.noNA, ], ncores = NCORES)
) # 92 min

snp_qq(gwas) + xlim(2, NA)
gwas.gc <- snp_gc(gwas)
snp_qq(gwas.gc) + xlim(2, NA)
plot(gwas.gc)
snp_manhattan(gwas.gc, CHR, POS, npoints = 20e3) +
  ggplot2::geom_hline(yintercept = -log10(5e-8), color = "red")

ind.noNA2 <- intersect(ind.noNA, which(pop0 == "British"))

system.time(
  gwas2 <- big_univLinReg(G, height[ind.noNA2], ind.train = ind.noNA2,
                         covar.train = obj.svd$u[ind.noNA2, ], ncores = NCORES)
) # 104 min

snp_qq(gwas2) + xlim(2, NA)
summary(abs(gwas2$estim))
gwas2.gc <- snp_gc(gwas2)
snp_qq(gwas2.gc) + xlim(2, NA)
plot(gwas2.gc)
snp_manhattan(gwas2.gc, CHR, POS, npoints = 20e3) +
  ggplot2::geom_hline(yintercept = -log10(5e-8), color = "red")

system.time(
  gwas3 <- big_univLinReg(G$copy(as.numeric(G$code256 > 1.5)),
                          height[ind.noNA2], ind.train = ind.noNA2,
                          covar.train = obj.svd$u[ind.noNA2, ], ncores = NCORES)
) # 104 min

snp_qq(gwas3) + xlim(2, NA)
summary(abs(gwas3$estim))
gwas3.gc <- snp_gc(gwas3)
snp_qq(gwas3.gc) + xlim(2, NA)
plot(gwas3.gc)
snp_manhattan(gwas3.gc, CHR, POS, npoints = 20e3) +
  ggplot2::geom_hline(yintercept = -log10(5e-8), color = "red")

qplot(-predict(gwas3.gc), -predict(gwas2.gc)) +
  theme_bigstatsr() +
  geom_abline(slope = 1, color = "red") +
  xlim(2, NA)


sex <- factor(df0[match(ukb$fam$sample.ID, df0$eid), "22001-0.0"],
              levels = c(0, 1), labels = c("Female", "Male"))
ind.noNA3 <- intersect(ind.noNA2, which(!is.na(sex) & height > 130))

ind.train <- sort(sample(ind.noNA3, size = 350e3))
ind.test <- setdiff(ind.noNA3, ind.train)
COV <- cbind(obj.svd$u)
base <- predict(lm(height ~ sex, subset = ind.train), sex)
system.time(
  prs <- big_spLinReg(G, height[ind.train], ind.train = ind.train,
                      covar.train = COV[ind.train, ],
                      base.train = base[ind.train],
                      alphas = c(1, 0.05, 0.001),
                      ncores = NCORES)
) # 21 h for alpha = 1 // 30h with grid-search
attr(prs, "alpha")  # 1
sum(rowSums(sapply(prs, function(x) x$beta.X) > 0) > 0) # 47013
sapply(prs, function(x) x$beta.covar)
qplot(COV[ind.test, 5], COV[ind.test, 10], col = sex[ind.test]) +
  theme_bigstatsr() +
  scale_color_viridis_d()
# no gradient of height...
# no gradient of sex...

pred <- base[ind.test] +
  predict(prs, G, ind.row = ind.test, covar.row = COV[ind.test, ])

qplot(pred, height[ind.test], color = COV[ind.test, 10]) +
  theme_bigstatsr() +
  geom_abline(slope = 1, color = "red") +
  coord_equal() +
  labs(y = "True height", x = "Predicted height", color = "Sex") +
  scale_color_viridis_c()

cor(pred, height[ind.test])
library(dplyr)
data_frame(pred = pred, true = height[ind.test], sex = sex[ind.test]) %>%
  group_by(sex) %>%
  summarize(cor = cor(pred, true))
# # A tibble: 2 x 2
#   sex      cor
#   <fct>  <dbl>
# 1 Female 0.617
# 2 Male   0.618

ind.trainF <- intersect(ind.train, which(sex == "Female"))
system.time(
  prsF <- big_spLinReg(G, height[ind.trainF], ind.train = ind.trainF,
                      covar.train = COV[ind.trainF, ],
                      alphas = c(1, 0.05, 0.001),
                      ncores = NCORES)
) # 19h
str(prsF)
sapply(prsF, function(x) x$beta.covar)

ind.testF <- intersect(ind.test, which(sex == "Female"))
predF <- predict(prsF, G, ind.row = ind.testF, covar.row = COV[ind.testF, ])

predF.train <- predict(prsF, G, ind.row = ind.trainF,
                       covar.row = COV[ind.trainF, ])
adj <- lm(y ~ pred, data = data.frame(y = height[ind.trainF], pred = predF.train))
summary(adj)
predF.adj <- predict(adj, data.frame(pred = predF))
qplot(predF.adj, height[ind.testF], color = COV[ind.testF, 5]) +
  theme_bigstatsr() +
  geom_abline(slope = 1, color = "red") +
  coord_equal() +
  labs(y = "True height", x = "Predicted height", color = "") +
  scale_color_viridis_c()

cor(predF.adj, height[ind.testF])  # 58.6

hist(predF.adj - height[ind.testF], breaks = 20) # +|- 10


ind <- sample(ind.trainF, 20e3)
qplot(predict(prsF, G, ind.row = ind, covar.row = COV[ind, ]),
      height[ind], color = COV[ind, 5]) +
  theme_bigstatsr() +
  geom_abline(slope = 1, color = "red") +
  coord_equal() +
  labs(y = "True height", x = "Predicted height", color = "") +
  scale_color_viridis_c()
