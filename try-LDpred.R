library(bigsnpr)
simu <- snp_attach("data/ukbb4simu.rds")
G <- simu$genotypes
CHR <- as.integer(simu$map$chromosome)
POS <- simu$map$physical.pos
INFO <- readRDS("data/ukbb4simu_info.rds")$info
SD <- readRDS("data/ukbb4simu_stats.rds")$scale
load("data/ukbb4simu_ind.RData")
NCORES <- nb_cores()

# train <- snp_attach("data/ukbb4simu_train.rds")
# G.train <- train$genotypes
# # G.test  <- snp_attach("data/ukbb4simu_test.rds")$genotypes
#
# beta <- runif(ncol(G.train))
# pval <- runif(ncol(G.train))

library(dplyr)
tmp <- tempfile(tmpdir = "res_simu")
file_sumstats <- paste0(tmp, ".txt")
mutate(train$map, beta = beta_gwas, pval = 10^-lpval,
       chromosome = as.integer(chromosome)) %>%
  bigreadr::fwrite2(file_sumstats, sep = "\t")
file_hdf5 <- paste0(tmp, ".hdf5")

reticulate::use_python("/opt/rh/rh-python36/root/usr/bin/python")
reticulate::py_config()
# system("python3 --version", intern = TRUE)
ldpred <- "../ldpred/LDpred.py"
# system(glue::glue("python3 {ldpred} coord --help"))
system(glue::glue(
  "python3 {ldpred} coord",
  " --gf data/ukbb4simu_train",
  " --ssf {file_sumstats}",
  " --maf 0 --skip-coordination",
  " --rs marker.ID --A1 allele2 --A2 allele1 --pos physical.pos --chr chromosome",
  " --pval pval --eff beta --beta",
  " --N 315609",
  " --out {file_hdf5}"
))

system.time(
  system(glue::glue(
    "python3 {ldpred} gibbs",
    " --cf {file_hdf5}",
    " --ldr 333",
    " --ldf res_simu/LDpred",
    " --h2 {h2}",
    " --N 315609",
    " --out {tmp}"
  ))
) # 3.2H / 1.7H

beta_ldpred <-
  list.files("res_simu", pattern = paste0("^", basename(tmp), "_LDpred.*\\.txt$"),
             full.names = TRUE) %>%
  sapply(function(file) {
    res_ldpred <- bigreadr::fread2(file, select = c(3, 7))
    beta_ldpred <- numeric(1e6)
    beta_ldpred[match(res_ldpred$sid, simu$map$marker.ID)] <- res_ldpred[[2]]
    beta_ldpred
  })

pred_train_ldpred <- big_prodMat(G.train, beta_ldpred)
auc_train_ldpred <- apply(pred_train_ldpred, 2, AUC, target = y[ind.train])

pred_ldpred <- big_prodVec(G.test, beta_ldpred[, which.max(auc_train_ldpred)])
AUCBoot(pred_ldpred, y[ind.test]) # 71.7 [70.1-73.3]

unlink(paste0(tmp, "*"))
