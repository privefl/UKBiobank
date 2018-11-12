sumstats_joined <- readRDS("breast-cancer/sumstats_joined.rds")[c("chr", "V1")]

library(bigsnpr)
system.time(
  test <- snp_readBGEN(
    bgenfiles = glue::glue("ukb_imp_chr{chr}_v3.bgen", chr = 1:3),
    backingfile = "test_bgen13",
    list_snp_id = split(sumstats_joined$V1, sumstats_joined$chr),
    ind_row = 1:5000,
    ncores = 3
  )
) # 1.8h // 46 min // 10 min ?????
# tmp <- readRDS("/tmp/RtmpAsN7Oy/file665f42353338.rds")

snp <- snp_attach(test)
G <- snp$genotypes
dim(G)
hist(G[, 1])
str(snp$map)

#### 1 core
#     user   system  elapsed
# 1624.951   35.137 3109.432
#     user   system  elapsed
# 1619.318   14.406 1636.606
#### 3 cores
#  user  system elapsed
# 3.711   1.959 859.998
#  user  system elapsed
# 3.544   1.889 643.408


system.time(
  test <- snp_readBGEN(
    bgenfiles = glue::glue("ukb_imp_chr{chr}_v3.bgen", chr = 1),
    backingfile = "test_bgen7",
    list_snp_id = list(sumstats_joined$V1[1:10e3]),
    ind_row = 1:5000
  )
) # 57 sec
