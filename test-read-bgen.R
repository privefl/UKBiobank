sumstats_joined <- readRDS("sumstats_all_chr.rds")

IDs <- with(sumstats_joined, paste(chr, position_b37, a0, a1, sep = ":"))

library(bigsnpr)
system.time(
  test <- snp_readBGEN(
    bgenfiles = glue::glue("ukb_imp_chr{chr}_v3.bgen", chr = 20:22),
    backingfile = "test_bgen16",
    list_snp_id = split(IDs, sumstats_joined$chr),
    ind_row = 1:50000,
    bgi_dir = "ukb_imp_bgi/",
    ncores = 3
  )
) # 1.8h // 46 min // 10 min ?????
# tmp <- readRDS("ukb_imp_bgi//ukb_imp_chr20_v3_not_found.rds")

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
