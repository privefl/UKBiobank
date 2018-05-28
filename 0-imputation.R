library(bigsnpr)
ukb <- snp_attach("UKB-merged.rds")
G <- ukb$genotypes
# system.time(
#   counts <- big_counts(G)
# ) # 17 min
# 
# counts[, 1:10]
# hist(counts[4, ])

CHR <- ukb$map$chromosome
infos <- snp_fastImpute(G, CHR, n.cor = 20e3, ncores = 6)

