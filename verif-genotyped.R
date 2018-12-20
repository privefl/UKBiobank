library(bigsnpr)

ukb <- snp_attach("data/UKB-merged.rds")
map <- ukb$map

file <- paste0("ukb_imp_mfi/ukb_mfi_chr", 1, "_v3.txt")
infos_chr <- bigreadr::fread2(file, showProgress = FALSE)
map_chr1 <- subset(map, chromosome == 1)
mean(infos_chr$V6 > 0.001 & infos_chr$V8 > 0.9)
test <- map_chr1$physical.pos %in% infos_chr$V3
mean(test)

ind <- match(map_chr1$physical.pos, infos_chr$V3)
infos_chr[head(ind), ]
head(map_chr1)
summary(infos_chr$V8[ind])


