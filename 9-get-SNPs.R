sumstats <- bigreadr::fread2("breast-cancer/oncoarray_bcac_public_release_oct17.txt",
                             na.strings = "NULL")
sumstats <- subset(sumstats, bcac_icogs2_P1df_Wald < 0.1)
dim(sumstats) # 1,513,288 x 66

library(foreach)
system.time({
  sumstats_joined <- foreach(chr = 1:22, .combine = "rbind") %do% {
    cat("Processing chromosome", chr, "..\n")
    file <- paste0("ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
    infos_chr <- bigreadr::fread2(file, showProgress = FALSE)
    ind <- which(sumstats$chr == chr)
    dplyr::inner_join(sumstats[ind, ], infos_chr,
                      by = c(position_b37 = "V3", a0 = "V4", a1 = "V5"),
                      na_matches = "never")
  }
}) # <8min

dim(sumstats_joined) # 1,448,906 x 71
str(sumstats_joined)
saveRDS(sumstats_joined, "breast-cancer/sumstats_all_chr.rds")

