sumstats <- bigreadr::fread2("breast-cancer/oncoarray_bcac_public_release_oct17.txt",
                             na.strings = "NULL")
sumstats <- subset(sumstats, bcac_icogs2_P1df_Wald < 0.1)
dim(sumstats) # 1,513,288 x 66

library(foreach)
system.time({
  sumstats_joined <- foreach(chr = 20:22, .combine = "rbind") %do% {
    cat("Processing chromosome", chr, "..\n")
    file <- paste0("ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
    infos_chr <- bigreadr::fread2(file, showProgress = FALSE)
    ind <- which(sumstats$chr == chr)
    dplyr::inner_join(sumstats[ind, ], infos_chr,
                      by = c(position_b37 = "V3", a0 = "V4", a1 = "V5"),
                      na_matches = "never")
  }
}) # <8min

dim(sumstats_joined) # 1,471,855 x 69  (Chr 1: 114,483)

hist(pmin(sumstats_joined$bcac_icogs2_r2, sumstats_joined$V8))

str(sumstats_joined)
saveRDS(sumstats_joined, "sumstats_all_chr.rds")


# ID <- bigsnpr::snp_attach("UKB-merged_sub1.rds")$fam$sample.ID

## with R package: super slow
data <- bigsnpr.bgen::bgen.load(
  "ukb_imp_chr1_v3.bgen",
  index.filename = "ukb_imp_bgi/ukb_imp_chr1_v3.bgen.bgi",
  rsids = sumstats_joined$V2[1:100 + 100]
) # 242 Mb for 10 // 1.4 GB for 100
str(data)

## with PLINK 2.0
plink2 <- path.expand("~/Bureau/plink2")
tmp <- tempfile()
write(sumstats_joined$V2, tmp2 <- paste0(tmp, "_SNPs.txt"))
bgen_file <- "ukb_imp_chr1_v3.bgen"
sample_file <- "ukb25589_imp_chr1_v3_s487327.sample"
# Write bim and fam files from the BGEN file
system(glue::glue("{plink2} --bgen {bgen_file} --sample {sample_file}",
                  # " --make-just-bim --make-just-fam",
                  " --extract {tmp2}",
                  " --make-bpgen",
                  # " --sort-vars",
                  " --out {tmp}"))
# Get the variants and individual information
fam <- data.table::fread(paste0(tmp, ".fam"), data.table = FALSE)
names(fam) <- bigsnpr:::NAMES.FAM
n <- nrow(fam)
bim <- data.table::fread(paste0(tmp, ".bim"), data.table = FALSE)
names(bim) <- bigsnpr:::NAMES.MAP
m <- nrow(bim)
