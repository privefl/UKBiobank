library(doParallel)
registerDoParallel(cl <- makeCluster(12))
system.time({
  list_snp_id <- foreach(chr = 1:22) %dopar% {
    # cat("Processing chromosome", chr, "..\n")
    mfi <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
    infos_chr <- bigreadr::fread2(mfi, showProgress = FALSE)
    infos_chr_sub <- subset(infos_chr, V6 > 0.01 & V8 > 0.99)
    with(infos_chr_sub, paste(chr, V3, V4, V5, sep = "_"))
  }
}) # 80 sec
stopCluster(cl)

sample <- bigreadr::fread2("ukb25589_imp_chr1_v3_s487327.sample")
str(sample)
(N <- readBin("data/ukb_imp_chr1_v3.bgen", what = 1L, size = 4, n = 4)[4])
sample <- sample[-1, ]
stopifnot(nrow(sample) == N)

df0 <- readRDS("pheno.rds")
rel <- bigreadr::fread2("ukb25589_rel_s488346.dat")
sum(ind.rm <- df0$eid %in% rel[["ID2"]])  ## 43,400
dim(df0[!ind.rm, ][])
df1 <- subset(df0, !eid %in% rel[["ID2"]] & is_caucasian)
set.seed(2)
ind.indiv <- sample(na.omit(match(df1$eid, sample$ID_2)), 50e3) # among 335,609

system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = c(6, 8)),
    list_snp_id = list_snp_id[c(6, 8)],
    backingfile = "simus/UKB_unrel_chr68",
    ind_row = ind.indiv,
    bgi_dir = "data/ukb_imp_bgi",
    ncores = 10
  )
) # 72 min


library(bigsnpr)
ukb <- snp_attach("simus/UKB_unrel_chr68.rds")
G <- ukb$genotypes$copy(code = round(CODE_DOSAGE))
CHR <- as.numeric(ukb$map$chromosome)
POS <- ukb$map$physical.pos
dim(G) # 50,000 686,066
file.size(G$backingfile) / 1024^3  # 32 GB


set.seed(1)
ind.gwas <- sort(sample(50e3, 40e3))
ind.nogwas <- setdiff(1:50e3, ind.gwas)
ind.train <- sort(sample(ind.nogwas, 8e3))
ind.test <- setdiff(ind.nogwas, ind.train)
M <- 1000
pheno <- pkg.paper.PRS::get_pheno(G, h2 = 0.5, M = M, K = 0.2)
saveRDS(pheno, glue::glue("simus/simu{M}.rds"))

gwas <- big_univLogReg(G, pheno$pheno[ind.gwas], ind.train = ind.gwas,
                       ncores = nb_cores())
snp_qq(gwas) + ggplot2::xlim(1, NA)
snp_manhattan(gwas, CHR, POS, ind.highlight = pheno$set, npoints = 200e3)

map <- ukb$map
map$beta <- gwas$estim
map$pval <- predict(gwas, log10 = FALSE)
str(map)
bigreadr::fwrite2(map[-2], "simus/sumstats.txt")
readLines("simus/sumstats.txt", n = 5)

ukb$fam <- data.frame(
  family.ID = 1:50e3,
  sample.ID = 1:50e3,
  paternal.ID = 0,
  maternal.ID = 0,
  sex = as.numeric(sample[ind.indiv, "sex"]),
  affection = pheno$pheno
)

ukb$map$genetic.dist <- 0
snp_writeBed(ukb, "simus/data_gwas.bed",  ind.row = ind.gwas)
snp_writeBed(ukb, "simus/data_train.bed", ind.row = ind.train)
snp_writeBed(ukb, "simus/data_test.bed",  ind.row = ind.test)

# big_write(G, "simus/data_train.txt", ind.row = ind.train,
#           every_nrow = 500, progress = TRUE)
# bigreadr::nlines("simus/data_train.txt")
# big_write(G, "simus/data_test.txt", ind.row = ind.test,
#           every_nrow = 500, progress = TRUE)
