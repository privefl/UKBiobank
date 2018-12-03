library(bigsnpr)
ukb_imp <- snp_attach("data/UKB_imp_BC.rds")
G <- ukb_imp$genotypes
CHR <- as.integer(ukb_imp$map$chromosome)
POS <- ukb_imp$map$physical.pos
sumstats_joined <- readRDS("breast-cancer/sumstats_all_chr.rds")
stopifnot(ncol(G) == nrow(sumstats_joined))

ind.row <- sort(sample(nrow(G), 50e3))
lpval <- -log10(sumstats_joined$bcac_icogs2_P1df_Wald)
summary(10^(-lpval))
hist(r2_impute <- pmin(sumstats_joined$bcac_icogs2_r2, sumstats_joined$V8))
hist(af <- sumstats_joined$bcac_icogs2_eaf_controls)
maf <- pmin(af, 1 - af)

excl <- (r2_impute < 0.8)
table(CHR[!excl])
table(CHR)
system.time(
  ind.keep <- snp_clumping(
    G, CHR, ind.row = ind.row, S = lpval,
    exclude = which(excl), thr.r2 = 0.2,
    infos.pos = POS, is.size.in.bp = TRUE,
    ncores = 4
  )
)
#### ncores = 6 -> must have had first ones in cache
# Chromosome 1 - init - 19 sec
# Chromosome 4 - init - 43 sec
# Chromosome 1 - clumping - 37 sec
# Chromosome 4 - clumping - 28 sec
# Chromosome 3 - init - 166 sec
# Chromosome 3 - clumping - 28 sec
# Chromosome 6 - init - 447 sec
# Chromosome 5 - init - 526 sec
# Chromosome 6 - clumping - 133 sec
# Chromosome 5 - clumping - 116 sec
# Chromosome 7 - init - 760 sec
# Chromosome 10 - init - 815 sec
# Chromosome 8 - init - 720 sec
# Chromosome 2 - init - 996 sec
# Chromosome 11 - init - 384 sec
# Chromosome 12 - init - 460 sec
# Chromosome 7 - clumping - 289 sec
# Chromosome 10 - clumping - 404 sec
# Chromosome 8 - clumping - 391 sec
# Chromosome 11 - clumping - 348 sec
# Chromosome 12 - clumping - 420 sec
# Chromosome 9 - init - 485 sec
# Chromosome 9 - clumping - 149 sec

#### ncores = 2
# Chromosome 2 - init - 474 sec
# Chromosome 1 - init - 478 sec
# Chromosome 1 - clumping - 31 sec
# Chromosome 2 - clumping - 39 sec
# Chromosome 5 - init - 371 sec
# Chromosome 4 - init - 396 sec
# Chromosome 5 - clumping - 40 sec
# Chromosome 4 - clumping - 25 sec
# Chromosome 6 - init - 306 sec
# Chromosome 3 - init - 329 sec
# Chromosome 6 - clumping - 44 sec
# Chromosome 3 - clumping - 29 sec


#### ncores = 4
# Chromosome 4 - init - 372 sec
# Chromosome 4 - clumping - 24 sec
# Chromosome 1 - init - 692 sec
# Chromosome 5 - init - 734 sec
# Chromosome 1 - clumping - 218 sec
# Chromosome 2 - init - 930 sec
# Chromosome 6 - init - 607 sec
# Chromosome 5 - clumping - 314 sec
# Chromosome 6 - clumping - 297 sec
# Chromosome 7 - init - 385 sec
# Chromosome 2 - clumping - 576 sec
# Chromosome 7 - clumping - 78 sec
# Chromosome 3 - init - 679 sec
# Chromosome 10 - init - 413 sec
# Chromosome 10 - clumping - 61 sec
# Chromosome 12 - init - 343 sec
# Chromosome 8 - init - 380 sec
# Chromosome 12 - clumping - 40 sec
# Chromosome 3 - clumping - 317 sec
# Chromosome 8 - clumping - 60 sec
# Chromosome 11 - init - 266 sec
# Chromosome 11 - clumping - 34 sec
# Chromosome 14 - init - 373 sec
# Chromosome 9 - init - 401 sec
# Chromosome 17 - init - 374 sec
# Chromosome 14 - clumping - 44 sec
# Chromosome 17 - clumping - 28 sec
# Chromosome 9 - clumping - 72 sec
# Chromosome 13 - init - 304 sec
# Chromosome 13 - clumping - 15 sec
# Chromosome 18 - init - 339 sec
# Chromosome 16 - init - 367 sec
# Chromosome 15 - init - 327 sec
# Chromosome 19 - init - 304 sec
# Chromosome 18 - clumping - 13 sec
# Chromosome 16 - clumping - 15 sec
# Chromosome 19 - clumping - 8 sec
# Chromosome 15 - clumping - 11 sec
# Chromosome 21 - init - 124 sec
# Chromosome 21 - clumping - 5 sec
# Chromosome 22 - init - 132 sec
# Chromosome 22 - clumping - 6 sec
# Chromosome 20 - init - 155 sec
# Chromosome 20 - clumping - 7 sec
