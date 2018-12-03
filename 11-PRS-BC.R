library(bigsnpr)
ukb_imp <- snp_attach("data/UKB_imp_BC.rds")
G <- ukb_imp$genotypes
CHR <- as.integer(ukb_imp$map$chromosome)
POS <- ukb_imp$map$physical.pos
sumstats_joined <- readRDS("breast-cancer/sumstats_all_chr.rds")
stopifnot(ncol(G) == nrow(sumstats_joined))

hist(lpval <- -log10(sumstats_joined$bcac_icogs2_P1df_Wald))
hist(r2_impute <- pmin(sumstats_joined$bcac_icogs2_r2, sumstats_joined$V8))
hist(betas <- sumstats_joined$bcac_icogs2_beta)

grid <- expand.grid(
  thr.imp = c(0.3, 0.6, 0.8, 0.9, 0.95),
  thr.clmp = c(0.05, 0.2, 0.5, 0.8),
  chr = unique(CHR)
)

n_thr_pval <- 20
scores <- FBM(nrow(G), nrow(grid) * n_thr_pval,
              backingfile = "breast-cancer/scores")$save()


set.seed(1)
ind.train <- sort(sample(nrow(G), 150e3))
ind.test <- setdiff(rows_along(G), ind.train)


library(doParallel)
registerDoParallel(cl <- makeCluster(12, outfile = ""))
pb <- txtProgressBar(0, nrow(grid), style = 3)
foreach(ic = rows_along(grid)) %dopar% {

  chr      <- grid[ic, "chr"]
  thr.imp  <- grid[ic, "thr.imp"]
  thr.clmp <- grid[ic, "thr.clmp"]

  ind.keep <- bigsnpr::snp_clumping(
    G, CHR, S = lpval, thr.r2 = thr.clmp,
    ind.row = sort(sample(ind.train, 50e3)),
    exclude = which(CHR != chr | r2_impute < thr.imp),
    infos.pos = POS, is.size.in.bp = TRUE
  )

  lp.thr <- exp(seq(log(0.1), log(0.99 * max(lpval[ind.keep])),
                    length.out = n_thr_pval))

  prs <- bigsnpr::snp_PRS(G, betas.keep = betas[ind.keep], ind.keep = ind.keep,
                          lpS.keep = lpval[ind.keep], thr.list = lp.thr)

  scores[, (ic - 1) * n_thr_pval + 1:n_thr_pval] <- prs

  setTxtProgressBar(pb, ic)
  NULL
}
close(pb)
stopCluster(cl)


# dim(prs)
# dim(G)
# y <- ukb_imp$fam$has_cancer
# y[y] <- NA  ## set to missing all types of cancer
# y[ukb_imp$fam$cancer_type == "breast cancer"] <- 1  ## keep only BC
# ind <- which(ukb_imp$fam$sex == "Female" & !is.na(y))  ## only women
#
# print(max(aucs <- apply(prs[ind, ], 2, AUC, target = y[ind])))  # 0.61 -> 0.618 -> 0.619
# plot(lp.thr, aucs, pch = 20, log = "x")
