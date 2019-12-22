lam_gc_chr <- by(gwas$score^2, CHR, FUN = function(x) median(x) / qchisq(0.5, df = 1))
lam_gc_chr[c(FALSE, TRUE)]

table(CHR[gwas$score^2 > 25])
snp_qq(gwas[!is_odd_chr, ]) + xlim(1, NA) +
  aes(color = as.factor(CHR[!is_odd_chr]))


chr <- 8
repeat {
  ind_g <- sample(which(CHR == chr), 1)
  g <- G[, ind_g]
  af <- mean(g) / 2
  maf <- min(af, 1 - af)
  if (maf > 0.3) break
}
system.time(
  gwas_g_all_chr <- lapply(setdiff(1:22, chr), function(chr) {
    big_univLinReg(G, g, ind.col = which(CHR == chr),
                   covar.train = cbind(PC, Reduce('+', poly[-chr])),
                   ncores = 12)
  })
)
gwas_g <- do.call("rbind", gwas_g_all_chr)
plot(gwas_g)
snp_qq(gwas_g) + xlim(1, NA)

res <- matrix(NA, 22, 2)
by(gwas_g$score^2, CHR[CHR != chr], FUN = function(x) median(x) / qchisq(0.5, df = 1))[]
res[as.integer(names(.Last.value)), 2] <- .Last.value

# gwas_g1 <- gwas_g
by(gwas_g1$score^2, CHR[CHR != 1], FUN = function(x) median(x) / qchisq(0.5, df = 1))[]
res[as.integer(names(.Last.value)), 1] <- .Last.value
plot(res)
summary(lm(I(res[, 1] - 1) ~ I(res[, 2] - 1)))
bigsnpr:::getLambdaGC(gwas_g1)
bigsnpr:::getLambdaGC(gwas_g)
