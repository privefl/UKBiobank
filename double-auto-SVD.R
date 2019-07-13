library(bigsnpr)

obj.bed <- bed.1000G <- bed("tmp-data/1000G_phase3_common_hapmap.bed")
map.1000G <- bigreadr::fread2(
  bed.1000G$bimfile,
  col.names = c("chr", "rsid", "osef", "pos", "a1", "a0")
)

info_snp <- snp_match(
  cbind(map.1000G[-3], beta = 1),
  setNames(bed("data/celiacQC_sub1.bed")$map[-3], c("chr", "rsid", "pos", "a1", "a0")),
  join_by_pos = FALSE
)
info_snp


fam <- bed.1000G$fam
ped <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped")
fam2 <- dplyr::left_join(fam[c(2, 5)], ped[c(1:5, 7)],
                         by = c("sample.ID" = "Individual ID", "sex" = "Gender"))
pop <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv")
fam2 <- dplyr::left_join(fam2, pop[1:3], by = c("Population" = "Population Code"))
str(fam2)


# AutoSVD algo
NCORES <- 10
nPC <- 30
ind.col <- info_snp$`_NUM_ID_.ss`
infos.chr <- obj.bed$map$chromosome
infos.pos <- obj.bed$map$physical.pos

# First clumping
ind.keep <- bed_clumping(obj.bed,
                         exclude = setdiff(cols_along(obj.bed), ind.col),
                         thr.r2 = 0.2,
                         size = 500,
                         ncores = NCORES)

svd.1000G <- obj.svd <- bed_randomSVD(bed.1000G, k = nPC,
                                      ind.col = ind.keep,
                                      ncores = NCORES)

library(ggplot2)
plot(svd.1000G) + scale_y_log10()
plot(svd.1000G, type = "scores")
system.time(
  p_scores <- plot(svd.1000G, type = "scores", scores = 1:nPC,
                   cols = 5, coeff = 0.7)
)
system.time(
  p_loadings <- plot(svd.1000G, type = "loadings", loadings = 22:nPC,
                     cols = 3, coeff = 0.4)
)

# -log10(p-values) of being an outlier
lpval <- -stats::predict(
  bed_pcadapt(obj.bed, obj.svd$u, ind.col = ind.keep, ncores = NCORES))
# Bonferroni-corrected threshold
lim <- -log10(0.05 / length(lpval))

tukey <- function(stat, coef = 1.5) {
  quantile(stat, probs = 0.75) + coef * IQR(stat)
}

# Roll mean to get only consecutive outliers
lpval2 <- bigsnpr:::rollMean(lpval, size = 50)

plot(lpval2, pch = 20, log = "y"); abline(h = lim, col = "red")
stat <- log(lpval2)
q <- tukey(stat, 2); abline(h = exp(q), col = "blue")
p <- pchisq(((stat - median(stat)) / mad(stat))^2, df = 1, lower.tail = FALSE)
q_p <- qchisq(0.05 / length(stat), df = 1, lower.tail = FALSE)
abline(h = exp(sqrt(q_p) * mad(stat) + median(stat)), col = "green")

hist(stat, breaks = 50) ## rollmean makes it more normal
q <- tukey(stat, 2); abline(v = q, col = "blue")
abline(v = sqrt(q_p) * mad(stat) + median(stat), col = "green")

which(lpval2 > lim)
which(stat > q)
which(keep <- p > (0.05 / length(stat)))
(ind.range <- bigsnpr:::getIntervals(which(!keep), n = 10))
cbind(
  infos.chr[ind.keep[ind.range[, 1]]],
  matrix(infos.pos[ind.keep[ind.range]], ncol = 2)
)
# [1,]    2 135208077 137155850
# [2,]    6  25370316  33045558
# [3,]    8   6643301  13066409

U <- obj.svd$u
maha <- robust::covRob(U, estim = "pairwiseGK", distance = FALSE, corr = FALSE)
Rcpp::sourceCpp('compute-dist2.cpp')
system.time({
  eigs <- eigen(unname(maha$cov))
  U2 <- t(U %*% sweep(eigs$vectors, 2, sqrt(eigs$values), '/'))
  ord <- order(U2[1, ])
  dist15 <- dist_kNN_all_clever(U2[, ord], 10) / 10
  dist15 <- log(dist15[match(seq_along(ord), ord)])
})
hist(dist15, breaks = 50)
q2 <- tukey(dist15, 1); abline(v = q2, col = "red")

cowplot::plot_grid(plotlist = lapply(1:15, function(k) {
  plot(svd.1000G, type = "scores", scores = c(2 * k - 1, 2 * k), coeff = 0.5) +
    aes(color = dist15) + scale_colour_viridis_c() +
    theme(legend.position = "none") + ggtitle(NULL)
}))

cowplot::plot_grid(plotlist = lapply(1:15, function(k) {
  plot(svd.1000G, type = "scores", scores = c(2 * k - 1, 2 * k), coeff = 0.5) +
    aes(color = (dist15 > q2)) + scale_colour_viridis_d() +
    theme(legend.position = "none") + ggtitle(NULL)
}))


#### ITERATION 2 ####

no_out_ind <- which(dist15 < q2)
stopifnot(length(stat) == length(ind.keep))
no_out_var <- ind.keep[keep]

svd.1000G <- obj.svd2 <- bed_randomSVD(bed.1000G, k = nPC,
                                       ind.row = no_out_ind,
                                       ind.col = no_out_var,
                                       ncores = NCORES)

library(ggplot2)
plot(svd.1000G) + scale_y_log10()
plot(svd.1000G, type = "scores")
system.time(
  p_scores <- plot(svd.1000G, type = "scores", scores = 1:nPC,
                   cols = 5, coeff = 0.7)
)
system.time(
  p_loadings <- plot(svd.1000G, type = "loadings", loadings = 22:nPC,
                     cols = 3, coeff = 0.4)
)

# -log10(p-values) of being an outlier
lpval <- -stats::predict(
  bed_pcadapt(obj.bed, obj.svd2$u, ind.row = no_out_ind, ind.col = no_out_var,
              ncores = NCORES))
plot(lpval, pch = 20)
# Bonferroni-corrected threshold
lim <- -log10(0.05 / length(lpval))

# Roll mean to get only consecutive outliers
lpval2 <- bigsnpr:::rollMean(lpval, size = 50)

plot(lpval2, pch = 20, log = "y"); abline(h = lim, col = "red")
stat <- log(lpval2)
q <- tukey(stat, 2); abline(h = exp(q), col = "blue")
p <- pchisq(((stat - median(stat)) / mad(stat))^2, df = 1, lower.tail = FALSE)
q_p <- qchisq(0.05 / length(stat), df = 1, lower.tail = FALSE)
abline(h = exp(sqrt(q_p) * mad(stat) + median(stat)), col = "green")

hist(stat, breaks = 50) ## rollmean makes it more normal
q <- tukey(stat, 2); abline(v = q, col = "blue")
abline(v = sqrt(q_p) * mad(stat) + median(stat), col = "green")

which(lpval2 > lim)
which(stat > q)
which(keep <- p > (0.05 / length(stat)))
sum(!keep)
(ind.range <- bigsnpr:::getIntervals(which(!keep), n = 10))
cbind(
  infos.chr[ind.keep[ind.range[, 1]]],
  matrix(infos.pos[ind.keep[ind.range]], ncol = 2)
)


U <- obj.svd2$u
maha <- robust::covRob(U, estim = "pairwiseGK", distance = FALSE, corr = FALSE)
Rcpp::sourceCpp('compute-dist2.cpp')
system.time({
  eigs <- eigen(unname(maha$cov))
  U2 <- t(U %*% sweep(eigs$vectors, 2, sqrt(eigs$values), '/'))
  ord <- order(U2[1, ])
  dist15 <- dist_kNN_all_clever(U2[, ord], 10) / 10
  dist15 <- log(dist15[match(seq_along(ord), ord)])
})
hist(dist15, breaks = 50)
q2 <- tukey(dist15, 1); abline(v = q2, col = "red")

cowplot::plot_grid(plotlist = lapply(1:15, function(k) {
  plot(svd.1000G, type = "scores", scores = c(2 * k - 1, 2 * k), coeff = 0.5) +
    aes(color = dist15) + scale_colour_viridis_c() +
    theme(legend.position = "none") + ggtitle(NULL)
}))

cowplot::plot_grid(plotlist = lapply(1:15, function(k) {
  plot(svd.1000G, type = "scores", scores = c(2 * k - 1, 2 * k), coeff = 0.5) +
    aes(color = (dist15 > q2)) + scale_colour_viridis_d() +
    theme(legend.position = "none") + ggtitle(NULL)
}))


#### ITERATION 3 ####

no_out_ind2 <- no_out_ind[dist15 < q2]
stopifnot(length(stat) == length(no_out_var))
no_out_var2 <- no_out_var[stat < q]

svd.1000G <- obj.svd3 <- bed_randomSVD(bed.1000G, k = nPC,
                                       ind.row = no_out_ind2,
                                       ind.col = no_out_var2,
                                       ncores = NCORES)

library(ggplot2)
plot(svd.1000G) + scale_y_log10()
plot(svd.1000G, type = "scores")
system.time(
  p_scores <- plot(svd.1000G, type = "scores", scores = 1:nPC,
                   cols = 5, coeff = 0.7)
)
system.time(
  p_loadings <- plot(svd.1000G, type = "loadings", loadings = 22:nPC,
                     cols = 3, coeff = 0.4)
)

# -log10(p-values) of being an outlier
lpval <- -stats::predict(
  bed_pcadapt(obj.bed, obj.svd3$u, ind.row = no_out_ind2, ind.col = no_out_var2,
              ncores = NCORES))
# Bonferroni-corrected threshold
lim <- -log10(0.05 / length(lpval))

# Roll mean to get only consecutive outliers
lpval2 <- bigsnpr:::rollMean(lpval, size = 50)

plot(lpval2, pch = 20, log = "y"); abline(h = lim, col = "red")
stat <- log(lpval2)
q <- tukey(stat, 2); abline(h = exp(q), col = "blue")
p <- pchisq(((stat - median(stat)) / mad(stat))^2, df = 1, lower.tail = FALSE)
q_p <- qchisq(0.05 / length(stat), df = 1, lower.tail = FALSE)
abline(h = exp(sqrt(q_p) * mad(stat) + median(stat)), col = "green")

hist(stat, breaks = 50) ## rollmean makes it more normal
q <- tukey(stat, 2); abline(v = q, col = "blue")
abline(v = sqrt(q_p) * mad(stat) + median(stat), col = "green")

which(lpval2 > lim)
which(stat > q)
which(keep <- p > (0.05 / length(stat)))
sum(!keep)
(ind.range <- bigsnpr:::getIntervals(which(!keep), n = 10))
cbind(
  infos.chr[ind.keep[ind.range[, 1]]],
  matrix(infos.pos[ind.keep[ind.range]], ncol = 2)
)

U <- obj.svd3$u
maha <- robust::covRob(U, estim = "pairwiseGK", distance = FALSE, corr = FALSE)
Rcpp::sourceCpp('compute-dist2.cpp')
system.time({
  eigs <- eigen(unname(maha$cov))
  U2 <- t(U %*% sweep(eigs$vectors, 2, sqrt(eigs$values), '/'))
  ord <- order(U2[1, ])
  dist15 <- dist_kNN_all_clever(U2[, ord], 10) / 10
  dist15 <- log(dist15[match(seq_along(ord), ord)])
})
hist(dist15, breaks = 50)
q2 <- tukey(dist15, 1); abline(v = q2, col = "red")

cowplot::plot_grid(plotlist = lapply(1:15, function(k) {
  plot(svd.1000G, type = "scores", scores = c(2 * k - 1, 2 * k), coeff = 0.5) +
    aes(color = dist15) + scale_colour_viridis_c() +
    theme(legend.position = "none") + ggtitle(NULL)
}))

cowplot::plot_grid(plotlist = lapply(1:15, function(k) {
  plot(svd.1000G, type = "scores", scores = c(2 * k - 1, 2 * k), coeff = 0.5) +
    aes(color = (dist15 > q2)) + scale_colour_viridis_d() +
    theme(legend.position = "none") + ggtitle(NULL)
}))


svd.1000G <- obj.svd3 <- bed_randomSVD(bed.1000G, k = 19,
                                       ind.row = no_out_ind[dist15 < q2],
                                       ind.col = no_out_var[stat < q & lpval < 5],
                                       ncores = NCORES)
plot(svd.1000G, type = "scores", scores = 1:19, coeff = 0.7)

# -log10(p-values) of being an outlier
system.time(
  gwas <- bed_pcadapt(obj.bed, obj.svd3$u, ind.row = no_out_ind2, ncores = NCORES)
)
snp_qq(gwas) + xlim(1, NA)
snp_manhattan(gwas, infos.chr, infos.pos, npoints = 20e3) +
  geom_hline(yintercept = -log10(5e-8), col = "red")

sum(signif <- (predict(gwas, log10 = FALSE) < 1e-12))
indep <- bed_clumping(obj.bed, ind.row = no_out_ind[dist15 < q2],
                      S = -predict(gwas, log10 = FALSE), thr.r2 = 0.05, size = 2000,
                      exclude = which(!signif), ncores = NCORES)
(rsid <- obj.bed$map$marker.ID[indep])
devtools::source_gist("42b41d771bbeae63245b8304ef283c70", filename = "get-genes.R")
genes <- snp_gene(rsid)
indep2 <- indep[!is.na(genes)]


library(dplyr)
test <- fam2 %>%
  mutate(id = rows_along(fam2)) %>%
  slice(no_out_ind2) %>%
  select(Population, `Super Population`, id) %>%
  group_by(`Super Population`) %>%
  # group_by(`Super Population`, Population) %>%
  summarize(stat = list(bigsnpr:::bed_stats(obj.bed, id, ind_col = indep2)))

AF <- setNames(as.data.frame(t(sapply(test$stat, function(s) s$sum / s$nb_nona_col))),
               rsid[!is.na(genes)])
test %>%
  select(-stat) %>%
  cbind(N = sapply(test$stat, function(s) s$nb_nona_col)[1, ], af = AF) %>%
  tibble::repair_names() %>%
  print() %>%
  mutate_at(-(1:2), ~ `if`(mean(. > 1), 2 - ., .)) %>%
  tidyr::gather("num_id", "af", -(1:2)) %>%
  ggplot() +
  theme_bigstatsr() +
  geom_col(aes(x = "", y = af, fill = `Super Population`), position = position_dodge(width = 0.9)) +
  facet_wrap(~ num_id, nrow = 3) +
  theme(legend.position = "top") +
  labs(x = "", y = "Allele frequency")
