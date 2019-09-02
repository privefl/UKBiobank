# # test <- genie::hclust2(objects = svd.1000G$u, )
# d <- dist(svd.1000G$u)
# h <- fastcluster::hclust(d)
# plot(h)
# cl <- as.factor(stats::cutree(h, k = 5))
#
# plot(svd.1000G, type = "scores", scores = 13:14) + aes(color = cl)
# plot(svd.1000G, type = "scores", scores = 19:20) + aes(color = cl)
# plot(svd.1000G, type = "scores", scores = 19:20 + 2) + aes(color = cl)
# plot(svd.1000G, type = "scores", scores = 19:20 + 4) + aes(color = cl)
#
# K_seq <- seq(2, 40)
# stat <- sapply(K_seq, function(k) {
#   mean(cluster::silhouette(cutree(h, k), d)[, 3])
# })
# plot(K_seq, stat, pch = 20)
#
# cl <- as.factor(stats::cutree(h, k = K_seq[which.min(stat)]))
#
# plot(svd.1000G, type = "scores", scores = 13:14) + aes(color = cl)
# plot(svd.1000G, type = "scores", scores = 19:20) + aes(color = cl)
# plot(svd.1000G, type = "scores", scores = 19:20 + 2) + aes(color = cl)
# plot(svd.1000G, type = "scores", scores = 19:20 + 4) + aes(color = cl)
# num_cl <- table(cl)
# plot(dist0, num_cl[cl], log = "xy", pch = 20)
# plot(dist0, num_cl[cl], log = "xy", pch = 20)
#
# dist7 <- log(exp(dist) / sqrt(as.vector(num_cl[cl])))
# hist(dist7, breaks = 20)

hist(dist <- log(robust::covRob(svd.1000G$u[, 1:30], estim = "pairwiseGK")$dist))
# hist(dist <- log(rowSums(sweep(svd.1000G$u, 2, 1:30, '*')^2)))

nb <- sapply(1:30, function(k) {
  pc1 <- svd.1000G$u[, k]
  sapply(pc1, function(x) sum(abs(x - pc1) < mad(pc1)))
})


# plot(rowMeans(nb), nb2 <- apply(nb, 1, sd))
# nb2 <- exp(rowMeans(log(nb)))
nb2 <- apply(nb, 1, quantile, probs = 0.25)
nb2 <- apply(nb, 1, function(x) exp(mean(log(x[x < median(x)]))))
nb3 <- colMeans(as.matrix(d))
# d2  <- as.matrix(d)
# diag(d2) <- 1
# nb3 <- colMeans(log(d2))
plot(nb3, nb2, pch = 20)
nb2 <- nb3
hist(log(nb2))
# abline(v = sd(pc1), col = "red")
# abline(v = 2 * sd(pc1), col = "blue")
summary(nb2)
plot(nb2, nb, pch = 20)

hist(dist7 <- log(nb2), breaks = 20)


# d <- dist(svd.1000G$u)
d <- dist(predict(svd.1000G))
# d <- dist(sweep(svd.1000G$u, 2, svd.1000G$d^2, '/'))
hist(dist7 <- log(apply(as.matrix(d), 2, function(x) crossprod(sort(x), seq_along(x)^2))))
hist(dist7.2 <- log(apply(as.matrix(d), 2, function(x) crossprod(sort(x^2), seq_along(x)^2))))
plot(dist7.2, dist7, pch = 20)
# hist(dist7 <- log(dist7))
dist2 <- log(robust::covRob(svd.1000G$u, estim = "pairwiseGK")$dist)
plot(dist2, dist7, pch = 20)


d <- dist(svd.1000G$u)
hist(dist7 <- log(apply(as.matrix(d), 2, function(x) mean(sort(x)[2:11]))))
hist(dist7.2 <- apply(as.matrix(d), 2, function(x) mean(log(sort(x)[2:11]))))
plot(dist7, dist2, pch = 20)
plot(dist7.2, dist7, pch = 20)  # same

# dist2 <- (dist2 - median(dist2)) / mad(dist2)
# dist7 <- (dist7 - median(dist7)) / mad(dist7)
# plot(dist2, dist7, pch = 20); abline(0, 1, col = "red")
# hist(dist7 <- dist7 + dist2, breaks = 50)

hist(dist7, breaks = 50)
print(q  <- quantile(dist7, 0.75) + 1   * IQR(dist7)); abline(v = q,  col = "blue")
print(q2 <- quantile(dist7, 0.75) + 1.5 * IQR(dist7)); abline(v = q2, col = "red")
print(q3 <- quantile(dist7, 0.75) + 2   * IQR(dist7)); abline(v = q3, col = "green")
print(q4 <- quantile(dist7, 0.75) + 3   * IQR(dist7)); abline(v = q4, col = "green")

sum(dist8 <- (dist7 > q))

cowplot::plot_grid(plotlist = lapply(1:15, function(k) {
  plot(svd.1000G, type = "scores", scores = c(2 * k - 1, 2 * k), coeff = 0.5) +
    aes(color = dist7) + scale_colour_viridis_c() +
    theme(legend.position = "none") + ggtitle(NULL)
}))


cowplot::plot_grid(plotlist = lapply(1:15, function(k) {
  plot(svd.1000G, type = "scores", scores = c(2 * k - 1, 2 * k), coeff = 0.5) +
    aes(color = dist8) + scale_colour_viridis_d() +
    theme(legend.position = "none") + ggtitle(NULL)
}))
table(fam2$Population[dist8])
# ACB ASW GBR GIH ITU LWK PEL PJL STU
#   3  12   2   4   2   1   1   1   8
# dput(unname(which(dist8)) + 0)
out <- c(19, 23, 675, 705, 938, 954, 1044, 1417, 1427, 1430, 1495, 1502,
         1516, 1517, 1548, 1565, 1567, 2108, 2249, 2252, 2264, 2275, 2277,
         2278, 2279, 2280, 2282, 2289, 2292, 2293, 2417, 2430, 2438, 2445)
out2 <- c(78, 80, 487, 1022, 1023, 1042)


ind <- sort(unique(unlist(
  lapply(1:30, function(k) {
    pc1 <- svd.1000G$u[, k]
    which(abs(pc1) > 6 * sd(pc1))
  })
)))
hist(dist7, breaks = 50)
abline(v = dist7[ind], col = "red")

cowplot::plot_grid(
  cowplot::plot_grid(plotlist = lapply(1:15, function(k) {
    plot(svd.1000G, type = "scores", scores = c(2 * k - 1, 2 * k), coeff = 0.5) +
      aes(color = rows_along(svd.1000G$u) %in% ind) + scale_colour_viridis_d() +
      theme(legend.position = "none") + ggtitle(NULL)
  }), nrow = 5),
  cowplot::plot_grid(plotlist = lapply(1:15, function(k) {
    plot(svd.1000G, type = "scores", scores = c(2 * k - 1, 2 * k), coeff = 0.5) +
      aes(color = dist8) + scale_colour_viridis_d() +
      theme(legend.position = "none") + ggtitle(NULL)
  }), nrow = 5),
  scale = 0.9, labels = c("A", "B"), label_size = 20
)

cowplot::plot_grid(plotlist = c(
  lapply(1:15, function(k) {
    plot(svd.1000G, type = "scores", scores = c(2 * k - 1, 2 * k), coeff = 0.5) +
      aes(color = dist7) + scale_colour_viridis_c() +
      theme(legend.position = "none") + ggtitle(NULL)
  }),
  list(
    plot(svd.1000G, type = "scores", scores = 19:20, coeff = 0.5) +
      aes(color = fam2$Population[-out][-out2]) +
      theme(legend.position = "none") + ggtitle(NULL)
  )
))


plot(svd.1000G, type = "scores", scores = 19:20) +
  aes(color = fam2$Population[-out][-out2]) +
  theme(legend.position = "none") + ggtitle(NULL)
