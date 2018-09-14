library(bigsnpr)
ukb <- snp_attach("UKB-merged_sub1.rds")
G <- ukb$genotypes
CHR <- ukb$map$chromosome
POS <- ukb$map$physical.pos
dim(G) # 449926 x 564148
# big_counts(G)[4, ] # OK
(NCORES <- nb_cores())

ind.rm <- snp_indLRLDR(CHR, POS)
system.time(
  ind.keep <- snp_clumping(G, CHR, thr.r2 = 0.1, exclude = ind.rm,
                           ncores = NCORES)
) # 195K
# utilisateur     système      écoulé
#   27488.856     220.668   29665.410
#    user   system  elapsed
# 402.985  203.360 6352.371
saveRDS(ind.keep, "keep_pruning.rds")

system.time(
  obj.svd <- big_randomSVD(G, snp_scaleBinom(), ind.col = ind.keep,
                           verbose = TRUE, ncores = NCORES)
)
#  user   system  elapsed
# 5.275   10.437 1700.229
# obj.svd <- snp_autoSVD(G, CHR, POS, thr.r2 = 0.1, ncores = 1)

saveRDS(obj.svd, "UKB_SVD.rds")
plot(obj.svd)  #  K = 5
# Compare with https://biobank.ctsu.ox.ac.uk/crystal/docs/genotyping_qc.pdf (p11)
plot(obj.svd, type = "scores")
plot(obj.svd, type = "scores", scores = 3:4)
plot(obj.svd, type = "scores", scores = 5:6)
plot(obj.svd, type = "loadings", loadings = 1:10, coeff = 0.5)
ggplot2::ggsave("tmp.png")

# loadings <- bigreadr::fread2("snp_pca_map.txt")
# match(loadings$V2, ukb$map$marker.ID)  ## not all in common

obj.pcadapt <- snp_pcadapt(G, obj.svd$u[, 1:5])
saveRDS(obj.pcadapt, "UKB_pcadapt.rds")
plot(obj.pcadapt)
ggplot2::ggsave("tmp.png")
snp_qq(obj.pcadapt)
ggplot2::ggsave("tmp.png")
obj.pcadapt.gc <- snp_gc(obj.pcadapt)
plot(obj.pcadapt.gc)
snp_manhattan(obj.pcadapt.gc, CHR, POS, npoints = 20e3) +
  ggplot2::geom_hline(yintercept = -log10(5e-8), color = "red")
ggplot2::ggsave("tmp.png")
