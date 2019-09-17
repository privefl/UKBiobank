loadings <- bigreadr::fread2("https://biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/snp_pca_map.txt")
for (k in 8:22) hist(loadings[[k]])

hist(S <- rowSums(loadings[8:22]^2))
median(S)
qchisq(0.5, df = 15)
S2 <- S * qchisq(0.5, df = 15) / median(S)
hist(p <- pchisq(S2, df = 15, lower.tail = FALSE))
plot(-log10(p), pch = 20, col = loadings$V1)
loadings$V2[p < 1e-40]
