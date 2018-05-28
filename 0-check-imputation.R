library(bigsnpr)

infos.imp <- big_attach("UKB-merged-infos-impute.rds")

cat("Imputed:", paste0(round(mean(!is.na(infos.imp[1, ])), 3) * 100, "%"))

# summary(infos.imp[2, ])
# N <- 449923
# error <- infos.imp[1, ] * infos.imp[2, ] * N
# summary(error)
# 
# # Percentage of missing values
# plot(infos.imp[1, ], pch = 20)

# Imputation error vs %missing
pvals <- c(0.01, 0.005, 0.002, 0.001); colvals <- 2:5
df <- data.frame(pNA = infos.imp[1, ], pError = infos.imp[2, ])
plot(subset(df, pNA > 0.01), pch = 20)
idc <- lapply(seq_along(pvals), function(i) {
  curve(pvals[i] / x, from = 0, lwd = 2,
        col = colvals[i], add = TRUE)
})
legend("topright", legend = pvals, title = "p(NA & Error)",
       col = colvals, lty = 1, lwd = 2)
