library(bigsnpr)

ukb_imp <- snp_attach("data/UKB_imp_BC_save.rds")
G <- ukb_imp$genotypes
CHR <- as.integer(ukb_imp$map$chromosome)
POS <- ukb_imp$map$physical.pos

grid <- expand.grid(
  thr.imp = c(0.3, 0.6, 0.8, 0.9, 0.95),
  thr.clmp = c(0.05, 0.2, 0.5, 0.8),
  chr = unique(CHR)
)

n_thr_pval <- 20
scores <- big_attach("breast-cancer/scores.rds")

hist(scores[, ncol(scores)])

y <- ukb_imp$fam$has_cancer
y[y] <- NA  ## set to missing all types of cancer
y[ukb_imp$fam$cancer_type == "breast cancer"] <- 1  ## keep only BC
ind <- which(ukb_imp$fam$sex == "Female" & !is.na(y))  ## only women

ind.test2 <- intersect(ind.test, ind)
aucs <- big_apply(scores, a.FUN = function(X, ind, ind.test, y.test) {
  apply(X[ind.test, ind], 2, bigstatsr::AUC, y.test)
}, a.combine = 'c', ind.test = ind.test2, y.test = y[ind.test2], ncores = 12)

ord <- unlist(lapply(order(grid$chr), function(ic) {
  (ic - 1) * n_thr_pval + 1:n_thr_pval
}))
plot(aucs[ord], col = sort(rep(grid$chr, n_thr_pval)), pch = 20)

ind.train2 <- intersect(ind.train, ind)
table(y[ind.train2])
system.time(
  prs_final <- big_spLogReg(scores, y[ind.train2], ind.train2,
                            ncores = 12, alphas = 10^(-(0:4)))
) # 17 min

str(prs_final)
coefs <- rowMeans(sapply(prs_final, function(mod) mod$beta.X))
sum(coefs < 0)  # 0

grid2 <- grid %>% mutate(num = list(seq_len(n_thr_pval))) %>% tidyr::unnest() %>%
  mutate(coef = coefs, used = coefs != 0)

library(ggplot2)
grid2 %>%
  group_by(thr.imp, thr.clmp) %>%
  summarise(m_coef = mean(coef), p_used = mean(used)) %>%
  ggplot(aes(thr.clmp, m_coef, color = as.factor(thr.imp))) +
  geom_line(size = 2) +
  theme_bigstatsr()

grid2 %>%
  group_by(thr.imp, thr.clmp) %>%
  summarise(m_coef = mean(coef), p_used = mean(used)) %>%
  ggplot(aes(thr.clmp, p_used, color = as.factor(thr.imp))) +
  geom_line(size = 2) +
  theme_bigstatsr()

preds_final <- predict(prs_final, scores, ind.test2)
AUC(preds_final, y[ind.test2])


size <- nrow(grid) / 22 * n_thr_pval
scores2 <- scores[ind.test2, 1:size]
for (i in 1:21) inplace::`%+<-%`(scores2, scores[ind.test2, 1:size + i * size])
aucs2 <- apply(scores2, 2, AUC, y[ind.test2])

grid %>% mutate(num = list(seq_len(n_thr_pval))) %>% tidyr::unnest() %>%
  tidyr::nest(chr) %>%
  mutate(auc = aucs2) %>%
  ggplot() +
  geom_line(aes(num, auc, color = as.factor(thr.imp)), size = 1) +
  facet_wrap(~ thr.clmp) +
  geom_hline(yintercept = AUC(preds_final, y[ind.test2]), linetype = 2) +
  theme_bigstatsr() +
  xlim(5, NA) +
  labs(x = "p-value threshold number", y = "AUC", color = "Imputation\nthreshold")

grid %>% mutate(num = list(seq_len(n_thr_pval))) %>% tidyr::unnest() %>%
  tidyr::nest(chr) %>%
  mutate(auc = aucs2) %>%
  ggplot() +
  geom_line(aes(num, auc, color = as.factor(thr.clmp)), size = 1) +
  facet_wrap(~ thr.imp) +
  geom_hline(yintercept = AUC(preds_final, y[ind.test2]), linetype = 2) +
  theme_bigstatsr() +
  xlim(5, NA) +
  labs(x = "p-value threshold number", y = "AUC", color = "Clumping\nthreshold") +
  theme(legend.position = c(0.85, 0.2))

library(ggplot2)
qplot(preds_final, fill = as.factor(y[ind.test2]), geom = "density", alpha = I(0.4))
