library(bigsnpr)
ukbb <- snp_attach("data/UKBB_T1D.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 359 GB

# remotes::install_github("privefl/paper2-PRS/pkg.paper.PRS")
hla <- snp_indLRLDR(as.integer(ukbb$map$chromosome), ukbb$map$physical.pos,
                    LD.wiki34[12,])
y.simu <- pkg.paper.PRS::get_pheno(G, h2 = 0.4, M = 5000, K = 0.1)$pheno

ind.train <- sort(sample(nrow(G), 250e3))
ind.test <- setdiff(rows_along(G), ind.train)

system.time(
  lin <- big_spLinReg(G, y.simu[ind.train], ind.train = ind.train, K = 5,
                      ncores = nb_cores())
)
lin_pred <- predict(lin, G, ind.test)
AUCBoot(lin_pred, y.simu[ind.test]) # 82.7 [82.5-82.9]

system.time(
  log <- big_spLogReg(G, y.simu[ind.train], ind.train = ind.train, K = 5,
                      ncores = nb_cores())
)
log_pred <- predict(log, G, ind.test, proba = FALSE)
AUCBoot(log_pred, y.simu[ind.test]) # 82.7 [82.5-83.0]

plot(lin_pred, log_pred, pch = 20); abline(MASS::rlm(log_pred ~ lin_pred), col = "red")
# Useful in the tails?

hist(lin_pred, "FD")
hist(log_pred, "FD")
lin_pred2 <- glm(y.simu[ind.test] ~ lin_pred, family = "binomial")$fitted.values
log_pred2 <- 1 / (1 + exp(-log_pred))
plot(lin_pred2, log_pred2, pch = 20, col = y.simu[ind.test] + 1); abline(0, 1, col = "blue")

prev <- function(y.test, pred.test) {
  if (calib) {
    by(y.test, cut(pred.test, 0:nq / nq), mean)
  } else {
    by(y.test, cut(pred.test, quantile(pred.test, 0:nq / nq)), mean)
  }
}

nq <- 10; calib <- TRUE
prev_lin <- prev(y.simu[ind.test], lin_pred2)
prev_log <- prev(y.simu[ind.test], log_pred2)
plot(0.5:nq/nq, prev_log); points(0.5:nq/nq, prev_lin, pch = 20, col = "red")
if (calib) abline(0, 1, lty = 2, col = "blue") else abline(h = mean(y.simu), lty = 2, col = "blue")

generalhoslem::logitgof(y.simu[ind.test], lin_pred2, g = nq)
generalhoslem::logitgof(y.simu[ind.test], log_pred2, g = nq)

library(ggplot2)
qplot(lin_pred, fill = factor(y.simu[ind.test])) +
  labs(fill = "Pheno") +
  theme_bigstatsr()
qplot(log_pred, fill = factor(y.simu[ind.test])) +
  labs(fill = "Pheno") +
  theme_bigstatsr()

pAUC <- function(pred, target, p = 0.1) {
  val.min <- min(target)
  q <- quantile(pred[target == val.min], probs = 1 - p)
  ind <- (target != val.min) | (pred > q)
  bigstatsr::AUC(pred[ind], target[ind]) * p
}
pAUC(lin_pred, y.simu[ind.test], p = 0.05)
pAUC(log_pred, y.simu[ind.test], p = 0.05)

