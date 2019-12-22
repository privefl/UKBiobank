# Implementation of LDSC reg in R -> should verify, but based on
# https://helda.helsinki.fi/bitstream/handle/10138/273501/MT_ldsc_hautakangas.pdf
LDSC <- function(chi2, ld, M, N, THR = 30) {

  WEIGHTS <- function(w_ld, pred) {
    het_w <- 1 / (2 * pred^2)    # heteroscedasticity weights
    oc_w <- 1 / pmax(w_ld, 1)    # overcounting weights
    het_w * oc_w                 # merged weights
  }

  WEIGHT <- function(w, x) {
    w1 <- sqrt(w)
    w_norm <- w1 / sum(w1)
    x * w_norm
  }

  #### step 1 ####

  x1 <- ld[chi2 < THR]
  y1 <- chi2[chi2 < THR]

  pred0 <- y1
  for (i in 1:100) {
    w1 <- WEIGHTS(pred0, x1)
    mylm <- lm(y1 ~ x1, weights = w1)
    pred <- predict(mylm)
    if (max(abs(pred - pred0)) < 1e-4) break
    pred0 <- pred
  }
  xw1 <- WEIGHT(w1, cbind(x1, 1))
  yw1 <- WEIGHT(w1, y1)

  blocks <- bigstatsr:::CutBySize(length(x1), nb = 200)
  n_blocks <- nrow(blocks)
  xtx_block_values <- xty_block_values <- list()
  for (i in 1:n_blocks) {
    s <- bigstatsr:::seq2(blocks[i, ])
    X <- xw1[s, , drop = FALSE]
    xtx_block_values[[i]] <- crossprod(X)
    xty_block_values[[i]] <- crossprod(X, yw1[s])
  }
  xtx <- Reduce('+', xtx_block_values)
  xty <- Reduce('+', xty_block_values)
  est <- solve(xtx, xty)

  delete_values <- matrix(NA_real_, nrow(blocks), 2)
  for (i in rows_along(blocks)) {
    delete_values[i, ] <-
      solve(xtx - xtx_block_values[[i]], xty - xty_block_values[[i]])
  }
  delete_values
  pseudovalues <- sweep((1 - n_blocks) * delete_values, 2, n_blocks * est, '+')

  jknife_cov <- cov(pseudovalues) / n_blocks
  jknife_var <- diag(jknife_cov)
  jknife_se  <- sqrt(jknife_var)
  jknife_est <- colMeans(pseudovalues)

  step1_int    <- jknife_est[2]
  step1_int_se <- jknife_se[2]

  #### step 2 ####

  x <- ld
  y <- chi2
  yp <- chi2 - step1_int

  pred0 <- y
  for (i in 1:100) {
    w2 <- WEIGHTS(pred0, x)
    mylm <- lm(yp ~ x + 0, weights = w2)
    pred <- predict(mylm) + step1_int
    if (max(abs(pred - pred0)) < 1e-4) break
    pred0 <- pred
  }
  xw2 <- WEIGHT(w2, x)
  yw2 <- WEIGHT(w2, yp)

  xtx_block_values <- xty_block_values <- list()
  for (i in 1:n_blocks) {
    s <- bigstatsr:::seq2(blocks[i, ])
    X <- xw2[s]
    xtx_block_values[[i]] <- crossprod(X)
    xty_block_values[[i]] <- crossprod(X, yw2[s])
  }
  xtx <- Reduce('+', xtx_block_values)
  xty <- Reduce('+', xty_block_values)
  est <- drop(solve(xtx, xty))

  delete_values <- rep(NA_real_, nrow(blocks))
  for (i in rows_along(blocks)) {
    delete_values[i] <-
      solve(xtx - xtx_block_values[[i]], xty - xty_block_values[[i]])
  }
  delete_values
  pseudovalues <- n_blocks * est - (n_blocks - 1) * delete_values

  jknife_cov <- cov(as.matrix(pseudovalues)) / n_blocks
  jknife_var <- diag(jknife_cov)
  jknife_se <- sqrt(jknife_var)
  jknife_est <- mean(pseudovalues)

  c(step1_int, step1_int_se, M / N * jknife_est, M / N * jknife_se)
}

################################################################################

# Another LD part (population structure)
LDSC2 <- function(chi2, ld, ld_part, M, N, THR = 30) {

  WEIGHTS <- function(w_ld, pred) {
    het_w <- 1 / (2 * pred^2)    # heteroscedasticity weights
    oc_w <- 1 / pmax(w_ld, 1)    # overcounting weights
    het_w * oc_w                 # merged weights
  }

  WEIGHT <- function(w, x) {
    w1 <- sqrt(w)
    w_norm <- w1 / sum(w1)
    x * w_norm
  }

  #### step 1 ####

  x1 <- ld[chi2 < THR]
  x1_part <- ld_part[chi2 < THR]
  y1 <- chi2[chi2 < THR]

  pred0 <- y1
  for (i in 1:100) {
    w1 <- WEIGHTS(pred0, x1)
    mylm <- lm(y1 ~ x1 + x1_part, weights = w1)
    pred <- predict(mylm)
    if (max(abs(pred - pred0)) < 1e-4) break
    pred0 <- pred
  }
  xw1 <- WEIGHT(w1, cbind(x1, 1, x1_part))
  yw1 <- WEIGHT(w1, y1)

  blocks <- bigstatsr:::CutBySize(length(x1), nb = 200)
  n_blocks <- nrow(blocks)
  xtx_block_values <- xty_block_values <- list()
  for (i in 1:n_blocks) {
    s <- bigstatsr:::seq2(blocks[i, ])
    X <- xw1[s, , drop = FALSE]
    xtx_block_values[[i]] <- crossprod(X)
    xty_block_values[[i]] <- crossprod(X, yw1[s])
  }
  xtx <- Reduce('+', xtx_block_values)
  xty <- Reduce('+', xty_block_values)
  est <- solve(xtx, xty)

  delete_values <- matrix(NA_real_, nrow(blocks), 3)
  for (i in rows_along(blocks)) {
    delete_values[i, ] <-
      solve(xtx - xtx_block_values[[i]], xty - xty_block_values[[i]])
  }
  delete_values
  pseudovalues <- sweep((1 - n_blocks) * delete_values, 2, n_blocks * est, '+')

  jknife_cov <- cov(pseudovalues) / n_blocks
  jknife_var <- diag(jknife_cov)
  jknife_se  <- sqrt(jknife_var)
  jknife_est <- colMeans(pseudovalues)

  step1_int    <- jknife_est[2]
  step1_int_se <- jknife_se[2]

  #### step 2 ####

  x <- ld
  x_part <- ld_part
  y <- chi2
  yp <- chi2 - step1_int
  pred <- yp

  pred0 <- y
  for (i in 1:100) {
    w2 <- WEIGHTS(pred0, x)
    mylm <- lm(yp ~ x + ld_part + 0, weights = w2)
    pred <- predict(mylm) + step1_int
    if (max(abs(pred - pred0)) < 1e-4) break
    pred0 <- pred
  }
  xw2 <- WEIGHT(w2, cbind(x, x_part))
  yw2 <- WEIGHT(w2, yp)

  xtx_block_values <- xty_block_values <- list()
  for (i in 1:n_blocks) {
    s <- bigstatsr:::seq2(blocks[i, ])
    X <- xw2[s, , drop = FALSE]
    xtx_block_values[[i]] <- crossprod(X)
    xty_block_values[[i]] <- crossprod(X, yw2[s])
  }
  xtx <- Reduce('+', xtx_block_values)
  xty <- Reduce('+', xty_block_values)
  est <- drop(solve(xtx, xty))

  delete_values <- matrix(NA_real_, nrow(blocks), 2)
  for (i in rows_along(blocks)) {
    delete_values[i, ] <-
      solve(xtx - xtx_block_values[[i]], xty - xty_block_values[[i]])
  }
  delete_values
  pseudovalues <- sweep((1 - n_blocks) * delete_values, 2, n_blocks * est, '+')

  jknife_cov <- cov(pseudovalues) / n_blocks
  jknife_var <- diag(jknife_cov)
  jknife_se  <- sqrt(jknife_var)
  jknife_est <- colMeans(pseudovalues)

  c(step1_int, step1_int_se, M / N * jknife_est[1], M / N * jknife_se[1])
}


################################################################################

# Another LD part (population structure)
LDSC3 <- function(chi2, ld, ld_part, M, N, THR = 30) {

  WEIGHTS <- function(w_ld, pred) {
    het_w <- 1 / (2 * pred^2)    # heteroscedasticity weights
    oc_w <- 1 / pmax(w_ld, 1)    # overcounting weights
    het_w * oc_w                 # merged weights
  }

  WEIGHT <- function(w, x) {
    w1 <- sqrt(w)
    w_norm <- w1 / sum(w1)
    x * w_norm
  }

  #### step 1 ####

  x1 <- ld[chi2 < THR]
  x1_part <- ld_part[chi2 < THR]
  y1 <- chi2[chi2 < THR]

  pred0 <- y1
  for (i in 1:100) {
    w1 <- WEIGHTS(pred0, x1)
    mylm <- lm(y1 ~ x1 + x1_part, weights = w1)
    pred <- predict(mylm)
    if (max(abs(pred - pred0)) < 1e-4) break
    pred0 <- pred
  }
  xw1 <- WEIGHT(w1, cbind(x1, 1, x1_part))
  yw1 <- WEIGHT(w1, y1)

  blocks <- bigstatsr:::CutBySize(length(x1), nb = 200)
  n_blocks <- nrow(blocks)
  xtx_block_values <- xty_block_values <- list()
  for (i in 1:n_blocks) {
    s <- bigstatsr:::seq2(blocks[i, ])
    X <- xw1[s, , drop = FALSE]
    xtx_block_values[[i]] <- crossprod(X)
    xty_block_values[[i]] <- crossprod(X, yw1[s])
  }
  xtx <- Reduce('+', xtx_block_values)
  xty <- Reduce('+', xty_block_values)
  est <- solve(xtx, xty)

  delete_values <- matrix(NA_real_, nrow(blocks), 3)
  for (i in rows_along(blocks)) {
    delete_values[i, ] <-
      solve(xtx - xtx_block_values[[i]], xty - xty_block_values[[i]])
  }
  delete_values
  pseudovalues <- sweep((1 - n_blocks) * delete_values, 2, n_blocks * est, '+')

  jknife_cov <- cov(pseudovalues) / n_blocks
  jknife_var <- diag(jknife_cov)
  jknife_se  <- sqrt(jknife_var)
  jknife_est <- colMeans(pseudovalues)

  c(jknife_est[2], jknife_se[2], M / N * jknife_est[1], M / N * jknife_se[1])
}
################################################################################

# Another LD part (population structure)
LDSC4 <- function(chi2, ld, ld_part, M, N, THR = 30) {

  WEIGHTS <- function(w_ld, pred) {
    het_w <- 1 / (2 * pred^2)    # heteroscedasticity weights
    oc_w <- 1 / pmax(w_ld, 1)    # overcounting weights
    het_w * oc_w                 # merged weights
  }

  WEIGHT <- function(w, x) {
    w1 <- sqrt(w)
    w_norm <- w1 / sum(w1)
    x * w_norm
  }

  #### step 1 ####

  x1 <- ld[chi2 < THR]
  x1_part <- ld_part[chi2 < THR]
  y1 <- chi2[chi2 < THR]

  pred0 <- y1
  for (i in 1:100) {
    w1 <- WEIGHTS(pred0, x1)
    mylm <- lm(y1 ~ x1 + x1_part, weights = w1)
    pred <- predict(mylm)
    if (max(abs(pred - pred0)) < 1e-4) break
    pred0 <- pred
  }
  xw1 <- WEIGHT(w1, cbind(x1, 1, x1_part))
  yw1 <- WEIGHT(w1, y1)

  blocks <- bigstatsr:::CutBySize(length(x1), nb = 200)
  n_blocks <- nrow(blocks)
  xtx_block_values <- xty_block_values <- list()
  for (i in 1:n_blocks) {
    s <- bigstatsr:::seq2(blocks[i, ])
    X <- xw1[s, , drop = FALSE]
    xtx_block_values[[i]] <- crossprod(X)
    xty_block_values[[i]] <- crossprod(X, yw1[s])
  }
  xtx <- Reduce('+', xtx_block_values)
  xty <- Reduce('+', xty_block_values)
  est <- solve(xtx, xty)

  delete_values <- matrix(NA_real_, nrow(blocks), 3)
  for (i in rows_along(blocks)) {
    delete_values[i, ] <-
      solve(xtx - xtx_block_values[[i]], xty - xty_block_values[[i]])
  }
  delete_values
  pseudovalues <- sweep((1 - n_blocks) * delete_values, 2, n_blocks * est, '+')

  jknife_cov <- cov(pseudovalues) / n_blocks
  jknife_var <- diag(jknife_cov)
  jknife_se  <- sqrt(jknife_var)
  jknife_est <- colMeans(pseudovalues)

  step1_int    <- jknife_est[2]
  step1_int_se <- jknife_se[2]
  step1_x1_part <- jknife_est[3]

  #### step 2 ####

  x <- ld
  x_part <- ld_part
  y <- chi2
  yp <- chi2 - (step1_int + x_part * step1_x1_part)
  pred <- yp

  pred0 <- y
  for (i in 1:100) {
    w2 <- WEIGHTS(pred0, x)
    mylm <- lm(yp ~ x + 0, weights = w2)
    pred <- predict(mylm) + step1_int + x_part * step1_x1_part
    if (max(abs(pred - pred0)) < 1e-4) break
    pred0 <- pred
  }
  xw2 <- WEIGHT(w2, x)
  yw2 <- WEIGHT(w2, yp)

  xtx_block_values <- xty_block_values <- list()
  for (i in 1:n_blocks) {
    s <- bigstatsr:::seq2(blocks[i, ])
    X <- xw2[s]
    xtx_block_values[[i]] <- crossprod(X)
    xty_block_values[[i]] <- crossprod(X, yw2[s])
  }
  xtx <- Reduce('+', xtx_block_values)
  xty <- Reduce('+', xty_block_values)
  est <- drop(solve(xtx, xty))

  delete_values <- rep(NA_real_, nrow(blocks))
  for (i in rows_along(blocks)) {
    delete_values[i] <-
      solve(xtx - xtx_block_values[[i]], xty - xty_block_values[[i]])
  }
  delete_values
  pseudovalues <- n_blocks * est - (n_blocks - 1) * delete_values

  jknife_cov <- cov(as.matrix(pseudovalues)) / n_blocks
  jknife_var <- diag(jknife_cov)
  jknife_se <- sqrt(jknife_var)
  jknife_est <- mean(pseudovalues)

  c(step1_int, step1_int_se, M / N * jknife_est, M / N * jknife_se)
}

