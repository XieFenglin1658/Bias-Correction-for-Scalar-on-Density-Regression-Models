#!/usr/bin/env Rscript
# run_simex_single.R
# Modification: No longer select the best extrapolation method, but save results for three extrapolation methods (without cheating)
# Modification: Use Tlqd(obs, method="spline") instead of KDE + log_quantile_transform
# Changed k=10 back to k=20 because real data is not suitable for k=10

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_simex_single.R <n> <sigma> <rep>")
}
n_val <- as.numeric(args[1])
sigma_val <- as.numeric(args[2])
rep_id <- as.integer(args[3])

seed <- 2026 + n_val * 1000 + round(sigma_val * 1000) * 10 + rep_id
set.seed(seed)

.libPaths(c("/insomnia001/home/fx2212/miniconda3/envs/fda-43/lib/R/library",
            "/insomnia001/home/fx2212/R/library"))

library(refund)
library(fdadensity)
library(mgcv)
library(fda)
source("utils.R")

load("task2_unif_results.RData")

q_grid <- common_grid
q_eval <- q_grid
w_unif <- 0.3
m0_values <- c(10, 20, 50, 100, 200, 500, 1000)
B_reps <- 100

true_beta <- function(q) {
  mean1 <- 0.3; sd1 <- 0.1; mean2 <- 0.7; sd2 <- 0.1
  amp1 <- 0.25; amp2 <- 0.25
  raw_beta <- amp1 * dnorm(q, mean1, sd1) - amp2 * dnorm(q, mean2, sd2)
  q_grid_mean <- seq(0, 1, length.out = 1000)
  raw_grid <- amp1 * dnorm(q_grid_mean, mean1, sd1) - amp2 * dnorm(q_grid_mean, mean2, sd2)
  mean_beta <- mean(raw_grid)
  return(raw_beta - mean_beta)
}
true_vals <- true_beta(q_eval)

generate_mixture_observations <- function(params, n_obs, w_unif) {
  alpha <- params$alpha; beta <- params$beta; mean_val <- params$mean
  sd_val <- params$sd; w_beta_orig <- params$weight_beta; w_normal_orig <- params$weight_normal
  w_beta <- w_beta_orig * (1 - w_unif)
  w_normal <- w_normal_orig * (1 - w_unif)
  component <- sample(c("beta", "normal", "unif"), n_obs, replace = TRUE,
                      prob = c(w_beta, w_normal, w_unif))
  observations <- numeric(n_obs)
  n_beta <- sum(component == "beta"); n_normal <- sum(component == "normal"); n_unif <- sum(component == "unif")
  if (n_beta > 0) observations[component == "beta"] <- rbeta(n_beta, alpha, beta)
  if (n_normal > 0) {
    normal_obs <- rnorm(n_normal, mean_val, sd_val)
    observations[component == "normal"] <- pmin(pmax(normal_obs, 0), 1)
  }
  if (n_unif > 0) observations[component == "unif"] <- runif(n_unif, 0, 1)
  return(observations)
}

run_one_dataset <- function(seed, n_fixed, sigma_val, 
                            params_pool, Y_true_pool, 
                            m0_values, B_reps) {
  set.seed(seed)
  n_large <- length(params_pool)
  selected_idx <- sample(1:n_large, n_fixed, replace = FALSE)
  Y_true <- Y_true_pool[selected_idx]
  
  original_obs_list <- vector("list", n_fixed)
  for (i in 1:n_fixed) {
    original_obs_list[[i]] <- generate_mixture_observations(params_pool[[selected_idx[i]]], 1000, w_unif)
  }
  
  beta_m0_matrix <- matrix(NA, nrow = length(q_grid), ncol = length(m0_values))
  colnames(beta_m0_matrix) <- paste0("m0=", m0_values)
  
  for (j in seq_along(m0_values)) {
    m0 <- m0_values[j]
    beta_reps <- matrix(NA, nrow = length(q_grid), ncol = B_reps)
    for (b in 1:B_reps) {
      set.seed(seed + m0 * 1000 + b)
      X_est <- matrix(0, nrow = n_fixed, ncol = length(q_grid))
      for (i in 1:n_fixed) {
        boot_idx <- sample(1:1000, size = m0, replace = TRUE)
        obs <- original_obs_list[[i]][boot_idx]
        trans <- Tlqd(obs, method = "spline")
        X_est[i, ] <- trans$y
      }
      Y_obs <- Y_true + rnorm(n_fixed, 0, sigma_val)
      pfr_data <- data.frame(Y = Y_obs, X = I(X_est))
      fit <- tryCatch(
        pfr(Y ~ lf(X, argvals = q_grid, k = 20, bs = "ps"), data = pfr_data),
        error = function(e) NULL
      )
      if (!is.null(fit)) {
        beta_est <- coef(fit, n = length(q_grid))
        beta_reps[, b] <- approx(beta_est$X.argvals, beta_est$value, xout = q_eval, rule = 2)$y
      }
    }
    beta_m0_matrix[, j] <- rowMeans(beta_reps, na.rm = TRUE)
  }
  
  mse_m0 <- apply(beta_m0_matrix, 2, function(curve) mean((curve - true_vals)^2, na.rm = TRUE))
  
  # ----- Three extrapolation methods (all saved) -----
  lambda <- 1 / sqrt(m0_values)
  
  # Linear
  beta_lin <- numeric(length(q_grid))
  for (i in 1:length(q_grid)) {
    y <- beta_m0_matrix[i, ]
    fit <- lm(y ~ lambda)
    beta_lin[i] <- predict(fit, newdata = data.frame(lambda = 0))
  }
  mse_lin <- mean((beta_lin - true_vals)^2)
  
  # Quadratic
  beta_quad <- numeric(length(q_grid))
  for (i in 1:length(q_grid)) {
    y <- beta_m0_matrix[i, ]
    fit <- tryCatch(
      lm(y ~ lambda + I(lambda^2)),
      error = function(e) lm(y ~ lambda)
    )
    if (any(is.na(coef(fit)))) fit <- lm(y ~ lambda)
    beta_quad[i] <- predict(fit, newdata = data.frame(lambda = 0))
  }
  mse_quad <- mean((beta_quad - true_vals)^2)
  
  # Nonlinear
  df_nonlin <- expand.grid(q = q_grid, lambda = lambda)
  df_nonlin$beta <- as.vector(beta_m0_matrix)
  gam_fit <- mgcv::gam(beta ~ te(q, lambda, k = c(10,5)), data = df_nonlin, method = "REML")
  beta_nonlin <- predict(gam_fit, newdata = data.frame(q = q_grid, lambda = 0))
  mse_nonlin <- mean((beta_nonlin - true_vals)^2)
  
  list(
    mse_m0 = mse_m0,
    mse_lin = mse_lin,
    mse_quad = mse_quad,
    mse_nonlin = mse_nonlin,
    beta_m0_matrix = beta_m0_matrix,
    beta_lin = beta_lin,
    beta_quad = beta_quad,
    beta_nonlin = beta_nonlin,
    q_grid = q_grid
  )
}

res <- run_one_dataset(seed, n_val, sigma_val, 
                       params_list_large, Y_true_large, 
                       m0_values, B_reps)

outfile <- sprintf("result_n%d_sigma%.3f_rep%d.rds", n_val, sigma_val, rep_id)
saveRDS(res, file = outfile)