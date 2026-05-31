#!/usr/bin/env Rscript
# Receive command line arguments: m_val, reps
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_unif_lf.R <m_val> <reps>")
}
m_val <- as.numeric(args[1])
reps <- as.numeric(args[2])

# Set library paths
.libPaths(c("/insomnia001/home/fx2212/miniconda3/envs/fda-43/lib/R/library",
            "/insomnia001/home/fx2212/R/library"))

# Load packages
library(refund)
library(fdadensity)
library(foreach)
library(doParallel)
library(doRNG)
library(fda)

# Load T0_Tlqd_function.R, which contains the Tlqd function
source("T0_Tlqd_function.R")

# Fixed parameters
n_fixed <- 100
sigma_fixed <- 0.01
use_true_Y <- FALSE
w_unif <- 0.3

# Load base data
load("task2_unif_results.RData")  # contains params_list_large, Y_true_large, common_grid, q_grid_integral,
# transformed_densities_true_large, etc.
common_grid <- common_grid
q_grid_integral <- q_grid_integral
q_eval_beta <- common_grid
n_large_available <- length(params_list_large)

# Define true beta function (consistent with data generation)
true_beta <- function(q) {
  mean1 <- 0.3; sd1 <- 0.1; mean2 <- 0.7; sd2 <- 0.1
  amp1 <- 0.25; amp2 <- 0.25
  raw_beta <- amp1 * dnorm(q, mean1, sd1) - amp2 * dnorm(q, mean2, sd2)
  q_grid_mean <- seq(0, 1, length.out = 1000)
  raw_grid <- amp1 * dnorm(q_grid_mean, mean1, sd1) - amp2 * dnorm(q_grid_mean, mean2, sd2)
  mean_beta <- mean(raw_grid)
  return(raw_beta - mean_beta)
}

# Function to generate mixture observations (including uniform component)
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

# Single replicate function
run_one_rep <- function(rep) {
  library(fda)
  set.seed(round(m_val * 1000 + sigma_fixed * 10 + rep * 10000))
  
  # Randomly sample n_fixed individuals from the large pool
  selected_idx <- sample(1:n_large_available, n_fixed, replace = FALSE)
  
  # Estimate density and transformation function for each individual (using Tlqd)
  X_est <- matrix(0, nrow = n_fixed, ncol = length(q_grid_integral))
  for (k in 1:n_fixed) {
    i <- selected_idx[k]
    obs <- generate_mixture_observations(params_list_large[[i]], m_val, w_unif)
    # Directly call Tlqd to obtain transformed function values
    trans <- Tlqd(obs)   # returns list(x, y), where y are the transformed values
    X_est[k, ] <- trans$y
  }
  
  # Obtain true transformation function values
  X_true <- matrix(0, nrow = n_fixed, ncol = length(q_grid_integral))
  for (k in 1:n_fixed) {
    i <- selected_idx[k]
    X_true[k, ] <- transformed_densities_true_large[[i]](q_grid_integral)
  }
  
  # Save X information to file (for later plotting)
  saveRDS(list(
    m = m_val,
    rep = rep,
    selected_idx = selected_idx,
    X_est = X_est,
    X_true = X_true,
    q_grid = q_grid_integral
  ), file = paste0("X_m", m_val, "_rep", rep, ".rds"))
  
  # Response variable (true Y plus noise)
  Y_true_sel <- Y_true_large[selected_idx]
  Y_obs <- Y_true_sel + rnorm(n_fixed, 0, sigma_fixed)
  
  # Fit pfr using lf()
  pfr_data <- data.frame(Y = Y_obs, X = I(X_est))
  fit <- tryCatch(
    pfr(Y ~ lf(X, argvals = q_grid_integral, k = 20, bs = "ps"),
        data = pfr_data),
    error = function(e) NULL
  )
  
  if (!is.null(fit)) {
    beta_est <- coef(fit, n = length(q_grid_integral))
    beta_est_vals <- approx(beta_est$X.argvals, beta_est$value,
                            xout = q_eval_beta, rule = 2)$y
    return(beta_est_vals)
  } else {
    return(rep(NA_real_, length(q_eval_beta)))
  }
}

# Parallel execution (same as original script)
num_cores <- 32
cl <- makeCluster(num_cores - 2)
registerDoParallel(cl)
registerDoRNG(123456)

cat(sprintf("Running for m = %d, n = %d, sigma = %.3f, reps = %d\n",
            m_val, n_fixed, sigma_fixed, reps))

rep_results <- foreach(rep = 1:reps, .combine = 'cbind',
                       .packages = c("refund", "fdadensity")) %dopar% {
                         run_one_rep(rep)
                       }

stopCluster(cl)

# Compute statistics
avg_curve <- rowMeans(rep_results, na.rm = TRUE)
true_beta_vals <- true_beta(q_eval_beta)
bias_curve <- avg_curve - true_beta_vals
squared_errors <- (rep_results - true_beta_vals)^2
mse_curve <- rowMeans(squared_errors, na.rm = TRUE)
variance_curve <- apply(rep_results, 1, var, na.rm = TRUE)

mse_value <- mean(mse_curve, na.rm = TRUE)
bias_sq_value <- mean(bias_curve^2, na.rm = TRUE)
variance_value <- mean(variance_curve, na.rm = TRUE)

# Save results
out_prefix <- paste0("unif_n", n_fixed, "_m", m_val, "_sigma", gsub("\\.", "_", sigma_fixed))
save(rep_results, avg_curve, bias_curve, mse_curve, variance_curve,
     mse_value, bias_sq_value, variance_value,
     file = paste0(out_prefix, "_results.RData"))
