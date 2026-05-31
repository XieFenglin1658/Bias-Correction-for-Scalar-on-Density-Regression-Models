#!/usr/bin/env Rscript
# ==================== Generate base data including a uniform component (weight 0.3) ====================
# Used as an individual pool for subsequent simulations, containing 10000 individuals
# Save file: task2_unif_results.RData

# Set library paths (if needed)
.libPaths(c("/insomnia001/home/fx2212/miniconda3/envs/r-env/lib/R/library",
            "/insomnia001/depts/5sigma/users/fx2212/rpackages/",
            .libPaths()))

# Load packages
library(refund)
library(fdadensity)

# Fixed parameters
common_grid <- seq(0, 1, length.out = 512)
q_grid_integral <- common_grid
intercept_true <- 1.0
w_unif <- 0.3  # fixed weight for uniform component

# Define true beta function (same as before)
true_beta <- function(q) {
  mean1 <- 0.3; sd1 <- 0.1; mean2 <- 0.7; sd2 <- 0.1
  amp1 <- 0.25; amp2 <- 0.25
  raw_beta <- amp1 * dnorm(q, mean1, sd1) - amp2 * dnorm(q, mean2, sd2)
  q_grid_mean <- seq(0, 1, length.out = 1000)
  raw_grid <- amp1 * dnorm(q_grid_mean, mean1, sd1) - amp2 * dnorm(q_grid_mean, mean2, sd2)
  mean_beta <- mean(raw_grid)
  return(raw_beta - mean_beta)
}

# Generate individual parameters (same as before, but weights will be adjusted later)
set.seed(2026)
n_large <- 10000
params_list_large <- vector("list", n_large)
for(i in 1:n_large) {
  w_beta <- runif(1, 0.4, 0.7)
  params_list_large[[i]] <- list(
    alpha = runif(1, 1.8, 2.5),
    beta = runif(1, 4, 6),
    mean = runif(1, 0.6, 0.8),
    sd = runif(1, 0.08, 0.15),
    weight_beta = w_beta,
    weight_normal = 1 - w_beta
  )
}

# Function to compute true density values (including uniform component)
compute_true_density_vals <- function(params, grid, w_unif) {
  alpha <- params$alpha
  beta <- params$beta
  mean_val <- params$mean
  sd_val <- params$sd
  w_beta_orig <- params$weight_beta
  w_normal_orig <- params$weight_normal
  
  w_beta <- w_beta_orig * (1 - w_unif)
  w_normal <- w_normal_orig * (1 - w_unif)
  
  # Unnormalized density
  dens_raw <- w_beta * dbeta(grid, alpha, beta) +
    w_normal * dnorm(grid, mean_val, sd_val) +
    w_unif * dunif(grid, 0, 1)
  
  # Numerical integration (trapezoidal rule) for normalization
  h <- grid[2] - grid[1]
  norm_const <- sum(dens_raw) * h
  dens_norm <- dens_raw / norm_const
  return(dens_norm)
}

# Log-density quantile transformation function
log_density_quantile_transform <- function(density_vals, dSup = common_grid) {
  tryCatch({
    lqd <- suppressWarnings(-dens2lqd(dens = density_vals, dSup = dSup, N = 512))
    q_grid <- seq(0, 1, length.out = length(lqd))
    lqd_func <- approxfun(q_grid, lqd, rule = 2,
                          yleft = lqd[1], yright = lqd[length(lqd)])
    return(lqd_func)
  }, error = function(e) {
    density_vals <- density_vals + 1e-6
    lqd <- suppressWarnings(-dens2lqd(dens = density_vals, dSup = dSup, N = 512))
    q_grid <- seq(0, 1, length.out = length(lqd))
    lqd_func <- approxfun(q_grid, lqd, rule = 2,
                          yleft = lqd[1], yright = lqd[length(lqd)])
    return(lqd_func)
  })
}

# Compute true density values and transformation functions for all individuals
cat("Computing true densities and transformation functions...\n")
true_density_vals_large <- vector("list", n_large)
transformed_densities_true_large <- vector("list", n_large)

for (i in 1:n_large) {
  if (i %% 1000 == 0) cat("Progress:", i, "/", n_large, "\n")
  dens_vals <- compute_true_density_vals(params_list_large[[i]], common_grid, w_unif)
  true_density_vals_large[[i]] <- dens_vals
  transformed_densities_true_large[[i]] <- log_density_quantile_transform(dens_vals)
}

# Construct design matrix (true X)
X_matrix_true <- matrix(0, nrow = n_large, ncol = length(q_grid_integral))
for (i in 1:n_large) {
  X_matrix_true[i, ] <- transformed_densities_true_large[[i]](q_grid_integral)
}

# Compute true Y (no noise)
Y_true_large <- numeric(n_large)
for (i in 1:n_large) {
  beta_vals <- true_beta(q_grid_integral)
  log_density_vals <- X_matrix_true[i, ]
  h <- q_grid_integral[2] - q_grid_integral[1]
  n_points <- length(q_grid_integral)
  integral <- (h/2) * (beta_vals[1] * log_density_vals[1] +
                         beta_vals[n_points] * log_density_vals[n_points] +
                         2 * sum(beta_vals[-c(1, n_points)] * log_density_vals[-c(1, n_points)]))
  Y_true_large[i] <- intercept_true + integral
}

# Save data
save(params_list_large, true_density_vals_large,
     transformed_densities_true_large, X_matrix_true, Y_true_large,
     common_grid, q_grid_integral, w_unif,
     file = "task2_unif_results.RData")

cat("Base data generated and saved as task2_unif_results.RData\n")