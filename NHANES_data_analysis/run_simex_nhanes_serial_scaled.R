#!/usr/bin/env Rscript
# run_simex_nhanes_serial_scaled.R
# Usage: Rscript run_simex_nhanes_serial_scaled.R <m0> <B_reps> <start_rep>
# Performs B_reps bootstrap simulations for a single m0, starting from start_rep (serial)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_simex_nhanes_serial_scaled.R <m0> <B_reps> <start_rep>")
}
m0 <- as.numeric(args[1])
B_reps <- as.numeric(args[2])
start_rep <- as.numeric(args[3])

.libPaths(c("/insomnia001/home/fx2212/miniconda3/envs/fda-43/lib/R/library",
            "/insomnia001/home/fx2212/R/library"))

library(refund)
library(fdadensity)
library(mgcv)
library(fda)

# Load scaled data
simex_data <- readRDS("nhanes_simex_data_scaled.rds")
original_conn <- simex_data$original_conn   # list of scaled MIMS vectors
covariates <- simex_data$covariates
global_max <- simex_data$global_max

# Convert categorical variables to numeric
covariates$gender <- as.numeric(covariates$gender == "Male")
covariates$CHD <- as.numeric(covariates$CHD == "Yes")
covariates$event <- as.numeric(covariates$event)

# Remove missing values
complete_idx <- complete.cases(covariates[, c("time", "event", "age", "BMI", "gender", "CHD")])
covariates <- covariates[complete_idx, ]
original_conn <- original_conn[complete_idx]

n_individuals <- nrow(covariates)
total_events <- sum(covariates$event)
cat("Total sample size:", n_individuals, ", number of events:", total_events, "\n")
cat("Global maximum:", global_max, "\n")
cat("m0 =", m0, ", B_reps =", B_reps, ", start_rep =", start_rep, "\n")

# Define T0 function (density estimation, fixed domain [-0.1, 1])
T0_fixed <- function(corr_i, nb = 19) {
  dens.eps <- 1e-6
  domain <- c(-0.1, 1)
  basis1 <- fda::create.bspline.basis(domain, nb)
  densf <- fda::density.fd(corr_i, fda::fd(matrix(0, nb, 1), basis1))
  x <- seq(domain[1], domain[2], length.out = 1024)
  y <- c(exp(fda::eval.fd(x, densf$Wfdobj)) / densf$C)
  ind1 <- min(which(y > dens.eps * 512))
  ind2 <- max(which(y > dens.eps * 512))
  list(x = x[ind1:ind2], y = y[ind1:ind2])
}

# Define Tlqd function (density + LQD transformation, fixed domain [-0.1, 1])
Tlqd_fixed <- function(corr_i, nb = 19) {
  dens.eps <- 1e-6
  domain <- c(-0.1, 1)
  basis1 <- fda::create.bspline.basis(domain, nb)
  densf <- fda::density.fd(corr_i, fda::fd(matrix(0, nb, 1), basis1))
  x <- seq(domain[1], domain[2], length.out = 1024)
  y <- c(exp(fda::eval.fd(x, densf$Wfdobj)) / densf$C)
  ind1 <- min(which(y > dens.eps * 512))
  ind2 <- max(which(y > dens.eps * 512))
  di <- list(x = x[ind1:ind2], y = y[ind1:ind2])
  tryCatch({
    y_lqd <- suppressWarnings(-fdadensity::dens2lqd(di$y, di$x, N = 512))
    x_lqd <- seq(0, 1, length.out = 512)
    return(list(x = x_lqd, y = y_lqd))
  }, error = function(e) {
    di$y <- di$y + dens.eps
    y_lqd <- suppressWarnings(-fdadensity::dens2lqd(di$y, di$x, N = 512))
    x_lqd <- seq(0, 1, length.out = 512)
    return(list(x = x_lqd, y = y_lqd))
  })
}

# Result matrix
beta_reps <- matrix(NA, nrow = 512, ncol = B_reps)

for (rep in 1:B_reps) {
  global_rep <- start_rep + rep - 1
  cat("Processing global rep", global_rep, "\n")
  set.seed(12345 + m0 * 1000 + global_rep)
  
  X_est <- matrix(0, n_individuals, 512)
  for (i in 1:n_individuals) {
    obs <- sample(original_conn[[i]], size = m0, replace = TRUE)
    trans <- Tlqd_fixed(obs)
    X_est[i, ] <- trans$y
  }
  
  if (any(is.na(X_est)) || any(is.infinite(X_est))) {
    cat("Rep", rep, "X_est contains NA/Inf\n")
    next
  }
  
  df <- data.frame(
    time = covariates$time,
    event = covariates$event,
    age = covariates$age,
    BMI = covariates$BMI,
    gender = covariates$gender,
    CHD = covariates$CHD,
    X = I(X_est)
  )
  
  n_events <- sum(df$event)
  if (n_events < 10) {
    cat("Rep", rep, "insufficient events:", n_events, "\n")
    next
  }
  
  fit <- tryCatch(
    pfr(time ~ age + BMI + gender + CHD + lf(X, bs = "ps", k = 20, argvals = seq(0,1,len=512)),
        method = "REML", data = df, weights = df$event, family = cox.ph),
    error = function(e) {
      cat("Rep", rep, "pfr error:", e$message, "\n")
      return(NULL)
    }
  )
  
  if (!is.null(fit)) {
    beta_est <- coef(fit, n = 512, select = 1)
    beta_reps[, rep] <- beta_est$value
    cat("Rep", rep, "successful\n")
  } else {
    cat("Rep", rep, "failed\n")
  }
}

success_count <- sum(apply(beta_reps, 2, function(x) !any(is.na(x))))
cat(sprintf("Number of successful replicates: %d / %d (%.2f%%)\n", success_count, B_reps, 100 * success_count / B_reps))

avg_beta <- rowMeans(beta_reps, na.rm = TRUE)
outfile <- paste0("simex_m0_", m0, "_B", B_reps, "_start", start_rep, ".rds")
saveRDS(list(m0 = m0, avg_beta = avg_beta, all_beta = beta_reps), file = outfile)