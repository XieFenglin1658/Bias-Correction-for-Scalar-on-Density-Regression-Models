# Generate scaled naive curve (using all data)
library(refund)
library(mgcv)
library(fda)
library(tidyverse)

# Load scaled data
setwd("/Users/fenglinxie/Desktop/simulation_lf_uniform_Tlqd_SIMEX_k20/NHANES_cluster_serial/data_scaled_Todd")

simex_data <- readRDS("nhanes_simex_data_scaled.rds")
original_conn <- simex_data$original_conn
covariates <- simex_data$covariates
covariates$gender <- as.numeric(covariates$gender == "Male")
covariates$CHD <- as.numeric(covariates$CHD == "Yes")
covariates$event <- as.numeric(covariates$event)

#Tlqd function
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


X_full <- matrix(0, nrow = nrow(covariates), ncol = 512)
for (i in 1:nrow(covariates)) {
  trans <- Tlqd_fixed(original_conn[[i]])
  X_full[i, ] <- trans$y
}
df_full <- data.frame(covariates, X = I(X_full))

fit_naive <- pfr(time ~ age + BMI + gender + CHD + lf(X, bs = "ps", k = 20, argvals = seq(0,1,len=512)),
                 method = "REML", data = df_full, weights = df_full$event, family = cox.ph)
beta_curve <- coef(fit_naive, n = 512, select = 1)
naive_data <- data.frame(q = beta_curve$X.argvals, beta = beta_curve$value)
saveRDS(naive_data, file = "nhanes_naive_beta_scaled.rds")
