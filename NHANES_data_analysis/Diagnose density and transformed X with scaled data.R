# ==================== Diagnose density and transformed X with scaled data ====================
# This script scales all MIMS values by the global maximum (so that values are in [0,1]).
# Then density estimation uses a fixed domain of [-0.1, 1] to avoid boundary issues.

library(ggplot2)
library(patchwork)

# Load data
setwd("/Users/fenglinxie/Desktop/simulation_lf_uniform_Tlqd_SIMEX_k20/Real_data_NHANES")

simex_data <- readRDS("nhanes_simex_data.rds")
original_conn <- simex_data$original_conn   # list of MIMS vectors
covariates <- simex_data$covariates

# Compute global maximum MIMS value across all individuals and all observations
all_mims <- unlist(original_conn)
global_max <- max(all_mims, na.rm = TRUE)
cat("Global maximum MIMS value:", global_max, "\n")

# Scale all MIMS vectors to [0,1] by dividing by global_max
scaled_conn <- lapply(original_conn, function(x) x / global_max)

# Select first N individuals to visualize
n_show <- 20
selected_idx <- 1:n_show

# Define T0 function (density only) with fixed domain [-0.1, 1]
T0_fixed <- function(corr_i, nb = 19) {
  dens.eps <- 1e-6
  # Use fixed domain instead of data range
  domain <- c(-0.1, 1)
  basis1 <- fda::create.bspline.basis(domain, nb)
  densf <- fda::density.fd(corr_i, fda::fd(matrix(0, nb, 1), basis1))
  x <- seq(domain[1], domain[2], length.out = 1024)
  y <- c(exp(fda::eval.fd(x, densf$Wfdobj)) / densf$C)
  # Truncate where density is very small (optional, but keep for stability)
  ind1 <- min(which(y > dens.eps * 512))
  ind2 <- max(which(y > dens.eps * 512))
  list(x = x[ind1:ind2], y = y[ind1:ind2])
}

# Define Tlqd function (density + LQD transform) with fixed domain [-0.1, 1]
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

# Function to process a single m0 (numeric or "full") using scaled data
process_m0_scaled <- function(m0, scaled_conn, selected_idx) {
  cat("\n========================================\n")
  cat("Processing m0 =", m0, "\n")
  
  density_ranges <- data.frame()
  lqd_ranges <- data.frame()
  density_plots <- list()
  lqd_plots <- list()
  
  for (i in selected_idx) {
    if (is.character(m0) && m0 == "full") {
      obs <- scaled_conn[[i]]
      sample_type <- "full data (no resampling)"
    } else {
      set.seed(12345 + as.numeric(m0) * 1000 + i)
      obs <- sample(scaled_conn[[i]], size = as.numeric(m0), replace = TRUE)
      sample_type <- paste0("bootstrap sample of size ", m0)
    }
    
    cat("\nIndividual", i, "(SEQN =", covariates$SEQN[i], ")\n")
    cat("  Sample type:", sample_type, "\n")
    cat("  Scaled MIMS range (original max", global_max, "):", range(obs), "\n")
    
    dens <- T0_fixed(obs)
    cat("  Density x range:", range(dens$x), "\n")
    cat("  Density y range:", range(dens$y), "\n")
    
    lqd_res <- Tlqd_fixed(obs)
    cat("  Transformed density range:", range(lqd_res$y), "\n")
    
    density_ranges <- rbind(density_ranges,
                            data.frame(individual = i, m0 = m0,
                                       x_min = min(dens$x), x_max = max(dens$x),
                                       y_min = min(dens$y), y_max = max(dens$y)))
    lqd_ranges <- rbind(lqd_ranges,
                        data.frame(individual = i, m0 = m0,
                                   y_min = min(lqd_res$y), y_max = max(lqd_res$y)))
    
    # Density plot
    df_dens <- data.frame(x = dens$x, y = dens$y)
    p_dens <- ggplot(df_dens, aes(x = x, y = y)) +
      geom_line(color = "blue") +
      labs(title = paste0("Ind ", i, ", m0 = ", m0), x = "Scaled MIMS", y = "Density") +
      theme_minimal()
    density_plots[[i]] <- p_dens
    
    # Transformed X plot
    df_lqd <- data.frame(q = lqd_res$x, X = lqd_res$y)
    p_lqd <- ggplot(df_lqd, aes(x = q, y = X)) +
      geom_line(color = "red") +
      labs(title = paste0("Ind ", i, ", m0 = ", m0), x = "q", y = "X(q)") +
      theme_minimal()
    lqd_plots[[i]] <- p_lqd
  }
  
  cat("\n--- Density summary (m0 =", m0, ") ---\n")
  print(density_ranges)
  cat("\n--- Transformed density range summary (m0 =", m0, ") ---\n")
  print(lqd_ranges)
  
  return(list(density_plots = density_plots, lqd_plots = lqd_plots,
              density_ranges = density_ranges, lqd_ranges = lqd_ranges))
}

# List of m0 values to test
m0_values <- c(500, 1000, "full")

all_density_ranges <- data.frame()
all_lqd_ranges <- data.frame()

for (m0 in m0_values) {
  res <- process_m0_scaled(m0, scaled_conn, selected_idx)
  all_density_ranges <- rbind(all_density_ranges, res$density_ranges)
  all_lqd_ranges <- rbind(all_lqd_ranges, res$lqd_ranges)
  
  # Save combined plots for first 4 individuals
  if (length(res$density_plots) >= 4) {
    combined_dens <- res$density_plots[[1]] + res$density_plots[[2]] + 
      res$density_plots[[3]] + res$density_plots[[4]] +
      plot_layout(ncol = 2)
    ggsave(paste0("diagnose_density_scaled_m0_", m0, ".png"), combined_dens, width = 10, height = 8)
    
    combined_lqd <- res$lqd_plots[[1]] + res$lqd_plots[[2]] + 
      res$lqd_plots[[3]] + res$lqd_plots[[4]] +
      plot_layout(ncol = 2)
    ggsave(paste0("diagnose_lqd_scaled_m0_", m0, ".png"), combined_lqd, width = 10, height = 8)
  }
}

# Print overall comparison
cat("\n\n========== Overall Comparison ==========\n")
cat("Density ranges:\n")
print(all_density_ranges)
cat("\nTransformed density ranges:\n")
print(all_lqd_ranges)

# Optional: Detailed inspection of individual 12 with m0=1000
ind <- 12
m0 <- 1000
cat("\n\n===== Detailed inspection of individual", ind, "(SEQN", covariates$SEQN[ind], ") =====\n")
print(covariates[ind, ])

# Original scaled MIMS histogram
obs_full_scaled <- scaled_conn[[ind]]
hist(obs_full_scaled, breaks = 50, 
     main = paste("Scaled MIMS distribution for SEQN", covariates$SEQN[ind]),
     xlab = "Scaled MIMS", col = "lightblue")

# Bootstrap sample
set.seed(12345 + m0 * 1000 + ind)
obs <- sample(scaled_conn[[ind]], size = m0, replace = TRUE)
dens <- T0_fixed(obs)
df_dens <- data.frame(x = dens$x, y = dens$y)
p_dens <- ggplot(df_dens, aes(x = x, y = y)) +
  geom_line(color = "blue") +
  labs(title = paste("Ind", ind, "- Density (m0=1000, scaled)"),
       x = "Scaled MIMS", y = "Density") +
  theme_minimal()
print(p_dens)

lqd_res <- Tlqd_fixed(obs)
df_lqd <- data.frame(q = lqd_res$x, X = lqd_res$y)
p_lqd <- ggplot(df_lqd, aes(x = q, y = X)) +
  geom_line(color = "red") +
  labs(title = paste("Ind", ind, "- Transformed X(q) (m0=1000, scaled)"),
       x = "q", y = "X(q)") +
  theme_minimal()
print(p_lqd)