# ==================== Set working directory ====================
# Please modify to your local path where result files are stored
work_dir <- "/Users/fenglinxie/Downloads/simulation_lf_n100_sigma0.01_more_mvalues_Tlqd" 
setwd(work_dir)

# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)  # for color palettes

# -------------------- 1. Define true beta function --------------------
true_beta <- function(q) {
  mean1 <- 0.3; sd1 <- 0.1; mean2 <- 0.7; sd2 <- 0.1
  amp1 <- 0.25; amp2 <- 0.25
  raw_beta <- amp1 * dnorm(q, mean1, sd1) - amp2 * dnorm(q, mean2, sd2)
  q_grid_mean <- seq(0, 1, length.out = 1000)
  raw_grid <- amp1 * dnorm(q_grid_mean, mean1, sd1) - amp2 * dnorm(q_grid_mean, mean2, sd2)
  mean_beta <- mean(raw_grid)
  return(raw_beta - mean_beta)
}

# Generate q grid (should match simulation, 512 points)
q_grid <- seq(0, 1, length.out = 512)
true_vals <- true_beta(q_grid)

# -------------------- 2. Read result files for all m --------------------
m_values <- c(50, 100, 200, 300, 500, 1000, 2000, 3000, 5000, 10000, 100000)
results_list <- list()
stats_df <- data.frame()

for (m in m_values) {
  file_name <- paste0("unif_n100_m", m, "_sigma0_01_results.RData")
  if (file.exists(file_name)) {
    load(file_name)  # should contain avg_curve, mse_value, bias_sq_value, variance_value
    results_list[[as.character(m)]] <- avg_curve
    stats_df <- rbind(stats_df, data.frame(
      m = m,
      MISE = mse_value,        # note: original variable named mse_value, but output label changed to MISE
      Bias2 = bias_sq_value,
      Variance = variance_value
    ))
  } else {
    warning(paste("File", file_name, "not found. Skipping."))
  }
}

# Check if any data was read
if (length(results_list) == 0) stop("No result files found.")

# -------------------- 3. Prepare spaghetti plot data --------------------
# Convert list to data frame with columns as m values and rows as q points
df_curves <- as.data.frame(results_list)
names(df_curves) <- names(results_list)
df_curves$q <- q_grid

# Convert to long format
df_long <- df_curves %>%
  pivot_longer(cols = -q, names_to = "m", values_to = "beta")

# Convert m to factor, sorted in decreasing order (so legend shows m decreasing top to bottom)
m_levels <- sort(m_values, decreasing = TRUE)  # from largest to smallest
df_long$m <- factor(df_long$m, levels = as.character(m_levels))

# Add true curve data frame, with m set to "True"
df_true <- data.frame(q = q_grid, beta = true_vals, m = "True")

# Combine the two data frames
df_all <- bind_rows(df_long, df_true)

# Convert m to factor, placing "True" at the end (or any position)
df_all$m <- factor(df_all$m, levels = c(as.character(m_levels), "True"))

# -------------------- 4. Draw spaghetti plot --------------------
# Assign colors to m values (use viridis palette, reverse direction so largest m is darkest)
n_m <- length(m_levels)
viridis_colors <- viridis(n_m, direction = -1)  # reverse direction
names(viridis_colors) <- as.character(m_levels)

# Add black for the true curve
all_colors <- c(viridis_colors, "True" = "black")

p_spaghetti <- ggplot(df_all, aes(x = q, y = beta, color = m)) +
  geom_line(size = 0.7) +
  scale_color_manual(values = all_colors, name = "m") +
  labs(title = expression(paste("Estimated ", beta, " curves for different m (", italic(N), " = 100, ", sigma, " = 0.01)")),
       x = "q", y = expression(beta(q))) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 24))

print(p_spaghetti)

# Save plot
ggsave("spaghetti_plot.png", p_spaghetti, width = 10, height = 6, dpi = 300)

# -------------------- 5. Plot MISE, Bias², Variance vs. m --------------------
stats_df <- stats_df[order(stats_df$m), ]

stats_long <- stats_df %>%
  pivot_longer(cols = c(MISE, Bias2, Variance), names_to = "metric", values_to = "value")

p_stats <- ggplot(stats_long, aes(x = m, y = value, color = metric, group = metric)) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_log10() +
  scale_color_manual(values = c("MISE" = "red", "Bias2" = "blue", "Variance" = "darkgreen")) +
  labs(title = expression(paste("Performance metrics vs. m (", italic(N), " = 100, ", sigma, " = 0.01)")),
       x = "m (log scale)", y = "Value") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 20))

print(p_stats)

ggsave("metrics_vs_m.png", p_stats, width = 8, height = 5, dpi = 300)

# Print statistics table
print(stats_df)