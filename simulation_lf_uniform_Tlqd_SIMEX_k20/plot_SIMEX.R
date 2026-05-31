# ==================== Summarize all results and plot (without selection bias) ====================
setwd("/Users/fenglinxie/Downloads/average_then_MSE_lf_uniform_Tlqd_k20")

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

true_beta <- function(q) {
  mean1 <- 0.3; sd1 <- 0.1; mean2 <- 0.7; sd2 <- 0.1
  amp1 <- 0.25; amp2 <- 0.25
  raw_beta <- amp1 * dnorm(q, mean1, sd1) - amp2 * dnorm(q, mean2, sd2)
  q_grid_mean <- seq(0, 1, length.out = 1000)
  raw_grid <- amp1 * dnorm(q_grid_mean, mean1, sd1) - amp2 * dnorm(q_grid_mean, mean2, sd2)
  mean_beta <- mean(raw_grid)
  return(raw_beta - mean_beta)
}

files <- list.files(pattern = "result_.*\\.rds")
cat("Found", length(files), "result files\n")

# Initialize lists
mse_m0_list <- list()
mse_lin_list <- list()
mse_quad_list <- list()
mse_nonlin_list <- list()
beta_m0_list <- list()
beta_lin_list <- list()
beta_quad_list <- list()
beta_nonlin_list <- list()

for (f in files) {
  parts <- strsplit(f, "_")[[1]]
  n <- as.numeric(gsub("n", "", parts[2]))
  sigma <- as.numeric(gsub("sigma", "", parts[3]))
  rep <- as.numeric(gsub("rep|\\.rds", "", parts[4]))
  
  res <- readRDS(f)
  if (!exists("q_grid")) q_grid <- res$q_grid
  
  m0_vals <- c(10,20,50,100,200,500,1000)
  for (i in seq_along(m0_vals)) {
    mse_m0_list[[length(mse_m0_list) + 1]] <- data.frame(
      n = n, sigma = sigma, rep = rep,
      m0 = m0_vals[i],
      mse = res$mse_m0[i]
    )
    beta_m0_list[[length(beta_m0_list) + 1]] <- data.frame(
      n = n, sigma = sigma, rep = rep,
      m0 = m0_vals[i],
      q = q_grid,
      beta = res$beta_m0_matrix[, i],
      type = "m0"
    )
  }
  
  mse_lin_list[[length(mse_lin_list) + 1]] <- data.frame(
    n = n, sigma = sigma, rep = rep,
    mse = res$mse_lin, method = "Linear"
  )
  mse_quad_list[[length(mse_quad_list) + 1]] <- data.frame(
    n = n, sigma = sigma, rep = rep,
    mse = res$mse_quad, method = "Quadratic"
  )
  mse_nonlin_list[[length(mse_nonlin_list) + 1]] <- data.frame(
    n = n, sigma = sigma, rep = rep,
    mse = res$mse_nonlin, method = "Nonlinear"
  )
  
  beta_lin_list[[length(beta_lin_list) + 1]] <- data.frame(
    n = n, sigma = sigma, rep = rep,
    q = q_grid,
    beta = res$beta_lin,
    method = "Linear"
  )
  beta_quad_list[[length(beta_quad_list) + 1]] <- data.frame(
    n = n, sigma = sigma, rep = rep,
    q = q_grid,
    beta = res$beta_quad,
    method = "Quadratic"
  )
  beta_nonlin_list[[length(beta_nonlin_list) + 1]] <- data.frame(
    n = n, sigma = sigma, rep = rep,
    q = q_grid,
    beta = res$beta_nonlin,
    method = "Nonlinear"
  )
}

# Combine data
mse_m0_df <- do.call(rbind, mse_m0_list)
mse_simex_df <- bind_rows(do.call(rbind, mse_lin_list),
                          do.call(rbind, mse_quad_list),
                          do.call(rbind, mse_nonlin_list))
beta_m0_df <- do.call(rbind, beta_m0_list)
beta_simex_df <- bind_rows(do.call(rbind, beta_lin_list),
                           do.call(rbind, beta_quad_list),
                           do.call(rbind, beta_nonlin_list))

# Compute average MSE for each (n, sigma, m0)
avg_mse_m0 <- mse_m0_df %>%
  group_by(n, sigma, m0) %>%
  summarise(mean_mse = mean(mse, na.rm = TRUE), .groups = "drop")

# Compute average MSE for the three extrapolation methods (each (n, sigma, method))
avg_mse_simex <- mse_simex_df %>%
  group_by(n, sigma, method) %>%
  summarise(mean_mse = mean(mse, na.rm = TRUE), .groups = "drop")

# Compute average beta curves (for m0 curves)
avg_beta_m0 <- beta_m0_df %>%
  group_by(n, sigma, m0, q) %>%
  summarise(mean_beta = mean(beta, na.rm = TRUE), .groups = "drop")

# Compute average beta curves for the three extrapolation methods
avg_beta_simex <- beta_simex_df %>%
  group_by(n, sigma, method, q) %>%
  summarise(mean_beta = mean(beta, na.rm = TRUE), .groups = "drop")

# True curve
true_vals <- true_beta(q_grid)
true_df <- data.frame(q = q_grid, beta = true_vals)
true_full <- crossing(true_df, n = unique(avg_beta_m0$n), sigma = unique(avg_beta_m0$sigma))

# ==================== Figure 1: MSE vs m0 ====================
p_mse <- ggplot(avg_mse_m0, aes(x = m0, y = mean_mse, color = as.factor(sigma))) +
  geom_line() + geom_point() +
  facet_wrap(~ n, labeller = label_both, scales = "free_y") +
  scale_x_log10() +
  labs(x = "m0 (log scale)", y = "Average MSE", color = "σ") +
  theme_minimal()
plot(p_mse)
ggsave("MSE_trends.png", p_mse, width = 10, height = 6, dpi = 300)

# ==================== Figure 2: Bar plot comparing average MSE of three extrapolation methods ====================
p_mse_simex <- ggplot(avg_mse_simex, aes(x = method, y = mean_mse, fill = method)) +
  geom_col(position = "dodge", show.legend = FALSE) +
  facet_grid(n ~ sigma, 
             labeller = labeller(n = function(x) paste0("N = ", x),
                                 sigma = function(x) paste0("σ = ", x))) +
  labs(x = "Extrapolation method", y = "Average MISE") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

print(p_mse_simex)
ggsave("MSE_simex_comparison.png", p_mse_simex, width = 10, height = 6, dpi = 300)

# ==================== Figure 3: Average MSE vs m0 + Nonlinear SIMEX point (m=2000, marked as ∞) ====================
# Extract average MSE for Nonlinear SIMEX
nonlin_mse <- avg_mse_simex %>% filter(method == "Nonlinear") %>%
  mutate(m0 = 2000, type = "Nonlinear SIMEX")

# Combine data
mse_plot <- bind_rows(
  avg_mse_m0 %>% mutate(type = "m0"),
  nonlin_mse %>% dplyr::select(n, sigma, m0, mean_mse, type)
)

p_mse_with_simex <- ggplot(mse_plot, aes(x = m0, y = mean_mse, color = as.factor(sigma))) +
  geom_line(data = subset(mse_plot, type == "m0"), aes(group = sigma)) +
  geom_point(data = subset(mse_plot, type == "m0"), size = 2) +
  geom_point(data = subset(mse_plot, type == "Nonlinear SIMEX"), 
             aes(shape = type), size = 3, stroke = 1) +
  facet_wrap(~ n, labeller = label_both, scales = "free_y") +
  scale_x_log10(breaks = c(10, 50, 200, 1000, 2000),
                labels = c("10", "50", "200", "1000", "∞")) +
  labs(x = "m0 (log scale)", y = "Average MSE", color = "σ", shape = "Curve type") +
  scale_shape_manual(values = c("Nonlinear SIMEX" = 18)) +
  theme_minimal()

plot(p_mse_with_simex)
ggsave("MSE_trends_with_SIMEX_point.png", p_mse_with_simex, width = 10, height = 6, dpi = 300)

# ==================== Figure 4: Average beta curves (m0 blue gradient, Nonlinear SIMEX red, True black) ====================
# Prepare data: combine m0 curves, Nonlinear SIMEX curve, and true curve
df_m0 <- avg_beta_m0 %>%
  mutate(curve_type = as.character(m0)) %>%
  rename(beta = mean_beta)   # unify column name to beta

df_simex <- avg_beta_simex %>%
  filter(method == "Nonlinear") %>%
  mutate(curve_type = "Nonlinear SIMEX") %>%
  rename(beta = mean_beta)

df_true <- true_full %>%
  mutate(curve_type = "True")   # true_full already has beta column

all_curves <- bind_rows(df_m0, df_simex, df_true)

# Define colors: m0 use blue gradient, Nonlinear SIMEX red, True black
m0_levels <- sort(unique(avg_beta_m0$m0))
colors_m0 <- brewer.pal(length(m0_levels), "Blues")
names(colors_m0) <- as.character(m0_levels)
all_colors <- c(colors_m0, "Nonlinear SIMEX" = "red", "True" = "black")

# Legend order: m0 from small to large, then Nonlinear SIMEX, then True
legend_order <- c(as.character(m0_levels), "Nonlinear SIMEX", "True")

# Convert curve_type to factor to control legend order
all_curves$curve_type <- factor(all_curves$curve_type, levels = legend_order)

# Get all parameter combinations
combos <- expand.grid(n = unique(avg_beta_m0$n), sigma = unique(avg_beta_m0$sigma))

for (i in 1:nrow(combos)) {
  n_val <- combos$n[i]
  sigma_val <- combos$sigma[i]
  
  sub <- all_curves %>% filter(n == n_val, sigma == sigma_val)
  
  p <- ggplot(sub, aes(x = q, y = beta, color = curve_type, group = curve_type)) +
    geom_line(size = 0.8) +
    scale_color_manual(name = "Curve", values = all_colors, 
                       breaks = legend_order,
                       labels = c(as.character(m0_levels), "Nonlinear SIMEX", "True")) +
    labs(title = paste0("n = ", n_val, ", σ = ", sigma_val),
         x = "q", y = "β(q)") +
    theme_minimal() +
    theme(legend.position = "right")
  
  print(p)
  filename <- paste0("beta_curves_m0_simex_true_n", n_val, "_sigma", gsub("\\.", "_", sigma_val), ".png")
  ggsave(filename, p, width = 8, height = 5, dpi = 300)
}

# ==================== Figure 5: Spaghetti plots (individual curves for naive and nonlinear SIMEX + average curves) ====================
# Extract naive individual curves (m0=1000)
naive_ind <- beta_m0_df %>% filter(m0 == 1000) %>% 
  mutate(curve_type = "Naive")
# Extract nonlinear SIMEX individual curves
simex_ind <- beta_simex_df %>% filter(method == "Nonlinear") %>% 
  mutate(curve_type = "Nonlinear SIMEX")
# True curve for legend (handled separately)
true_for_plot <- true_df %>% mutate(curve_type = "True", rep = NA)  # add rep column for bind_rows, but handled separately

combined_ind <- bind_rows(naive_ind, simex_ind)

for (n_val in unique(combined_ind$n)) {
  for (sigma_val in unique(combined_ind$sigma)) {
    sub <- combined_ind %>% filter(n == n_val, sigma == sigma_val)
    # Compute average curves (by curve_type and q)
    avg_sub <- sub %>% group_by(curve_type, q) %>% 
      summarise(beta_avg = mean(beta), .groups = "drop")
    
    # True curve data (same parameter combination)
    true_sub <- true_df %>% mutate(n = n_val, sigma = sigma_val)
    
    p <- ggplot() +
      # Individual curves (semi-transparent)
      geom_line(data = sub, 
                aes(x = q, y = beta, group = interaction(rep, curve_type), color = curve_type),
                size = 0.3, alpha = 0.3) +
      # Average curves (solid lines, darker colors)
      geom_line(data = avg_sub, 
                aes(x = q, y = beta_avg, color = curve_type), size = 1) +
      # True curve (black, mapped via aes(color="True"))
      geom_line(data = true_sub, aes(x = q, y = beta, color = "True"), size = 1) +
      scale_color_manual(name = "Curve type",
                         values = c("Naive" = "blue", "Nonlinear SIMEX" = "red", "True" = "black"),
                         breaks = c("Naive", "Nonlinear SIMEX", "True")) +
      labs(title = paste0("n = ", n_val, ", σ = ", sigma_val), 
           x = "q", y = "β(q)") +
      theme_minimal()
    print(p)
    ggsave(paste0("spaghetti_naive_simex_n", n_val, "_sigma", gsub("\\.", "_", sigma_val), ".png"),
           p, width = 8, height = 5, dpi = 300)
  }
}

# ==================== Figure 6: Variance curves (naive vs nonlinear SIMEX) and summary table ====================
# Extract naive individual curves (m0=1000) and compute variance
naive_var <- beta_m0_df %>% filter(m0 == 1000) %>%
  group_by(n, sigma, q) %>%
  summarise(variance = var(beta, na.rm = TRUE), .groups = "drop") %>%
  mutate(type = "Naive")

# Extract nonlinear SIMEX individual curves and compute variance
simex_var <- beta_simex_df %>% filter(method == "Nonlinear") %>%
  group_by(n, sigma, q) %>%
  summarise(variance = var(beta, na.rm = TRUE), .groups = "drop") %>%
  mutate(type = "Nonlinear SIMEX")

variance_df <- bind_rows(naive_var, simex_var)

p_var <- ggplot(variance_df, aes(x = q, y = variance, color = type)) +
  geom_line(size = 0.8) +
  facet_grid(n ~ sigma, labeller = label_both) +
  labs(x = "q", y = "Variance", color = "Curve type") +
  scale_color_manual(values = c("Naive" = "blue", "Nonlinear SIMEX" = "red")) +
  theme_minimal(base_size = 18)

plot(p_var)
ggsave("variance_curves.png", p_var, width = 10, height = 6, dpi = 300)


####### Analyze variance and bias of MSE ##########
# Compute squared bias for naive average curve
bias_naive <- avg_beta_m0 %>%
  filter(m0 == 1000) %>%
  left_join(true_df, by = "q") %>%
  group_by(n, sigma) %>%
  summarise(
    bias_sq_naive = mean((mean_beta - beta)^2, na.rm = TRUE),
    .groups = "drop"
  )

# Compute squared bias for nonlinear SIMEX average curve
bias_simex <- avg_beta_simex %>%
  filter(method == "Nonlinear") %>%
  left_join(true_df, by = "q") %>%
  group_by(n, sigma) %>%
  summarise(
    bias_sq_simex = mean((mean_beta - beta)^2, na.rm = TRUE),
    .groups = "drop"
  )

# Combine variance and squared bias, compute total MSE
avg_var <- variance_df %>%
  group_by(n, sigma, type) %>%
  summarise(mean_variance = mean(variance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = type, values_from = mean_variance)

comparison <- bias_naive %>%
  full_join(bias_simex, by = c("n", "sigma")) %>%
  left_join(avg_var, by = c("n", "sigma")) %>%
  mutate(
    MSE_naive = bias_sq_naive + Naive,
    MSE_simex = bias_sq_simex + `Nonlinear SIMEX`
  ) %>%
  dplyr::select(n, sigma, 
                bias_sq_naive, bias_sq_simex,
                variance_naive = Naive, variance_simex = `Nonlinear SIMEX`,
                MSE_naive, MSE_simex)

print(comparison)
# write.csv(comparison, "bias_variance_mse_comparison.csv", row.names = FALSE)



# ==================== Modified figure: Average root-MISE vs m0 + Nonlinear SIMEX point ====================
nonlin_mse <- avg_mse_simex %>% filter(method == "Nonlinear") %>%
  mutate(m0 = 2000, type = "Nonlinear SIMEX")
nonlin_mse <- nonlin_mse %>% mutate(root_mise = sqrt(mean_mse))

avg_mse_m0 <- avg_mse_m0 %>% mutate(root_mise = sqrt(mean_mse))

mse_plot <- bind_rows(
  avg_mse_m0 %>% mutate(type = "unadjusted"),
  nonlin_mse %>% dplyr::select(n, sigma, m0, root_mise, type)
) %>%
  mutate(shape_group = case_when(
    type == "unadjusted" & m0 == 1000 ~ "m0 = 1000 (naive)",
    type == "unadjusted" & m0 != 1000 ~ "m0 < 1000",
    type == "Nonlinear SIMEX" ~ "SIMEX (Nonlinear)"
  ))

p_mse_with_simex <- ggplot() +
  geom_line(data = subset(mse_plot, type == "unadjusted" & m0 != 1000), 
            aes(x = m0, y = root_mise, color = as.factor(sigma), group = sigma), size = 0.8) +
  geom_point(data = mse_plot, 
             aes(x = m0, y = root_mise, color = as.factor(sigma), shape = shape_group), size = 2) +
  facet_wrap(~ n, labeller = labeller(n = function(x) paste0("N = ", x)), scales = "fixed") +
  scale_x_log10(breaks = c(10, 50, 200, 1000, 2000),
                labels = c("10", "50", "200", "1000", "∞")) +
  coord_cartesian(ylim = c(0.18, 0.6)) +
  labs(x = "Bootstrap sample size (log scale)", 
       y = "Average RMISE", 
       color = expression(sigma), 
       shape = "Curve type") +
  scale_color_manual(values = c("0.01" = "blue", "0.05" = "red"),
                     labels = c(expression(sigma == 0.01), expression(sigma == 0.05))) +
  scale_shape_manual(values = c("m0 < 1000" = 16, "m0 = 1000 (naive)" = 18, "SIMEX (Nonlinear)" = 17)) +
  theme_minimal(base_size = 16) +
  guides(shape = guide_legend(override.aes = list(linetype = 0)))

print(p_mse_with_simex)
ggsave("RMISE_trends_with_SIMEX_point.png", p_mse_with_simex, width = 10, height = 8, dpi = 300)


# ==================== Modified Figure 3: Beta curves for all parameter combinations (combined into one plot) ====================
# Prepare data: combine m0 curves, Nonlinear SIMEX curve, and true curve
df_m0 <- avg_beta_m0 %>%
  mutate(curve_type = as.character(m0)) %>%
  rename(beta = mean_beta)

df_simex <- avg_beta_simex %>%
  filter(method == "Nonlinear") %>%
  mutate(curve_type = "Nonlinear SIMEX") %>%
  rename(beta = mean_beta)

df_true <- true_full %>%
  mutate(curve_type = "True")

all_curves <- bind_rows(df_m0, df_simex, df_true)

# Define colors: m0 use blue gradient, Nonlinear SIMEX red, True black
m0_levels <- sort(unique(avg_beta_m0$m0))
colors_m0 <- brewer.pal(length(m0_levels), "Blues")
names(colors_m0) <- as.character(m0_levels)
all_colors <- c(colors_m0, "Nonlinear SIMEX" = "red", "True" = "black")

# Legend order: m0 from small to large, then Nonlinear SIMEX, then True
legend_order <- c(as.character(m0_levels), "Nonlinear SIMEX", "True")
all_curves$curve_type <- factor(all_curves$curve_type, levels = legend_order)

# Create facet label columns (to display N and σ)
all_curves$n_label <- paste0("N = ", all_curves$n)
all_curves$sigma_label <- paste0("σ = ", all_curves$sigma)

# Combined plot
p_combined <- ggplot(all_curves, aes(x = q, y = beta, color = curve_type, group = curve_type)) +
  geom_line(size = 0.8) +
  facet_grid(n_label ~ sigma_label, scales = "free_y") +
  scale_color_manual(
    name = "Curve",
    values = all_colors,
    breaks = legend_order,
    labels = c(as.character(m0_levels), "Nonlinear SIMEX", "True")
  ) +
  labs(x = "q", y = expression(beta(q))) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 18),
    legend.key.width = unit(1.5, "cm")
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

print(p_combined)
ggsave("beta_curves_combined.png", p_combined, width = 14, height = 8, dpi = 300)



# ==================== Modified Figure 4: Spaghetti plot (combined into one figure) ====================
# Extract naive individual curves (m0=1000)
naive_ind <- beta_m0_df %>% filter(m0 == 1000) %>% 
  mutate(curve_type = "Naive")

# Extract nonlinear SIMEX individual curves
simex_ind <- beta_simex_df %>% filter(method == "Nonlinear") %>% 
  mutate(curve_type = "Nonlinear SIMEX")

# True curve (shared across all combinations)
true_for_plot <- true_df %>% 
  crossing(n = unique(naive_ind$n), sigma = unique(naive_ind$sigma)) %>%
  mutate(curve_type = "True")

# Combine individual curves (for semi-transparent thin lines)
combined_ind <- bind_rows(naive_ind, simex_ind)

# Compute average curves (by curve_type, n, sigma, q)
avg_curves <- combined_ind %>%
  group_by(n, sigma, q, curve_type) %>%
  summarise(beta_avg = mean(beta, na.rm = TRUE), .groups = "drop")

# Create facet label columns
combined_ind$n_label <- paste0("N = ", combined_ind$n)
combined_ind$sigma_label <- paste0("σ = ", combined_ind$sigma)
avg_curves$n_label <- paste0("N = ", avg_curves$n)
avg_curves$sigma_label <- paste0("σ = ", avg_curves$sigma)
true_for_plot$n_label <- paste0("N = ", true_for_plot$n)
true_for_plot$sigma_label <- paste0("σ = ", true_for_plot$sigma)

# Set factor order: top to bottom N=50, N=100, N=200
n_levels <- c("N = 50", "N = 100", "N = 200")
combined_ind$n_label <- factor(combined_ind$n_label, levels = n_levels)
avg_curves$n_label <- factor(avg_curves$n_label, levels = n_levels)
true_for_plot$n_label <- factor(true_for_plot$n_label, levels = n_levels)

# Set sigma order (left to right, optional)
sigma_levels <- c("σ = 0.01", "σ = 0.05")
combined_ind$sigma_label <- factor(combined_ind$sigma_label, levels = sigma_levels)
avg_curves$sigma_label <- factor(avg_curves$sigma_label, levels = sigma_levels)
true_for_plot$sigma_label <- factor(true_for_plot$sigma_label, levels = sigma_levels)

# Plot
p_spaghetti <- ggplot() +
  # Individual curves (semi-transparent thin lines)
  geom_line(data = combined_ind, 
            aes(x = q, y = beta, group = interaction(rep, curve_type), color = curve_type),
            size = 0.3, alpha = 0.3) +
  # Average curves (solid lines, darker colors)
  geom_line(data = avg_curves, 
            aes(x = q, y = beta_avg, color = curve_type), size = 1) +
  # True curve (black)
  geom_line(data = true_for_plot, 
            aes(x = q, y = beta, color = "True"), size = 1) +
  facet_grid(n_label ~ sigma_label, scales = "free_y") +
  scale_color_manual(
    name = "Curve type",
    values = c("Naive" = "blue", "Nonlinear SIMEX" = "red", "True" = "black"),
    breaks = c("Naive", "Nonlinear SIMEX", "True")
  ) +
  labs(x = "q", y = expression(beta(q))) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 18),
    legend.key.width = unit(1, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1), nrow = 1))

print(p_spaghetti)
ggsave("spaghetti_combined_all.png", p_spaghetti, width = 14, height = 8, dpi = 300)