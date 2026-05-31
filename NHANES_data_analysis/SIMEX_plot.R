# ==================== Plot: Naive + m0 curves + SIMEX extrapolations (scaled data) ====================
library(ggplot2)
library(dplyr)
library(mgcv)

setwd("/Users/fenglinxie/Desktop/simulation_lf_uniform_Tlqd_SIMEX_k20/NHANES_cluster_serial/data_scaled_Todd")

# 1. Read naive curve
naive <- readRDS("nhanes_naive_beta_scaled.rds")
q_grid <- naive$q
beta_naive <- naive$beta

# 2. Read average curves for each m0 (combined files)
setwd("/Users/fenglinxie/Desktop/simulation_lf_uniform_Tlqd_SIMEX_k20/NHANES_cluster_serial/data_scaled_Todd/data")

# 2. Read all combined SIMEX result files (auto-detect m0 and number of replicates)
files <- list.files(pattern = "simex_m0_([0-9]+)_B([0-9]+)_combined\\.rds")
beta_m0 <- list()
m0_values <- c()

for (f in files) {
  parts <- strsplit(f, "_")[[1]]
  m0 <- as.numeric(parts[3])
  nrep <- as.numeric(gsub("B", "", parts[4]))
  res <- readRDS(f)
  beta_m0[[as.character(m0)]] <- res$avg_beta
  m0_values <- c(m0_values, m0)
  cat("Loaded", f, "with", nrep, "repetitions\n")
}

# Sort by m0 value (optional)
m0_values <- sort(unique(m0_values))
beta_m0 <- beta_m0[as.character(m0_values)]

# Build matrix (rows: q points, columns: m0)
beta_matrix <- do.call(cbind, beta_m0)
colnames(beta_matrix) <- m0_values
lambda <- 1 / sqrt(m0_values)

# 3. Extrapolation
n_q <- length(q_grid)
# Linear
beta_lin <- numeric(n_q)
for (i in 1:n_q) {
  y <- beta_matrix[i, ]
  fit <- lm(y ~ lambda)
  beta_lin[i] <- predict(fit, newdata = data.frame(lambda = 0))
}
# Quadratic
beta_quad <- numeric(n_q)
for (i in 1:n_q) {
  y <- beta_matrix[i, ]
  fit <- lm(y ~ lambda + I(lambda^2))
  beta_quad[i] <- predict(fit, newdata = data.frame(lambda = 0))
}
# Nonlinear (mgcv::gam)
df_simex <- expand.grid(q = q_grid, lambda = lambda)
df_simex$beta <- as.vector(beta_matrix)
df_simex <- df_simex[!is.na(df_simex$beta), ]
if (nrow(df_simex) > 10) {
  gam_fit <- gam(beta ~ te(q, lambda, k = c(10,5)), data = df_simex, method = "REML")
  beta_nonlin <- predict(gam_fit, newdata = data.frame(q = q_grid, lambda = 0))
} else {
  warning("Not enough data for nonlinear extrapolation; using quadratic instead.")
  beta_nonlin <- beta_quad
}

# 4. Prepare plotting data
df_m0 <- data.frame()
for (m0 in m0_values) {
  df_m0 <- rbind(df_m0,
                 data.frame(q = q_grid, beta = beta_m0[[as.character(m0)]], group = paste0("m0 = ", m0)))
}
df_simex_curves <- data.frame(
  q = rep(q_grid, 3),
  beta = c(beta_lin, beta_quad, beta_nonlin),
  method = rep(c("Linear", "Quadratic", "Nonlinear"), each = n_q)
)
df_naive <- data.frame(q = q_grid, beta = beta_naive, method = "Naive")

# 5. Plot
p <- ggplot() +
  # m0 curves (gray thin lines)
  geom_line(data = df_m0, aes(x = q, y = beta, group = group),
            color = "gray70", size = 0.5, alpha = 0.7) +
  # SIMEX extrapolation curves (colored thick lines)
  geom_line(data = df_simex_curves, aes(x = q, y = beta, color = method), size = 1.2) +
  # Naive curve (black thick line)
  geom_line(data = df_naive, aes(x = q, y = beta), color = "black", size = 1.2) +
  scale_color_manual(values = c("Linear" = "darkgreen", "Quadratic" = "blue", "Nonlinear" = "red")) +
  labs(x = "q", y = expression(beta(q)),
       title = "NHANES: Naive, m0 curves, and SIMEX extrapolations (scaled data)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(p)
ggsave("nhanes_simex_three.png", p, width = 10, height = 6, dpi = 300)



# ==================== Plot: Only nonlinear SIMEX curve ====================
# 4. Prepare plotting data
# Convert m0 curves to long format, add m0 numeric column for color mapping
df_m0 <- data.frame()
for (m0 in m0_values) {
  df_m0 <- rbind(df_m0,
                 data.frame(q = q_grid, beta = beta_m0[[as.character(m0)]], m0 = m0))
}
# Note: m0 is numeric and will be used for color gradient

# Nonlinear SIMEX curve
df_simex_curve <- data.frame(q = q_grid, beta = beta_nonlin, type = "Nonlinear SIMEX")

# Naive curve
df_naive <- data.frame(q = q_grid, beta = beta_naive, type = "Naive")

# 5. Set colors: m0 using blue gradient, Naive black, SIMEX dark green
# Define m0 color mapping: from light blue to dark blue (smallest m0 lightest, largest darkest)
colors_m0 <- scales::seq_gradient_pal(low = "#DEEBF5", high = "#08519C")(seq(0, 1, length.out = length(m0_values)))
names(colors_m0) <- m0_values

# 6. Plot
p <- ggplot() +
  # m0 curves (using color = as.factor(m0) mapped to gradient)
  geom_line(data = df_m0, aes(x = q, y = beta, color = as.factor(m0), group = m0), size = 0.6) +
  # Nonlinear SIMEX curve (green)
  geom_line(data = df_simex_curve, aes(x = q, y = beta, color = type), size = 1.2) +
  # Naive curve (black)
  geom_line(data = df_naive, aes(x = q, y = beta, color = type), size = 1.2) +
  scale_color_manual(
    name = "Curve",
    values = c(colors_m0, "Nonlinear SIMEX" = "red", "Naive" = "black"),
    breaks = c(as.character(m0_values), "Nonlinear SIMEX", "Naive"),
    labels = c(paste0("m0 = ", m0_values), "Nonlinear SIMEX", "Naive")
  ) +
  labs(x = "q", y = expression(beta(q)),
       title = "NHANES: Naive, m0 curves (colored by m0), and Nonlinear SIMEX") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 16))

print(p)
ggsave("nhanes_simex_nonlinear_gradient.png", p, width = 10, height = 6, dpi = 300)

# Save extrapolated curve data
# saveRDS(list(q = q_grid, nonlinear = beta_nonlin), file = "simex_nonlinear_curve_scaled.rds")