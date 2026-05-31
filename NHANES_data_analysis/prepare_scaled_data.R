# prepare_scaled_data.R
# Purpose: Divide all MIMS values by the global maximum to generate scaled data file

library(tidyverse)

# Load raw data
simex_data <- readRDS("nhanes_simex_data.rds")
original_conn <- simex_data$original_conn
covariates <- simex_data$covariates

# Calculate global maximum
all_mims <- unlist(original_conn)
global_max <- max(all_mims, na.rm = TRUE)
cat("Global maximum MIMS value:", global_max, "\n")

# Scale all MIMS vectors
scaled_conn <- lapply(original_conn, function(x) x / global_max)

# Save scaled data
saveRDS(list(original_conn = scaled_conn, 
             covariates = covariates,
             global_max = global_max), 
        file = "nhanes_simex_data_scaled.rds")