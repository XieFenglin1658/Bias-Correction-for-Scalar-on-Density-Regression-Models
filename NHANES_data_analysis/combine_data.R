setwd("/Users/fenglinxie/Desktop/simulation_lf_uniform_Tlqd_SIMEX_k20/NHANES_cluster_serial/data_scaled_Todd/data")

m0_values <- c(10, 20, 50, 100, 200, 500, 1000)
blocks <- 10
B_per_block <- 10

for (m0 in m0_values) {
  all_beta <- NULL
  for (block in 1:blocks) {
    start_rep <- (block - 1) * B_per_block + 1
    fname <- paste0("simex_m0_", m0, "_B10_start", start_rep, ".rds")
    if (file.exists(fname)) {
      res <- readRDS(fname)
      if (is.null(all_beta)) {
        all_beta <- res$all_beta
      } else {
        all_beta <- cbind(all_beta, res$all_beta)
      }
    }
  }
  if (!is.null(all_beta)) {
    avg_beta <- rowMeans(all_beta, na.rm = TRUE)
    outfile <- paste0("simex_m0_", m0, "_B", ncol(all_beta), "_combined.rds")
    saveRDS(list(m0 = m0, avg_beta = avg_beta, all_beta = all_beta), file = outfile)
    cat("Combined:", outfile, "contains", ncol(all_beta), "replicates\n")
  }
}