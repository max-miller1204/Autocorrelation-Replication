###############################################################################
# Table 5.1: m-dependent Gaussian product simulation
# Null rejection probabilities for rho(1) = 0
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
test_mode <- "--test" %in% args

source("/tmp/core_functions.R")
source("/tmp/data_generators.R")

if (test_mode) {
  n_sims <- 100
  n_perm <- 200
  sample_sizes <- c(50, 100)
  cat("Running in TEST mode (n_sims=100, n_perm=200)\n")
} else {
  n_sims <- 10000
  n_perm <- 2000
  sample_sizes <- c(10, 20, 50, 80, 100, 500, 1000)
}

alpha <- 0.05
m_values <- 0:3
test_names <- c("Stud. Perm.", "Unst. Perm.", "Ljung-Box", "Box-Pierce")

results <- data.frame()

for (m in m_values) {
  cat(sprintf("\n=== m = %d ===\n", m))
  for (n in sample_sizes) {
    cat(sprintf("  n = %d ... ", n))
    b_n <- floor(n^(1/3)) + 1

    rejections <- matrix(0, nrow = 4, ncol = 1)

    set.seed(42)
    for (sim in 1:n_sims) {
      x <- generate_m_dependent(n, m)

      # Studentized permutation test
      sp <- studentized_permutation_test(x, k = 1, n_perm = n_perm,
                                         b_n = b_n, epsilon = 1e-6)
      if (sp$p_value <= alpha) rejections[1] <- rejections[1] + 1

      # Unstudentized permutation test
      up <- unstudentized_permutation_test(x, k = 1, n_perm = n_perm)
      if (up$p_value <= alpha) rejections[2] <- rejections[2] + 1

      # Ljung-Box
      lb <- ljung_box_test(x, k = 1)
      if (lb$p_value <= alpha) rejections[3] <- rejections[3] + 1

      # Box-Pierce
      bp <- box_pierce_test(x, k = 1)
      if (bp$p_value <= alpha) rejections[4] <- rejections[4] + 1

      if (sim %% 1000 == 0) cat(sprintf("%d ", sim))
    }

    rej_probs <- rejections / n_sims
    cat(sprintf("done: SP=%.4f UP=%.4f LB=%.4f BP=%.4f\n",
                rej_probs[1], rej_probs[2], rej_probs[3], rej_probs[4]))

    for (t in 1:4) {
      results <- rbind(results, data.frame(
        m = m, test = test_names[t], n = n,
        rejection_prob = rej_probs[t]
      ))
    }
  }
}

# Reshape to wide format matching Table 5.1
library(tidyr)
results_wide <- pivot_wider(results, names_from = n, values_from = rejection_prob)

cat("\n\nTable 5.1: Null rejection probabilities (m-dependent Gaussian products)\n")
cat("======================================================================\n")
print(as.data.frame(results_wide), row.names = FALSE)

# Save output
output_dir <- if (test_mode) "/tmp" else "/tmp/output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
write.csv(results_wide, file.path(output_dir, "table_5_1.csv"), row.names = FALSE)
cat(sprintf("\nSaved to %s/table_5_1.csv\n", output_dir))
