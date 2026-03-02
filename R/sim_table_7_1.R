###############################################################################
# Table 7.1: Autocovariance test (Section 7)
# Rejection probabilities for tests of c_1 = 0
# Uses test statistic sqrt(n)*c_hat_n/T_hat_n instead of sqrt(n)*rho_hat_n/gamma_hat_n
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
test_mode <- "--test" %in% args

source("/tmp/core_functions.R")
source("/tmp/data_generators.R")

if (test_mode) {
  n_sims <- 100
  n_perm <- 200
  sample_sizes <- c(50, 100)
  cat("Running in TEST mode\n")
} else {
  n_sims <- 10000
  n_perm <- 2000
  sample_sizes <- c(10, 20, 50, 80, 100, 500, 1000)
}

alpha <- 0.05

# All scenarios: m-dependent (m=0,1,2,3) and alpha-mixing AR(2)
scenarios <- list(
  list(name = "m = 0", gen = function(n) generate_m_dependent(n, 0)),
  list(name = "m = 1", gen = function(n) generate_m_dependent(n, 1)),
  list(name = "m = 2", gen = function(n) generate_m_dependent(n, 2)),
  list(name = "m = 3", gen = function(n) generate_m_dependent(n, 3)),
  list(name = "AR(2), N(0,1) innov.",
       gen = function(n) generate_ar2(n, 0.5, "gaussian")),
  list(name = "AR(2) Prod., N(0,1) innov.",
       gen = function(n) generate_ar2_product(n, 0.5)),
  list(name = "AR(2), U[-1,1] innov.",
       gen = function(n) generate_ar2(n, 0.5, "uniform")),
  list(name = "AR(2), t_9.5 innov.",
       gen = function(n) generate_ar2(n, 0.5, "t", t_df = 9.5))
)

results <- data.frame()

for (scen in scenarios) {
  cat(sprintf("\n=== %s ===\n", scen$name))
  for (n in sample_sizes) {
    cat(sprintf("  n = %d ... ", n))
    b_n <- floor(n^(1/3)) + 1

    rejections <- 0
    set.seed(42)
    for (sim in 1:n_sims) {
      x <- scen$gen(n)
      sc <- studentized_autocovariance_test(x, k = 1, n_perm = n_perm,
                                            b_n = b_n, epsilon = 1e-6)
      if (sc$p_value <= alpha) rejections <- rejections + 1
    }
    rej_prob <- rejections / n_sims
    cat(sprintf("rej_prob = %.4f\n", rej_prob))

    results <- rbind(results, data.frame(
      Scenario = scen$name, n = n, rejection_prob = rej_prob
    ))
  }
}

library(tidyr)
results_wide <- pivot_wider(results, names_from = n, values_from = rejection_prob)

cat("\n\nTable 7.1: Rejection probabilities for c_1 = 0\n")
cat("================================================\n")
print(as.data.frame(results_wide), row.names = FALSE)

output_dir <- if (test_mode) "/tmp" else "/tmp/output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
write.csv(results_wide, file.path(output_dir, "table_7_1.csv"), row.names = FALSE)
cat(sprintf("\nSaved to %s/table_7_1.csv\n", output_dir))
