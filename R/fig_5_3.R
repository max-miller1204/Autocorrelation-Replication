###############################################################################
# Figure 5.3: Power curves under local alternatives
# AR(1): X_i = rho_n * X_{i-1} + epsilon_i, rho_n = h/sqrt(n)
# (Example 3.3, Eq. 3.44)
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
test_mode <- "--test" %in% args

source("/tmp/core_functions.R")
source("/tmp/data_generators.R")
library(ggplot2)

if (test_mode) {
  n_sims <- 100
  n_perm <- 200
  sample_sizes <- c(50, 100)
  h_values <- seq(0, 3, by = 0.5)
  cat("Running in TEST mode\n")
} else {
  n_sims <- 10000
  n_perm <- 2000
  sample_sizes <- c(10, 20, 50, 80, 100, 500, 1000)
  h_values <- seq(0, 3, by = 0.1)
}

alpha <- 0.05

# Generate AR(1) with local alternative rho_n = h/sqrt(n)
generate_local_ar1 <- function(n, h) {
  rho_n <- h / sqrt(n)
  eps <- rnorm(n)
  X <- numeric(n)
  X[1] <- eps[1]
  for (t in 2:n) {
    X[t] <- rho_n * X[t - 1] + eps[t]
  }
  X
}

results <- data.frame()

for (n in sample_sizes) {
  cat(sprintf("\n=== n = %d ===\n", n))
  b_n <- floor(n^(1/3)) + 1

  for (h in h_values) {
    cat(sprintf("  h = %.1f ... ", h))

    rejections <- 0
    set.seed(42)
    for (sim in 1:n_sims) {
      x <- generate_local_ar1(n, h)
      sp <- studentized_permutation_test(x, k = 1, n_perm = n_perm,
                                         b_n = b_n, epsilon = 1e-6)
      if (sp$p_value <= alpha) rejections <- rejections + 1
    }
    rej_prob <- rejections / n_sims
    cat(sprintf("rej_prob = %.4f\n", rej_prob))

    results <- rbind(results, data.frame(
      n = n, h = h, rejection_prob = rej_prob
    ))
  }
}

# Theoretical limiting power curve
# Power = 1 - Phi(z_{1-alpha} - h) where z_{1-alpha} = qnorm(1-alpha)
# For one-sided test: P(reject) = 1 - Phi(z_{1-alpha} - h)
z_alpha <- qnorm(1 - alpha)
theoretical <- data.frame(
  h = h_values,
  rejection_prob = 1 - pnorm(z_alpha - h_values)
)

# Plot
results$n_label <- factor(paste0("n = ", results$n),
                           levels = paste0("n = ", sample_sizes))

p <- ggplot() +
  geom_line(data = results,
            aes(x = h, y = rejection_prob, color = n_label),
            linewidth = 0.8) +
  geom_line(data = theoretical,
            aes(x = h, y = rejection_prob),
            linetype = "dashed", linewidth = 1, color = "black") +
  labs(x = "h", y = "Rejection probabilities",
       color = "") +
  scale_color_manual(
    values = c("n = 10" = "green3", "n = 20" = "blue",
               "n = 50" = "red", "n = 80" = "green4",
               "n = 100" = "orange", "n = 500" = "gray50",
               "n = 1000" = "gray30"),
    labels = c(paste0("n = ", sample_sizes), "Limiting probabilities")
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  annotate("text", x = 2.5, y = 0.15, label = "Limiting probabilities",
           fontface = "italic", size = 3)

output_dir <- if (test_mode) "/tmp" else "/tmp/output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
ggsave(file.path(output_dir, "figure_5_3.pdf"), p, width = 8, height = 6)
cat(sprintf("\nSaved figure_5_3.pdf to %s\n", output_dir))
