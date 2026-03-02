###############################################################################
# Figures 5.1 & 5.2: Density and QQ plots for m-dependent setting
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
test_mode <- "--test" %in% args

source("/tmp/core_functions.R")
source("/tmp/data_generators.R")
library(ggplot2)

if (test_mode) {
  n_sims <- 100
  n_perm <- 200
  cat("Running in TEST mode\n")
} else {
  n_sims <- 10000
  n_perm <- 2000
}

alpha <- 0.05
m_values <- 0:3

###############################################################################
# Figure 5.1: Kernel density estimates of test statistic vs permutation
# distribution for n=1000
###############################################################################
cat("\n=== Figure 5.1: Kernel density estimates (n=1000) ===\n")
n_fig1 <- if (test_mode) 100 else 1000

density_data <- data.frame()

for (m in m_values) {
  cat(sprintf("  m = %d\n", m))
  b_n <- floor(n_fig1^(1/3)) + 1

  test_stats <- numeric(n_sims)
  all_perm_stats <- numeric(0)

  set.seed(42)
  for (sim in 1:n_sims) {
    x <- generate_m_dependent(n_fig1, m)

    rho_obs <- sample_autocorrelation(x, 1)
    gamma_sq_obs <- gamma_hat_squared(x, 1, b_n)
    test_stats[sim] <- sqrt(n_fig1) * rho_obs / sqrt(gamma_sq_obs)

    # Collect permutation stats (pool across sims)
    x_perm <- sample(x)
    rho_perm <- sample_autocorrelation(x_perm, 1)
    gamma_sq_perm <- gamma_hat_squared(x_perm, 1, b_n)
    all_perm_stats <- c(all_perm_stats,
                        sqrt(n_fig1) * rho_perm / sqrt(gamma_sq_perm))

    if (sim %% 1000 == 0) cat(sprintf("    %d ", sim))
  }
  cat("\n")

  density_data <- rbind(density_data,
    data.frame(m = m, type = "Test Statistic", value = test_stats),
    data.frame(m = m, type = "Permutation Distribution", value = all_perm_stats)
  )
}

density_data$panel <- paste0("m = ", density_data$m, ", n = ", n_fig1)

p1 <- ggplot(density_data, aes(x = value, fill = type, color = type)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~ panel, ncol = 2, scales = "free") +
  labs(x = "x", y = "Density") +
  scale_fill_manual(values = c("Permutation Distribution" = "salmon",
                                "Test Statistic" = "cyan3"),
                    name = "") +
  scale_color_manual(values = c("Permutation Distribution" = "red",
                                 "Test Statistic" = "cyan4"),
                     name = "") +
  theme_minimal() +
  theme(legend.position = "right")

output_dir <- if (test_mode) "/tmp" else "/tmp/output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
ggsave(file.path(output_dir, "figure_5_1.pdf"), p1, width = 10, height = 8)
cat(sprintf("Saved figure_5_1.pdf to %s\n", output_dir))

###############################################################################
# Figure 5.2: QQ plots of p-values vs U[0,1] for n=50 and n=1000
###############################################################################
cat("\n=== Figure 5.2: QQ plots of p-values ===\n")
n_values_fig2 <- if (test_mode) c(50) else c(50, 1000)

qq_data <- data.frame()

for (m in m_values) {
  for (n in n_values_fig2) {
    cat(sprintf("  m = %d, n = %d\n", m, n))
    b_n <- floor(n^(1/3)) + 1
    p_values <- numeric(n_sims)

    set.seed(42)
    for (sim in 1:n_sims) {
      x <- generate_m_dependent(n, m)
      sp <- studentized_permutation_test(x, k = 1, n_perm = n_perm,
                                         b_n = b_n, epsilon = 1e-6)
      p_values[sim] <- sp$p_value
      if (sim %% 1000 == 0) cat(sprintf("    %d ", sim))
    }
    cat("\n")

    p_sorted <- sort(p_values)
    theoretical <- (1:n_sims) / (n_sims + 1)

    qq_data <- rbind(qq_data, data.frame(
      m = m, n = n,
      theoretical = theoretical,
      sample = p_sorted,
      panel = paste0("m = ", m, ", n = ", n)
    ))
  }
}

p2 <- ggplot(qq_data, aes(x = theoretical, y = sample)) +
  geom_point(size = 0.5, color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~ panel, ncol = 2, scales = "fixed") +
  labs(x = "U[0, 1] Quantiles", y = "Sample Quantiles") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_minimal()

ggsave(file.path(output_dir, "figure_5_2.pdf"), p2, width = 10, height = 12)
cat(sprintf("Saved figure_5_2.pdf to %s\n", output_dir))
