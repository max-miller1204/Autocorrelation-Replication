###############################################################################
# Section 6: Financial Application
# S&P 500 and Apple daily closing prices (2010-2019)
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
test_mode <- "--test" %in% args

source("/tmp/core_functions.R")

if (!requireNamespace("quantmod", quietly = TRUE)) {
  install.packages("quantmod", repos = "https://cloud.r-project.org")
}
library(quantmod)
library(ggplot2)

if (test_mode) {
  n_perm <- 200
  cat("Running in TEST mode (n_perm=200)\n")
} else {
  n_perm <- 2000
}

output_dir <- if (test_mode) "/tmp" else "/tmp/output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

###############################################################################
# Download data
###############################################################################
cat("Downloading S&P 500 data...\n")
getSymbols("^GSPC", src = "yahoo", from = "2010-01-01", to = "2019-12-31",
           auto.assign = TRUE)
cat("Downloading Apple data...\n")
getSymbols("AAPL", src = "yahoo", from = "2010-01-01", to = "2019-12-31",
           auto.assign = TRUE)

# Compute daily log-returns
spx_prices <- Cl(GSPC)
aapl_prices <- Cl(AAPL)

spx_returns <- diff(log(spx_prices))[-1]
aapl_returns <- diff(log(aapl_prices))[-1]

spx <- as.numeric(spx_returns)
aapl <- as.numeric(aapl_returns)

cat(sprintf("SPX: %d observations\n", length(spx)))
cat(sprintf("AAPL: %d observations\n", length(aapl)))

###############################################################################
# Table 6.1: Marginal p-values for lags k=1..10
###############################################################################
cat("\n=== Table 6.1: Marginal p-values ===\n")
max_lag <- 10

spx_pvalues <- numeric(max_lag)
aapl_pvalues <- numeric(max_lag)

for (k in 1:max_lag) {
  cat(sprintf("  Lag k = %d ... ", k))

  set.seed(42)
  sp_spx <- studentized_permutation_test(spx, k = k, n_perm = n_perm)
  spx_pvalues[k] <- sp_spx$p_value

  set.seed(42)
  sp_aapl <- studentized_permutation_test(aapl, k = k, n_perm = n_perm)
  aapl_pvalues[k] <- sp_aapl$p_value

  cat(sprintf("SPX=%.4f, AAPL=%.4f\n", spx_pvalues[k], aapl_pvalues[k]))
}

table_6_1 <- data.frame(
  k = 1:max_lag,
  SPX = round(spx_pvalues, 4),
  AAPL = round(aapl_pvalues, 4)
)

cat("\nTable 6.1: Marginal p-values (studentized permutation test)\n")
print(table_6_1, row.names = FALSE)

# Bonferroni correction
cat(sprintf("\nBonferroni threshold: alpha/r = 0.05/10 = 0.005\n"))
cat(sprintf("SPX min p-value: %.4f (reject: %s)\n",
            min(spx_pvalues), min(spx_pvalues) < 0.005))
cat(sprintf("AAPL min p-value: %.4f (reject: %s)\n",
            min(aapl_pvalues), min(aapl_pvalues) < 0.005))

# Ljung-Box comparison
lb_spx <- ljung_box_test_multi(spx, r = 10)
lb_aapl <- ljung_box_test_multi(aapl, r = 10)
cat(sprintf("\nLjung-Box p-values (r=10): SPX=%.4f, AAPL=%.4f\n",
            lb_spx$p_value, lb_aapl$p_value))

write.csv(table_6_1, file.path(output_dir, "table_6_1.csv"), row.names = FALSE)

###############################################################################
# Figure 6.1: Log-return time series plots
###############################################################################
cat("\n=== Figure 6.1: Log-return plots ===\n")

df_returns <- data.frame(
  t = rep(1:length(spx), 2),
  R_t = c(spx, aapl),
  series = rep(c("S&P 500", "Apple"), each = length(spx))
)
# Handle different lengths
if (length(spx) != length(aapl)) {
  df_returns <- rbind(
    data.frame(t = 1:length(spx), R_t = spx, series = "S&P 500"),
    data.frame(t = 1:length(aapl), R_t = aapl, series = "Apple")
  )
}

p1 <- ggplot(df_returns, aes(x = t, y = R_t)) +
  geom_line(color = "blue", linewidth = 0.3) +
  facet_wrap(~ series, ncol = 1, scales = "free_y") +
  labs(x = "t", y = expression(R[t])) +
  theme_minimal()

ggsave(file.path(output_dir, "figure_6_1.pdf"), p1, width = 8, height = 8)
cat("Saved figure_6_1.pdf\n")

###############################################################################
# Figure 6.2: Sample ACF plots with confidence bands
###############################################################################
cat("\n=== Figure 6.2: ACF plots ===\n")

acf_data <- data.frame()
for (series_name in c("S&P 500", "Apple")) {
  x <- if (series_name == "S&P 500") spx else aapl
  n <- length(x)
  acf_vals <- numeric(max_lag)
  for (k in 1:max_lag) {
    acf_vals[k] <- sample_autocorrelation(x, k)
  }
  ci_bound <- qnorm(0.975) / sqrt(n)
  acf_data <- rbind(acf_data, data.frame(
    lag = 1:max_lag,
    acf = acf_vals,
    ci_upper = ci_bound,
    ci_lower = -ci_bound,
    series = series_name
  ))
}

p2 <- ggplot(acf_data, aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_segment(aes(x = lag, xend = lag, y = 0, yend = acf), linewidth = 0.8) +
  geom_hline(aes(yintercept = ci_upper), linetype = "dotted", color = "blue") +
  geom_hline(aes(yintercept = ci_lower), linetype = "dotted", color = "blue") +
  facet_wrap(~ series, ncol = 1) +
  labs(x = "Lag", y = "ACF") +
  scale_x_continuous(breaks = 1:10) +
  theme_minimal()

ggsave(file.path(output_dir, "figure_6_2.pdf"), p2, width = 8, height = 8)
cat("Saved figure_6_2.pdf\n")
