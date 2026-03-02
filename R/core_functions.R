###############################################################################
# Core Statistical Functions for Romano & Tirlea (2020) Replication
# "Permutation Testing for Dependence in Time Series"
###############################################################################

# Sample autocorrelation at lag k (Eq. 1.6)
# rho_hat_n(k) = (1/(n-k)) * sum_{i=1}^{n-k} (X_i - X_bar)(X_{i+k} - X_bar) / sigma_hat^2
sample_autocorrelation <- function(x, k = 1) {
  n <- length(x)
  x_bar <- mean(x)
  sigma_sq <- mean((x - x_bar)^2)  # Eq. 1.7
  if (sigma_sq == 0) return(0)
  numerator <- sum((x[1:(n - k)] - x_bar) * (x[(1 + k):n] - x_bar)) / (n - k)
  numerator / sigma_sq
}

# Sample variance (Eq. 1.7)
sample_variance <- function(x) {
  mean((x - mean(x))^2)
}

# Sample autocovariance at lag k
# c_hat_n(k) = (1/(n-k)) * sum_{i=1}^{n-k} (X_i - X_bar)(X_{i+k} - X_bar)
sample_autocovariance <- function(x, k = 1) {
  n <- length(x)
  x_bar <- mean(x)
  sum((x[1:(n - k)] - x_bar) * (x[(1 + k):n] - x_bar)) / (n - k)
}

# Variance estimator components (Eq. 2.11)
# Y_i = (X_i - X_bar)(X_{i+1} - X_bar), Z_i = (X_i - X_bar)^2
# Generalized for lag k: Y_i = (X_i - X_bar)(X_{i+k} - X_bar)
variance_estimator_components <- function(x, k = 1, b_n = NULL) {
  n <- length(x)
  if (is.null(b_n)) b_n <- floor(n^(1/3)) + 1

  x_bar <- mean(x)
  xc <- x - x_bar

  # Y_i = (X_i - X_bar)(X_{i+k} - X_bar) for i = 1, ..., n-k
  Y <- xc[1:(n - k)] * xc[(1 + k):n]
  n_Y <- length(Y)
  Y_bar <- mean(Y)

  # Z_i = (X_i - X_bar)^2 for i = 1, ..., n
  Z <- xc^2
  Z_bar <- mean(Z)

  # K_hat^2 (variance estimator for Z)
  K_sq <- (1 / n) * sum((Z - Z_bar)^2)
  for (j in 1:b_n) {
    cross_sum <- sum((Z[1:(n - j)] - Z_bar) * (Z[(1 + j):n] - Z_bar))
    K_sq <- K_sq + (2 / n) * cross_sum
  }

  # T_hat^2 (variance estimator for Y)
  T_sq <- (1 / n) * sum((Y - Y_bar)^2)
  for (j in 1:b_n) {
    idx_max <- n_Y - j
    if (idx_max < 1) break
    cross_sum <- sum((Y[1:idx_max] - Y_bar) * (Y[(1 + j):(idx_max + j)] - Y_bar))
    T_sq <- T_sq + (2 / n) * cross_sum
  }

  # nu_hat (cross-covariance estimator between Y and Z)
  # First term: (1/n) * sum_{i=1}^{n-k} (Y_i - Y_bar)(Z_i - Z_bar)
  n_common <- min(n_Y, n)
  nu <- (1 / n) * sum((Y[1:n_common] - Y_bar) * (Z[1:n_common] - Z_bar))

  # Lagged cross terms
  for (j in 1:b_n) {
    # sum (Z_i - Z_bar)(Y_{i+j} - Y_bar) for valid indices
    idx_max_zy <- min(n - j, n_Y)
    if (idx_max_zy >= 1 && (1 + j) <= n_Y) {
      end_z <- min(n - j, n_Y - j)
      if (end_z >= 1) {
        nu <- nu + (1 / n) * sum((Z[1:end_z] - Z_bar) * (Y[(1 + j):(end_z + j)] - Y_bar))
      }
    }

    # sum (Y_i - Y_bar)(Z_{i+j} - Z_bar) for valid indices
    idx_max_yz <- min(n_Y, n - j)
    if (idx_max_yz >= 1 && (1 + j) <= n) {
      end_y <- min(n_Y, n - j)
      if (end_y >= 1) {
        nu <- nu + (1 / n) * sum((Y[1:end_y] - Y_bar) * (Z[(1 + j):(end_y + j)] - Z_bar))
      }
    }
  }

  list(K_sq = K_sq, T_sq = T_sq, nu = nu)
}

# Combined gamma_hat^2 (Eq. 2.12 / Eq. 3.36 for general lag k)
# gamma_hat^2 = (1/sigma_hat^4) * [T_hat^2 - rho_hat * nu_hat + rho_hat^2 * K_hat^2]
gamma_hat_squared <- function(x, k = 1, b_n = NULL, epsilon = 1e-6) {
  sigma_sq <- sample_variance(x)
  if (sigma_sq == 0) return(epsilon)
  rho_hat <- sample_autocorrelation(x, k)
  comp <- variance_estimator_components(x, k, b_n)
  result <- (1 / sigma_sq^2) * (comp$T_sq - rho_hat * comp$nu + rho_hat^2 * comp$K_sq)
  max(result, epsilon)  # truncate at epsilon
}

# Studentized permutation test
# Test statistic: sqrt(n) * rho_hat / gamma_hat
# One-sided test: reject if test stat is too large (positive autocorrelation)
# Returns p-value
studentized_permutation_test <- function(x, k = 1, n_perm = 2000,
                                         b_n = NULL, epsilon = 1e-6) {
  n <- length(x)
  if (is.null(b_n)) b_n <- floor(n^(1/3)) + 1

  rho_obs <- sample_autocorrelation(x, k)
  gamma_sq_obs <- gamma_hat_squared(x, k, b_n, epsilon)
  gamma_obs <- sqrt(gamma_sq_obs)
  T_obs <- sqrt(n) * rho_obs / gamma_obs

  # Permutation distribution
  T_perm <- numeric(n_perm)
  for (p in 1:n_perm) {
    x_perm <- sample(x)
    rho_perm <- sample_autocorrelation(x_perm, k)
    gamma_sq_perm <- gamma_hat_squared(x_perm, k, b_n, epsilon)
    gamma_perm <- sqrt(gamma_sq_perm)
    T_perm[p] <- sqrt(n) * rho_perm / gamma_perm
  }

  # One-sided p-value (upper tail)
  p_value <- mean(T_perm >= T_obs)
  list(p_value = p_value, T_obs = T_obs, T_perm = T_perm)
}

# Unstudentized permutation test
# Test statistic: sqrt(n) * rho_hat
unstudentized_permutation_test <- function(x, k = 1, n_perm = 2000) {
  n <- length(x)
  rho_obs <- sample_autocorrelation(x, k)
  T_obs <- sqrt(n) * rho_obs

  T_perm <- numeric(n_perm)
  for (p in 1:n_perm) {
    x_perm <- sample(x)
    T_perm[p] <- sqrt(n) * sample_autocorrelation(x_perm, k)
  }

  p_value <- mean(T_perm >= T_obs)
  list(p_value = p_value, T_obs = T_obs, T_perm = T_perm)
}

# Ljung-Box test (Eq. 5.1) for r=1
# Q_LB = n(n+2) * rho_hat^2 / (n-1)
# Compare to chi-squared(1) upper tail
ljung_box_test <- function(x, k = 1) {
  n <- length(x)
  rho_hat <- sample_autocorrelation(x, k)
  Q <- n * (n + 2) * rho_hat^2 / (n - k)
  p_value <- 1 - pchisq(Q, df = 1)
  list(p_value = p_value, Q = Q, rho_hat = rho_hat)
}

# Box-Pierce test (Eq. 5.2) for r=1
# Q_BP = n * rho_hat^2
# Compare to chi-squared(1) upper tail
box_pierce_test <- function(x, k = 1) {
  n <- length(x)
  rho_hat <- sample_autocorrelation(x, k)
  Q <- n * rho_hat^2
  p_value <- 1 - pchisq(Q, df = 1)
  list(p_value = p_value, Q = Q, rho_hat = rho_hat)
}

# Studentized autocovariance permutation test (Section 7)
# Test statistic: sqrt(n) * c_hat / T_hat (where T_hat is from Eq. 2.11)
studentized_autocovariance_test <- function(x, k = 1, n_perm = 2000,
                                            b_n = NULL, epsilon = 1e-6) {
  n <- length(x)
  if (is.null(b_n)) b_n <- floor(n^(1/3)) + 1

  c_obs <- sample_autocovariance(x, k)
  comp_obs <- variance_estimator_components(x, k, b_n)
  T_hat_obs <- sqrt(max(comp_obs$T_sq, epsilon))
  stat_obs <- sqrt(n) * c_obs / T_hat_obs

  stat_perm <- numeric(n_perm)
  for (p in 1:n_perm) {
    x_perm <- sample(x)
    c_perm <- sample_autocovariance(x_perm, k)
    comp_perm <- variance_estimator_components(x_perm, k, b_n)
    T_hat_perm <- sqrt(max(comp_perm$T_sq, epsilon))
    stat_perm[p] <- sqrt(n) * c_perm / T_hat_perm
  }

  p_value <- mean(stat_perm >= stat_obs)
  list(p_value = p_value, stat_obs = stat_obs, stat_perm = stat_perm)
}

# Ljung-Box test for multiple lags (portmanteau)
ljung_box_test_multi <- function(x, r = 10) {
  n <- length(x)
  rho_sq <- numeric(r)
  for (k in 1:r) {
    rho_sq[k] <- sample_autocorrelation(x, k)^2 / (n - k)
  }
  Q <- n * (n + 2) * sum(rho_sq)
  p_value <- 1 - pchisq(Q, df = r)
  list(p_value = p_value, Q = Q)
}

###############################################################################
# Self-test
###############################################################################
if (sys.nframe() == 0) {
  cat("Running self-tests for core_functions.R\n")

  set.seed(42)
  x_iid <- rnorm(200)

  # Test sample autocorrelation: should be near 0 for i.i.d. data
  rho1 <- sample_autocorrelation(x_iid, 1)
  cat(sprintf("  rho_hat(1) for i.i.d. N(0,1): %.4f (expect ~0)\n", rho1))
  stopifnot(abs(rho1) < 0.2)

  # Test sample variance
  sv <- sample_variance(x_iid)
  cat(sprintf("  sigma_hat^2 for N(0,1): %.4f (expect ~1)\n", sv))
  stopifnot(abs(sv - 1) < 0.5)

  # Test gamma_hat_squared
  g2 <- gamma_hat_squared(x_iid, 1)
  cat(sprintf("  gamma_hat^2 for i.i.d.: %.4f (expect ~1)\n", g2))

  # Test AR(1) with known autocorrelation
  set.seed(42)
  phi <- 0.5
  n_ar <- 500
  x_ar <- numeric(n_ar)
  x_ar[1] <- rnorm(1)
  for (i in 2:n_ar) x_ar[i] <- phi * x_ar[i - 1] + rnorm(1)
  rho_ar <- sample_autocorrelation(x_ar, 1)
  cat(sprintf("  rho_hat(1) for AR(1) phi=0.5: %.4f (expect ~0.5)\n", rho_ar))
  stopifnot(abs(rho_ar - 0.5) < 0.15)

  # Test Ljung-Box on i.i.d.: should not reject
  lb <- ljung_box_test(x_iid, 1)
  cat(sprintf("  Ljung-Box p-value for i.i.d.: %.4f (expect >0.05)\n", lb$p_value))

  # Test Box-Pierce on i.i.d.
  bp <- box_pierce_test(x_iid, 1)
  cat(sprintf("  Box-Pierce p-value for i.i.d.: %.4f (expect >0.05)\n", bp$p_value))

  # Quick permutation test (small n_perm for speed)
  sp <- studentized_permutation_test(x_iid[1:50], 1, n_perm = 200)
  cat(sprintf("  Stud. perm. p-value for i.i.d. (n=50): %.4f\n", sp$p_value))

  up <- unstudentized_permutation_test(x_iid[1:50], 1, n_perm = 200)
  cat(sprintf("  Unstud. perm. p-value for i.i.d. (n=50): %.4f\n", up$p_value))

  cat("All self-tests passed!\n")
}
