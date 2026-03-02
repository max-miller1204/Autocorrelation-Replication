###############################################################################
# Data Generation Functions for Romano & Tirlea (2020) Replication
###############################################################################

# (a) m-dependent Gaussian product process (Example 2.1, Eq. 2.17)
# X_i = prod(Z_j, j=i..i+m-1) where Z_j are i.i.d. N(0,1)
# For m=0, X_i = Z_i (i.i.d. standard normal)
generate_m_dependent <- function(n, m = 0) {
  if (m == 0) {
    return(rnorm(n))
  }
  # Need Z_1, ..., Z_{n+m-1}
  Z <- rnorm(n + m - 1)
  X <- numeric(n)
  for (i in 1:n) {
    X[i] <- prod(Z[i:(i + m - 1)])
  }
  X
}

# (b) AR(2) process with rho(1)=0 (Example 3.2, Eq. 3.27)
# X_t = phi*X_{t-1} + rho*X_{t-2} + epsilon_t
# For rho(1)=0: phi/(1-rho) = 0, so phi=0
# Thus X_t = rho*X_{t-2} + epsilon_t with rho=0.5
# innovation_type: "gaussian", "uniform", "t"
generate_ar2 <- function(n, rho_param = 0.5, innovation_type = "gaussian",
                         t_df = 9.5, burn_in = 500) {
  n_total <- n + burn_in

  # Generate innovations
  eps <- switch(innovation_type,
    "gaussian" = rnorm(n_total),
    "uniform"  = runif(n_total, -1, 1),
    "t"        = rt(n_total, df = t_df),
    stop("Unknown innovation type")
  )

  # phi=0, so X_t = rho*X_{t-2} + epsilon_t
  X <- numeric(n_total)
  X[1] <- eps[1]
  X[2] <- eps[2]
  for (t in 3:n_total) {
    X[t] <- rho_param * X[t - 2] + eps[t]
  }

  # Discard burn-in
  X[(burn_in + 1):n_total]
}

# (c) AR(2) product process (Eq. 5.4)
# {X_{2i}} and {X_{2i-1}} are independent stationary AR(1) processes
# with parameter rho, and X_{2i} = Y_{2i} * Y_{2(i+1)}
# More precisely: the even-indexed and odd-indexed subsequences are
# independent AR(1)(rho) processes. Then form product:
# X_{2i} = Y_{2i} * Y_{2(i+1)} where Y is stationary AR(1)(rho)
# with standard Gaussian innovations
generate_ar2_product <- function(n, rho_param = 0.5, burn_in = 500) {
  # We need two independent AR(1) processes
  n_half <- ceiling(n / 2) + 1  # extra for product
  n_total <- n_half + burn_in

  generate_ar1 <- function(nn, rho) {
    eps <- rnorm(nn)
    Y <- numeric(nn)
    Y[1] <- rnorm(1) / sqrt(1 - rho^2)
    for (t in 2:nn) {
      Y[t] <- rho * Y[t - 1] + eps[t]
    }
    Y
  }

  Y_even <- generate_ar1(n_total, rho_param)[(burn_in + 1):n_total]
  Y_odd <- generate_ar1(n_total, rho_param)[(burn_in + 1):n_total]

  # Construct the interleaved process
  X <- numeric(n)
  for (i in seq_len(ceiling(n / 2))) {
    idx_odd <- 2 * i - 1
    idx_even <- 2 * i
    if (idx_odd <= n) X[idx_odd] <- Y_odd[i]
    if (idx_even <= n) X[idx_even] <- Y_even[i] * Y_even[i + 1]
  }
  X
}

# (d) AR(2) with uniform innovations U[-1,1]: use generate_ar2 with innovation_type="uniform"
# (e) AR(2) with t_9.5 innovations: use generate_ar2 with innovation_type="t"

###############################################################################
# Self-test
###############################################################################
if (sys.nframe() == 0) {
  cat("Running self-tests for data_generators.R\n")

  set.seed(42)

  # Test m-dependent
  for (m in 0:3) {
    x <- generate_m_dependent(10000, m)
    rho1 <- cor(x[1:(length(x) - 1)], x[2:length(x)])
    cat(sprintf("  m=%d: mean=%.3f, sd=%.3f, rho(1)~=%.4f (expect ~0)\n",
                m, mean(x), sd(x), rho1))
  }

  # Test AR(2) with rho=0.5, phi=0
  x_ar2 <- generate_ar2(10000, 0.5, "gaussian")
  rho1_ar2 <- cor(x_ar2[1:(length(x_ar2) - 1)], x_ar2[2:length(x_ar2)])
  rho2_ar2 <- cor(x_ar2[1:(length(x_ar2) - 2)], x_ar2[3:length(x_ar2)])
  cat(sprintf("  AR(2) Gaussian: rho(1)~=%.4f (expect ~0), rho(2)~=%.4f (expect ~0.5)\n",
              rho1_ar2, rho2_ar2))
  stopifnot(abs(rho1_ar2) < 0.1)
  stopifnot(abs(rho2_ar2 - 0.5) < 0.15)

  # Test AR(2) uniform
  x_unif <- generate_ar2(5000, 0.5, "uniform")
  rho1_u <- cor(x_unif[1:(length(x_unif) - 1)], x_unif[2:length(x_unif)])
  cat(sprintf("  AR(2) Uniform: rho(1)~=%.4f (expect ~0)\n", rho1_u))

  # Test AR(2) t
  x_t <- generate_ar2(5000, 0.5, "t", t_df = 9.5)
  rho1_t <- cor(x_t[1:(length(x_t) - 1)], x_t[2:length(x_t)])
  cat(sprintf("  AR(2) t_9.5: rho(1)~=%.4f (expect ~0)\n", rho1_t))

  # Test AR(2) product
  x_prod <- generate_ar2_product(5000, 0.5)
  rho1_p <- cor(x_prod[1:(length(x_prod) - 1)], x_prod[2:length(x_prod)])
  cat(sprintf("  AR(2) Product: rho(1)~=%.4f (expect ~0)\n", rho1_p))

  cat("All self-tests passed!\n")
}
