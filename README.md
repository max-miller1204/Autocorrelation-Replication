# Replication: Permutation Testing for Dependence in Time Series

Replication of Romano & Tirlea (2020), "Permutation Testing for Dependence in Time Series" (Stanford University).

The paper proposes a **studentized permutation test** for autocorrelation that is exact under i.i.d. and asymptotically valid under weak dependence, outperforming classical Ljung-Box and Box-Pierce tests in non-i.i.d. settings.

## Project Structure

```
R/
  core_functions.R        # Core statistical functions (autocorrelation, variance
                          #   estimators, permutation tests, Ljung-Box, Box-Pierce)
  data_generators.R       # Data generation (m-dependent, AR(2) variants)
  sim_table_5_1.R         # Table 5.1: m-dependent Gaussian product simulation
  sim_table_5_2.R         # Table 5.2: alpha-mixing AR(2) simulation
  fig_5_1_5_2.R           # Figures 5.1 (density plots) and 5.2 (QQ plots)
  fig_5_3.R               # Figure 5.3: power curves under local alternatives
  financial_application.R # Section 6: S&P 500 and Apple stock analysis
  sim_table_7_1.R         # Table 7.1: autocovariance test
output/
  table_5_1.csv           # Null rejection probabilities (m-dependent)
  table_5_2.csv           # Null rejection probabilities (alpha-mixing AR(2))
  table_6_1.csv           # Financial data marginal p-values
  table_7_1.csv           # Autocovariance test rejection probabilities
  figure_5_1.pdf          # Kernel density: test statistic vs permutation distribution
  figure_5_2.pdf          # QQ plots of p-values vs U[0,1]
  figure_5_3.pdf          # Power curves under local alternatives
  figure_6_1.pdf          # Daily log-return time series (S&P 500, Apple)
  figure_6_2.pdf          # Sample ACF plots with confidence bands
```

## Requirements

- R 4.4.3 (provided via devcontainer)
- Packages: `tidyverse`, `ggplot2`, `quantmod`

## Running the Simulations

All scripts run inside the Docker devcontainer (`peaceful_lamport`). Copy the source files and execute:

```bash
# Copy core libraries
docker cp R/core_functions.R peaceful_lamport:/tmp/
docker cp R/data_generators.R peaceful_lamport:/tmp/

# Copy and run a simulation script
docker cp R/sim_table_5_1.R peaceful_lamport:/tmp/
docker exec peaceful_lamport Rscript /tmp/sim_table_5_1.R
```

### Test Mode

Every script accepts a `--test` flag for quick validation (100 sims, 200 permutations, n={50,100}):

```bash
docker exec peaceful_lamport Rscript /tmp/sim_table_5_1.R --test
```

### Full Run Parameters

| Parameter      | Test Mode | Full Run  |
|----------------|-----------|-----------|
| Simulations    | 100       | 10,000    |
| Permutations   | 200       | 2,000     |
| Sample sizes   | 50, 100   | 10, 20, 50, 80, 100, 500, 1000 |

Full runs take several hours per table/figure.

## Implemented Equations

| Component | Paper Reference |
|-----------|----------------|
| Sample autocorrelation | Eq. (1.6) |
| Sample variance | Eq. (1.7) |
| Variance estimator components (K-hat, T-hat, nu-hat) | Eq. (2.11) |
| Combined gamma-hat-squared | Eq. (2.12) |
| Generalized lag-k test | Eqs. (3.35)-(3.36) |
| Ljung-Box statistic | Eq. (5.1) |
| Box-Pierce statistic | Eq. (5.2) |
| m-dependent Gaussian products | Example 2.1, Eq. (2.17) |
| AR(2) with rho(1)=0 | Example 3.2, Eq. (3.27) |
| AR(2) product process | Eq. (5.4) |
| Autocovariance test | Theorem 7.1, Eqs. (7.5)-(7.6) |
| Local alternatives power | Example 3.3, Eq. (3.44) |

## Key Results

The replication confirms the paper's findings:

- **Studentized permutation test** maintains rejection probability near the nominal 0.05 across all dependence settings.
- **Unstudentized permutation test** over-rejects substantially under m-dependence (up to ~0.33 for m=3).
- **Ljung-Box and Box-Pierce** over-reject dramatically in non-i.i.d. settings (up to ~0.60 for m=3).
- **Financial application**: Ljung-Box rejects for S&P 500 (p=0.001) but the studentized permutation test does not, suggesting the Ljung-Box rejection is a false positive driven by non-i.i.d. dependence structure in returns.

## Reference

Romano, J. P. & Tirlea, M. A. (2020). Permutation Testing for Dependence in Time Series. *Stanford University Working Paper*.
