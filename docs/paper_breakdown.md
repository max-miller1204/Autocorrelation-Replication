# Paper Breakdown: "Permutation Testing for Dependence in Time Series"

**Authors:** Joseph P. Romano & Marius A. Tirlea (Stanford University, 2020)

---

## One-Sentence Summary

This paper shows that the standard way of using permutation tests to check for autocorrelation in time series data can give badly wrong answers, and proposes a fix ("studentization") that makes the test reliable for a broad class of dependent data.

---

## 1. What Problem Does This Paper Solve?

### The Setup
Imagine you have a time series -- a sequence of data points over time (e.g., daily stock returns). A fundamental question is: **are consecutive observations correlated with each other?** This is called *autocorrelation*.

- **Autocorrelation = 0** means knowing today's value tells you nothing about tomorrow's.
- **Autocorrelation != 0** means there's a pattern -- past values carry information about future ones.

### The Existing Tools
There are classic tests for this (Ljung-Box, Box-Pierce), but they **assume the data follows a specific model** (like an ARMA process). If your data doesn't fit that model, these tests can be unreliable.

**Permutation tests** are appealing because they're *nonparametric* -- they don't assume a specific data-generating model. The idea: shuffle your data randomly many times, compute the autocorrelation each time, and see if the real autocorrelation is unusually large compared to the shuffled versions.

### The Problem
Permutation tests work perfectly when data points are **independent and identically distributed (i.i.d.)** -- that is, each observation is drawn independently from the same distribution. But here's the catch:

> **"Zero autocorrelation" and "independence" are NOT the same thing.**

A sequence can have zero autocorrelation (no linear relationship between consecutive terms) but still have dependence in other ways (e.g., the *variance* of the process might change over time). When this happens, a naive permutation test can:

1. **Type 1 error problems**: Reject the null hypothesis way too often (or too rarely) -- meaning the test's stated significance level (e.g., 5%) is a lie.
2. **Type 3 / directional errors**: Conclude that autocorrelation is *positive* when it's actually *negative* (or vice versa). This is perhaps even worse than a simple false positive.

---

## 2. The Key Insight: Why Naive Permutation Tests Fail

When you shuffle (permute) a dependent time series, the shuffled version behaves like i.i.d. data. The permutation distribution of the test statistic converges to a normal distribution with **variance 1**.

But the *actual* test statistic (computed from the un-shuffled, dependent data) converges to a normal distribution with a **different variance**, called gamma-1-squared. This variance depends on the dependence structure of the data.

**If these two variances don't match, the permutation test is invalid.** The test is comparing your statistic against the wrong reference distribution.

The only case where they automatically match is when the data truly is i.i.d. (where gamma-1-squared = 1 by definition).

---

## 3. The Fix: Studentization

The solution is elegant: **divide the test statistic by an estimate of its standard deviation** (this is called "studentization" -- think of it like computing a t-statistic instead of a z-statistic).

Specifically, instead of using the raw test statistic sqrt(n) * rho-hat (where rho-hat is the sample autocorrelation), they use:

> **sqrt(n) * rho-hat / gamma-hat**

where gamma-hat is an estimator of the standard deviation that accounts for the dependence structure.

### Why This Works
After studentization:
- The **actual test statistic** converges to N(0, 1) -- standard normal.
- The **permutation distribution** also converges to N(0, 1).

Since both converge to the same distribution, the permutation test becomes asymptotically valid -- its rejection probability converges to the nominal level alpha.

### The Best of Both Worlds
- **When data is truly i.i.d.**: The permutation test is *exactly* level alpha (not just approximately). This is a property that Ljung-Box and Box-Pierce don't have.
- **When data is dependent but uncorrelated**: The test is *asymptotically* valid (approaches the correct level as sample size grows).

---

## 4. Two Classes of Dependence Considered

The paper proves results for two progressively broader classes of time series:

### a) m-dependent sequences (Section 2)
- **What it means**: Observations more than m time steps apart are independent. For example, if m=2, then X_1 and X_4 are independent, but X_1 and X_3 might not be.
- **Intuition**: The dependence has a finite "memory."
- **Technical requirement**: Finite 8th moments (for general m-dependence) or finite 4th moments (for i.i.d.).
- **Proof technique**: Uses Stein's method to show that sums of locally dependent random variables are approximately normal.

### b) Alpha-mixing sequences (Section 3)
- **What it means**: Observations far apart in time are *approximately* independent, with the dependence dying off. This covers a much larger class, including all stationary ARMA processes.
- **Intuition**: The dependence fades gradually rather than cutting off sharply.
- **Technical requirement**: Finite (8 + 4*delta) moments and summable mixing coefficients.
- **Proof technique**: Uses a truncation argument and the central limit theorem of Wald and Wolfowitz (1943).

---

## 5. Constructing the Studentized Estimator (gamma-hat-squared)

The variance estimator gamma-hat-squared is built from three components:

1. **T-hat-squared**: Estimates the variance of the product terms (X_i - X-bar)(X_{i+1} - X-bar), using a HAC-style (Heteroskedasticity and Autocorrelation Consistent) windowed sum.
2. **nu-hat**: Estimates a cross-covariance term between the product terms and the squared terms.
3. **K-hat-squared**: Estimates the variance of the squared terms (X_i - X-bar)^2, again using a windowed sum.

These are combined as:

> gamma-hat-squared = (1 / sigma-hat^4) * [T-hat-squared - rho-hat * nu-hat + rho-hat-squared * K-hat-squared]

The window width b_n must grow to infinity but slower than sqrt(n). In practice, b_n = floor(n^{1/3}) + 1 works well.

---

## 6. Multiple Testing (Section 4)

Often you want to test autocorrelation at multiple lags simultaneously (e.g., is rho(1) = rho(2) = ... = rho(10) = 0?). The paper shows:

- You can run separate permutation tests for each lag k and combine the p-values using **Bonferroni correction** (reject if any p-value < alpha/r, where r is the number of lags tested).
- This is **not overly conservative** because the sample autocorrelations at different lags are asymptotically independent (in the i.i.d. case) or nearly so.
- The Sidak correction can also be used and is only marginally more powerful.

---

## 7. Simulation Results (Section 5)

The paper validates the theory with Monte Carlo simulations comparing four tests:
- **Studentized permutation test** (their proposal)
- **Unstudentized permutation test** (the naive version)
- **Ljung-Box test**
- **Box-Pierce test**

### Key findings (at nominal level alpha = 0.05):

| Scenario | Stud. Perm. | Unstud. Perm. | Ljung-Box | Box-Pierce |
|---|---|---|---|---|
| i.i.d. (m=0) | ~0.05 | ~0.05 | ~0.05 | ~0.04 |
| m-dependent (m=1,2,3) | ~0.05-0.07 | 0.10-0.33 | 0.09-0.60 | 0.05-0.60 |
| AR(2) processes | ~0.04-0.05 | 0.06-0.17 | 0.10-0.27 | 0.08-0.26 |

**Takeaways:**
- The **studentized permutation test** stays near 0.05 across all settings.
- The **unstudentized permutation test** gets progressively worse as dependence increases.
- **Ljung-Box and Box-Pierce** perform terribly in non-i.i.d. settings -- rejection rates of 20-60% when they should be 5%.

---

## 8. Application to Financial Data (Section 6)

The test is applied to daily log-returns of the **S&P 500 index** and **Apple stock** (2010-2019) to test the Efficient Market Hypothesis (which implies returns should be uncorrelated).

- **Permutation test**: After Bonferroni correction for 10 lags, no significant autocorrelation is found for either series. Conclusion: no evidence against market efficiency.
- **Ljung-Box test**: Rejects the null for the S&P 500 (p = 0.001), suggesting autocorrelation exists. But given the simulation evidence that Ljung-Box over-rejects in non-i.i.d. settings, this is likely a false positive.

This is a compelling real-world example of why the choice of test matters.

---

## 9. Testing Autocovariance (Section 7)

The paper extends all results to testing whether the first-order **autocovariance** (not autocorrelation) equals zero. Since the sample variance is permutation invariant, the autocovariance-based test produces nearly identical rejection probabilities to the autocorrelation-based test. This is confirmed by additional simulations.

---

## 10. Bottom Line / What to Tell Your Advisor

**The paper's main contribution in plain language:**

1. Standard permutation tests for autocorrelation are **broken** when data is dependent -- they can reject far too often and even get the direction of the correlation wrong.

2. The fix is to **studentize** the test statistic (divide by an appropriate standard deviation estimate). This makes the permutation distribution match the actual distribution of the test statistic.

3. The resulting test has a **unique advantage**: it's *exact* (perfectly calibrated) when data is i.i.d., and *asymptotically valid* (approaches perfect calibration) for a broad class of dependent data including ARMA processes. No other test has both properties.

4. Existing tests like Ljung-Box and Box-Pierce **fail dramatically** in non-i.i.d. settings, with false positive rates reaching 60% at a nominal 5% level. The studentized permutation test stays near 5%.

---

## Glossary of Key Terms

| Term | Plain English |
|---|---|
| **Autocorrelation** | How much a time series is correlated with a lagged version of itself |
| **Permutation test** | A test that shuffles the data to build a reference distribution, rather than assuming a parametric model |
| **Studentization** | Dividing a statistic by an estimate of its standard deviation (like going from Z-score to t-score) |
| **m-dependent** | Observations more than m steps apart are independent |
| **Alpha-mixing** | Observations far apart are approximately independent, with dependence decaying |
| **Type 1 error** | Rejecting the null when it's actually true (false positive) |
| **Type 3 / directional error** | Getting the *direction* of the effect wrong (saying positive when it's negative) |
| **Asymptotically valid** | The test approaches the correct error rate as sample size goes to infinity |
| **Exact level alpha** | The test has precisely the stated error rate, not just approximately |
| **ARMA process** | AutoRegressive Moving Average -- a common parametric time series model |
| **Exchangeable** | The joint distribution doesn't change when you reorder the observations |
| **HAC estimator** | Heteroskedasticity and Autocorrelation Consistent -- a variance estimator designed for dependent data |
