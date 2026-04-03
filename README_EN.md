# kw-approx: Higher Order Asymptotic Approximations for Kruskal-Wallis Statistics

A Python package implementing higher-order asymptotic approximations to the null distribution and critical values of the Kruskal-Wallis test statistic.

## Paper

> **Lee, J.-S., Murakami, H., & Ha, H.-T. (2026).** *Higher Order Asymptotic Approximations of Kruskal-Wallis Statistics Based on Skewness and Kurtosis.* Preprint submitted to Journals.

## Background

The Kruskal-Wallis test is a rank-based, distribution-free alternative to one-way ANOVA for testing the equality of location parameters across $k$ treatment groups.

### Kruskal-Wallis H Statistic

$$H = \frac{12}{N(N+1)} \sum_{i=1}^{k} \frac{R_i^2}{n_i} - 3(N+1)$$

- $N = \sum_{i=1}^{k} n_i$: total sample size
- $n_i$: sample size of group $i$
- $R_i$: rank sum of group $i$

The conventional $\chi^2_{k-1}$ approximation can be inaccurate for small or unbalanced samples. This package provides saddlepoint, Edgeworth, Gram-Charlier, and polynomially adjusted gamma approximations that incorporate exact skewness and kurtosis corrections.

## Project Structure

```
KS-test/
|
|-- kw_approx/                    # Core package
|   |-- __init__.py               #   Package init (v0.3.0)
|   |-- kruskal_wallis.py         #   H statistic computation
|   |-- moments.py                #   Moments/cumulants (exact/simulation)
|   |-- saddlepoint.py            #   Saddlepoint approx (ER1, ER2, Wang, KT + L-R)
|   |-- edgeworth.py              #   Edgeworth expansion (chi-sq + Laguerre)
|   |-- gram_charlier.py          #   Gram-Charlier Type A (Hermite)
|   |-- pam.py                    #   PAG(d) polynomially adjusted gamma
|   |-- exact.py                  #   Exact distribution (recursive enumeration)
|   |-- simulation.py             #   Monte Carlo simulation
|   +-- approximator.py           #   Unified interface (KWApproximator)
|
|-- examples/
|   +-- reproduce_paper_tables.py #   Reproduce paper Tables 4.1-4.7 + extensions
|
|-- tests/
|   +-- test_approximations.py    #   36 test cases
|
|-- result/                       # Outputs
|   |-- section2_3_exact_moments_proof.tex  # Section 2.3 proof supplement
|   |-- errata_formulas.tex                 # Formula errata
|   +-- paper_tables.tex                    # LaTeX tables
|
|-- README.md                     # Korean documentation
|-- README_EN.md                  # English documentation
|-- CLAUDE.md
|-- pyproject.toml
+-- Kruskal_Wallis_Test.pdf       # Paper PDF
```

## Core Pipeline

```
                         +------------------+
                         |   sample_sizes   |
                         | (n1, n2, ..., nk)|
                         +--------+---------+
                                  |
                                  v
                     +------------+-------------+
                     |        KWMoments          |
                     |  (exact enum / MC sim)    |
                     |                           |
                     |  k<=3, N<=15: exact       |
                     |  k>=4, N<=13: exact       |
                     |  otherwise:   simulation  |
                     +------------+-------------+
                                  |
                                  v
                       +--------------------+
                       | Cumulants k1...k4  |
                       | (exact finite-     |
                       |  sample, NOT       |
                       |  asymptotic chi-sq)|
                       +----+----+----+----+
                            |    |    |    |
            +---------------+    |    |    +----------------+
            |           +--------+    +---------+           |
            v           v                       v           v
     +-----------+ +-----------+         +-----------+ +-----------+
     |Saddlepoint| |Edgeworth  |         |Gram-      | |   PAG(d)  |
     |  ER1/ER2  | |(chi-sq +  |         |Charlier   | | (gamma x  |
     | Wang/KT   | | Laguerre) |         |Type A     | | polynomial|
     +-----------+ +-----------+         |(Hermite)  | |  d=4,6)   |
     |  +- CC    |                       +-----------+ +-----------+
     |  variants |
     +-----+-----+
           |
           v
    +-------------+
    |Lugannani-   |
    |Rice tail    |---------->  P(H >= v)
    |probability  |
    +-------------+
```

## Installation

```bash
pip install numpy scipy
pip install -e .
```

## Quick Start

```python
from kw_approx import KWApproximator

# Three groups, 3 observations each
approx = KWApproximator([3, 3, 3])

# Tail probability P(H >= 4.62)
p_value = approx.tail_probability(4.62, method='ER1')
print(f"ER1 P-value: {p_value:.6f}")

# Compare all methods
results = approx.compare_methods(4.62)
for method, p in results.items():
    print(f"{method}: {p:.6f}")
```

## Approximation Methods

### Method Naming (Paper <-> Code)

```
+----------+----------------+---------------------+-----------------+
| Paper    | Code method    | CGF                 | Note            |
+----------+----------------+---------------------+-----------------+
| SD1      | 'ER1'          | Easton-Ronchetti 1  | Saddlepoint     |
| SD2      | 'Wang'         | Wang damped         | Saddlepoint     |
| SDC1     | 'ER1_cc'       | ER1 + CC            | Cont. corrected |
| SDC2     | 'Wang_cc'      | Wang + CC           | Cont. corrected |
| CHI      | 'chi_square'   | -                   | Baseline        |
| ED       | 'edgeworth'    | -                   | Laguerre-based  |
| GC-A     | 'gram_charlier'| -                   | Hermite-based   |
| PAG(4)   | 'pam'          | -                   | Gamma x poly    |
| PAG(6)   | 'pam6'         | -                   | Gamma x poly    |
| -        | 'exact'        | -                   | Small N only    |
| -        | 'simulation'   | -                   | Monte Carlo     |
+----------+----------------+---------------------+-----------------+
```

### Cumulants -- Exact Finite-Sample

Following Section 4 of the paper, all approximation methods use **exact finite-sample cumulants**, not asymptotic chi-square cumulants:

- Small N (k<=3, N<=15 or k>=4, N<=13): computed from exact distribution
- Large N: Monte Carlo simulation (50,000 iterations, deterministic seed)

> **Note:** Asymptotic chi-sq cumulants show large errors for small samples.
> E.g., for (3,3,3): kappa_2: exact 2.72 vs asymptotic 4.0

### Saddlepoint Approximation -- Lugannani-Rice

Daniels (1954) saddlepoint density approximation:

$$f_{SP}(x) = \left(2\pi K_H^{(2)}(\hat{t})\right)^{-1/2} \exp(K_H(\hat{t}) - x\hat{t})$$

**Lugannani-Rice tail probability:**

$$\Pr(H \geq v) \approx 1 - \Phi(\hat{w}) + \phi(\hat{w})\left(\frac{1}{\hat{u}} - \frac{1}{\hat{w}}\right)$$

where:
- $\hat{w} = \text{sgn}(\hat{t})\sqrt{2(\hat{t}v - K_H(\hat{t}))}$
- $\hat{u} = \hat{t}\sqrt{K_H^{(2)}(\hat{t})}$

### CGF Approximation Methods

**ER1 (Easton-Ronchetti 1st):**
$$K_H(t) \approx \sum_{i=1}^{4} \frac{\kappa_i t^i}{i!}$$

**ER2 (Easton-Ronchetti 2nd):**
$$K_H(t) \approx \kappa_1 t + \frac{\kappa_2}{2}t^2 + \log\left(1 + \frac{\kappa_3}{6}t^3 + \frac{3\kappa_4}{72}t^4 + \frac{\kappa_3^2}{72}t^6\right)$$

**Wang (damped):**
$$K_H(t) \approx \kappa_1 t + \frac{\kappa_2}{2}t^2 + \left(\frac{\kappa_3}{6}t^3 + \frac{\kappa_4}{24}t^4\right)\eta_p(t)$$

where $\eta_p(t) = \exp(-\kappa_2 p^2 t^2 / 2)$, and $p$ is the minimum damping parameter ensuring $K''_W(t;p) \geq 0$.

### Edgeworth Expansion

Chi-square baseline with generalized Laguerre polynomial corrections for skewness and kurtosis:

$$F_{ED}(x) = G_{\nu}(x) - g_{\nu}(x)\left[\frac{\gamma_1}{6}L_3^{(\nu/2-1)}(x/2) + \frac{\gamma_2}{24}L_4^{(\nu/2-1)}(x/2) + \frac{\gamma_1^2}{72}L_6^{(\nu/2-1)}(x/2)\right]$$

### Gram-Charlier Type A

Normal baseline with Hermite polynomial corrections ($H_3, H_4, H_6$).

### PAG(d) -- Polynomially Adjusted Gamma

Gamma density multiplied by a degree-$d$ polynomial:

$$f_{PAG}(x; d) = \psi(x) \sum_{i=0}^{d} \xi_i x^i$$

Coefficients $\xi_0, \ldots, \xi_d$ are determined by matching the first $d+1$ moments via moment matrix inversion. Variants $d=4$ and $d=6$ are provided.

## Moments and Cumulants

Under $H_0$:

- **Mean**: $E(H) = k - 1$
- **Variance**: Wallace (1959) exact formula:

$$\text{Var}(H) = 2(k-1) - \frac{2A_W}{5\,N(N+1)} - \frac{6}{5}\sum_{i=1}^{k}\frac{1}{n_i}$$

where $A_W = 3k(k-2) + N(2k^2 - 6k + 1)$.

**Cumulants** (exact finite-sample):
- $\kappa_1 = E(H) = k - 1$
- $\kappa_2 = \text{Var}(H)$
- $\kappa_3 = \mu_3 - 3\mu_2\mu_1 + 2\mu_1^3$
- $\kappa_4 = \mu_4 - 4\mu_3\mu_1 - 3\mu_2^2 + 12\mu_2\mu_1^2 - 6\mu_1^4$

## Detailed Usage

### Individual Approximation Classes

```python
from kw_approx import SaddlepointApproximation, ExactDistribution, MonteCarloSimulation

sample_sizes = [3, 3, 3]

# Saddlepoint with different CGF methods
sp_er1 = SaddlepointApproximation(sample_sizes, cgf_method='ER1')
sp_wang = SaddlepointApproximation(sample_sizes, cgf_method='Wang')

print(f"ER1:  {sp_er1.tail_probability_lr(4.62):.6f}")
print(f"Wang: {sp_wang.tail_probability_lr(4.62):.6f}")

# Exact (small samples only)
exact = ExactDistribution(sample_sizes)
print(f"Exact: {exact.tail_probability(4.62):.6f}")

# Monte Carlo (large samples)
sim = MonteCarloSimulation(sample_sizes, n_simulations=10000, seed=42)
print(f"Simulation: {sim.tail_probability(4.62):.6f}")
```

### Critical Value Computation

```python
from kw_approx import KWApproximator

approx = KWApproximator([5, 5, 5])

for alpha in [0.10, 0.05, 0.01]:
    cv_exact = approx.critical_value(alpha, method='exact')
    cv_chi2 = approx.critical_value(alpha, method='chi_square')
    cv_er1 = approx.critical_value(alpha, method='ER1')
    print(f"alpha={alpha}: Exact={cv_exact:.4f}, Chi2={cv_chi2:.4f}, ER1={cv_er1:.4f}")
```

## Recommended Methods by Sample Size

```
+---------------+----------+---------------+---------------------------+
| Sample size N | Groups k | Recommended   | Note                      |
+---------------+----------+---------------+---------------------------+
| N <= 15       | k <= 3   | exact         | Exact distribution        |
| N <= 13       | k = 4    | exact         | Exact distribution        |
| 15 < N < 100  | any      | ER1, Wang     | Saddlepoint + L-R         |
| N >= 100      | any      | chi_square    | Asymptotic sufficient     |
| any           | any      | simulation    | MC reference (slow)       |
+---------------+----------+---------------+---------------------------+
```

## Testing

```bash
# Run all tests (36 cases)
python -m pytest tests/test_approximations.py -v

# Run a single test class
python -m pytest tests/test_approximations.py::TestKWMoments -v

# Reproduce paper tables
uv run python examples/reproduce_paper_tables.py
```

## Dependencies

- Python >= 3.13
- NumPy >= 2.4.1
- SciPy >= 1.17.0

## References

1. Kruskal, W. H. & Wallis, A. (1952). Use of ranks in one-criterion variance analysis. *JASA*, 47, 583-621.
2. Wallace, D. L. (1959). Simplified beta-approximations to the Kruskal-Wallis H test. *JASA*, 54, 225-230.
3. Daniels, H. E. (1954). Saddlepoint approximations in statistics. *Annals of Mathematical Statistics*, 25, 631-650.
4. Iman, R. L., Quade, D., & Alexander, D. A. (1975). Exact probability levels for the Kruskal-Wallis test.
5. Lugannani, R. & Rice, S. O. (1980). Saddlepoint approximation for the distribution of the sum of independent random variables. *Advances in Applied Probability*, 12, 475-490.
6. Easton, G. S. & Ronchetti, E. (1986). General saddlepoint approximations with applications to L statistics. *JASA*, 81, 420-430.
7. Wang, S. (1992). General saddlepoint approximations in the bootstrap. *Statistics & Probability Letters*, 13, 61-66.
8. Ha, H.-T. & Provost, S. B. (2007). A viable alternative to resorting to statistical tables. *Communications in Statistics*, 36, 1135-1151.
9. Hall, P. (1992). *The Bootstrap and Edgeworth Expansion*. Springer.

## License

MIT License
