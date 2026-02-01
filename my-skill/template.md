# Kruskal-Wallis Exact Distribution Verification Report

## Overview

This report verifies the implementation of exact distribution computation for the Kruskal-Wallis H statistic against the reference paper:

**Murakami, Lee & Ha (2009)**: "Higher Order Asymptotic Approximations of Kruskal-Wallis Statistics Based on Skewness and Kurtosis"

---

## 1. Exact Distribution Verification

### 1.1 Brute Force vs Recursive Algorithm

For small sample sizes (N < 10), we verify our recursive algorithm against brute force enumeration of all permutations.

| Configuration | N | Distinct H Values | Brute Force Match |
|--------------|---|-------------------|-------------------|
| {{CONFIG_1}} | {{N_1}} | {{N_VALUES_1}} | {{MATCH_1}} |
| {{CONFIG_2}} | {{N_2}} | {{N_VALUES_2}} | {{MATCH_2}} |
| {{CONFIG_3}} | {{N_3}} | {{N_VALUES_3}} | {{MATCH_3}} |

### 1.2 Paper Table 4.1 Verification

**Configuration**: (3, 3, 3), N = 9, alpha = 0.10

| Metric | Paper Value | Our Value | Status |
|--------|-------------|-----------|--------|
| Critical Value (E-Q) | 4.622222 | {{OUR_CV_4_1}} | {{STATUS_CV_4_1}} |
| Exact Probability (E-P) | 0.100000 | {{OUR_EP_4_1}} | {{STATUS_EP_4_1}} |

### 1.3 Paper Table 4.2 Verification (Three Groups, alpha = 0.10)

| Configuration | Paper CV | Our CV | Paper EP | Our EP | CV Match | EP Match |
|--------------|----------|--------|----------|--------|----------|----------|
| 3,3,3 | 4.622222 | {{CV_333}} | 0.100000 | {{EP_333}} | {{M_CV_333}} | {{M_EP_333}} |
| 4,4,4 | 4.500000 | {{CV_444}} | 0.104242 | {{EP_444}} | {{M_CV_444}} | {{M_EP_444}} |
| 5,5,5 | 4.500000 | {{CV_555}} | 0.101502 | {{EP_555}} | {{M_CV_555}} | {{M_EP_555}} |

### 1.4 Paper Table 4.6 Verification (Four Groups, alpha = 0.10)

| Configuration | Paper CV | Our CV | Paper EP | Our EP | CV Match | EP Match |
|--------------|----------|--------|----------|--------|----------|----------|
| 2,2,2,2 | 5.500000 | {{CV_2222}} | 0.114286 | {{EP_2222}} | {{M_CV_2222}} | {{M_EP_2222}} |
| 3,2,2,3 | 5.727273 | {{CV_3223}} | 0.100159 | {{EP_3223}} | {{M_CV_3223}} | {{M_EP_3223}} |
| 3,3,3,3 | 5.974359 | {{CV_3333}} | 0.102727 | {{EP_3333}} | {{M_CV_3333}} | {{M_EP_3333}} |

---

## 2. Simulation Verification

### 2.1 Simulation vs Exact (Small N)

Using 200,000 simulations with fixed seed (42) for reproducibility.

| Configuration | H Value | Exact P | Sim P | 95% CI | In CI |
|--------------|---------|---------|-------|--------|-------|
| {{SIM_CONFIG_1}} | {{SIM_H_1}} | {{SIM_EXACT_1}} | {{SIM_P_1}} | {{SIM_CI_1}} | {{SIM_INCI_1}} |
| {{SIM_CONFIG_2}} | {{SIM_H_2}} | {{SIM_EXACT_2}} | {{SIM_P_2}} | {{SIM_CI_2}} | {{SIM_INCI_2}} |
| {{SIM_CONFIG_3}} | {{SIM_H_3}} | {{SIM_EXACT_3}} | {{SIM_P_3}} | {{SIM_CI_3}} | {{SIM_INCI_3}} |

### 2.2 Convergence Analysis

| N Simulations | Estimate | Error | 95% CI Width |
|---------------|----------|-------|--------------|
| 1,000 | {{EST_1K}} | {{ERR_1K}} | {{CIW_1K}} |
| 10,000 | {{EST_10K}} | {{ERR_10K}} | {{CIW_10K}} |
| 100,000 | {{EST_100K}} | {{ERR_100K}} | {{CIW_100K}} |
| 200,000 | {{EST_200K}} | {{ERR_200K}} | {{CIW_200K}} |

---

## 3. Cross-Validation

### 3.1 SciPy Consistency

H statistic computation matches SciPy's `scipy.stats.kruskal`:

| Configuration | Our H | SciPy H | Match |
|--------------|-------|---------|-------|
| {{SCIPY_CONFIG_1}} | {{SCIPY_OUR_1}} | {{SCIPY_SCIPY_1}} | {{SCIPY_MATCH_1}} |
| {{SCIPY_CONFIG_2}} | {{SCIPY_OUR_2}} | {{SCIPY_SCIPY_2}} | {{SCIPY_MATCH_2}} |

### 3.2 R Comparison (Optional)

Using R's `coin` package for exact Kruskal-Wallis test:

| Configuration | Python CV | R CV | Python EP | R EP |
|--------------|-----------|------|-----------|------|
| {{R_CONFIG_1}} | {{R_PY_CV_1}} | {{R_R_CV_1}} | {{R_PY_EP_1}} | {{R_R_EP_1}} |

---

## 4. Definition Analysis

### 4.1 Critical Value Definition

The paper uses **P(H >= cv) = alpha** convention (E-Q method).

For configuration (3,3,3):

| H value | P(H=h) | P(H>=h) | P(H>h) |
|---------|--------|---------|--------|
| {{DEF_H_1}} | {{DEF_PMF_1}} | {{DEF_SF_GE_1}} | {{DEF_SF_GT_1}} |
| {{DEF_H_2}} | {{DEF_PMF_2}} | {{DEF_SF_GE_2}} | {{DEF_SF_GT_2}} |
| {{DEF_H_3}} | {{DEF_PMF_3}} | {{DEF_SF_GE_3}} | {{DEF_SF_GT_3}} |

### 4.2 Potential Discrepancy Sources

1. **>= vs > definition**: Paper uses P(H >= cv), some implementations use P(H > cv)
2. **Tie handling**: Different methods for handling ties in ranks
3. **Quantile rule**: Type of quantile interpolation used
4. **Floating point precision**: Rounding of H values

---

## 5. Summary

| Test Category | Status | Notes |
|---------------|--------|-------|
| Brute Force Validation | {{STATUS_BF}} | Ground truth for small N |
| Paper Table 4.1 | {{STATUS_T41}} | (3,3,3) configuration |
| Paper Table 4.2 | {{STATUS_T42}} | Three groups, exact cases |
| Paper Table 4.6 | {{STATUS_T46}} | Four groups, exact cases |
| Simulation Accuracy | {{STATUS_SIM}} | 200k iterations |
| SciPy Consistency | {{STATUS_SCIPY}} | H statistic computation |
| R Cross-Validation | {{STATUS_R}} | Optional verification |

---

## 6. Conclusions

{{CONCLUSIONS}}

---

*Generated by: Kruskal-Wallis Validation Suite*
*Date: {{DATE}}*
