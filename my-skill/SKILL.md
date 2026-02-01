# Kruskal-Wallis Verification Skill

## Purpose

This skill provides comprehensive verification of Kruskal-Wallis exact distribution computations against:

1. **Brute Force Ground Truth**: For small N, enumerate all permutations
2. **Paper Reference Values**: Murakami, Lee & Ha (2009) tables
3. **Independent Implementations**: SciPy and R cross-validation
4. **High-Precision Simulation**: 200k iterations with fixed seed

## Files

```
my-skill/
├── SKILL.md              # This file
├── template.md           # Report template with placeholders
├── examples/
│   └── sample.md         # Example output showing expected format
└── scripts/
    ├── validate.sh       # Main validation entry point
    ├── verify_exact.py   # Exact distribution verification
    ├── verify_simulation.py  # Simulation verification
    └── verify_r_comparison.R # R cross-validation script
```

## Usage

### Quick Test
```bash
./scripts/validate.sh --quick
```

### Full Verification
```bash
./scripts/validate.sh --all
```

### Individual Components
```bash
./scripts/validate.sh --exact   # Only exact distribution tests
./scripts/validate.sh --sim     # Only simulation tests
./scripts/validate.sh --r       # Only R cross-validation
```

## Verification Strategy

### 1. Exact Distribution (Small N)

For N < 10, we use brute force enumeration as ground truth:
- Enumerate all N! permutations
- Compute H statistic for each
- Compare probability distribution with recursive algorithm

### 2. Paper Table Comparison

Compare against paper tables:
- **Table 4.1**: (3,3,3), alpha = 0.10
- **Table 4.2**: Three groups with varying sizes
- **Table 4.6**: Four groups with varying sizes

### 3. Simulation Verification

High-precision Monte Carlo:
- 200,000 iterations (vs paper's 10,000)
- Fixed seed (42) for reproducibility
- Wilson score 95% confidence intervals
- Convergence analysis

### 4. Cross-Validation

Independent verification:
- **SciPy**: `scipy.stats.kruskal` for H statistic
- **R**: `coin` package for exact distribution

## Expected Outcomes

### If All Tests Pass
- Implementation is correct
- Paper values are reproducible
- Definition matches paper's P(H >= cv) = alpha

### If Tests Fail

Check for:
1. **Definition mismatch**: >= vs > in tail probability
2. **Tie handling**: Different tie-breaking methods
3. **Quantile rules**: Different interpolation methods
4. **Paper typos**: Rare but possible

## Output Format

See `examples/sample.md` for the expected output format. The report includes:

1. Brute force validation results
2. Paper table comparisons
3. Simulation accuracy with confidence intervals
4. Cross-validation results
5. Definition analysis
6. Summary and conclusions

## Dependencies

### Python
- numpy
- scipy
- kw_approx package (this project)

### R (Optional)
- coin package
- gtools package
- jsonlite package

## Reference

Murakami, H., Lee, W., & Ha, H. (2009). Higher Order Asymptotic Approximations of Kruskal-Wallis Statistics Based on Skewness and Kurtosis.
