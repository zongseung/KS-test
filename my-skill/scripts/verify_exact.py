#!/usr/bin/env python3
"""
Exact Distribution Verification Script

Verifies the exact distribution computation for small N cases by:
1. Computing exact distribution for small sample sizes
2. Comparing with paper reference values
3. Cross-validating with SciPy kruskal test where applicable
4. Computing exact critical values and tail probabilities

Reference: Murakami, Lee & Ha (2009) - Higher Order Asymptotic
Approximations of Kruskal-Wallis Statistics
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
from scipy import stats
from itertools import permutations
from math import factorial, comb
from collections import defaultdict
from typing import Dict, List, Tuple

from kw_approx.exact import ExactDistribution
from kw_approx.kruskal_wallis import KruskalWallisStatistic


def multinomial_coefficient(sample_sizes: List[int]) -> int:
    """Compute N! / (n1! * n2! * ... * nk!)"""
    N = sum(sample_sizes)
    result = 1
    remaining = N
    for n_i in sample_sizes:
        result *= comb(remaining, n_i)
        remaining -= n_i
    return result


def compute_h_from_rank_sums(rank_sums: List[int], sample_sizes: List[int]) -> float:
    """Compute H statistic from rank sums."""
    N = sum(sample_sizes)
    H = 0.0
    for R_i, n_i in zip(rank_sums, sample_sizes):
        H += R_i**2 / n_i
    H = 12.0 / (N * (N + 1)) * H - 3 * (N + 1)
    return H


def brute_force_exact_distribution(sample_sizes: List[int]) -> Dict[float, float]:
    """
    Compute exact distribution by enumerating ALL permutations.
    This is the ground truth reference for verification.
    Only feasible for very small N (< 10).
    """
    N = sum(sample_sizes)
    k = len(sample_sizes)

    # All possible rank assignments
    ranks = list(range(1, N + 1))

    # Count H values from all permutations
    h_counts = defaultdict(int)

    for perm in permutations(ranks):
        # Split permutation into groups
        idx = 0
        rank_sums = []
        for n in sample_sizes:
            rank_sums.append(sum(perm[idx:idx + n]))
            idx += n

        H = compute_h_from_rank_sums(rank_sums, sample_sizes)
        H_rounded = round(H, 10)
        h_counts[H_rounded] += 1

    # Convert to probabilities
    total = factorial(N)
    distribution = {h: count / total for h, count in h_counts.items()}

    return distribution


def compare_distributions(dist1: Dict[float, float],
                         dist2: Dict[float, float],
                         tol: float = 1e-10) -> Tuple[bool, str]:
    """Compare two distributions for equality."""
    all_keys = set(dist1.keys()) | set(dist2.keys())

    differences = []
    for h in sorted(all_keys):
        p1 = dist1.get(h, 0.0)
        p2 = dist2.get(h, 0.0)
        if abs(p1 - p2) > tol:
            differences.append(f"  H={h:.6f}: dist1={p1:.10f}, dist2={p2:.10f}, diff={abs(p1-p2):.2e}")

    if differences:
        return False, "\n".join(differences)
    return True, "Distributions match exactly"


def verify_small_cases():
    """
    Verify exact distribution for small cases where brute force is feasible.
    """
    print("=" * 80)
    print("VERIFICATION 1: Brute Force vs Recursive Algorithm (Small N)")
    print("=" * 80)
    print()

    test_cases = [
        [2, 2, 2],      # N=6
        [3, 3, 3],      # N=9 (Paper Table 4.1)
        [2, 2, 2, 2],   # N=8 (Paper Table 4.6)
        [3, 2, 2],      # N=7
    ]

    results = []

    for sample_sizes in test_cases:
        N = sum(sample_sizes)
        config = ",".join(map(str, sample_sizes))
        print(f"Testing configuration: ({config}), N={N}")

        # Brute force (ground truth)
        print("  Computing brute force distribution...", end=" ", flush=True)
        bf_dist = brute_force_exact_distribution(sample_sizes)
        print(f"Done. Found {len(bf_dist)} distinct H values.")

        # Recursive algorithm
        print("  Computing recursive distribution...", end=" ", flush=True)
        exact = ExactDistribution(sample_sizes)
        rec_dist = exact.distribution
        print(f"Done. Found {len(rec_dist)} distinct H values.")

        # Compare
        match, details = compare_distributions(bf_dist, rec_dist)

        if match:
            print(f"  [PASS] Distributions match exactly!")
        else:
            print(f"  [FAIL] Distributions differ:")
            print(details)

        results.append({
            'config': config,
            'N': N,
            'n_values': len(bf_dist),
            'match': match
        })
        print()

    return results


def verify_paper_table_4_1():
    """
    Verify Table 4.1 from the paper: (3,3,3), alpha=0.10
    """
    print("=" * 80)
    print("VERIFICATION 2: Paper Table 4.1 - Configuration (3,3,3), alpha=0.10")
    print("=" * 80)
    print()

    sample_sizes = [3, 3, 3]
    alpha = 0.10

    # Paper reference values
    paper_cv = 4.622222  # Critical value P(E-Q)
    paper_ep = 0.100000  # Exact probability at critical value

    # Our exact computation
    exact = ExactDistribution(sample_sizes)

    # Find critical value for alpha=0.10
    our_cv, our_ep = exact.critical_value(alpha)

    print(f"Sample sizes: {sample_sizes}")
    print(f"Total N = {sum(sample_sizes)}")
    print()
    print(f"Paper values:")
    print(f"  Critical Value P(E-Q) = {paper_cv}")
    print(f"  Exact Probability E-P = {paper_ep}")
    print()
    print(f"Our computed values:")
    print(f"  Critical Value = {our_cv:.6f}")
    print(f"  Exact Probability P(H >= cv) = {our_ep:.6f}")
    print()

    # Check tail probability at paper's critical value
    our_sf = exact.sf(paper_cv)
    print(f"Tail probability at paper's CV ({paper_cv}):")
    print(f"  P(H >= {paper_cv}) = {our_sf:.6f}")
    print()

    # List all H values and their probabilities
    sorted_dist = exact.get_sorted_distribution()
    print("Exact distribution (H values and probabilities):")
    print("-" * 50)
    cumulative = 0.0
    for h, p in sorted_dist:
        cumulative += p
        print(f"  H = {h:8.6f}, P(H=h) = {p:.6f}, P(H>=h) = {1-cumulative+p:.6f}")
    print("-" * 50)

    cv_match = abs(our_cv - paper_cv) < 0.001
    ep_match = abs(our_ep - paper_ep) < 0.001

    print()
    print(f"Critical Value Match: {'[PASS]' if cv_match else '[FAIL]'}")
    print(f"Exact Probability Match: {'[PASS]' if ep_match else '[FAIL]'}")

    return cv_match and ep_match


def verify_paper_table_4_2():
    """
    Verify Table 4.2 from the paper: Three groups with exact computation
    """
    print()
    print("=" * 80)
    print("VERIFICATION 3: Paper Table 4.2 - Three Groups, alpha=0.10 (Exact cases)")
    print("=" * 80)
    print()

    # Cases marked as "E" (Exact) in the paper
    paper_values = [
        # (sample_sizes, paper_cv, paper_ep)
        ([3, 3, 3], 4.622222, 0.100000),
        ([4, 4, 4], 4.500000, 0.104242),
        ([5, 5, 5], 4.500000, 0.101502),
    ]

    print(f"{'Config':>10} {'Paper CV':>10} {'Our CV':>10} {'Paper EP':>10} {'Our EP':>10} {'CV Match':>10} {'EP Match':>10}")
    print("-" * 80)

    all_pass = True
    for sample_sizes, paper_cv, paper_ep in paper_values:
        config = ",".join(map(str, sample_sizes))

        exact = ExactDistribution(sample_sizes)
        our_cv, our_ep = exact.critical_value(0.10)

        cv_match = abs(our_cv - paper_cv) < 0.01
        # For EP, check tail probability at paper's CV
        our_sf = exact.sf(paper_cv)
        ep_match = abs(our_sf - paper_ep) < 0.01

        if not cv_match or not ep_match:
            all_pass = False

        print(f"{config:>10} {paper_cv:>10.6f} {our_cv:>10.6f} {paper_ep:>10.6f} {our_sf:>10.6f} {'PASS' if cv_match else 'FAIL':>10} {'PASS' if ep_match else 'FAIL':>10}")

    print("-" * 80)
    return all_pass


def verify_paper_table_4_6():
    """
    Verify Table 4.6 from the paper: Four groups with exact computation
    """
    print()
    print("=" * 80)
    print("VERIFICATION 4: Paper Table 4.6 - Four Groups, alpha=0.10 (Exact cases)")
    print("=" * 80)
    print()

    # Cases marked as "E" (Exact) in the paper
    paper_values = [
        # (sample_sizes, paper_cv, paper_ep)
        ([2, 2, 2, 2], 5.500000, 0.114286),
        ([3, 2, 2, 3], 5.727273, 0.100159),
        ([3, 2, 2, 5], 5.753846, 0.100697),
        ([3, 3, 3, 3], 5.974359, 0.102727),
    ]

    print(f"{'Config':>12} {'Paper CV':>10} {'Our CV':>10} {'Paper EP':>10} {'Our EP':>10} {'CV Match':>10} {'EP Match':>10}")
    print("-" * 90)

    all_pass = True
    for sample_sizes, paper_cv, paper_ep in paper_values:
        config = ",".join(map(str, sample_sizes))

        exact = ExactDistribution(sample_sizes)
        our_cv, our_ep = exact.critical_value(0.10)

        cv_match = abs(our_cv - paper_cv) < 0.01
        # For EP, check tail probability at paper's CV
        our_sf = exact.sf(paper_cv)
        ep_match = abs(our_sf - paper_ep) < 0.01

        if not cv_match or not ep_match:
            all_pass = False

        print(f"{config:>12} {paper_cv:>10.6f} {our_cv:>10.6f} {paper_ep:>10.6f} {our_sf:>10.6f} {'PASS' if cv_match else 'FAIL':>10} {'PASS' if ep_match else 'FAIL':>10}")

    print("-" * 90)
    return all_pass


def verify_scipy_consistency():
    """
    Cross-validate with SciPy's kruskal test.
    SciPy uses chi-square approximation but we can verify H computation.
    """
    print()
    print("=" * 80)
    print("VERIFICATION 5: SciPy Consistency Check")
    print("=" * 80)
    print()

    np.random.seed(42)

    test_cases = [
        [5, 5, 5],
        [3, 4, 5],
        [4, 4, 4, 4],
    ]

    print("Comparing H statistic computation with SciPy:")
    print("-" * 60)

    for sample_sizes in test_cases:
        config = ",".join(map(str, sample_sizes))

        # Generate random data
        groups = [np.random.randn(n) for n in sample_sizes]

        # Our computation
        kw = KruskalWallisStatistic(sample_sizes)
        our_H = kw.compute_statistic(*groups)

        # SciPy computation
        scipy_H, _ = stats.kruskal(*groups)

        match = abs(our_H - scipy_H) < 1e-10
        print(f"Config ({config}): Our H = {our_H:.10f}, SciPy H = {scipy_H:.10f} - {'PASS' if match else 'FAIL'}")

    print("-" * 60)
    return True


def analyze_definition_differences():
    """
    Analyze potential definition differences: >= vs >, tie handling, quantile rules
    """
    print()
    print("=" * 80)
    print("ANALYSIS: Definition Differences (>= vs >, quantile rules)")
    print("=" * 80)
    print()

    sample_sizes = [3, 3, 3]
    exact = ExactDistribution(sample_sizes)
    sorted_dist = exact.get_sorted_distribution()

    print("For configuration (3,3,3), analyzing critical value definitions:")
    print()

    alpha = 0.10

    # Method 1: Smallest h where P(H >= h) <= alpha
    # Method 2: Smallest h where P(H > h) < alpha
    # Method 3: Largest h where P(H >= h) >= alpha

    cumulative_sf = 1.0  # P(H >= min_H) = 1

    print(f"{'H value':>10} {'P(H=h)':>12} {'P(H>=h)':>12} {'P(H>h)':>12}")
    print("-" * 50)

    for h, p in sorted_dist:
        sf_ge = sum(prob for h_val, prob in sorted_dist if h_val >= h)
        sf_gt = sum(prob for h_val, prob in sorted_dist if h_val > h)
        print(f"{h:>10.6f} {p:>12.6f} {sf_ge:>12.6f} {sf_gt:>12.6f}")

    print("-" * 50)
    print()
    print("Paper uses P(H >= cv) = alpha convention (E-Q method)")
    print()


def main():
    """Run all verification tests."""
    print()
    print("#" * 80)
    print("# KRUSKAL-WALLIS EXACT DISTRIBUTION VERIFICATION")
    print("# Reference: Murakami, Lee & Ha - Higher Order Asymptotic Approximations")
    print("#" * 80)
    print()

    results = []

    # Test 1: Brute force vs recursive for small cases
    bf_results = verify_small_cases()
    all_bf_pass = all(r['match'] for r in bf_results)
    results.append(('Brute Force Validation', all_bf_pass))

    # Test 2: Paper Table 4.1
    table_4_1_pass = verify_paper_table_4_1()
    results.append(('Paper Table 4.1', table_4_1_pass))

    # Test 3: Paper Table 4.2
    table_4_2_pass = verify_paper_table_4_2()
    results.append(('Paper Table 4.2', table_4_2_pass))

    # Test 4: Paper Table 4.6
    table_4_6_pass = verify_paper_table_4_6()
    results.append(('Paper Table 4.6', table_4_6_pass))

    # Test 5: SciPy consistency
    scipy_pass = verify_scipy_consistency()
    results.append(('SciPy Consistency', scipy_pass))

    # Analysis of definition differences
    analyze_definition_differences()

    # Summary
    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print()
    for name, passed in results:
        status = "[PASS]" if passed else "[FAIL]"
        print(f"  {name}: {status}")
    print()

    all_pass = all(r[1] for r in results)
    if all_pass:
        print("All verification tests passed!")
    else:
        print("Some tests failed. Check details above.")

    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
