#!/usr/bin/env python3
"""
Monte Carlo Simulation Verification Script

Verifies simulation results with:
1. Fixed seed for reproducibility
2. Large iteration counts (200k) for narrow confidence intervals
3. Comparison with exact values (where available) and paper values

Reference: Murakami, Lee & Ha (2009) - Higher Order Asymptotic
Approximations of Kruskal-Wallis Statistics
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
from scipy import stats
from typing import List, Tuple, Dict
import warnings

from kw_approx.simulation import MonteCarloSimulation, simulate_tail_probability
from kw_approx.exact import ExactDistribution


def wilson_confidence_interval(p_hat: float, n: int, z: float = 1.96) -> Tuple[float, float]:
    """
    Wilson score confidence interval for a proportion.
    More accurate than normal approximation for extreme probabilities.
    """
    denominator = 1 + z**2 / n
    center = (p_hat + z**2 / (2 * n)) / denominator
    margin = z * np.sqrt(p_hat * (1 - p_hat) / n + z**2 / (4 * n**2)) / denominator
    return max(0, center - margin), min(1, center + margin)


def simulation_with_confidence(sample_sizes: List[int], h: float,
                                n_simulations: int = 200000,
                                seed: int = 42) -> Tuple[float, float, float]:
    """
    Compute tail probability with 95% confidence interval.

    Returns
    -------
    (point_estimate, ci_lower, ci_upper)
    """
    sim = MonteCarloSimulation(sample_sizes, n_simulations, seed)
    p_hat = sim.tail_probability(h)
    ci_lower, ci_upper = wilson_confidence_interval(p_hat, n_simulations)
    return p_hat, ci_lower, ci_upper


def verify_simulation_vs_exact():
    """
    Compare simulation results with exact distribution for small N.
    """
    print("=" * 90)
    print("VERIFICATION 1: Simulation vs Exact Distribution (Small N)")
    print("=" * 90)
    print()
    print("Using 200,000 simulations with fixed seed for reproducibility")
    print()

    test_cases = [
        # (sample_sizes, alpha, paper_cv)
        ([3, 3, 3], 0.10, 4.622222),
        ([4, 4, 4], 0.10, 4.500000),
        ([5, 5, 5], 0.10, 4.500000),
        ([2, 2, 2, 2], 0.10, 5.500000),
        ([3, 3, 3, 3], 0.10, 5.974359),
    ]

    n_simulations = 200000
    seed = 42

    print(f"{'Config':>12} {'H value':>10} {'Exact P':>12} {'Sim P':>12} {'95% CI':>20} {'Status':>10}")
    print("-" * 90)

    for sample_sizes, alpha, h in test_cases:
        config = ",".join(map(str, sample_sizes))

        # Exact computation
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            exact = ExactDistribution(sample_sizes)
            exact_p = exact.sf(h)

        # Simulation
        sim_p, ci_lower, ci_upper = simulation_with_confidence(sample_sizes, h, n_simulations, seed)

        # Check if exact value falls within CI
        in_ci = ci_lower <= exact_p <= ci_upper
        status = "PASS" if in_ci else "FAIL"

        ci_str = f"[{ci_lower:.6f}, {ci_upper:.6f}]"
        print(f"{config:>12} {h:>10.6f} {exact_p:>12.6f} {sim_p:>12.6f} {ci_str:>20} {status:>10}")

    print("-" * 90)
    print()


def verify_paper_simulation_values():
    """
    Verify simulation-based values from the paper (larger N cases).
    """
    print("=" * 90)
    print("VERIFICATION 2: Paper Simulation-Based Values (Larger N)")
    print("=" * 90)
    print()
    print("Paper used 10,000 simulations. We use 200,000 for better precision.")
    print()

    # Cases marked as "S" (Simulation) in the paper
    # These are larger N cases where exact is not feasible
    paper_values = [
        # (sample_sizes, alpha, paper_cv, paper_ep)
        ([6, 6, 6], 0.10, 4.526316, 0.101700),
        ([7, 7, 7], 0.10, 4.460111, 0.101300),
        ([8, 8, 8], 0.10, 4.595000, 0.100600),
        ([9, 9, 9], 0.10, 4.490300, 0.101600),
        ([10, 10, 10], 0.10, 4.560000, 0.100100),
    ]

    n_simulations = 200000
    seed = 42

    print(f"{'Config':>12} {'Paper CV':>10} {'Paper EP':>10} {'Sim P(>=CV)':>12} {'95% CI':>22} {'In CI':>8}")
    print("-" * 90)

    for sample_sizes, alpha, paper_cv, paper_ep in paper_values:
        config = ",".join(map(str, sample_sizes))

        # Our simulation
        sim_p, ci_lower, ci_upper = simulation_with_confidence(sample_sizes, paper_cv, n_simulations, seed)

        # Check if paper value falls within CI
        in_ci = ci_lower <= paper_ep <= ci_upper
        status = "YES" if in_ci else "NO"

        ci_str = f"[{ci_lower:.6f}, {ci_upper:.6f}]"
        print(f"{config:>12} {paper_cv:>10.6f} {paper_ep:>10.6f} {sim_p:>12.6f} {ci_str:>22} {status:>8}")

    print("-" * 90)
    print()
    print("Note: 'In CI' indicates whether the paper's value falls within our 95% CI.")
    print("If 'NO', the paper value may have been computed with different seed or method.")
    print()


def verify_critical_value_estimation():
    """
    Verify critical value estimation via simulation.
    """
    print("=" * 90)
    print("VERIFICATION 3: Critical Value Estimation")
    print("=" * 90)
    print()

    test_cases = [
        ([3, 3, 3], 0.10),
        ([5, 5, 5], 0.10),
        ([3, 3, 3], 0.05),
        ([5, 5, 5], 0.05),
    ]

    n_simulations = 200000
    seed = 42

    print(f"{'Config':>12} {'Alpha':>8} {'Exact CV':>12} {'Sim CV':>12} {'Exact EP':>12} {'Sim EP':>12}")
    print("-" * 80)

    for sample_sizes, alpha in test_cases:
        config = ",".join(map(str, sample_sizes))

        # Exact
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            exact = ExactDistribution(sample_sizes)
            exact_cv, exact_ep = exact.critical_value(alpha)

        # Simulation
        sim = MonteCarloSimulation(sample_sizes, n_simulations, seed)
        sim_cv, sim_ep = sim.critical_value(alpha)

        print(f"{config:>12} {alpha:>8.2f} {exact_cv:>12.6f} {sim_cv:>12.6f} {exact_ep:>12.6f} {sim_ep:>12.6f}")

    print("-" * 80)
    print()


def multiple_seed_analysis():
    """
    Analyze variation across different random seeds.
    """
    print("=" * 90)
    print("ANALYSIS: Variation Across Different Seeds")
    print("=" * 90)
    print()

    sample_sizes = [7, 7, 7]
    h = 4.460111  # Paper's CV for (7,7,7)
    n_simulations = 200000

    seeds = [42, 123, 456, 789, 1234, 5678, 9999, 11111, 22222, 33333]

    print(f"Configuration: ({','.join(map(str, sample_sizes))})")
    print(f"Testing H = {h} with {n_simulations} simulations per seed")
    print()

    results = []
    for seed in seeds:
        sim = MonteCarloSimulation(sample_sizes, n_simulations, seed)
        p = sim.tail_probability(h)
        results.append(p)
        print(f"  Seed {seed:>5}: P(H >= {h}) = {p:.6f}")

    print()
    print(f"Mean: {np.mean(results):.6f}")
    print(f"Std:  {np.std(results):.6f}")
    print(f"Range: [{min(results):.6f}, {max(results):.6f}]")
    print()


def convergence_analysis():
    """
    Analyze how simulation estimates converge as N increases.
    """
    print("=" * 90)
    print("ANALYSIS: Simulation Convergence")
    print("=" * 90)
    print()

    sample_sizes = [3, 3, 3]
    h = 4.622222

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        exact = ExactDistribution(sample_sizes)
        exact_p = exact.sf(h)

    print(f"Configuration: ({','.join(map(str, sample_sizes))})")
    print(f"H = {h}, Exact P = {exact_p:.6f}")
    print()

    n_values = [1000, 5000, 10000, 50000, 100000, 200000, 500000]
    seed = 42

    print(f"{'N Simulations':>15} {'Estimate':>12} {'Error':>12} {'95% CI Width':>15}")
    print("-" * 60)

    for n_sim in n_values:
        sim_p, ci_lower, ci_upper = simulation_with_confidence(sample_sizes, h, n_sim, seed)
        error = abs(sim_p - exact_p)
        ci_width = ci_upper - ci_lower

        print(f"{n_sim:>15} {sim_p:>12.6f} {error:>12.6f} {ci_width:>15.6f}")

    print("-" * 60)
    print()
    print("Note: Error should decrease approximately as 1/sqrt(N)")
    print()


def main():
    """Run all simulation verification tests."""
    print()
    print("#" * 90)
    print("# KRUSKAL-WALLIS SIMULATION VERIFICATION")
    print("# Reference: Murakami, Lee & Ha - Higher Order Asymptotic Approximations")
    print("#" * 90)
    print()

    # Test 1: Simulation vs Exact
    verify_simulation_vs_exact()

    # Test 2: Paper simulation values
    verify_paper_simulation_values()

    # Test 3: Critical value estimation
    verify_critical_value_estimation()

    # Analysis: Multiple seeds
    multiple_seed_analysis()

    # Analysis: Convergence
    convergence_analysis()

    print()
    print("=" * 90)
    print("VERIFICATION COMPLETE")
    print("=" * 90)
    print()

    return 0


if __name__ == "__main__":
    sys.exit(main())
