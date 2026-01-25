

"""
Example: Reproducing Tables from Murakami & Ha Paper -> just example paper

This script demonstrates how to use the kw_approx package to reproduce
the numerical results from the paper:

"Higher Order Asymptotic Approximations of Kruskal-Wallis Statistics"
by Hidetoshi Murakami and Hyung-Tae Ha

Usage:
    python examples/reproduce_paper_tables.py
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from kw_approx import KWApproximator, ExactDistribution


def print_header(title: str):
    """Print a formatted header."""
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


def reproduce_table_4_1():
    """
    Reproduce Table 4.1: Kruskal-Wallis Statistic with 3 Groups
    10% significance level
    
    n1=3, n2=3, n3=3, P=4.62222
    """
    print_header("Table 4.1: Kruskal-Wallis Statistic with 3 Groups (10% significance)")
    
    sample_sizes = [3, 3, 3]
    h = 4.62222
    alpha = 0.10
    
    print(f"\nSample sizes: n1={sample_sizes[0]}, n2={sample_sizes[1]}, n3={sample_sizes[2]}")
    print(f"Test statistic H = {h}")
    print(f"Significance level α = {alpha}")
    
    approx = KWApproximator(sample_sizes)
    
    # Compute all approximations
    print(f"\n{'Method':<25} {'P(H ≥ h)':<15} {'Notes'}")
    print("-" * 60)
    
    methods_info = [
        ('exact', 'Exact (P)'),
        ('chi_square', 'Chi-square (E-Q)'),
        ('saddlepoint', 'Saddlepoint (SD1)'),
        ('saddlepoint_cc', 'Saddlepoint CC (SDC1)'),
        ('gram_charlier', 'Gram-Charlier (GC-A)'),
        ('pam', 'PAM (degree 4)'),
        ('pam6', 'PAM (degree 6)'),
    ]
    
    for method, name in methods_info:
        try:
            p = approx.tail_probability(h, method)
            print(f"{name:<25} {p:<15.6f}")
        except Exception as e:
            print(f"{name:<25} {'Error':<15} {str(e)[:30]}")
    
    # Paper reference values (from Table 4.1)
    print("\n" + "-" * 60)
    print("Reference values from paper:")
    print(f"  E-P (Exact):      0.0842833")
    print(f"  SD1:              0.0789061")
    print(f"  SDC1:             0.124486")
    print(f"  GC-A:             0.396732")
    print(f"  PAM:              0.098122")
    print(f"  PAM(6):           0.0933829")


def reproduce_table_4_4():
    """
    Reproduce Table 4.4: Kruskal-Wallis Statistic with 4 Groups
    10% significance level
    
    n1=3, n2=2, n3=2, n4=5
    """
    print_header("Table 4.4: Kruskal-Wallis Statistic with 4 Groups (10% significance)")
    
    sample_sizes = [3, 2, 2, 5]
    h = 5.587179  # From simulation
    alpha = 0.10
    
    print(f"\nSample sizes: n1={sample_sizes[0]}, n2={sample_sizes[1]}, "
          f"n3={sample_sizes[2]}, n4={sample_sizes[3]}")
    print(f"Test statistic H = {h}")
    print(f"Significance level α = {alpha}")
    
    approx = KWApproximator(sample_sizes)
    
    # Compute all approximations
    print(f"\n{'Method':<25} {'P(H ≥ h)':<15}")
    print("-" * 45)
    
    methods = ['exact', 'chi_square', 'saddlepoint', 'saddlepoint_cc', 
               'gram_charlier', 'pam', 'pam6']
    
    for method in methods:
        try:
            p = approx.tail_probability(h, method)
            print(f"{method:<25} {p:<15.6f}")
        except Exception as e:
            print(f"{method:<25} {'Error':<15}")
    
    # Paper reference values (from Table 4.4)
    print("\n" + "-" * 45)
    print("Reference values from paper:")
    print(f"  SD1:              0.105404")
    print(f"  SD2:              0.112762")
    print(f"  SDC1:             0.148341")
    print(f"  Asymptotic:       0.133516")
    print(f"  PAM:              0.109421")
    print(f"  PAM(6):           0.111407")


def compare_critical_values():
    """Compare critical values at different significance levels."""
    print_header("Critical Values Comparison")
    
    sample_sizes = [3, 3, 3]
    approx = KWApproximator(sample_sizes)
    
    print(f"\nSample sizes: {sample_sizes}, k={len(sample_sizes)}, N={sum(sample_sizes)}")
    
    alpha_values = [0.10, 0.05, 0.01]
    
    print(f"\n{'Method':<20}", end="")
    for alpha in alpha_values:
        print(f"{'α=' + str(alpha):<12}", end="")
    print()
    print("-" * 56)
    
    methods = ['chi_square', 'saddlepoint', 'pam', 'pam6', 'gram_charlier', 'exact']
    
    for method in methods:
        print(f"{method:<20}", end="")
        for alpha in alpha_values:
            try:
                cv = approx.critical_value(alpha, method)
                print(f"{cv:<12.4f}", end="")
            except:
                print(f"{'N/A':<12}", end="")
        print()


def examine_exact_distribution():
    """Examine the exact distribution for small samples."""
    print_header("Exact Distribution Analysis")
    
    sample_sizes = [2, 2, 2]
    print(f"\nSample sizes: {sample_sizes}")
    
    exact = ExactDistribution(sample_sizes)
    
    # Get distribution
    dist = exact.get_sorted_distribution()
    
    print(f"\nExact null distribution of H:")
    print(f"{'H value':<12} {'P(H = h)':<15} {'P(H ≥ h)':<15}")
    print("-" * 42)
    
    for h, prob in dist:
        sf = exact.sf(h)
        print(f"{h:<12.6f} {prob:<15.6f} {sf:<15.6f}")
    
    # Summary statistics
    summary = exact.summary()
    print(f"\nSummary Statistics:")
    print(f"  Mean:           {summary['mean']:.6f}")
    print(f"  Variance:       {summary['variance']:.6f}")
    print(f"  Skewness:       {summary['skewness']:.6f}")
    print(f"  Excess Kurtosis:{summary['excess_kurtosis']:.6f}")
    print(f"  # of values:    {summary['n_distinct_values']}")


def demonstrate_method_selection():
    """Demonstrate automatic method selection based on sample size."""
    print_header("Method Recommendation by Sample Size")
    
    test_cases = [
        [3, 3, 3],           # Small: exact
        [5, 5, 5],           # Small-medium: exact or pam6
        [10, 10, 10],        # Medium: pam6
        [20, 20, 20],        # Medium-large: pam
        [50, 50, 50],        # Large: saddlepoint
    ]
    
    print(f"\n{'Sample Sizes':<20} {'N':<8} {'Recommended':<15}")
    print("-" * 45)
    
    for sizes in test_cases:
        approx = KWApproximator(sizes)
        recommended = approx.recommend_method()
        print(f"{str(sizes):<20} {sum(sizes):<8} {recommended:<15}")


def practical_example():
    """Practical example with real data simulation."""
    print_header("Practical Example: Comparing 3 Treatments")
    
    np.random.seed(42)
    
    # Simulate data from 3 groups with different means
    group1 = np.random.normal(10, 2, 5)  # Control
    group2 = np.random.normal(12, 2, 5)  # Treatment A
    group3 = np.random.normal(11, 2, 5)  # Treatment B
    
    print("\nSimulated data:")
    print(f"  Group 1 (Control):     {group1.round(2)}")
    print(f"  Group 2 (Treatment A): {group2.round(2)}")
    print(f"  Group 3 (Treatment B): {group3.round(2)}")
    
    # Compute test statistic
    from scipy import stats
    scipy_result = stats.kruskal(group1, group2, group3)
    
    print(f"\nSciPy Kruskal-Wallis test:")
    print(f"  H statistic: {scipy_result.statistic:.4f}")
    print(f"  P-value:     {scipy_result.pvalue:.4f}")
    
    # Our approximations
    sample_sizes = [len(group1), len(group2), len(group3)]
    approx = KWApproximator(sample_sizes)
    
    H = scipy_result.statistic
    
    print(f"\nOur approximations for H = {H:.4f}:")
    results = approx.compare_methods(H)
    for method, p in sorted(results.items()):
        if p is not None:
            print(f"  {method:<20}: {p:.4f}")
    
    print(f"\nRecommended method for N={sum(sample_sizes)}: {approx.recommend_method()}")


if __name__ == '__main__':
    print("\n" + "#" * 70)
    print("# Reproducing Results from Murakami & Ha Paper")
    print("# Higher Order Asymptotic Approximations of Kruskal-Wallis Statistics")
    print("#" * 70)
    
    reproduce_table_4_1()
    reproduce_table_4_4()
    compare_critical_values()
    examine_exact_distribution()
    demonstrate_method_selection()
    practical_example()
    
    print("\n" + "=" * 70)
    print("Examples completed!")
    print("=" * 70)
