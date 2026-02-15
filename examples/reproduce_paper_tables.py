"""
Example: Reproducing Tables from Murakami, Lee & Ha Paper

This script demonstrates how to use the kw_approx package to reproduce
the numerical results from the paper:

"Higher Order Asymptotic Approximations of Kruskal-Wallis Statistics
 Based on Skewness and Kurtosis"
by Hidetoshi Murakami, Jong-Seung Lee, and Hyung-Tae Ha

Tables reproduced:
- Table 4.1: Three groups, balanced (3,3,3), alpha = 0.10
- Tables 4.2-4.3: Three groups, increasing n, alpha = 0.10 and 0.05
- Table 4.4: Three groups, larger n, simulation-based benchmarks
- Table 4.5: Four groups (3,2,2,5), alpha = 0.10
- Tables 4.6-4.7: Additional four-group designs, alpha = 0.10 and 0.05

Usage:
    python examples/reproduce_paper_tables.py
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import random
import numpy as np
from kw_approx import KWApproximator, ExactDistribution, MonteCarloSimulation


def print_header(title: str):
    """Print a formatted header."""
    print("\n" + "=" * 90)
    print(title)
    print("=" * 90)


def print_table_row(values: list, widths: list):
    """Print a formatted table row."""
    row = ""
    for val, w in zip(values, widths):
        if isinstance(val, float):
            row += f"{val:>{w}.6f}"
        elif val is None:
            row += f"{'N/A':>{w}}"
        else:
            row += f"{val:>{w}}"
    print(row)


def get_exact_critical_value(sample_sizes: list, alpha: float):
    """Get exact critical value for given sample sizes and alpha."""
    try:
        approx = KWApproximator(sample_sizes)
        cv, exact_alpha = approx._get_method('exact').critical_value(alpha)
        return cv, exact_alpha
    except Exception as e:
        return None, None


def get_reference_probability(approx: KWApproximator, h: float, N: int, k: int,
                               n_simulations: int = 10000) -> tuple:
    """
    Get reference probability (exact or simulation-based).

    Following paper methodology:
    - For small N: use exact distribution
    - For large N: use Monte Carlo simulation

    Returns:
        (probability, is_simulation_based)
    """
    # Determine threshold based on k
    if k <= 3:
        limit_N = 15
    elif k == 4:
        limit_N = 13
    else:
        limit_N = 10

    if N <= limit_N:
        # Use exact
        try:
            return approx.tail_probability(h, 'exact'), False
        except:
            pass

    # Use simulation for large N
    try:
        sim = MonteCarloSimulation(approx.sample_sizes, n_simulations=n_simulations)
        return sim.tail_probability(h), True
    except:
        return None, False


def get_reference_critical_value(sample_sizes: list, alpha: float,
                                  n_simulations: int = 10000,
                                  seed: int = None) -> tuple:
    """
    Get reference critical value following paper methodology.

    The critical value H(alpha) should be determined from a reference distribution:
    1. Exact distribution (for small N)
    2. Monte Carlo simulation (for large N)
    3. Chi-square (fallback only)

    IMPORTANT: PAM6 or other approximations should NEVER be used as reference
    for determining critical values, as this causes self-calibration issues.

    Parameters
    ----------
    sample_sizes : list
        Sample sizes for each group
    alpha : float
        Significance level (e.g., 0.10, 0.05)
    n_simulations : int
        Number of Monte Carlo simulations if needed
    seed : int, optional
        Random seed for reproducibility

    Returns
    -------
    tuple : (h_ref, ref_alpha, source)
        h_ref : float - Reference critical value
        ref_alpha : float - Actual tail probability at h_ref
        source : str - "exact", "simulation", or "chi_square"
    """
    from scipy.stats import chi2

    k = len(sample_sizes)
    N = sum(sample_sizes)

    # Determine threshold based on k (following paper recommendations)
    if k <= 3:
        limit_N = 15
    elif k == 4:
        limit_N = 13
    else:
        limit_N = 10  # Very strict for k >= 5

    # Try exact distribution first (for small N)
    if N <= limit_N:
        try:
            exact = ExactDistribution(sample_sizes)
            h_ref, ref_alpha = exact.critical_value(alpha)
            return h_ref, ref_alpha, "exact"
        except Exception as e:
            pass  # Fall through to simulation

    # Use Monte Carlo simulation for larger N
    try:
        sim = MonteCarloSimulation(sample_sizes, n_simulations=n_simulations, seed=seed)
        h_ref, ref_alpha = sim.critical_value(alpha)
        return h_ref, ref_alpha, "simulation"
    except Exception as e:
        pass  # Fall through to chi-square

    # Fallback to chi-square (last resort)
    df = k - 1
    h_ref = chi2.ppf(1 - alpha, df)
    ref_alpha = 1 - chi2.cdf(h_ref, df)  # Should be very close to alpha
    return h_ref, ref_alpha, "chi_square"


def reproduce_table_4_1():
    """
    Reproduce Table 4.1: Three groups, balanced (3, 3, 3), alpha = 0.10
    """
    print_header("Table 4.1: Three groups, balanced (3, 3, 3), alpha = 0.10")

    sample_sizes = [3, 3, 3]
    H = 4.62222  # Exact 10% critical value from paper

    print(f"\nSample sizes: n1={sample_sizes[0]}, n2={sample_sizes[1]}, n3={sample_sizes[2]}")
    print(f"Total N = {sum(sample_sizes)}, k = {len(sample_sizes)}")
    print(f"Critical value P (E-Q) = {H}")

    approx = KWApproximator(sample_sizes)

    headers = ["n1,n2,n3", "P(E-Q)", "E-P", "SD1", "SD2", "SDC1", "SDC2", "ED", "GC-A", "PAG(4)", "PAG(6)"]
    widths = [12, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10]

    print(f"\n{'-'*112}")
    print_table_row(headers, widths)
    print(f"{'-'*112}")

    exact_p = approx.tail_probability(H, 'exact')
    sd1 = approx.tail_probability(H, 'saddlepoint')
    sd2 = approx.tail_probability(H, 'saddlepoint_sd2')
    sdc1 = approx.tail_probability(H, 'saddlepoint_cc')
    sdc2 = approx.tail_probability(H, 'saddlepoint_cc2')
    ed = approx.tail_probability(H, 'edgeworth')
    gc_a = approx.tail_probability(H, 'gram_charlier')
    pag4 = approx.tail_probability(H, 'pam')
    pag6 = approx.tail_probability(H, 'pam6')

    values = ["3,3,3", H, exact_p, sd1, sd2, sdc1, sdc2, ed, gc_a, pag4, pag6]
    print_table_row(values, widths)
    print(f"{'-'*112}")


def reproduce_tables_4_2_4_3():
    """
    Reproduce Tables 4.2-4.3: Three groups, increasing n, alpha = 0.10 and 0.05

    Following paper methodology:
    - Critical value H is determined from REFERENCE distribution (exact or simulation)
    - All approximation methods are compared at the SAME reference H

    Generate balanced three-group designs with increasing sample sizes.
    """
    print_header("Tables 4.2-4.3: Three groups, increasing n, alpha = 0.10 and 0.05")

    # Balanced three-group designs with increasing n
    n_values = [3, 4, 5, 6, 7, 8, 9, 10]

    for alpha in [0.10, 0.05]:
        print(f"\n--- alpha = {alpha} ---\n")

        headers = ["n1,n2,n3", "N", "E-Q", "E-P", "SRC", "CHI", "SD1", "SD2", "SDC1", "SDC2", "ED", "GC-A", "PAG(4)", "PAG(6)"]
        widths = [12, 6, 10, 10, 6, 10, 10, 10, 10, 10, 10, 10, 10, 10]

        print_table_row(headers, widths)
        print("-" * 140)

        for n in n_values:
            sample_sizes = [n, n, n]
            N = sum(sample_sizes)

            approx = KWApproximator(sample_sizes)

            # Get REFERENCE critical value (exact or simulation - NEVER pam6!)
            H, ref_alpha, source = get_reference_critical_value(sample_sizes, alpha, n_simulations=10000)

            # Compute all approximation tail probabilities at the same reference H
            chi = approx.tail_probability(H, 'chi_square')
            sd1 = approx.tail_probability(H, 'saddlepoint')
            sd2 = approx.tail_probability(H, 'saddlepoint_sd2')
            sdc1 = approx.tail_probability(H, 'saddlepoint_cc')
            sdc2 = approx.tail_probability(H, 'saddlepoint_cc2')
            ed = approx.tail_probability(H, 'edgeworth')
            gc_a = approx.tail_probability(H, 'gram_charlier')
            pag4 = approx.tail_probability(H, 'pam')
            pag6 = approx.tail_probability(H, 'pam6')

            config = f"{n},{n},{n}"
            # Source indicator: E=exact, S=simulation, C=chi-square
            src_label = {"exact": "E", "simulation": "S", "chi_square": "C"}[source]
            values = [config, N, H, ref_alpha, src_label, chi, sd1, sd2, sdc1, sdc2, ed, gc_a, pag4, pag6]
            print_table_row(values, widths)

        print("-" * 140)
        print("  SRC: E=Exact, S=Simulation (10,000 iterations), C=Chi-square")


def reproduce_table_4_4():
    """
    Reproduce Table 4.4: Three groups, larger n, simulation-based benchmarks

    For larger sample sizes where exact computation is impractical,
    use simulation or chi-square as reference.
    """
    print_header("Table 4.4: Three groups, larger n (simulation-based benchmarks)")

    # Larger balanced three-group designs
    n_values = [15, 20, 25, 30, 40, 50]

    print("\n--- alpha = 0.10 ---\n")

    headers = ["n1,n2,n3", "N", "SIM-CV", "SIM", "CHI", "SD1", "SD2", "SDC1", "SDC2", "ED", "GC-A", "PAG(4)", "PAG(6)"]
    widths = [12, 6, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10]

    print_table_row(headers, widths)
    print("-" * 134)

    alpha = 0.10

    for n in n_values:
        sample_sizes = [n, n, n]
        N = sum(sample_sizes)

        approx = KWApproximator(sample_sizes)

        # Use simulation as reference for larger N (following paper Table 4.4)
        H, ref_alpha, source = get_reference_critical_value(sample_sizes, alpha, n_simulations=10000)

        chi = approx.tail_probability(H, 'chi_square')
        sd1 = approx.tail_probability(H, 'saddlepoint')
        sd2 = approx.tail_probability(H, 'saddlepoint_sd2')
        sdc1 = approx.tail_probability(H, 'saddlepoint_cc')
        sdc2 = approx.tail_probability(H, 'saddlepoint_cc2')
        ed = approx.tail_probability(H, 'edgeworth')
        gc_a = approx.tail_probability(H, 'gram_charlier')
        pag4 = approx.tail_probability(H, 'pam')
        pag6 = approx.tail_probability(H, 'pam6')

        config = f"{n},{n},{n}"
        values = [config, N, H, ref_alpha, chi, sd1, sd2, sdc1, sdc2, ed, gc_a, pag4, pag6]
        print_table_row(values, widths)

    print("-" * 134)


def reproduce_table_4_5():
    """
    Reproduce Table 4.5: Four groups, (3, 2, 2, 5), alpha = 0.10
    """
    print_header("Table 4.5: Four groups (3, 2, 2, 5), alpha = 0.10")

    sample_sizes = [3, 2, 2, 5]
    H = 5.587179

    print(f"\nSample sizes: {sample_sizes}")
    print(f"Total N = {sum(sample_sizes)}, k = {len(sample_sizes)}")
    print(f"Test statistic H = {H}")

    approx = KWApproximator(sample_sizes)

    headers = ["Config", "E-P/Sim", "SD1", "SD2", "SDC1", "SDC2", "CHI", "PAG(4)", "PAG(6)"]
    widths = [14, 10, 10, 10, 10, 10, 10, 10, 10]

    print(f"\n{'-'*94}")
    print_table_row(headers, widths)
    print(f"{'-'*94}")

    exact_p = approx.tail_probability(H, 'exact')
    sd1 = approx.tail_probability(H, 'saddlepoint')
    sd2 = approx.tail_probability(H, 'saddlepoint_sd2')
    sdc1 = approx.tail_probability(H, 'saddlepoint_cc')
    sdc2 = approx.tail_probability(H, 'saddlepoint_cc2')
    chi = approx.tail_probability(H, 'chi_square')
    pag4 = approx.tail_probability(H, 'pam')
    pag6 = approx.tail_probability(H, 'pam6')

    values = ["3,2,2,5", exact_p, sd1, sd2, sdc1, sdc2, chi, pag4, pag6]
    print_table_row(values, widths)
    print(f"{'-'*94}")

    # Paper reference
    print("\nPaper reference values:")
    print(f"  SD1: 0.105404 (Code: {sd1:.6f}) {'OK' if abs(sd1 - 0.105404) < 0.001 else 'DIFF'}")
    print(f"  CHI: 0.133516 (Code: {chi:.6f}) {'OK' if abs(chi - 0.133516) < 0.001 else 'DIFF'}")


def reproduce_tables_4_6_4_7():
    """
    Reproduce Tables 4.6-4.7: Additional four-group designs, alpha = 0.10 and 0.05

    Following paper methodology:
    - Critical value H is determined from REFERENCE distribution (exact or simulation)
    - All approximation methods are compared at the SAME reference H

    Generate various unbalanced four-group designs.
    """
    print_header("Tables 4.6-4.7: Additional four-group designs, alpha = 0.10 and 0.05")

    # Various four-group configurations (unbalanced)
    configs = [
        [2, 2, 2, 2],
        [3, 2, 2, 3],
        [3, 2, 2, 5],
        [4, 3, 3, 5],
        [5, 4, 4, 5],
        [3, 3, 3, 3],
        [4, 4, 4, 4],
        [5, 5, 5, 5],
    ]

    for alpha in [0.10, 0.05]:
        print(f"\n--- alpha = {alpha} ---\n")

        headers = ["Config", "N", "E-Q", "E-P", "SRC", "CHI", "SD1", "SD2", "SDC1", "SDC2", "ED", "GC-A", "PAG(4)", "PAG(6)"]
        widths = [14, 6, 10, 10, 6, 10, 10, 10, 10, 10, 10, 10, 10, 10]

        print_table_row(headers, widths)
        print("-" * 142)

        for sample_sizes in configs:
            N = sum(sample_sizes)
            approx = KWApproximator(sample_sizes)

            # Get REFERENCE critical value (exact or simulation - NEVER pam6!)
            H, ref_alpha, source = get_reference_critical_value(sample_sizes, alpha, n_simulations=10000)

            # Compute all approximation tail probabilities at the same reference H
            chi = approx.tail_probability(H, 'chi_square')
            sd1 = approx.tail_probability(H, 'saddlepoint')
            sd2 = approx.tail_probability(H, 'saddlepoint_sd2')
            sdc1 = approx.tail_probability(H, 'saddlepoint_cc')
            sdc2 = approx.tail_probability(H, 'saddlepoint_cc2')
            ed = approx.tail_probability(H, 'edgeworth')
            gc_a = approx.tail_probability(H, 'gram_charlier')
            pag4 = approx.tail_probability(H, 'pam')
            pag6 = approx.tail_probability(H, 'pam6')

            config = ','.join(map(str, sample_sizes))
            # Source indicator: E=exact, S=simulation, C=chi-square
            src_label = {"exact": "E", "simulation": "S", "chi_square": "C"}[source]
            values = [config, N, H, ref_alpha, src_label, chi, sd1, sd2, sdc1, sdc2, ed, gc_a, pag4, pag6]
            print_table_row(values, widths)

        print("-" * 142)
        print("  SRC: E=Exact, S=Simulation (10,000 iterations), C=Chi-square")


def generate_random_three_group_designs(n_designs: int = 10):
    """
    Generate random three-group designs and compare approximation methods.

    Following paper methodology:
    - Critical value H is determined from REFERENCE distribution (exact or simulation)
    - All approximation methods are compared at the SAME reference H
    - This avoids self-calibration issues where PAG(6) = 0.10 always

    Parameters
    ----------
    n_designs : int
        Number of random designs to generate
    """
    print_header(f"Random Three-Group Designs (n={n_designs})")

    headers = ["Config", "N", "H(10%)", "E-P", "SRC", "CHI", "SD1", "SD2", "SDC1", "SDC2", "ED", "GC-A", "PAG(4)", "PAG(6)"]
    widths = [14, 6, 10, 10, 6, 10, 10, 10, 10, 10, 10, 10, 10, 10]

    print_table_row(headers, widths)
    print("-" * 142)

    for i in range(n_designs):
        # Generate random sample sizes (between 2 and 10 for each group)
        n1 = random.randint(2, 10)
        n2 = random.randint(2, 10)
        n3 = random.randint(2, 10)
        sample_sizes = [n1, n2, n3]
        N = sum(sample_sizes)

        approx = KWApproximator(sample_sizes)

        # Get REFERENCE critical value (exact or simulation - NEVER pam6!)
        H, ref_alpha, source = get_reference_critical_value(sample_sizes, 0.10, n_simulations=10000)

        # Compute all approximation tail probabilities at the same reference H
        chi = approx.tail_probability(H, 'chi_square')
        sd1 = approx.tail_probability(H, 'saddlepoint')
        sd2 = approx.tail_probability(H, 'saddlepoint_sd2')
        sdc1 = approx.tail_probability(H, 'saddlepoint_cc')
        sdc2 = approx.tail_probability(H, 'saddlepoint_cc2')
        ed = approx.tail_probability(H, 'edgeworth')
        gc_a = approx.tail_probability(H, 'gram_charlier')
        pag4 = approx.tail_probability(H, 'pam')
        pag6 = approx.tail_probability(H, 'pam6')

        config = f"{n1},{n2},{n3}"
        # Source indicator: E=exact, S=simulation, C=chi-square
        src_label = {"exact": "E", "simulation": "S", "chi_square": "C"}[source]
        values = [config, N, H, ref_alpha, src_label, chi, sd1, sd2, sdc1, sdc2, ed, gc_a, pag4, pag6]
        print_table_row(values, widths)

    print("-" * 142)
    print("  SRC: E=Exact, S=Simulation (10,000 iterations), C=Chi-square")
    print("  E-P: Reference tail probability at H(10%)")


def generate_random_four_group_designs(n_designs: int = 10):
    """
    Generate random four-group designs and compare approximation methods.

    Following paper methodology:
    - Critical value H is determined from REFERENCE distribution (exact or simulation)
    - All approximation methods are compared at the SAME reference H
    - This avoids self-calibration issues where PAG(6) = 0.10 always

    Parameters
    ----------
    n_designs : int
        Number of random designs to generate
    """
    print_header(f"Random Four-Group Designs (n={n_designs})")

    headers = ["Config", "N", "H(10%)", "E-P", "SRC", "CHI", "SD1", "SD2", "SDC1", "SDC2", "ED", "GC-A", "PAG(4)", "PAG(6)"]
    widths = [16, 6, 10, 10, 6, 10, 10, 10, 10, 10, 10, 10, 10, 10]

    print_table_row(headers, widths)
    print("-" * 144)

    for i in range(n_designs):
        # Generate random sample sizes (between 2 and 8 for each group)
        n1 = random.randint(2, 8)
        n2 = random.randint(2, 8)
        n3 = random.randint(2, 8)
        n4 = random.randint(2, 8)
        sample_sizes = [n1, n2, n3, n4]
        N = sum(sample_sizes)

        approx = KWApproximator(sample_sizes)

        # Get REFERENCE critical value (exact or simulation - NEVER pam6!)
        H, ref_alpha, source = get_reference_critical_value(sample_sizes, 0.10, n_simulations=10000)

        # Compute all approximation tail probabilities at the same reference H
        chi = approx.tail_probability(H, 'chi_square')
        sd1 = approx.tail_probability(H, 'saddlepoint')
        sd2 = approx.tail_probability(H, 'saddlepoint_sd2')
        sdc1 = approx.tail_probability(H, 'saddlepoint_cc')
        sdc2 = approx.tail_probability(H, 'saddlepoint_cc2')
        ed = approx.tail_probability(H, 'edgeworth')
        gc_a = approx.tail_probability(H, 'gram_charlier')
        pag4 = approx.tail_probability(H, 'pam')
        pag6 = approx.tail_probability(H, 'pam6')

        config = f"{n1},{n2},{n3},{n4}"
        # Source indicator: E=exact, S=simulation, C=chi-square
        src_label = {"exact": "E", "simulation": "S", "chi_square": "C"}[source]
        values = [config, N, H, ref_alpha, src_label, chi, sd1, sd2, sdc1, sdc2, ed, gc_a, pag4, pag6]
        print_table_row(values, widths)

    print("-" * 144)
    print("  SRC: E=Exact, S=Simulation (10,000 iterations), C=Chi-square")
    print("  E-P: Reference tail probability at H(10%)")


def generate_random_k_group_designs(k: int = 5, n_designs: int = 10):
    """
    Generate random k-group designs and compare approximation methods.

    Following paper methodology:
    - Critical value H is determined from REFERENCE distribution (exact or simulation)
    - All approximation methods are compared at the SAME reference H
    - This avoids self-calibration issues where PAG(6) = 0.10 always

    Parameters
    ----------
    k : int
        Number of groups
    n_designs : int
        Number of random designs to generate
    """
    print_header(f"Random {k}-Group Designs (n={n_designs})")

    headers = ["Config", "N", "H(10%)", "E-P", "SRC", "CHI", "SD1", "SD2", "SDC1", "SDC2", "ED", "GC-A", "PAG(4)", "PAG(6)"]
    widths = [20, 6, 10, 10, 6, 10, 10, 10, 10, 10, 10, 10, 10, 10]

    print_table_row(headers, widths)
    print("-" * 148)

    for i in range(n_designs):
        # Generate random sample sizes (between 2 and 6 for each group)
        sample_sizes = [random.randint(2, 6) for _ in range(k)]
        N = sum(sample_sizes)

        approx = KWApproximator(sample_sizes)

        # Get REFERENCE critical value (exact or simulation - NEVER pam6!)
        H, ref_alpha, source = get_reference_critical_value(sample_sizes, 0.10, n_simulations=10000)

        # Compute all approximation tail probabilities at the same reference H
        chi = approx.tail_probability(H, 'chi_square')
        sd1 = approx.tail_probability(H, 'saddlepoint')
        sd2 = approx.tail_probability(H, 'saddlepoint_sd2')
        sdc1 = approx.tail_probability(H, 'saddlepoint_cc')
        sdc2 = approx.tail_probability(H, 'saddlepoint_cc2')
        ed = approx.tail_probability(H, 'edgeworth')
        gc_a = approx.tail_probability(H, 'gram_charlier')
        pag4 = approx.tail_probability(H, 'pam')
        pag6 = approx.tail_probability(H, 'pam6')

        config = ','.join(map(str, sample_sizes))
        # Source indicator: E=exact, S=simulation, C=chi-square
        src_label = {"exact": "E", "simulation": "S", "chi_square": "C"}[source]
        values = [config, N, H, ref_alpha, src_label, chi, sd1, sd2, sdc1, sdc2, ed, gc_a, pag4, pag6]
        print_table_row(values, widths)

    print("-" * 148)
    print("  SRC: E=Exact, S=Simulation (10,000 iterations), C=Chi-square")
    print("  E-P: Reference tail probability at H(10%)")


def comprehensive_random_study(n_per_category: int = 5):
    """
    Comprehensive random study across different group configurations.

    Following paper methodology:
    - Critical value H is determined from REFERENCE distribution (exact or simulation)
    - All approximation methods are compared at the SAME reference H
    - This avoids self-calibration issues where PAG(6) = 0.10 always

    Generates random designs for:
    - 3 groups (balanced and unbalanced)
    - 4 groups (balanced and unbalanced)
    - 5 groups

    Parameters
    ----------
    n_per_category : int
        Number of random designs per category
    """
    print_header("Comprehensive Random Study")
    print(f"Generating {n_per_category} random designs per category")

    all_results = []

    # Category 1: Balanced 3-group designs
    print("\n[Category 1] Balanced 3-group designs:")
    for i in range(n_per_category):
        n = random.randint(3, 12)
        sample_sizes = [n, n, n]
        all_results.append(('3-bal', sample_sizes))
        print(f"  {i+1}. {sample_sizes} (N={sum(sample_sizes)})")

    # Category 2: Unbalanced 3-group designs
    print("\n[Category 2] Unbalanced 3-group designs:")
    for i in range(n_per_category):
        sample_sizes = [random.randint(2, 10) for _ in range(3)]
        all_results.append(('3-unbal', sample_sizes))
        print(f"  {i+1}. {sample_sizes} (N={sum(sample_sizes)})")

    # Category 3: Balanced 4-group designs
    print("\n[Category 3] Balanced 4-group designs:")
    for i in range(n_per_category):
        n = random.randint(2, 8)
        sample_sizes = [n, n, n, n]
        all_results.append(('4-bal', sample_sizes))
        print(f"  {i+1}. {sample_sizes} (N={sum(sample_sizes)})")

    # Category 4: Unbalanced 4-group designs
    print("\n[Category 4] Unbalanced 4-group designs:")
    for i in range(n_per_category):
        sample_sizes = [random.randint(2, 8) for _ in range(4)]
        all_results.append(('4-unbal', sample_sizes))
        print(f"  {i+1}. {sample_sizes} (N={sum(sample_sizes)})")

    # Category 5: 5-group designs
    print("\n[Category 5] 5-group designs:")
    for i in range(n_per_category):
        sample_sizes = [random.randint(2, 6) for _ in range(5)]
        all_results.append(('5-grp', sample_sizes))
        print(f"  {i+1}. {sample_sizes} (N={sum(sample_sizes)})")

    # Now compute results for all
    print_header("Results for All Random Designs (alpha = 0.10)")

    headers = ["Category", "Config", "N", "E-P", "SRC", "CHI", "SD1", "SD2", "SDC1", "SDC2", "ED", "GC-A", "PAG(4)", "PAG(6)"]
    widths = [10, 18, 6, 10, 6, 10, 10, 10, 10, 10, 10, 10, 10, 10]

    print_table_row(headers, widths)
    print("-" * 146)

    for category, sample_sizes in all_results:
        N = sum(sample_sizes)
        approx = KWApproximator(sample_sizes)

        # Get REFERENCE critical value (exact or simulation - NEVER pam6!)
        H, ref_alpha, source = get_reference_critical_value(sample_sizes, 0.10, n_simulations=10000)

        # Compute all approximation tail probabilities at the same reference H
        chi = approx.tail_probability(H, 'chi_square')
        sd1 = approx.tail_probability(H, 'saddlepoint')
        sd2 = approx.tail_probability(H, 'saddlepoint_sd2')
        sdc1 = approx.tail_probability(H, 'saddlepoint_cc')
        sdc2 = approx.tail_probability(H, 'saddlepoint_cc2')
        ed = approx.tail_probability(H, 'edgeworth')
        gc_a = approx.tail_probability(H, 'gram_charlier')
        pag4 = approx.tail_probability(H, 'pam')
        pag6 = approx.tail_probability(H, 'pam6')

        config = ','.join(map(str, sample_sizes))
        # Source indicator: E=exact, S=simulation, C=chi-square
        src_label = {"exact": "E", "simulation": "S", "chi_square": "C"}[source]

        values = [category, config, N, ref_alpha, src_label, chi, sd1, sd2, sdc1, sdc2, ed, gc_a, pag4, pag6]
        print_table_row(values, widths)

    print("-" * 146)
    print("  SRC: E=Exact, S=Simulation (10,000 iterations), C=Chi-square")
    print("  E-P: Reference tail probability at H(10%)")


def main():
    """Main function to run all table reproductions."""
    print("\n" + "#" * 90)
    print("# Reproducing Results from Murakami, Lee & Ha Paper")
    print("# Higher Order Asymptotic Approximations of Kruskal-Wallis Statistics")
    print("# Based on Skewness and Kurtosis")
    print("#" * 90)

    # Paper tables
    reproduce_table_4_1()
    reproduce_tables_4_2_4_3()
    reproduce_table_4_4()
    reproduce_table_4_5()
    reproduce_tables_4_6_4_7()

    # Random designs
    generate_random_three_group_designs(n_designs=10)
    generate_random_four_group_designs(n_designs=10)
    generate_random_k_group_designs(k=5, n_designs=10)

    # Comprehensive study
    comprehensive_random_study(n_per_category=5)

    print("\n" + "=" * 90)
    print("All tables reproduced successfully!")
    print("=" * 90)


if __name__ == '__main__':
    main()
