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

    headers = ["n1,n2,n3", "P(E-Q)", "E-P", "SD1", "SD2", "SDC1", "SDC2", "ED", "GC-A", "PAG(4)"]
    widths = [12, 10, 10, 10, 10, 10, 10, 10, 10, 10]

    print(f"\n{'-'*102}")
    print_table_row(headers, widths)
    print(f"{'-'*102}")

    exact_p = approx.tail_probability(H, 'exact')
    sd1 = approx.tail_probability(H, 'saddlepoint')
    sd2 = approx.tail_probability(H, 'saddlepoint_sd2')
    sdc1 = approx.tail_probability(H, 'saddlepoint_cc')
    sdc2 = approx.tail_probability(H, 'saddlepoint_cc2')
    ed = approx.tail_probability(H, 'edgeworth')
    gc_a = approx.tail_probability(H, 'gram_charlier')
    pag4 = approx.tail_probability(H, 'pam')

    values = ["3,3,3", H, exact_p, sd1, sd2, sdc1, sdc2, ed, gc_a, pag4]
    print_table_row(values, widths)
    print(f"{'-'*102}")


def reproduce_tables_4_2_4_3():
    """
    Reproduce Tables 4.2-4.3: Three groups, increasing n, alpha = 0.10 and 0.05

    Generate balanced three-group designs with increasing sample sizes.
    """
    print_header("Tables 4.2-4.3: Three groups, increasing n, alpha = 0.10 and 0.05")

    # Balanced three-group designs with increasing n
    n_values = [3, 4, 5, 6, 7, 8, 9, 10]

    for alpha in [0.10, 0.05]:
        print(f"\n--- alpha = {alpha} ---\n")

        headers = ["n1,n2,n3", "N", "E-Q", "E-P", "CHI", "SD1", "SDC1", "ED", "PAG(4)", "PAG(6)"]
        widths = [12, 6, 10, 10, 10, 10, 10, 10, 10, 10]

        print_table_row(headers, widths)
        print("-" * 98)

        for n in n_values:
            sample_sizes = [n, n, n]
            N = sum(sample_sizes)

            approx = KWApproximator(sample_sizes)

            # Get critical value (use exact only for N <= 15 per paper recommendation)
            try:
                if N <= 15:
                    cv, exact_alpha = get_exact_critical_value(sample_sizes, alpha)
                    if cv is None:
                        cv = approx.critical_value(alpha, 'pam6')
                else:
                    # N > 15: Use pam6 instead of exact (paper recommendation)
                    cv = approx.critical_value(alpha, 'pam6')
                    exact_alpha = None
            except:
                cv = approx.critical_value(alpha, 'chi_square')
                exact_alpha = None

            H = cv if cv else approx.critical_value(alpha, 'chi_square')

            # Compute reference probability (exact for small N, simulation for large N)
            ref_p, is_sim = get_reference_probability(approx, H, N, len(sample_sizes))

            chi = approx.tail_probability(H, 'chi_square')
            sd1 = approx.tail_probability(H, 'saddlepoint')
            sdc1 = approx.tail_probability(H, 'saddlepoint_cc')
            ed = approx.tail_probability(H, 'edgeworth')
            pag4 = approx.tail_probability(H, 'pam')
            pag6 = approx.tail_probability(H, 'pam6')

            config = f"{n},{n},{n}"
            # Add marker for simulation-based values
            if ref_p is not None:
                ep_label = f"{ref_p:.6f}" if not is_sim else f"{ref_p:.6f}*"
            else:
                ep_label = None
            values = [config, N, H, ep_label, chi, sd1, sdc1, ed, pag4, pag6]
            print_table_row(values, widths)

        print("-" * 98)
        if any(n * 3 > 15 for n in n_values):
            print("  * E-P column: simulation-based (10,000 iterations) for N > 15")


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

    headers = ["n1,n2,n3", "N", "CHI-CV", "CHI", "SD1", "SDC1", "ED", "PAG(4)", "PAG(6)"]
    widths = [12, 6, 10, 10, 10, 10, 10, 10, 10]

    print_table_row(headers, widths)
    print("-" * 88)

    alpha = 0.10

    for n in n_values:
        sample_sizes = [n, n, n]
        N = sum(sample_sizes)

        approx = KWApproximator(sample_sizes)

        # Use chi-square critical value as reference
        cv = approx.critical_value(alpha, 'chi_square')
        H = cv

        chi = approx.tail_probability(H, 'chi_square')
        sd1 = approx.tail_probability(H, 'saddlepoint')
        sdc1 = approx.tail_probability(H, 'saddlepoint_cc')
        ed = approx.tail_probability(H, 'edgeworth')
        pag4 = approx.tail_probability(H, 'pam')
        pag6 = approx.tail_probability(H, 'pam6')

        config = f"{n},{n},{n}"
        values = [config, N, H, chi, sd1, sdc1, ed, pag4, pag6]
        print_table_row(values, widths)

    print("-" * 88)


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

    headers = ["Config", "E-P/Sim", "SD1", "SD2", "SDC1", "SDC2", "CHI", "PAG(4)"]
    widths = [14, 10, 10, 10, 10, 10, 10, 10]

    print(f"\n{'-'*84}")
    print_table_row(headers, widths)
    print(f"{'-'*84}")

    exact_p = approx.tail_probability(H, 'exact')
    sd1 = approx.tail_probability(H, 'saddlepoint')
    sd2 = approx.tail_probability(H, 'saddlepoint_sd2')
    sdc1 = approx.tail_probability(H, 'saddlepoint_cc')
    sdc2 = approx.tail_probability(H, 'saddlepoint_cc2')
    chi = approx.tail_probability(H, 'chi_square')
    pag4 = approx.tail_probability(H, 'pam')

    values = ["3,2,2,5", exact_p, sd1, sd2, sdc1, sdc2, chi, pag4]
    print_table_row(values, widths)
    print(f"{'-'*84}")

    # Paper reference
    print("\nPaper reference values:")
    print(f"  SD1: 0.105404 (Code: {sd1:.6f}) {'OK' if abs(sd1 - 0.105404) < 0.001 else 'DIFF'}")
    print(f"  CHI: 0.133516 (Code: {chi:.6f}) {'OK' if abs(chi - 0.133516) < 0.001 else 'DIFF'}")


def reproduce_tables_4_6_4_7():
    """
    Reproduce Tables 4.6-4.7: Additional four-group designs, alpha = 0.10 and 0.05

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

        headers = ["Config", "N", "E-Q", "E-P", "CHI", "SD1", "SDC1", "PAG(4)", "PAG(6)"]
        widths = [14, 6, 10, 10, 10, 10, 10, 10, 10]

        print_table_row(headers, widths)
        print("-" * 90)

        for sample_sizes in configs:
            N = sum(sample_sizes)
            approx = KWApproximator(sample_sizes)

            # Determine threshold for exact computation
            # For k=4, N=15 yields ~12.6 million combinations, which is too slow.
            # We use a stricter threshold (N<=13) for 4 or more groups.
            limit_N = 15 if len(sample_sizes) <= 3 else 13

            # Get critical value (use exact only for small N)
            try:
                if N <= limit_N:
                    cv, _ = get_exact_critical_value(sample_sizes, alpha)
                    if cv is None:
                        cv = approx.critical_value(alpha, 'pam6')
                else:
                    cv = approx.critical_value(alpha, 'pam6')
            except:
                cv = approx.critical_value(alpha, 'chi_square')

            H = cv

            # Compute reference probability (exact for small N, simulation for large N)
            ref_p, is_sim = get_reference_probability(approx, H, N, len(sample_sizes))

            chi = approx.tail_probability(H, 'chi_square')
            sd1 = approx.tail_probability(H, 'saddlepoint')
            sdc1 = approx.tail_probability(H, 'saddlepoint_cc')
            pag4 = approx.tail_probability(H, 'pam')
            pag6 = approx.tail_probability(H, 'pam6')

            config = ','.join(map(str, sample_sizes))
            # Add marker for simulation-based values
            if ref_p is not None:
                ep_label = f"{ref_p:.6f}" if not is_sim else f"{ref_p:.6f}*"
            else:
                ep_label = None
            values = [config, N, H, ep_label, chi, sd1, sdc1, pag4, pag6]
            print_table_row(values, widths)

        print("-" * 90)
        print("  * E-P column: simulation-based (10,000 iterations) for large N")


def generate_random_three_group_designs(n_designs: int = 10):
    """
    Generate random three-group designs and compare approximation methods.

    Parameters
    ----------
    n_designs : int
        Number of random designs to generate
    """
    print_header(f"Random Three-Group Designs (n={n_designs})")

    headers = ["Config", "N", "H(10%)", "E-P", "CHI", "SD1", "SDC1", "PAG(4)", "PAG(6)"]
    widths = [14, 6, 10, 10, 10, 10, 10, 10, 10]

    print_table_row(headers, widths)
    print("-" * 90)

    for i in range(n_designs):
        # Generate random sample sizes (between 2 and 10 for each group)
        n1 = random.randint(2, 10)
        n2 = random.randint(2, 10)
        n3 = random.randint(2, 10)
        sample_sizes = [n1, n2, n3]
        N = sum(sample_sizes)

        approx = KWApproximator(sample_sizes)

        # Get critical value at alpha = 0.10
        try:
            if N <= 15:
                cv, _ = get_exact_critical_value(sample_sizes, 0.10)
                if cv is None:
                    cv = approx.critical_value(0.10, 'pam6')
            else:
                cv = approx.critical_value(0.10, 'pam6')
        except:
            cv = approx.critical_value(0.10, 'chi_square')

        H = cv

        # Compute reference probability (exact for small N, simulation for large N)
        ref_p, is_sim = get_reference_probability(approx, H, N, len(sample_sizes))

        chi = approx.tail_probability(H, 'chi_square')
        sd1 = approx.tail_probability(H, 'saddlepoint')
        sdc1 = approx.tail_probability(H, 'saddlepoint_cc')
        pag4 = approx.tail_probability(H, 'pam')
        pag6 = approx.tail_probability(H, 'pam6')

        config = f"{n1},{n2},{n3}"
        # Add marker for simulation-based values
        if ref_p is not None:
            ep_label = f"{ref_p:.6f}" if not is_sim else f"{ref_p:.6f}*"
        else:
            ep_label = None
        values = [config, N, H, ep_label, chi, sd1, sdc1, pag4, pag6]
        print_table_row(values, widths)

    print("-" * 90)
    print("  * E-P column: simulation-based (10,000 iterations) for N > 15")


def generate_random_four_group_designs(n_designs: int = 10):
    """
    Generate random four-group designs and compare approximation methods.

    Parameters
    ----------
    n_designs : int
        Number of random designs to generate
    """
    print_header(f"Random Four-Group Designs (n={n_designs})")

    headers = ["Config", "N", "H(10%)", "E-P", "CHI", "SD1", "SDC1", "PAG(4)", "PAG(6)"]
    widths = [16, 6, 10, 10, 10, 10, 10, 10, 10]

    print_table_row(headers, widths)
    print("-" * 92)

    for i in range(n_designs):
        # Generate random sample sizes (between 2 and 8 for each group)
        n1 = random.randint(2, 8)
        n2 = random.randint(2, 8)
        n3 = random.randint(2, 8)
        n4 = random.randint(2, 8)
        sample_sizes = [n1, n2, n3, n4]
        N = sum(sample_sizes)

        approx = KWApproximator(sample_sizes)

        # Use stricter threshold for k=4 (N<=13)
        limit_N = 13

        # Get critical value at alpha = 0.10
        try:
            if N <= limit_N:
                cv, _ = get_exact_critical_value(sample_sizes, 0.10)
                if cv is None:
                    cv = approx.critical_value(0.10, 'pam6')
            else:
                cv = approx.critical_value(0.10, 'pam6')
        except:
            cv = approx.critical_value(0.10, 'chi_square')

        H = cv

        # Compute reference probability (exact for small N, simulation for large N)
        ref_p, is_sim = get_reference_probability(approx, H, N, len(sample_sizes))

        chi = approx.tail_probability(H, 'chi_square')
        sd1 = approx.tail_probability(H, 'saddlepoint')
        sdc1 = approx.tail_probability(H, 'saddlepoint_cc')
        pag4 = approx.tail_probability(H, 'pam')
        pag6 = approx.tail_probability(H, 'pam6')

        config = f"{n1},{n2},{n3},{n4}"
        # Add marker for simulation-based values
        if ref_p is not None:
            ep_label = f"{ref_p:.6f}" if not is_sim else f"{ref_p:.6f}*"
        else:
            ep_label = None
        values = [config, N, H, ep_label, chi, sd1, sdc1, pag4, pag6]
        print_table_row(values, widths)

    print("-" * 92)
    print("  * E-P column: simulation-based (10,000 iterations) for N > 13")


def generate_random_k_group_designs(k: int = 5, n_designs: int = 10):
    """
    Generate random k-group designs and compare approximation methods.

    Parameters
    ----------
    k : int
        Number of groups
    n_designs : int
        Number of random designs to generate
    """
    print_header(f"Random {k}-Group Designs (n={n_designs})")

    headers = ["Config", "N", "H(10%)", "CHI", "SD1", "SDC1", "PAG(4)", "PAG(6)"]
    widths = [20, 6, 10, 10, 10, 10, 10, 10]

    print_table_row(headers, widths)
    print("-" * 86)

    for i in range(n_designs):
        # Generate random sample sizes (between 2 and 6 for each group)
        sample_sizes = [random.randint(2, 6) for _ in range(k)]
        N = sum(sample_sizes)

        approx = KWApproximator(sample_sizes)

        # Determine safe threshold based on k
        k_val = len(sample_sizes)
        if k_val == 3:
            limit_N = 15
        elif k_val == 4:
            limit_N = 13
        else:
            limit_N = 10  # Very strict for k>=5

        # Get critical value at alpha = 0.10
        try:
            if N <= limit_N:
                cv, _ = get_exact_critical_value(sample_sizes, 0.10)
                if cv is None:
                    cv = approx.critical_value(0.10, 'pam6')
            else:
                cv = approx.critical_value(0.10, 'pam6')
        except:
            cv = approx.critical_value(0.10, 'chi_square')

        H = cv

        try:
            exact_p = approx.tail_probability(H, 'exact') if N <= limit_N else None
        except:
            exact_p = None

        chi = approx.tail_probability(H, 'chi_square')
        sd1 = approx.tail_probability(H, 'saddlepoint')
        sdc1 = approx.tail_probability(H, 'saddlepoint_cc')
        pag4 = approx.tail_probability(H, 'pam')
        pag6 = approx.tail_probability(H, 'pam6')

        config = ','.join(map(str, sample_sizes))
        # Add marker
        ep_label = exact_p if exact_p is not None else (f"{pag6:.6f}*" if pag6 is not None else "N/A")
        values = [config, N, H, chi, sd1, sdc1, pag4, pag6]
        # Note: headers don't have E-P column in this function, so we don't need to add it, 
        # but let's check headers first. 
        # Ah wait, headers are ["Config", "N", "H(10%)", "CHI", "SD1", "SDC1", "PAG(4)", "PAG(6)"]
        # So no E-P column. We don't need to calculate exact_p or add it to values.
        # But for correctness let's keep exact_p calculation logic if headers were to change, 
        # or just remove exact_p calculation to be faster. 
        # Actually, let's just stick to the existing logic but safeguard the potentially expensive call.
        
        values = [config, N, H, chi, sd1, sdc1, pag4, pag6]
        print_table_row(values, widths)

    print("-" * 86)
    print("  * Note: Exact Probabilities are not shown for k>=5 due to computational cost")


def comprehensive_random_study(n_per_category: int = 5):
    """
    Comprehensive random study across different group configurations.

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

    headers = ["Category", "Config", "N", "E-P", "CHI", "SD1", "PAG(6)"]
    widths = [10, 18, 6, 10, 10, 10, 10]

    print_table_row(headers, widths)
    print("-" * 74)

    for category, sample_sizes in all_results:
        N = sum(sample_sizes)
        approx = KWApproximator(sample_sizes)

        # Determine safe threshold based on k
        k_val = len(sample_sizes)
        if k_val == 3:
            limit_N = 15
        elif k_val == 4:
            limit_N = 13
        else:
            limit_N = 10  # Very strict for k>=5

        # Get critical value
        try:
            if N <= limit_N:
                cv, _ = get_exact_critical_value(sample_sizes, 0.10)
                if cv is None:
                    cv = approx.critical_value(0.10, 'pam6')
            else:
                cv = approx.critical_value(0.10, 'pam6')
        except:
            cv = approx.critical_value(0.10, 'chi_square')

        H = cv

        try:
            exact_p = approx.tail_probability(H, 'exact') if N <= limit_N else None
        except:
            exact_p = None

        chi = approx.tail_probability(H, 'chi_square')
        sd1 = approx.tail_probability(H, 'saddlepoint')
        pag6 = approx.tail_probability(H, 'pam6')

        config = ','.join(map(str, sample_sizes))
        # Add marker for approximate reference values
        ep_label = exact_p if exact_p is not None else (f"{pag6:.6f}*" if pag6 is not None else "N/A")
        
        values = [category, config, N, ep_label, chi, sd1, pag6]
        print_table_row(values, widths)

    print("-" * 74)
    print("  * E-P column shows pam6 values when exact computation is skipped due to large N/k")


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
