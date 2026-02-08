"""
Unified Approximator Interface for Kruskal-Wallis Statistics

Provides a single interface to access all approximation methods and
compare their results.

References:
- Murakami & Ha: Higher Order Asymptotic Approximations of Kruskal-Wallis Statistics
"""

import numpy as np
from typing import List, Dict, Optional, Tuple, Union
from scipy import stats
import warnings

from .kruskal_wallis import KruskalWallisStatistic
from .moments import KWMoments
from .saddlepoint import SaddlepointApproximation
from .pam import PolynomialAdjustedGamma
from .gram_charlier import GramCharlierApproximation
from .edgeworth import EdgeworthApproximation
from .exact import ExactDistribution
from .simulation import MonteCarloSimulation


class KWApproximator:
    """
    Unified interface for all Kruskal-Wallis approximation methods.
    
    This class provides a convenient way to:
    1. Compute approximations using multiple methods
    2. Compare results across methods
    3. Select the most appropriate method for given sample sizes
    
    Parameters
    ----------
    sample_sizes : List[int]
        Sample sizes for each group [n1, n2, ..., nk]
    
    Attributes
    ----------
    kw : KruskalWallisStatistic
        Basic K-W statistic calculator
    moments : KWMoments
        Moment calculator
    methods : Dict
        Dictionary of initialized approximation objects
    
    Examples
    --------
    >>> approx = KWApproximator([3, 3, 3])
    >>> approx.tail_probability(4.62, method='pam')
    0.0933829
    >>> approx.compare_methods(4.62)
    {'chi_square': 0.0993, 'saddlepoint': 0.0789, ...}
    """
    
    AVAILABLE_METHODS = [
        'chi_square',      # Traditional chi-square approximation (CHI)
        'saddlepoint',     # Saddlepoint SD1 (ER1 CGF) - Lugannani-Rice
        'saddlepoint_sd2', # Saddlepoint SD2 (Wang CGF) - Lugannani-Rice
        'saddlepoint_cc',  # Saddlepoint with continuity correction (SDC1)
        'saddlepoint_cc2', # Saddlepoint SD2 (Wang) with continuity correction (SDC2)
        'saddlepoint_gamma',  # Gamma-based saddlepoint (Wood et al., 1993)
        'pam',             # Polynomial adjusted gamma PAG (degree 4)
        'pam6',            # Polynomial adjusted gamma PAG (degree 6)
        'gram_charlier',   # Gram-Charlier Type A (GC-A)
        'edgeworth',       # Edgeworth expansion (ED)
        'exact',           # Exact combinatorial (small samples only)
        'simulation'       # Monte Carlo simulation (for larger samples)
    ]
    
    def __init__(self, sample_sizes: List[int]):
        self.sample_sizes = list(sample_sizes)
        self.k = len(sample_sizes)
        self.N = sum(sample_sizes)
        
        # Initialize basic calculators
        self.kw = KruskalWallisStatistic(sample_sizes)
        self.moments = KWMoments(sample_sizes, max_moment=6)
        
        # Lazy initialization of approximation methods
        self._methods = {}
        self._exact = None
        self._simulation = None
        self._n_simulations = 10000  # Default simulation count
    
    def _get_method(self, method: str):
        """Get or initialize an approximation method."""
        if method in self._methods:
            return self._methods[method]

        if method == 'saddlepoint':
            self._methods[method] = SaddlepointApproximation(
                self.sample_sizes, cgf_method='ER1'
            )
        elif method == 'saddlepoint_sd2':
            self._methods[method] = SaddlepointApproximation(
                self.sample_sizes, cgf_method='Wang'
            )
        elif method == 'saddlepoint_cc':
            self._methods[method] = SaddlepointApproximation(
                self.sample_sizes, cgf_method='ER1'
            )
        elif method == 'saddlepoint_cc2':
            self._methods[method] = SaddlepointApproximation(
                self.sample_sizes, cgf_method='Wang'
            )
        elif method == 'saddlepoint_gamma':
            self._methods[method] = SaddlepointApproximation(
                self.sample_sizes, cgf_method='ER1'
            )
        elif method == 'pam':
            self._methods[method] = PolynomialAdjustedGamma(
                self.sample_sizes, degree=4
            )
        elif method == 'pam6':
            self._methods[method] = PolynomialAdjustedGamma(
                self.sample_sizes, degree=6
            )
        elif method == 'gram_charlier':
            self._methods[method] = GramCharlierApproximation(self.sample_sizes)
        elif method == 'edgeworth':
            self._methods[method] = EdgeworthApproximation(self.sample_sizes)
        elif method == 'exact':
            if self._exact is None:
                self._exact = ExactDistribution(self.sample_sizes)
            self._methods[method] = self._exact
        elif method == 'chi_square':
            # Chi-square doesn't need a special object
            self._methods[method] = None
        elif method == 'simulation':
            if self._simulation is None:
                self._simulation = MonteCarloSimulation(
                    self.sample_sizes, n_simulations=self._n_simulations
                )
            self._methods[method] = self._simulation
        else:
            raise ValueError(f"Unknown method: {method}. "
                           f"Available: {self.AVAILABLE_METHODS}")

        return self._methods[method]
    
    def tail_probability(self, h: float, method: str = 'pam') -> float:
        """
        Compute tail probability P(H >= h) using specified method.

        Parameters
        ----------
        h : float
            H statistic value
        method : str
            Approximation method to use

        Returns
        -------
        float
            Approximate P(H >= h)
        """
        if method == 'chi_square':
            return self.kw.chi_square_approximation(h)

        elif method == 'saddlepoint':
            sp = self._get_method(method)
            return sp.tail_probability_lr(h, continuity_correction=False)

        elif method == 'saddlepoint_sd2':
            sp = self._get_method(method)
            return sp.tail_probability_lr(h, continuity_correction=False)

        elif method == 'saddlepoint_cc':
            sp = self._get_method('saddlepoint')
            return sp.tail_probability_lr(h, continuity_correction=True)

        elif method == 'saddlepoint_cc2':
            sp = self._get_method('saddlepoint_sd2')
            return sp.tail_probability_lr(h, continuity_correction=True)

        elif method == 'saddlepoint_gamma':
            sp = self._get_method('saddlepoint')
            return sp.tail_probability_gamma_based(h)

        elif method in ['pam', 'pam6']:
            pam = self._get_method(method)
            return pam.tail_probability(h)

        elif method == 'gram_charlier':
            gc = self._get_method(method)
            return gc.tail_probability(h)

        elif method == 'edgeworth':
            ed = self._get_method(method)
            return ed.tail_probability(h)

        elif method == 'exact':
            exact = self._get_method(method)
            return exact.tail_probability(h)

        elif method == 'simulation':
            sim = self._get_method(method)
            return sim.tail_probability(h)

        else:
            raise ValueError(f"Unknown method: {method}")
    
    def cdf(self, h: float, method: str = 'pam') -> float:
        """
        Compute CDF P(H <= h) using specified method.
        
        Parameters
        ----------
        h : float
            H statistic value
        method : str
            Approximation method to use
            
        Returns
        -------
        float
            Approximate P(H <= h)
        """
        return 1.0 - self.tail_probability(h, method)
    
    def pdf(self, h: float, method: str = 'pam') -> float:
        """
        Compute approximate PDF at h using specified method.

        Parameters
        ----------
        h : float
            H statistic value
        method : str
            Approximation method to use

        Returns
        -------
        float
            Approximate density f_H(h)
        """
        if method == 'chi_square':
            return stats.chi2.pdf(h, self.k - 1)

        elif method in ['saddlepoint', 'saddlepoint_cc', 'saddlepoint_sd2', 'saddlepoint_cc2']:
            if 'sd2' in method:
                sp = self._get_method('saddlepoint_sd2')
            else:
                sp = self._get_method('saddlepoint')
            cc = 'cc' in method
            return sp.density_approximation(h, continuity_correction=cc)

        elif method in ['pam', 'pam6']:
            pam = self._get_method(method)
            return pam.pdf(h)

        elif method == 'gram_charlier':
            gc = self._get_method(method)
            return gc.pdf(h)

        elif method == 'edgeworth':
            ed = self._get_method(method)
            return ed.pdf(h)

        elif method == 'exact':
            exact = self._get_method(method)
            return exact.pmf(h)

        else:
            raise ValueError(f"Unknown method: {method}")
    
    def critical_value(self, alpha: float, method: str = 'pam') -> float:
        """
        Compute critical value for significance level alpha.

        Finds c such that P(H >= c) = alpha.

        Parameters
        ----------
        alpha : float
            Significance level (e.g., 0.05, 0.10)
        method : str
            Approximation method to use

        Returns
        -------
        float
            Critical value c
        """
        if method == 'chi_square':
            return stats.chi2.ppf(1 - alpha, self.k - 1)

        elif method in ['saddlepoint', 'saddlepoint_cc', 'saddlepoint_sd2', 'saddlepoint_cc2']:
            if 'sd2' in method:
                sp = self._get_method('saddlepoint_sd2')
            else:
                sp = self._get_method('saddlepoint')
            cc = 'cc' in method
            return sp.critical_value(alpha, continuity_correction=cc)

        elif method in ['pam', 'pam6']:
            pam = self._get_method(method)
            return pam.critical_value(alpha)

        elif method == 'gram_charlier':
            gc = self._get_method(method)
            return gc.critical_value(alpha)

        elif method == 'edgeworth':
            ed = self._get_method(method)
            return ed.critical_value(alpha)

        elif method == 'exact':
            exact = self._get_method(method)
            cv, exact_alpha = exact.critical_value(alpha)
            return cv

        else:
            raise ValueError(f"Unknown method: {method}")
    
    def compare_methods(self, h: float,
                       methods: Optional[List[str]] = None) -> Dict[str, float]:
        """
        Compare tail probabilities across multiple methods.

        Parameters
        ----------
        h : float
            H statistic value
        methods : List[str], optional
            Methods to compare. If None, uses all available methods
            (excluding 'exact' for large samples)

        Returns
        -------
        Dict[str, float]
            Dictionary mapping method names to tail probabilities
        """
        if methods is None:
            methods = ['chi_square', 'saddlepoint', 'saddlepoint_sd2',
                      'saddlepoint_cc', 'saddlepoint_cc2',
                      'pam', 'pam6', 'gram_charlier', 'edgeworth']
            # Only include exact for small samples
            if self.N <= 20:
                methods.append('exact')

        results = {}
        for method in methods:
            try:
                results[method] = self.tail_probability(h, method)
            except Exception as e:
                results[method] = None
                warnings.warn(f"Method {method} failed: {e}")

        return results
    
    def compare_critical_values(self, alpha: float,
                               methods: Optional[List[str]] = None) -> Dict[str, float]:
        """
        Compare critical values across multiple methods.

        Parameters
        ----------
        alpha : float
            Significance level
        methods : List[str], optional
            Methods to compare

        Returns
        -------
        Dict[str, float]
            Dictionary mapping method names to critical values
        """
        if methods is None:
            methods = ['chi_square', 'saddlepoint', 'saddlepoint_sd2',
                      'saddlepoint_cc', 'pam', 'pam6', 'gram_charlier', 'edgeworth']
            if self.N <= 20:
                methods.append('exact')

        results = {}
        for method in methods:
            try:
                results[method] = self.critical_value(alpha, method)
            except Exception as e:
                results[method] = None
                warnings.warn(f"Method {method} failed: {e}")

        return results
    
    def recommend_method(self) -> str:
        """
        Recommend the best approximation method based on sample sizes.
        
        Returns
        -------
        str
            Recommended method name
        """
        if self.N <= 15:
            return 'exact'
        elif self.N <= 30:
            return 'pam6'
        elif self.N <= 100:
            return 'pam'
        else:
            return 'saddlepoint'
    
    def generate_table(self, h_values: Optional[List[float]] = None,
                      methods: Optional[List[str]] = None,
                      alpha_values: Optional[List[float]] = None) -> Dict:
        """
        Generate comparison tables similar to those in the paper.

        Parameters
        ----------
        h_values : List[float], optional
            H statistic values to evaluate
        methods : List[str], optional
            Methods to include
        alpha_values : List[float], optional
            Significance levels for critical values

        Returns
        -------
        Dict
            Dictionary with 'probabilities' and 'critical_values' tables
        """
        if methods is None:
            methods = ['chi_square', 'saddlepoint', 'saddlepoint_sd2',
                      'saddlepoint_cc', 'pam', 'pam6', 'gram_charlier', 'edgeworth']
            if self.N <= 20:
                methods.append('exact')
        
        if alpha_values is None:
            alpha_values = [0.10, 0.05, 0.01]
        
        # Generate probability comparison table
        prob_table = {}
        if h_values is not None:
            for h in h_values:
                prob_table[h] = self.compare_methods(h, methods)
        
        # Generate critical value comparison table
        cv_table = {}
        for alpha in alpha_values:
            cv_table[alpha] = self.compare_critical_values(alpha, methods)
        
        return {
            'probabilities': prob_table,
            'critical_values': cv_table,
            'sample_sizes': self.sample_sizes,
            'methods': methods
        }
    
    def summary(self) -> Dict:
        """
        Get summary information about the approximation setup.
        
        Returns
        -------
        Dict
            Summary statistics and parameters
        """
        return {
            'sample_sizes': self.sample_sizes,
            'k': self.k,
            'N': self.N,
            'df': self.k - 1,
            'mean_H': self.moments.get_mean(),
            'var_H': self.moments.get_variance(),
            'std_H': self.moments.get_std(),
            'skewness': self.moments.get_skewness(),
            'kurtosis': self.moments.get_kurtosis(),
            'gamma_params': self.moments.get_gamma_params(),
            'recommended_method': self.recommend_method()
        }
    
    def test_statistic(self, *groups: np.ndarray) -> Tuple[float, Dict[str, float]]:
        """
        Compute test statistic from data and p-values using multiple methods.
        
        Parameters
        ----------
        *groups : np.ndarray
            Data arrays for each group
            
        Returns
        -------
        Tuple[float, Dict[str, float]]
            (H_statistic, dict of p-values by method)
        """
        H = self.kw.compute_statistic(*groups)
        p_values = self.compare_methods(H)
        return H, p_values
    
    def __repr__(self) -> str:
        return (f"KWApproximator(k={self.k}, N={self.N}, "
                f"sample_sizes={self.sample_sizes})")


def quick_test(sample_sizes: List[int], h: float, alpha: float = 0.10) -> None:
    """
    Quick test function to demonstrate all approximation methods.
    
    Parameters
    ----------
    sample_sizes : List[int]
        Sample sizes for each group
    h : float
        H statistic value to evaluate
    alpha : float
        Significance level for critical values
    """
    N = sum(sample_sizes)
    k = len(sample_sizes)
    
    print(f"=" * 60)
    print(f"Kruskal-Wallis Approximation Comparison")
    print(f"=" * 60)
    print(f"Sample sizes: {sample_sizes}")
    print(f"k = {k}, N = {N}")
    print(f"H value: {h}")
    print(f"Significance level: {alpha}")
    print(f"-" * 60)
    
    approx = KWApproximator(sample_sizes)
    recommended = approx.recommend_method()
    
    # Method selection based on paper recommendations
    print(f"\n[Paper-Recommended Method Selection]")
    if N <= 15:
        print(f"  N={N} <= 15: Using exact (exact distribution) for computation.")
    elif N <= 30:
        print(f"  N={N} > 15: Using pam6 (polynomial adjusted gamma, degree=6) for computation.")
    elif N <= 100:
        print(f"  N={N} > 30: Using pam (polynomial adjusted gamma, degree=4) for computation.")
    else:
        print(f"  N={N} > 100: Using saddlepoint for computation.")
    
    # Summary
    summary = approx.summary()
    print(f"\nDistribution Parameters:")
    print(f"  Mean(H):     {summary['mean_H']:.6f}")
    print(f"  Var(H):      {summary['var_H']:.6f}")
    print(f"  Skewness:    {summary['skewness']:.6f}")
    print(f"  Kurtosis:    {summary['kurtosis']:.6f}")
    print(f"  Recommended: {recommended}")
    
    # Tail probabilities
    print(f"\nTail Probabilities P(H >= {h}):")
    print(f"  {'Method':<20} {'P-value':>12}")
    print(f"  {'-'*20} {'-'*12}")
    
    probs = approx.compare_methods(h)
    for method, prob in probs.items():
        if prob is not None:
            print(f"  {method:<20} {prob:>12.6f}")
        else:
            print(f"  {method:<20} {'N/A':>12}")
    
    # Critical values
    print(f"\nCritical Values (alpha = {alpha}):")
    print(f"  {'Method':<20} {'Critical Value':>15}")
    print(f"  {'-'*20} {'-'*15}")
    
    cvs = approx.compare_critical_values(alpha)
    for method, cv in cvs.items():
        if cv is not None:
            print(f"  {method:<20} {cv:>15.6f}")
        else:
            print(f"  {method:<20} {'N/A':>15}")
    
    print(f"=" * 60)