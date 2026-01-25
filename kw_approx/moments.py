"""
Moments Calculation Module for Kruskal-Wallis Statistic (Improved Version)

Computes raw moments, central moments, cumulants, and related statistics
required for higher-order asymptotic approximations.

Key improvement: For small samples, uses exact distribution to compute
accurate moments rather than asymptotic formulas.

References:
- Murakami & Ha: Higher Order Asymptotic Approximations of Kruskal-Wallis Statistics
"""

import numpy as np
from typing import List, Tuple, Dict, Optional
from scipy.special import comb, factorial
from functools import lru_cache


class KWMoments:
    """
    Computes moments and cumulants of the Kruskal-Wallis H statistic.
    
    Under the null hypothesis, the distribution of H depends only on
    the sample sizes. This class computes exact moments using:
    1. Exact distribution enumeration (for small N <= 20)
    2. Asymptotic formulas (for large N)
    
    Parameters
    ----------
    sample_sizes : List[int]
        Sample sizes for each group [n1, n2, ..., nk]
    max_moment : int
        Maximum moment order to compute (default: 6)
    use_exact : bool or None
        If True, use exact distribution for moments.
        If False, use asymptotic formulas.
        If None (default), auto-select based on sample size.
        
    Attributes
    ----------
    k : int
        Number of groups
    N : int
        Total sample size
    raw_moments : Dict[int, float]
        Dictionary of raw moments mu_H(h) for h = 1, ..., max_moment
    cumulants : Dict[int, float]
        Dictionary of cumulants kappa_h for h = 1, ..., max_moment
    """
    
    def __init__(self, sample_sizes: List[int], max_moment: int = 6, 
                 use_exact: Optional[bool] = None):
        self.sample_sizes = np.array(sample_sizes, dtype=int)
        self.k = len(sample_sizes)
        self.N = np.sum(self.sample_sizes)
        self.max_moment = max_moment
        
        # Auto-select exact vs asymptotic
        if use_exact is None:
            self.use_exact = self.N <= 20
        else:
            self.use_exact = use_exact
        
        # Compute and cache moments
        self._raw_moments = {}
        self._central_moments = {}
        self._cumulants = {}
        
        self._compute_moments()
    
    def _compute_moments(self):
        """Compute all moments up to max_moment order."""
        if self.use_exact:
            self._compute_moments_from_exact()
        else:
            self._compute_moments_asymptotic()
        
        # Compute central moments from raw moments
        self._compute_central_moments()
        
        # Compute cumulants from central moments
        self._compute_cumulants()
    
    def _compute_moments_from_exact(self):
        """Compute moments from exact distribution (accurate for small N)."""
        # Import here to avoid circular dependency
        from .exact import ExactDistribution
        
        exact = ExactDistribution(list(self.sample_sizes))
        dist = exact.distribution
        
        # Compute raw moments directly from exact distribution
        for h in range(self.max_moment + 1):
            if h == 0:
                self._raw_moments[0] = 1.0
            else:
                self._raw_moments[h] = sum(
                    H_val**h * prob for H_val, prob in dist.items()
                )
    
    def _compute_moments_asymptotic(self):
        """Compute moments using asymptotic formulas (for large N)."""
        # Mean of H under H0: E[H] = k - 1
        self._raw_moments[0] = 1.0
        self._raw_moments[1] = float(self.k - 1)
        
        # Variance using exact formula
        var = self._compute_variance_exact()
        mu = self._raw_moments[1]
        self._raw_moments[2] = var + mu**2
        
        # Pre-compute central moments 1 and 2 for higher moment calculations
        self._central_moments[1] = 0.0
        self._central_moments[2] = var
        
        # Higher moments using asymptotic approximations
        for h in range(3, self.max_moment + 1):
            self._raw_moments[h] = self._compute_raw_moment_asymptotic(h)
    
    def _compute_variance_exact(self) -> float:
        """
        Compute exact variance of H under H0.
        
        Var(H) = 2(k-1) * (1 - correction terms for finite samples)
        """
        N = self.N
        k = self.k
        n = self.sample_sizes
        
        if k <= 1:
            return 0.0
        
        # From Kruskal (1952) and subsequent refinements
        # The exact variance formula involves sample sizes
        
        # Sum of 1/n_i
        sum_inv_n = np.sum(1.0 / n)
        
        # Finite sample correction
        # Var(H) = 2(k-1) * (N+1)/(N-1) * [1 - (sum(1/n_i) - k/N) / ((N+1)*(k-1))]
        
        if N > 1:
            factor1 = (N + 1) / (N - 1)
            factor2 = 1 - (sum_inv_n - k / N) / ((N + 1) * (k - 1))
            var = 2 * (k - 1) * factor1 * factor2
        else:
            var = 2 * (k - 1)
        
        return max(var, 1e-10)
    
    def _compute_raw_moment_asymptotic(self, h: int) -> float:
        """Compute E[H^h] for h >= 3 using asymptotic approximation."""
        mu = self._raw_moments[1]
        var = self._central_moments[2]
        
        # Use chi-square distribution as reference
        # Chi-square(k-1) has moments that can be computed exactly
        df = self.k - 1
        
        # For chi-square: E[X^h] can be computed from the MGF
        # E[X^h] = 2^h * Γ(df/2 + h) / Γ(df/2)
        from scipy.special import gamma as gamma_func
        chi2_moment_h = 2**h * gamma_func(df/2 + h) / gamma_func(df/2)
        
        # Scale to match our variance
        # If chi-square variance = 2*df, and our variance = var
        # Then scale factor is sqrt(var / (2*df))
        chi2_var = 2 * df
        if chi2_var > 0 and var > 0:
            scale = np.sqrt(var / chi2_var)
        else:
            scale = 1.0
        
        # Approximate raw moment by scaling chi-square moment
        # and adjusting mean
        if h == 3:
            m3 = self._compute_central_moment_3_asymptotic()
            return m3 + 3 * mu * var + mu**3
        elif h == 4:
            m3 = self._compute_central_moment_3_asymptotic()
            m4 = self._compute_central_moment_4_asymptotic()
            return m4 + 4 * mu * m3 + 6 * mu**2 * var + mu**4
        elif h == 5:
            m3 = self._compute_central_moment_3_asymptotic()
            m4 = self._compute_central_moment_4_asymptotic()
            m5 = self._compute_central_moment_5_asymptotic()
            return m5 + 5 * mu * m4 + 10 * mu**2 * m3 + 10 * mu**3 * var + mu**5
        elif h == 6:
            m3 = self._compute_central_moment_3_asymptotic()
            m4 = self._compute_central_moment_4_asymptotic()
            m5 = self._compute_central_moment_5_asymptotic()
            m6 = self._compute_central_moment_6_asymptotic()
            return (m6 + 6 * mu * m5 + 15 * mu**2 * m4 + 
                    20 * mu**3 * m3 + 15 * mu**4 * var + mu**6)
        else:
            # For very high moments, use scaled chi-square
            return scale**h * chi2_moment_h
    
    def _compute_central_moment_3_asymptotic(self) -> float:
        """Third central moment using asymptotic formula."""
        var = self._central_moments[2]
        df = self.k - 1
        
        # Chi-square skewness = sqrt(8/df), so m3 = skew * var^(3/2)
        skew = np.sqrt(8.0 / df) if df > 0 else 0.0
        return skew * var**(3/2)
    
    def _compute_central_moment_4_asymptotic(self) -> float:
        """Fourth central moment using asymptotic formula."""
        var = self._central_moments[2]
        df = self.k - 1
        
        # Chi-square excess kurtosis = 12/df
        kurt = 12.0 / df if df > 0 else 0.0
        return (3 + kurt) * var**2
    
    def _compute_central_moment_5_asymptotic(self) -> float:
        """Fifth central moment using asymptotic formula."""
        var = self._central_moments[2]
        m3 = self._compute_central_moment_3_asymptotic()
        m4 = self._compute_central_moment_4_asymptotic()
        
        skew = m3 / var**(3/2) if var > 0 else 0
        kurt = m4 / var**2 - 3 if var > 0 else 0
        
        # Approximate using relationship from gamma distribution
        return skew * (kurt + 6) * var**(5/2) * 0.4
    
    def _compute_central_moment_6_asymptotic(self) -> float:
        """Sixth central moment using asymptotic formula."""
        var = self._central_moments[2]
        m4 = self._compute_central_moment_4_asymptotic()
        
        kurt = m4 / var**2 - 3 if var > 0 else 0
        return (15 + 10 * kurt + kurt**2) * var**3
    
    def _compute_central_moments(self):
        """Convert raw moments to central moments."""
        mu = self._raw_moments[1]
        
        self._central_moments[1] = 0.0
        self._central_moments[2] = self._raw_moments[2] - mu**2
        
        if self.max_moment >= 3 and 3 in self._raw_moments:
            # m3 = E[(X-mu)^3] = E[X^3] - 3*mu*E[X^2] + 2*mu^3
            self._central_moments[3] = (self._raw_moments[3] - 
                                        3 * mu * self._raw_moments[2] + 
                                        2 * mu**3)
        
        if self.max_moment >= 4 and 4 in self._raw_moments:
            # m4 = E[X^4] - 4*mu*E[X^3] + 6*mu^2*E[X^2] - 3*mu^4
            self._central_moments[4] = (self._raw_moments[4] - 
                                        4 * mu * self._raw_moments[3] + 
                                        6 * mu**2 * self._raw_moments[2] - 
                                        3 * mu**4)
        
        if self.max_moment >= 5 and 5 in self._raw_moments:
            self._central_moments[5] = (self._raw_moments[5] - 
                                        5 * mu * self._raw_moments[4] + 
                                        10 * mu**2 * self._raw_moments[3] - 
                                        10 * mu**3 * self._raw_moments[2] + 
                                        4 * mu**5)
        
        if self.max_moment >= 6 and 6 in self._raw_moments:
            self._central_moments[6] = (self._raw_moments[6] - 
                                        6 * mu * self._raw_moments[5] + 
                                        15 * mu**2 * self._raw_moments[4] - 
                                        20 * mu**3 * self._raw_moments[3] + 
                                        15 * mu**4 * self._raw_moments[2] - 
                                        5 * mu**6)
    
    def _compute_cumulants(self):
        """
        Convert moments to cumulants.
        
        Cumulants are related to central moments by:
        - kappa_1 = mu (mean)
        - kappa_2 = sigma^2 (variance)
        - kappa_3 = m3 (third central moment)
        - kappa_4 = m4 - 3*sigma^4
        - etc.
        """
        self._cumulants[1] = self._raw_moments[1]
        self._cumulants[2] = self._central_moments[2]
        
        if 3 in self._central_moments:
            self._cumulants[3] = self._central_moments[3]
        
        if 4 in self._central_moments:
            self._cumulants[4] = (self._central_moments[4] - 
                                  3 * self._central_moments[2]**2)
        
        if 5 in self._central_moments:
            self._cumulants[5] = (self._central_moments[5] - 
                                  10 * self._central_moments[3] * 
                                  self._central_moments[2])
        
        if 6 in self._central_moments:
            self._cumulants[6] = (self._central_moments[6] - 
                                  15 * self._central_moments[4] * self._central_moments[2] -
                                  10 * self._central_moments[3]**2 + 
                                  30 * self._central_moments[2]**3)
    
    @property
    def raw_moments(self) -> Dict[int, float]:
        """Get dictionary of raw moments."""
        return self._raw_moments.copy()
    
    @property
    def central_moments(self) -> Dict[int, float]:
        """Get dictionary of central moments."""
        return self._central_moments.copy()
    
    @property
    def cumulants(self) -> Dict[int, float]:
        """Get dictionary of cumulants."""
        return self._cumulants.copy()
    
    def get_mean(self) -> float:
        """Get mean of H."""
        return self._raw_moments[1]
    
    def get_variance(self) -> float:
        """Get variance of H."""
        return self._central_moments[2]
    
    def get_std(self) -> float:
        """Get standard deviation of H."""
        return np.sqrt(self._central_moments[2])
    
    def get_skewness(self) -> float:
        """Get skewness of H."""
        if 3 not in self._central_moments:
            raise ValueError("Need max_moment >= 3 for skewness")
        var = self._central_moments[2]
        if var <= 0:
            return 0.0
        return self._central_moments[3] / var**(3/2)
    
    def get_kurtosis(self) -> float:
        """Get excess kurtosis of H."""
        if 4 not in self._central_moments:
            raise ValueError("Need max_moment >= 4 for kurtosis")
        var = self._central_moments[2]
        if var <= 0:
            return 0.0
        return self._central_moments[4] / var**2 - 3
    
    def get_gamma_params(self) -> Tuple[float, float]:
        """
        Get gamma distribution parameters matched to first two moments.
        
        If X ~ Gamma(alpha, beta), then:
        - E[X] = alpha * beta
        - Var[X] = alpha * beta^2
        
        Returns
        -------
        Tuple[float, float]
            (alpha, beta) parameters
        """
        mu = self._raw_moments[1]
        var = self._central_moments[2]
        
        if var <= 0 or mu <= 0:
            # Fallback to chi-square parameters
            df = self.k - 1
            return df / 2, 2.0
        
        beta = var / mu
        alpha = mu**2 / var
        
        return alpha, beta
    
    def __repr__(self) -> str:
        mode = "exact" if self.use_exact else "asymptotic"
        return (f"KWMoments(k={self.k}, N={self.N}, mode={mode}, "
                f"mean={self.get_mean():.4f}, var={self.get_variance():.4f})")
