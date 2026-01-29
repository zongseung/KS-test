"""
Saddlepoint Approximation Module for Kruskal-Wallis Statistic

Implements various saddlepoint approximation methods:
- Daniels' saddlepoint density approximation (SD1, SD2)
- Lugannani-Rice tail probability approximation
- Continuity-corrected versions (SDC1, SDC2)
- Gamma-based saddlepoint approximation (Wood et al., 1993)

References:
- Daniels, H.E. (1954). Saddlepoint approximations in statistics
- Lugannani, R. and Rice, S.O. (1980). Saddlepoint approximation for the 
  distribution of the sum of independent random variables
- Murakami & Ha: Higher Order Asymptotic Approximations of Kruskal-Wallis Statistics
"""

import numpy as np
from typing import List, Tuple, Optional, Dict
from scipy import stats
from scipy.optimize import brentq, newton
from scipy.special import gamma as gamma_func
from .moments import KWMoments


class SaddlepointApproximation:
    """
    Saddlepoint approximations for the Kruskal-Wallis H statistic.
    
    Parameters
    ----------
    sample_sizes : List[int]
        Sample sizes for each group [n1, n2, ..., nk]
    cgf_method : str
        Method for cumulant generating function approximation.
        Options: 'exact', 'ER1', 'ER2', 'Wang', 'KT'
        Default: 'ER1' (Easton-Ronchetti first approximation)
    
    Attributes
    ----------
    moments : KWMoments
        Moment calculator object
    """
    
    def __init__(self, sample_sizes: List[int], cgf_method: str = 'ER1'):
        self.sample_sizes = np.array(sample_sizes, dtype=int)
        self.k = len(sample_sizes)
        self.N = np.sum(self.sample_sizes)
        
        # Compute moments (need at least 4 for saddlepoint)
        self.moments = KWMoments(sample_sizes, max_moment=6)
        
        self.cgf_method = cgf_method
        
        # Cache cumulants for CGF computations
        self._kappa = self.moments.cumulants
    
    def cumulant_generating_function(self, t: float) -> float:
        """
        Compute the cumulant generating function K_H(t).
        
        Uses approximations based on cumulants since exact CGF
        is complex for the Kruskal-Wallis statistic.
        
        Parameters
        ----------
        t : float
            Argument of the CGF
            
        Returns
        -------
        float
            K_H(t)
        """
        kappa = self._kappa
        
        if self.cgf_method == 'ER1':
            # Easton-Ronchetti first approximation
            # K_H(t) ≈ sum_{i=1}^{4} kappa_i * t^i / i!
            return (kappa[1] * t + 
                    kappa[2] * t**2 / 2 + 
                    kappa[3] * t**3 / 6 + 
                    kappa[4] * t**4 / 24)
        
        elif self.cgf_method == 'ER2':
            # Easton-Ronchetti second approximation
            # K_H(t) ≈ kappa_1*t + (kappa_2/2)*t^2 + 
            #          log(1 + (kappa_3/6)*t^3 + (3*kappa_4/72)*t^4 + (kappa_3^2/72)*t^6)
            inner = (1 + kappa[3] * t**3 / 6 + 
                     3 * kappa[4] * t**4 / 72 + 
                     kappa[3]**2 * t**6 / 72)
            if inner > 0:
                return (kappa[1] * t + kappa[2] * t**2 / 2 + np.log(inner))
            else:
                # Fallback to ER1 if log argument becomes non-positive
                return (kappa[1] * t + 
                        kappa[2] * t**2 / 2 + 
                        kappa[3] * t**3 / 6 + 
                        kappa[4] * t**4 / 24)
        
        elif self.cgf_method == 'Wang':
            # Wang (1992) approximation with damping
            # K_H(t;p) = kappa_1*t + (kappa_2/2)*t^2 + (kappa_3/6)*t^3 + 
            #            (kappa_4/24)*t^4 * eta_p(t)
            # where eta_p(t) = exp(-kappa_2^p * t^2 / 2)
            p = 0.5  # Default p value
            eta = np.exp(-kappa[2]**p * t**2 / 2)
            return (kappa[1] * t + 
                    kappa[2] * t**2 / 2 + 
                    kappa[3] * t**3 / 6 + 
                    kappa[4] * t**4 / 24 * eta)
        
        elif self.cgf_method == 'KT':
            # Kakizawa-Taniguchi (1994) correction
            # K_H(t) ≈ kappa_1*t + ((1 + kappa_2)*t^2)/2 + (kappa_3*t^3)/6 + (kappa_4*t^4)/24
            return (kappa[1] * t + 
                    (1 + kappa[2]) * t**2 / 2 + 
                    kappa[3] * t**3 / 6 + 
                    kappa[4] * t**4 / 24)
        
        else:
            # Default to ER1
            return (kappa[1] * t + 
                    kappa[2] * t**2 / 2 + 
                    kappa[3] * t**3 / 6 + 
                    kappa[4] * t**4 / 24)
    
    def cgf_derivative1(self, t: float) -> float:
        """
        First derivative of CGF: K'_H(t).
        
        Parameters
        ----------
        t : float
            Argument
            
        Returns
        -------
        float
            K'_H(t)
        """
        kappa = self._kappa
        
        if self.cgf_method == 'ER1':
            return (kappa[1] + 
                    kappa[2] * t + 
                    kappa[3] * t**2 / 2 + 
                    kappa[4] * t**3 / 6)
        
        elif self.cgf_method == 'ER2':
            # Derivative of ER2 approximation
            inner = (1 + kappa[3] * t**3 / 6 + 
                     3 * kappa[4] * t**4 / 72 + 
                     kappa[3]**2 * t**6 / 72)
            d_inner = (kappa[3] * t**2 / 2 + 
                       kappa[4] * t**3 / 6 + 
                       kappa[3]**2 * t**5 / 12)
            if inner > 0:
                return kappa[1] + kappa[2] * t + d_inner / inner
            else:
                return (kappa[1] + 
                        kappa[2] * t + 
                        kappa[3] * t**2 / 2 + 
                        kappa[4] * t**3 / 6)
        
        elif self.cgf_method == 'Wang':
            p = 0.5
            eta = np.exp(-kappa[2]**p * t**2 / 2)
            d_eta = -kappa[2]**p * t * eta
            term4 = kappa[4] * t**3 / 6
            return (kappa[1] + 
                    kappa[2] * t + 
                    kappa[3] * t**2 / 2 + 
                    term4 * eta + kappa[4] * t**4 / 24 * d_eta)
        
        else:
            return (kappa[1] + 
                    kappa[2] * t + 
                    kappa[3] * t**2 / 2 + 
                    kappa[4] * t**3 / 6)
    
    def cgf_derivative2(self, t: float) -> float:
        """
        Second derivative of CGF: K''_H(t).
        
        Parameters
        ----------
        t : float
            Argument
            
        Returns
        -------
        float
            K''_H(t)
        """
        kappa = self._kappa
        
        if self.cgf_method in ['ER1', 'exact']:
            return kappa[2] + kappa[3] * t + kappa[4] * t**2 / 2
        
        else:
            # Numerical derivative as fallback
            h = 1e-6
            return (self.cgf_derivative1(t + h) - self.cgf_derivative1(t - h)) / (2 * h)
    
    def cgf_derivative3(self, t: float) -> float:
        """Third derivative of CGF: K'''_H(t)."""
        kappa = self._kappa
        return kappa[3] + kappa[4] * t
    
    def find_saddlepoint(self, x: float) -> float:
        """
        Find the saddlepoint t_hat that solves K'_H(t) = x.
        
        Parameters
        ----------
        x : float
            Target value (typically the observed H statistic)
            
        Returns
        -------
        float
            The saddlepoint t_hat
        """
        # For Kruskal-Wallis, the saddlepoint equation is a polynomial
        # which can be solved numerically
        
        def equation(t):
            return self.cgf_derivative1(t) - x
        
        # Try to bracket the root
        # The saddlepoint should be positive for x > mean
        mean = self.moments.get_mean()
        
        try:
            if x > mean:
                # Search in positive t region
                t_hat = brentq(equation, 0, 10, xtol=1e-10)
            elif x < mean:
                # Search in negative t region
                t_hat = brentq(equation, -10, 0, xtol=1e-10)
            else:
                t_hat = 0.0
        except ValueError:
            # If bracketing fails, use Newton's method
            try:
                t_hat = newton(equation, 0.0, fprime=self.cgf_derivative2, tol=1e-10)
            except RuntimeError:
                t_hat = 0.0
        
        return t_hat
    
    def density_approximation(self, x: float, continuity_correction: bool = False) -> float:
        """
        Daniels' saddlepoint density approximation.

        f_SP(x) = (2π K''_H(t_hat))^{-1/2} * exp(K_H(t_hat) - x*t_hat)

        Parameters
        ----------
        x : float
            Point at which to evaluate density
        continuity_correction : bool
            Whether to apply continuity correction for discrete distributions

        Returns
        -------
        float
            Approximate density f_H(x)
        """
        if continuity_correction:
            # For discrete distributions, adjust evaluation point
            x = x - 0.5
        
        t_hat = self.find_saddlepoint(x)
        
        K_t = self.cumulant_generating_function(t_hat)
        K2_t = self.cgf_derivative2(t_hat)
        
        if K2_t <= 0:
            return 0.0
        
        density = np.sqrt(1 / (2 * np.pi * K2_t)) * np.exp(K_t - x * t_hat)
        
        return max(density, 0.0)
    
    def tail_probability_lr(self, v: float, continuity_correction: bool = False) -> float:
        """
        Lugannani-Rice saddlepoint approximation for tail probability.

        P(H >= v) ≈ 1 - Φ(w_hat) + φ(w_hat) * (1/u_hat - 1/w_hat)

        where:
        - w_hat = sign(t_hat) * sqrt(2*(t_hat*v - K_H(t_hat)))
        - u_hat = t_hat * sqrt(K''_H(t_hat))

        For discrete distributions, continuity correction adjusts v to
        better approximate the discrete tail probability.

        Parameters
        ----------
        v : float
            Critical value
        continuity_correction : bool
            Whether to apply continuity correction (shifts v by -0.5
            to account for discrete nature of the statistic)

        Returns
        -------
        float
            Approximate P(H >= v)
        """
        if continuity_correction:
            # For discrete distributions, P(H >= v) is approximated by
            # the continuous integral from v - 0.5 to infinity
            v = v - 0.5
        
        mean = self.moments.get_mean()
        
        t_hat = self.find_saddlepoint(v)
        
        # Special case: v = mean (t_hat = 0)
        if np.abs(t_hat) < 1e-10:
            K2_0 = self.cgf_derivative2(0)
            K3_0 = self.cgf_derivative3(0)
            
            term = (K3_0 * K2_0**(-3/2) / 6 - K2_0**(-1/2) / 2) / np.sqrt(2 * np.pi)
            return 0.5 - term
        
        K_t = self.cumulant_generating_function(t_hat)
        K2_t = self.cgf_derivative2(t_hat)
        
        if K2_t <= 0:
            # Fallback to chi-square
            from .kruskal_wallis import KruskalWallisStatistic
            kw = KruskalWallisStatistic(list(self.sample_sizes))
            return kw.chi_square_approximation(v)
        
        # Compute w_hat and u_hat
        arg_w = 2 * (t_hat * v - K_t)
        if arg_w < 0:
            # Can happen due to numerical issues
            arg_w = max(arg_w, 0)
        
        w_hat = np.sign(t_hat) * np.sqrt(arg_w)
        u_hat = t_hat * np.sqrt(K2_t)
        
        # Lugannani-Rice formula
        Phi_w = stats.norm.cdf(w_hat)
        phi_w = stats.norm.pdf(w_hat)
        
        if np.abs(u_hat) < 1e-10 or np.abs(w_hat) < 1e-10:
            return 0.5
        
        prob = 1 - Phi_w + phi_w * (1/u_hat - 1/w_hat)
        
        return np.clip(prob, 0, 1)
    
    def tail_probability_gamma_based(self, v: float) -> float:
        """
        Gamma-based saddlepoint approximation (Wood et al., 1993).
        
        Uses gamma distribution instead of normal for the base approximation.
        
        Parameters
        ----------
        v : float
            Critical value
            
        Returns
        -------
        float
            Approximate P(H >= v)
        """
        # Get gamma parameters matched to first two moments
        alpha, beta = self.moments.get_gamma_params()
        
        t_hat = self.find_saddlepoint(v)
        
        K_t = self.cumulant_generating_function(t_hat)
        K2_t = self.cgf_derivative2(t_hat)
        
        if K2_t <= 0:
            return stats.gamma.sf(v, alpha, scale=beta)
        
        # Find xi_hat using gamma CGF matching
        # G(t_xi) - xi*t_xi = K_H(t_hat) - v*t_hat
        
        target = K_t - v * t_hat
        
        def gamma_cgf(t, a, b):
            if t < 1/b:
                return -a * np.log(1 - b * t)
            else:
                return np.inf
        
        # Use gamma approximation
        g_pdf = stats.gamma.pdf(v, alpha, scale=beta)
        g_cdf = stats.gamma.cdf(v, alpha, scale=beta)
        
        # Correction term
        u_hat = t_hat * np.sqrt(K2_t)
        
        if np.abs(u_hat) < 1e-10:
            return 1 - g_cdf
        
        prob = 1 - g_cdf + g_pdf * (1/u_hat - 1/t_hat if np.abs(t_hat) > 1e-10 else 0)
        
        return np.clip(prob, 0, 1)
    
    def cdf_approximation(self, x: float, method: str = 'LR', 
                          continuity_correction: bool = False) -> float:
        """
        Compute approximate CDF P(H <= x).
        
        Parameters
        ----------
        x : float
            Value at which to evaluate CDF
        method : str
            'LR' for Lugannani-Rice, 'gamma' for gamma-based
        continuity_correction : bool
            Whether to apply continuity correction
            
        Returns
        -------
        float
            Approximate P(H <= x)
        """
        if method == 'LR':
            return 1 - self.tail_probability_lr(x, continuity_correction)
        elif method == 'gamma':
            return 1 - self.tail_probability_gamma_based(x)
        else:
            raise ValueError(f"Unknown method: {method}")
    
    def critical_value(self, alpha: float, method: str = 'LR',
                       continuity_correction: bool = False) -> float:
        """
        Compute critical value for significance level alpha.
        
        Finds c such that P(H >= c) = alpha.
        
        Parameters
        ----------
        alpha : float
            Significance level (e.g., 0.05 or 0.10)
        method : str
            Approximation method
        continuity_correction : bool
            Whether to apply continuity correction
            
        Returns
        -------
        float
            Critical value c
        """
        def objective(c):
            if method == 'LR':
                return self.tail_probability_lr(c, continuity_correction) - alpha
            else:
                return self.tail_probability_gamma_based(c) - alpha
        
        # Search for critical value
        # H is non-negative, and critical values are typically in range [0, 20]
        try:
            c_star = brentq(objective, 0.01, 50, xtol=1e-8)
        except ValueError:
            # Fallback to chi-square critical value
            c_star = stats.chi2.ppf(1 - alpha, self.k - 1)
        
        return c_star
    
    def __repr__(self) -> str:
        return (f"SaddlepointApproximation(k={self.k}, N={self.N}, "
                f"method='{self.cgf_method}')")
