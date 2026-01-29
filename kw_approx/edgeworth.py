"""
Edgeworth Expansion Approximation for Kruskal-Wallis Statistic

Implements the Edgeworth series expansion for approximating the null distribution
of the Kruskal-Wallis H statistic, as described in Section 3.3 of:

"Higher Order Asymptotic Approximations of Kruskal-Wallis Statistics"
by Murakami, Lee, and Ha

The Edgeworth expansion provides corrections to the chi-square approximation
based on higher-order cumulants (skewness and kurtosis).

For a statistic that converges to chi-square(k-1), the expansion uses
Laguerre polynomials rather than Hermite polynomials.

References:
- Murakami, Lee & Ha: Higher Order Asymptotic Approximations of Kruskal-Wallis Statistics
- Berry-Esseen bounds for rate of convergence
"""

import numpy as np
from typing import List, Tuple, Optional
from scipy import stats
from scipy.special import gamma as gamma_func
from .moments import KWMoments


class EdgeworthApproximation:
    """
    Edgeworth expansion approximation for the Kruskal-Wallis H statistic.

    The Edgeworth expansion corrects the limiting chi-square(k-1) distribution
    using higher-order cumulants. For chi-square type statistics, this involves
    Laguerre polynomials.

    For standardized statistic T_N with CDF F_N(x), mean 0, variance 1:

        F_N(x) ≈ Φ(x) - (λ₃/(6√N))(x² - 1)φ(x) - (λ₄/(24N))(x³ - 3x)φ(x)
                 - (λ₃²/(72N))(x⁵ - 10x³ + 15x)φ(x)

    where λ₃ = E[T³] and λ₄ = E[T⁴] - 3 are standardized cumulants.

    Parameters
    ----------
    sample_sizes : List[int]
        Sample sizes for each group [n1, n2, ..., nk]

    Attributes
    ----------
    mu : float
        Mean of H under H0 (= k - 1)
    sigma : float
        Standard deviation of H under H0
    lambda3 : float
        Standardized third cumulant (skewness)
    lambda4 : float
        Standardized fourth cumulant (excess kurtosis)
    """

    def __init__(self, sample_sizes: List[int]):
        self.sample_sizes = np.array(sample_sizes, dtype=int)
        self.k = len(sample_sizes)
        self.N = np.sum(self.sample_sizes)
        self.df = self.k - 1  # degrees of freedom

        # Compute moments
        self.moments = KWMoments(sample_sizes, max_moment=6)

        # Extract key statistics
        self.mu = self.moments.get_mean()
        self.sigma = self.moments.get_std()
        self.lambda3 = self.moments.get_skewness()   # Standardized skewness
        self.lambda4 = self.moments.get_kurtosis()   # Excess kurtosis

        # For chi-square reference
        self._chi2_skewness = np.sqrt(8.0 / self.df) if self.df > 0 else 0.0
        self._chi2_kurtosis = 12.0 / self.df if self.df > 0 else 0.0

    def standardize(self, h: float) -> float:
        """
        Standardize h to z = (h - μ) / σ.

        Parameters
        ----------
        h : float
            Original H statistic value

        Returns
        -------
        float
            Standardized value z
        """
        if self.sigma <= 0:
            return 0.0
        return (h - self.mu) / self.sigma

    def pdf(self, h: float) -> float:
        """
        Edgeworth expansion approximation to the PDF of H.

        f_ED(h) = (1/σ) * φ(z) * [1 + correction terms]

        Parameters
        ----------
        h : float
            Point at which to evaluate density

        Returns
        -------
        float
            Approximate density value
        """
        if h <= 0:
            return 0.0

        z = self.standardize(h)

        # Standard normal PDF
        phi_z = stats.norm.pdf(z)

        # Edgeworth correction terms
        # Using Hermite polynomials for normal-based expansion
        He3 = z**3 - 3*z           # H_3(z)
        He4 = z**4 - 6*z**2 + 3    # H_4(z)
        He6 = z**6 - 15*z**4 + 45*z**2 - 15  # H_6(z)

        # Correction factor
        # Note: For density, we differentiate the CDF expansion
        correction = (1 +
                     self.lambda3 / 6 * He3 +
                     self.lambda4 / 24 * He4 +
                     self.lambda3**2 / 72 * He6)

        density = phi_z * correction / self.sigma

        return max(density, 0.0)

    def cdf(self, h: float) -> float:
        """
        Edgeworth expansion approximation to the CDF of H.

        F_ED(h) ≈ Φ(z) - φ(z) * [correction terms]

        From the paper (Section 3.3):
        F_N(x) ≈ Φ(x) - (λ₃/6√N)(x² - 1)φ(x) - (λ₄/24N)(x³ - 3x)φ(x)
                 - (λ₃²/72N)(x⁵ - 10x³ + 15x)φ(x)

        Parameters
        ----------
        h : float
            Point at which to evaluate CDF

        Returns
        -------
        float
            Approximate P(H <= h)
        """
        if h <= 0:
            return 0.0

        z = self.standardize(h)

        # Standard normal CDF and PDF
        Phi_z = stats.norm.cdf(z)
        phi_z = stats.norm.pdf(z)

        # Hermite polynomials (one degree lower for CDF)
        He2 = z**2 - 1             # H_2(z)
        He3 = z**3 - 3*z           # H_3(z)
        He5 = z**5 - 10*z**3 + 15*z  # H_5(z)

        # Edgeworth correction
        # Note: The correction involves 1/√N factors which are absorbed into λ
        correction = phi_z * (
            self.lambda3 / 6 * He2 +
            self.lambda4 / 24 * He3 +
            self.lambda3**2 / 72 * He5
        )

        cdf_value = Phi_z - correction

        return np.clip(cdf_value, 0.0, 1.0)

    def cdf_chi2_based(self, h: float) -> float:
        """
        Chi-square based Edgeworth expansion (alternative formulation).

        This uses the chi-square distribution as the base instead of normal,
        which is more natural for the Kruskal-Wallis statistic.

        The expansion uses Laguerre polynomials instead of Hermite polynomials.

        Parameters
        ----------
        h : float
            Point at which to evaluate CDF

        Returns
        -------
        float
            Approximate P(H <= h)
        """
        if h <= 0:
            return 0.0

        # Base chi-square CDF
        chi2_cdf = stats.chi2.cdf(h, self.df)
        chi2_pdf = stats.chi2.pdf(h, self.df)

        # Compute deviation from chi-square cumulants
        delta_skew = self.lambda3 - self._chi2_skewness
        delta_kurt = self.lambda4 - self._chi2_kurtosis

        # Laguerre-based correction (simplified)
        # L_0(x) = 1, L_1(x) = 1-x, L_2(x) = (x^2 - 4x + 2)/2, etc.
        x = h / 2  # Scale for chi-square

        L1 = 1 - x
        L2 = (x**2 - 4*x + 2) / 2
        L3 = (-x**3 + 9*x**2 - 18*x + 6) / 6

        # Correction based on cumulant differences
        correction = chi2_pdf * (
            delta_skew / 6 * L2 +
            delta_kurt / 24 * L3
        )

        cdf_value = chi2_cdf - correction

        return np.clip(cdf_value, 0.0, 1.0)

    def sf(self, h: float) -> float:
        """
        Survival function (tail probability) P(H >= h).

        Parameters
        ----------
        h : float
            Critical value

        Returns
        -------
        float
            Approximate P(H >= h)
        """
        return 1.0 - self.cdf(h)

    def tail_probability(self, h: float) -> float:
        """
        Alias for sf(). Compute P(H >= h).

        Parameters
        ----------
        h : float
            Critical value

        Returns
        -------
        float
            Approximate P(H >= h)
        """
        return self.sf(h)

    def critical_value(self, alpha: float) -> float:
        """
        Compute critical value for significance level alpha.

        Uses numerical root finding to solve P(H >= c) = alpha.

        Parameters
        ----------
        alpha : float
            Significance level (e.g., 0.05 or 0.10)

        Returns
        -------
        float
            Critical value c
        """
        from scipy.optimize import brentq

        def objective(c):
            return self.sf(c) - alpha

        try:
            # Use chi-square quantile as starting point
            c_chi2 = stats.chi2.ppf(1 - alpha, self.df)
            c_star = brentq(objective, max(0.01, c_chi2 - 5), c_chi2 + 10, xtol=1e-8)
        except ValueError:
            # Fallback to chi-square
            c_star = stats.chi2.ppf(1 - alpha, self.df)

        return c_star

    def ppf(self, q: float) -> float:
        """
        Percent point function (inverse CDF).

        Parameters
        ----------
        q : float
            Probability (between 0 and 1)

        Returns
        -------
        float
            Value h such that P(H <= h) = q
        """
        return self.critical_value(1 - q)

    def error_bound(self) -> float:
        """
        Estimate Berry-Esseen type error bound for the approximation.

        The classical Berry-Esseen bound is O(N^{-1/2}).
        With Edgeworth corrections, the error is O(N^{-1}).

        Returns
        -------
        float
            Estimated error bound
        """
        # Approximate error bound based on sample size
        # For Edgeworth, error is O(N^{-1}) in smooth case
        # but can be O(N^{-1/2}) for lattice distributions

        return 1.0 / self.N  # Optimistic bound

    def get_parameters(self) -> dict:
        """
        Get all approximation parameters.

        Returns
        -------
        dict
            Dictionary containing mu, sigma, lambda3, lambda4, etc.
        """
        return {
            'mu': self.mu,
            'sigma': self.sigma,
            'lambda3': self.lambda3,
            'lambda4': self.lambda4,
            'chi2_skewness': self._chi2_skewness,
            'chi2_kurtosis': self._chi2_kurtosis,
            'df': self.df,
            'N': self.N
        }

    def is_stable(self) -> bool:
        """
        Check if the Edgeworth approximation is likely stable.

        The approximation can produce oscillations (sawtooth error)
        for lattice-valued random variables and may give negative
        densities when skewness/kurtosis are extreme.

        Returns
        -------
        bool
            True if approximation is likely stable
        """
        # Approximation tends to be unstable for extreme skewness/kurtosis
        # or very small sample sizes
        return (abs(self.lambda3) < 2.0 and
                abs(self.lambda4) < 6.0 and
                self.N >= 9)

    def __repr__(self) -> str:
        return (f"EdgeworthApproximation(k={self.k}, N={self.N}, "
                f"df={self.df}, mu={self.mu:.4f}, sigma={self.sigma:.4f}, "
                f"lambda3={self.lambda3:.4f}, lambda4={self.lambda4:.4f})")
