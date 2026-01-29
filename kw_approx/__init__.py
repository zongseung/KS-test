"""
Kruskal-Wallis Higher Order Asymptotic Approximations Package

This package implements higher-order asymptotic approximation methods for the
Kruskal-Wallis statistic as described in:

Murakami, H., Lee, J.-S. and Ha, H.-T. "Higher Order Asymptotic Approximations
of Kruskal-Wallis Statistics Based on Skewness and Kurtosis"

Methods implemented:
- Saddlepoint Approximations (SD1, SD2, with continuity corrections SDC1, SDC2)
- Polynomially Adjusted Gamma Method (PAM/PAG)
- Gram-Charlier Series Approximation (GC-A)
- Edgeworth Expansion (ED)
- Exact combinatorial computation (for small samples)

References:
- Kruskal, W.H. & Wallis, A. (1952). Use of ranks in one-criterion variance analysis.
- Daniels, H.E. (1954). Saddlepoint approximations in statistics.
- Lugannani, R. & Rice, S.O. (1980). Saddlepoint approximation for distributions.
- Ha, H-T. & Provost, S.B. (2007). A viable alternative to statistical tables.
- Iman, R.L., Quade, D. & Alexander, D.A. (1975). Exact probability levels.
"""

from .kruskal_wallis import KruskalWallisStatistic
from .moments import KWMoments
from .saddlepoint import SaddlepointApproximation
from .pam import PolynomialAdjustedGamma
from .gram_charlier import GramCharlierApproximation
from .edgeworth import EdgeworthApproximation
from .exact import ExactDistribution
from .simulation import MonteCarloSimulation
from .approximator import KWApproximator

__version__ = "0.2.0"
__author__ = "Based on Murakami, Lee & Ha (2026)"

__all__ = [
    "KruskalWallisStatistic",
    "KWMoments",
    "SaddlepointApproximation",
    "PolynomialAdjustedGamma",
    "GramCharlierApproximation",
    "EdgeworthApproximation",
    "ExactDistribution",
    "MonteCarloSimulation",
    "KWApproximator",
]
