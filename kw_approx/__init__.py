"""
Kruskal-Wallis Higher Order Asymptotic Approximations Package

This package implements higher-order asymptotic approximation methods for the 
Kruskal-Wallis statistic as described in:

Murakami, H. and Ha, H.-T. "Higher Order Asymptotic Approximations of 
Kruskal-Wallis Statistics"

Methods implemented:
- Saddlepoint Approximations (SD1, SD2, with continuity corrections)
- Polynomially Adjusted Gamma Method (PAM)
- Gram-Charlier Series Approximation
- Exact combinatorial computation (for small samples)
"""

from .kruskal_wallis import KruskalWallisStatistic
from .moments import KWMoments
from .saddlepoint import SaddlepointApproximation
from .pam import PolynomialAdjustedGamma
from .gram_charlier import GramCharlierApproximation
from .exact import ExactDistribution
from .approximator import KWApproximator

__version__ = "0.1.0"
__author__ = "Based on Murakami & Ha (2024)"

__all__ = [
    "KruskalWallisStatistic",
    "KWMoments", 
    "SaddlepointApproximation",
    "PolynomialAdjustedGamma",
    "GramCharlierApproximation",
    "ExactDistribution",
    "KWApproximator",
]
