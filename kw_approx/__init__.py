"""
Kruskal-Wallis Higher Order Asymptotic Approximations Package

This package implements higher-order saddlepoint approximation methods for the
Kruskal-Wallis statistic as described in:

Murakami, H., Lee, J.-S. and Ha, H.-T. "Higher Order Asymptotic Approximations
of Kruskal-Wallis Statistics Based on Skewness and Kurtosis"

CGF methods: ER1, ER2, Wang, K-T (all with Lugannani-Rice tail probability)
All methods use exact finite-sample cumulants (Section 4 of the paper).

References:
- Kruskal, W.H. & Wallis, A. (1952). Use of ranks in one-criterion variance analysis.
- Daniels, H.E. (1954). Saddlepoint approximations in statistics.
- Lugannani, R. & Rice, S.O. (1980). Saddlepoint approximation for distributions.
- Iman, R.L., Quade, D. & Alexander, D.A. (1975). Exact probability levels.
"""

from .kruskal_wallis import KruskalWallisStatistic
from .moments import KWMoments
from .saddlepoint import SaddlepointApproximation
from .exact import ExactDistribution
from .simulation import MonteCarloSimulation
from .approximator import KWApproximator
from .pam import PolynomialAdjustedGamma
from .gram_charlier import GramCharlierApproximation
from .edgeworth import EdgeworthApproximation

__version__ = "0.3.0"
__author__ = "Based on Murakami, Lee & Ha (2026)"

__all__ = [
    "KruskalWallisStatistic",
    "KWMoments",
    "SaddlepointApproximation",
    "ExactDistribution",
    "MonteCarloSimulation",
    "KWApproximator",
    "PolynomialAdjustedGamma",
    "GramCharlierApproximation",
    "EdgeworthApproximation",
]
