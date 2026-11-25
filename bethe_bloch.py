"""Bethe-Bloch mean energy loss in matter."""

import numpy as np
from numpy import power
from numpy.typing import ArrayLike

try:
    from .constants import M_E
except ImportError:
    from constants import M_E

# Bethe bloch constants
MEAN_EX_ENERGY_BE = 63.7e-6
K_FACTOR = 0.307  # MeV cm^2 / g for A = 1 g/mol
DENSITY_BERYLLIUM = 1.848

BERYLLIUM = {"name": "Be", "Z": 4, "A": 9.012, "rho": DENSITY_BERYLLIUM,
             "I": MEAN_EX_ENERGY_BE}
ALUMINUM = {"name": "Al", "Z": 13, "A": 26.98, "rho": 2.699, "I": 166.0e-6}


def t_max(beta: ArrayLike, gamma: ArrayLike, M: float) -> ArrayLike:
    """Largest energy transfer to an electron, in MeV.

    Args:
        beta: Velocity p/E of the particle.
        gamma: Lorentz factor E/m of the particle.
        M: Mass of the particle in MeV.
    """
    return ((2*M_E*power(beta*gamma, 2))
            / (1 + 2*gamma*M_E/M + power(M_E/M, 2)))


def dEdx(beta: ArrayLike, gamma: ArrayLike, M: float, Z: float, A: float,
         rho: float = DENSITY_BERYLLIUM,
         Imin: float = MEAN_EX_ENERGY_BE) -> ArrayLike:
    """Bethe-Bloch mean energy loss in MeV/cm (positive).

    Args:
        beta: Velocity p/E of the particle.
        gamma: Lorentz factor E/m of the particle.
        M: Mass of the particle in MeV.
        Z: Atomic number of the material.
        A: Atomic weight of the material in g/mol.
        rho: Density of the material in g/cm^3.
        Imin: Mean excitation energy of the material in MeV.
    """
    prefactor = rho*K_FACTOR*(Z/A)*(1/beta**2)
    log_argument = ((2*M_E*power(beta*gamma, 2)*t_max(beta, gamma, M))
                    / power(Imin, 2))
    return prefactor * (0.5*np.log(log_argument) - power(beta, 2))
