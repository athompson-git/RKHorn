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
             "I": MEAN_EX_ENERGY_BE, "a": 0.8039, "k": 2.4339, "x0": 0.0592,
             "x1": 1.6922, "cbar": 2.7847, "delta0": 0.14}
ALUMINUM = {"name": "Al", "Z": 13, "A": 26.98, "rho": 2.699, "I": 166.0e-6,
            "a": 0.0802, "k": 3.6345, "x0": 0.1708, "x1": 3.0127,
            "cbar": 4.2395, "delta0": 0.12}


def t_max(beta: ArrayLike, gamma: ArrayLike, M: float) -> ArrayLike:
    """Largest energy transfer to an electron, in MeV.

    Args:
        beta: Velocity p/E of the particle.
        gamma: Lorentz factor E/m of the particle.
        M: Mass of the particle in MeV.
    """
    return ((2*M_E*power(beta*gamma, 2))
            / (1 + 2*gamma*M_E/M + power(M_E/M, 2)))


def density_effect(beta_gamma: ArrayLike, material: dict) -> ArrayLike:
    """Sternheimer density-effect correction delta.

    Args:
        beta_gamma: Momentum over mass, p/m, of the particle.
        material: Material dictionary with Sternheimer parameters.
    """
    x = np.log10(beta_gamma)
    high = 2*np.log(10.0)*x - material["cbar"]
    middle = high + material["a"]*power(
        np.clip(material["x1"] - x, 0.0, None), material["k"])
    low = material["delta0"]*power(10.0, 2*(x - material["x0"]))
    return np.where(x >= material["x1"], high,
                    np.where(x >= material["x0"], middle, low))


def dEdx(beta: ArrayLike, gamma: ArrayLike, M: float, Z: float, A: float,
         rho: float = DENSITY_BERYLLIUM,
         Imin: float = MEAN_EX_ENERGY_BE,
         delta: ArrayLike = 0.0) -> ArrayLike:
    """Bethe-Bloch mean energy loss in MeV/cm (positive).

    Args:
        beta: Velocity p/E of the particle.
        gamma: Lorentz factor E/m of the particle.
        M: Mass of the particle in MeV.
        Z: Atomic number of the material.
        A: Atomic weight of the material in g/mol.
        rho: Density of the material in g/cm^3.
        Imin: Mean excitation energy of the material in MeV.
        delta: Density-effect correction, from density_effect.
    """
    prefactor = rho*K_FACTOR*(Z/A)*(1/beta**2)
    log_argument = ((2*M_E*power(beta*gamma, 2)*t_max(beta, gamma, M))
                    / power(Imin, 2))
    return prefactor * (0.5*np.log(log_argument) - power(beta, 2)
                        - 0.5*delta)
