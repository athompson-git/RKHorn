"""
Bethe bloch formula
"""

from constants import *

import numpy as np
from numpy import pi, sqrt, power, sin, cos, arccos, arctan2


# Bethe bloch constants
MEAN_EX_ENERGY_BE = 48.0e-6  # 48 eV for Beryllium
K_FACTOR = 0.307  # MeV cm^2 / g for A = 1 g/mol
DENSITY_BERYLLIUM = 1.848

def Tmax(beta, gamma, M):
    # M: mass of propagating particle
    return (2*M_E*power(beta*gamma,2))/(1 + 2*gamma*M_E/M + power(M_E/M, 2))


def dEdx(beta, gamma, M, Z, A, rho=DENSITY_BERYLLIUM, Imin=MEAN_EX_ENERGY_BE):
    # M: mass of propagating particle
    # gamma: gamma factor E/m
    # beta: velocity p/E
    # Z: atomic number of absorbing material
    # A: atomic weight in g/mol
    # rho: density of material in g/cm^3
    # Imin: mean excitation energy of material
    # returns energy loss in MeV/cm

    prefactor = rho*K_FACTOR*(Z/A)*(1/beta**2)

    return prefactor * (0.5*np.log((2*M_E*power(beta*gamma, 2)*Tmax(beta, gamma, M))/power(Imin, 2)) \
                        - power(beta, 2))
