from constants import *

import numpy as np
from numpy import log, exp, pi, sqrt, power, sin, cos, arccos, heaviside
from scipy.stats import norm, expon


def p_decay(p, m, tau, l):
    # Probability that a particle will decay before reaching a distance l
    # momentum in lab frame p
    # lifetime tau in seconds
    # mass of decaying particle m
    # distance from source l in meters
    energy = sqrt(p**2 + m**2)
    boost = energy / m
    v = p / energy
    prob = exp(-l/(METER_BY_MEV*v*boost*tau/HBAR))
    return (1 - prob)




def pion_lifetime(p):
    energy = sqrt(p**2 + M_PI**2)
    boost = energy / M_PI
    v = p / energy
    return v*boost*PION_LIFETIME




def kaon_lifetime(p):
    energy = sqrt(p**2 + M_K**2)
    boost = energy / M_K
    v = p / energy
    return v*boost*KAON_LIFETIME




# Sanford-Wang and Feynman scaling diff xs
def meson_production_d2SdpdOmega(p, theta, p_proton, meson_type="pi_plus", sw_only=False):
    pB = p_proton
    mt = M_P
    # Sanford-Wang Parameterization
    if meson_type == "pi_plus":
        c1 = 220.7
        c2 = 1.080
        c3 = 1.0
        c4 = 1.978
        c5 = 1.32
        c6 = 5.572
        c7 = 0.0868
        c8 = 9.686
        c9 = 1.0
        #c1, c2, c3, c4, c5, c6, c7, c8, c9 = 1.20245,1.08, 2.15, 2.31,1.98,5.73,0.137,24.1, 1.0
        prefactor = c1 * power(p, c2) * (1 - p/(pB - c9))
        exponential = exp(-c3*power(p,c4)/power(pB,c5) - c6*theta*(p-c7*pB*power(cos(theta),c8)))
        return prefactor * exponential
    
    elif meson_type == "pi_minus":
        c1 = 213.7
        c2 = 0.9379
        c3 = 5.454
        c4 = 1.210
        c5 = 1.284
        c6 = 4.781
        c7 = 0.07338
        c8 = 8.329
        c9 = 1.0
        prefactor = c1 * power(p, c2) * (1 - p/(pB - c9))
        exponential = exp(-c3*power(p,c4)/power(pB,c5) - c6*theta*(p-c7*pB*power(cos(theta),c8)))
        return prefactor * exponential
    
    elif meson_type == "k_plus":
        if sw_only:
            # modified SW from 1110.0417
            c1 = 14.89
            c2 = 0.91
            c3 = 12.80
            c4 = 2.08
            c5 = 2.65
            c6 = 4.61
            c7 = 0.26
            c8 = 10.63
            c9 = 2.04
            #c1, c2, c3, c4, c5, c6, c7, c8, c9 = 1.20245,1.08, 2.15, 2.31,1.98,5.73,0.137,24.1, 1.0
            prefactor = c1 * power(p, c2) * (1 - p/(pB - c9))
            exponential = exp(-c3*power(p,c4)/power(pB,c5) - c6*theta*(p-c7*pB*power(cos(theta),c8)))
            return prefactor * exponential
    
        pT = p*sin(theta)
        pL = p*cos(theta)
        beta = pB / (mt*1e-3 + sqrt(pB**2 + (M_P*1e-3)**2))
        gamma = power(1-beta**2, -0.5)
        pLstar = gamma*(pL - sqrt(pL**2 + pT**2 + (M_K*1e-3)**2)*beta)
        s = (M_P*1e-3)**2 + (mt*1e-3)**2 + 2*sqrt(pB**2 + (M_P*1e-3)**2)*mt*1e-3
        xF = abs(2*pLstar/sqrt(s))

        c1 = 11.70
        c2 = 0.88
        c3 = 4.77
        c4 = 1.51
        c5 = 2.21
        c6 = 2.17
        c7 = 1.51
        prefactor = c1 * p**2 / sqrt(p**2 + (M_K*1e-3)**2)
        return prefactor * (1 - xF) * exp(-c2*pT - c3*power(xF, c4) - c5*pT**2 - c7*power(pT*xF, c6))
    
    elif meson_type == "k0S":
        c1 = 15.130
        c2 = 1.975
        c3 = 4.084
        c4 = 0.928
        c5 = 0.731
        c6 = 4.362
        c7 = 0.048
        c8 = 13.300
        c9 = 1.278
        prefactor = c1 * power(p, c2) * (1 - p/(pB - c9))
        exponential = exp(-c3*power(p,c4)/power(pB,c5) - c6*theta*(p-c7*pB*power(cos(theta),c8)))
        return prefactor * exponential
