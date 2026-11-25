try:
    from .constants import *
except ImportError:
    from constants import *

import numpy as np
from numpy import log, exp, pi, sqrt, power, sin, cos, arccos, heaviside
from scipy.stats import norm, expon


def p_decay(p, m, tau, distance):
    # Probability that a particle will decay before reaching a distance l
    # momentum in lab frame p
    # lifetime tau in seconds
    # mass of decaying particle m
    # distance from source l in meters
    energy = sqrt(p**2 + m**2)
    boost = energy / m
    v = p / energy
    prob = exp(-distance/(METER_BY_MEV*v*boost*tau/HBAR))
    return (1 - prob)




def lab_lifetime(p, m, tau):
    return sqrt(p**2 + m**2) / m * tau




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




def proton_be_xs(p, kind="tot"):
    if kind == "prod":
        return proton_be_xs(p, "ine") - proton_be_xs(p, "qel")
    a, b, n, c, d = {"tot": (307.8, 0.897, 0.003, -2.598, -4.973),
                     "ine": (186.7, 104.3, -1.039, 10.38, -15.83),
                     "qel": (164.8, -40.09, 0.408, 21.40, -61.45)}[kind]
    return a + b*power(p, n) + c*log(p)**2 + d*log(p)




PION_XS_PARAMS = {
    ("Be", +1, "ine"): (0.400, 5.142, 162.3, -99.79, -2.407, -0.423, 850.3, 1.201, 0.375),
    ("Be", +1, "qel"): (0.635, 3.784, -2.38, -81.84, -2.702, 3.173, 379.9, 1.201, 0.558),
    ("Al", +1, "ine"): (0.295, 2.307, 1537.4, -1109.4, 0.057, 14.40, 510.7, 1.189, 0.185),
    ("Al", +1, "qel"): (0.698, 2.134, 40.38, 89.20, -1.575, 0.335, 229.4, 1.189, 0.187),
    ("Be", -1, "ine"): (0.600, 2.874, 92.66, 112.2, -0.486, 7.500, 371.5, 1.201, 0.233),
    ("Be", -1, "qel"): (0.626, 2.504, -1.559, 46.41, -0.633, 1.874, 189.0, 1.201, 0.185),
    ("Al", -1, "ine"): (0.706, 1.685, 997.8, -457.8, 0.611, 233.4, 446.8, 1.189, 0.305),
    ("Al", -1, "qel"): (0.633, 2.199, 32.52, 85.15, -1.225, 1.383, 129.1, 1.189, 0.305),
}

def pion_nucleus_xs(p, charge=1, material="Be", kind="ine"):
    p = np.maximum(p, 0.15)
    if kind == "prod":
        return np.maximum(pion_nucleus_xs(p, charge, material, "ine") - pion_nucleus_xs(p, charge, material, "qel"), 0.0)
    theta0, thetas, a, b, n, c, n_r, m_r, g_r = PION_XS_PARAMS[(material, int(np.sign(charge)), kind)]
    m2 = (M_PI**2 + M_P**2 + 2*M_P*sqrt((1.0e3*p)**2 + M_PI**2)) * 1.0e-6
    resonance = n_r * m2*g_r**2 / ((m_r**2 - m2)**2 + m2*g_r**2)
    return np.maximum(resonance + 0.5*(1 + np.tanh(thetas*(p - theta0))) * (a + b*power(p, n) + c*log(p)**2), 0.0)




# Sanford-Wang and Feynman scaling diff xs
def meson_production_d2SdpdOmega(p, theta, p_proton, meson_type="pi_plus", sw_only=False, params=None):
    p_beam = p_proton
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
        if params is not None:
            c1, c2, c3, c4, c5, c6, c7, c8, c9 = params
        #c1, c2, c3, c4, c5, c6, c7, c8, c9 = 1.20245,1.08, 2.15, 2.31,1.98,5.73,0.137,24.1, 1.0
        prefactor = c1 * power(p, c2) * np.clip(1 - p/(p_beam - c9), 0.0, None)
        exponential = exp(-c3*power(p,c4)/power(p_beam,c5) - c6*theta*(p-c7*p_beam*power(np.clip(cos(theta), 0.0, None),c8)))
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
        if params is not None:
            c1, c2, c3, c4, c5, c6, c7, c8, c9 = params
        prefactor = c1 * power(p, c2) * np.clip(1 - p/(p_beam - c9), 0.0, None)
        exponential = exp(-c3*power(p,c4)/power(p_beam,c5) - c6*theta*(p-c7*p_beam*power(np.clip(cos(theta), 0.0, None),c8)))
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
            if params is not None:
                c1, c2, c3, c4, c5, c6, c7, c8, c9 = params
            #c1, c2, c3, c4, c5, c6, c7, c8, c9 = 1.20245,1.08, 2.15, 2.31,1.98,5.73,0.137,24.1, 1.0
            prefactor = c1 * power(p, c2) * np.clip(1 - p/(p_beam - c9), 0.0, None)
            exponential = exp(-c3*power(p,c4)/power(p_beam,c5) - c6*theta*(p-c7*p_beam*power(np.clip(cos(theta), 0.0, None),c8)))
            return prefactor * exponential
    
        p_t = p*sin(theta)
        p_l = p*cos(theta)
        beta = p_beam / (mt*1e-3 + sqrt(p_beam**2 + (M_P*1e-3)**2))
        gamma = power(1-beta**2, -0.5)
        p_l_star = gamma*(p_l - sqrt(p_l**2 + p_t**2 + (M_K*1e-3)**2)*beta)
        s = (M_P*1e-3)**2 + (mt*1e-3)**2 + 2*sqrt(p_beam**2 + (M_P*1e-3)**2)*mt*1e-3
        m_recoil, m_kaon = 2.053, M_K*1e-3
        p_l_star_max = sqrt((s - (m_kaon + m_recoil)**2) * (s - (m_kaon - m_recoil)**2)) / (2*sqrt(s))
        x_f = abs(p_l_star) / p_l_star_max

        c1 = 11.70
        c2 = 0.88
        c3 = 4.77
        c4 = 1.51
        c5 = 2.21
        c6 = 2.17
        c7 = 1.51
        if params is not None:
            c1, c2, c3, c4, c5, c6, c7 = params
        prefactor = c1 * p**2 / sqrt(p**2 + (M_K*1e-3)**2)
        return prefactor * (x_f <= 1.0) * exp(-c2*p_t - c3*power(x_f, c4) - c5*p_t**2 - c7*power(p_t*x_f, c6))
    
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
        if params is not None:
            c1, c2, c3, c4, c5, c6, c7, c8, c9 = params
        prefactor = c1 * power(p, c2) * np.clip(1 - p/(p_beam - c9), 0.0, None)
        exponential = exp(-c3*power(p,c4)/power(p_beam,c5) - c6*theta*(p-c7*p_beam*power(np.clip(cos(theta), 0.0, None),c8)))
        return prefactor * exponential
