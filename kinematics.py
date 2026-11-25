"""
Lorentz vector and vector classes, monte carlo methods for decays, boost method
"""

import numpy as np
from numpy import pi, sqrt, power, sin, cos, arccos, arctan2, log, arctan

try:
    from .constants import *
except ImportError:
    from constants import *




class Vector3:
    """
    Vector class for 3-dimensional vectors
    """
    def __init__(self, v1, v2, v3):
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        self.vec = np.array([v1, v2, v3])
    
    def __str__(self):
        return "({0},{1},{2})".format(self.v1, self.v2, self.v3)
    
    def __add__(self, other):
        v1_new = self.v1 + other.v1
        v2_new = self.v2 + other.v2
        v3_new = self.v3 + other.v3
        return Vector3(v1_new, v2_new, v3_new)
    
    def __sub__(self, other):
        v1_new = self.v1 - other.v1
        v2_new = self.v2 - other.v2
        v3_new = self.v3 - other.v3
        return Vector3(v1_new, v2_new, v3_new)
    
    def __neg__(self):
        return Vector3(-self.v1, -self.v2, -self.v3)
    
    def __mul__(self, other):
        if isinstance(other, (int, float)):  # Scalar multiplication
            return Vector3(self.v1 * other, self.v2 * other, self.v3 * other)
        elif isinstance(other, Vector3):  # Dot product
            return np.dot(self.vec, other.vec)
        else:
            raise TypeError("Unsupported operand type for *: 'Vector' and '{}'".format(type(other)))
    
    def __rmul__(self, other):
        return self.__mul__(other)  # Reuse __mul__
    
    def unit_vec(self, eps=1e-40):
        v = self.mag()
        if v < eps:
            return Vector3(0.0, 0.0, 0.0)
        return Vector3(self.v1/v, self.v2/v, self.v3/v)
    
    def mag2(self):
        return np.dot(self.vec, self.vec)
    
    def mag(self):
        return np.sqrt(np.dot(self.vec, self.vec))

    def phi(self):
        return arctan2(self.v2, self.v1)
    
    def theta(self):
        return arccos(self.v3 / self.mag())
    
    def set_v3(self, v1, v2, v3):
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        self.vec = np.array([v1, v2, v3])




class LorentzVector:
    """
    Lorentz vector class for 4-dimensional spacetime vectors with metric (1, -1, -1, -1)
    """
    def __init__(self, p0=0.0, p1=0.0, p2=0.0, p3=0.0):
        self.p0 = p0
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.pmu = np.array([p0, p1, p2, p3])
        self.mt = np.array([1, -1, -1, -1])
        self.momentum3 = Vector3(self.p1, self.p2, self.p3)
    
    def __str__(self):
        return "({0},{1},{2},{3})".format(self.p0, self.p1, self.p2, self.p3)
    
    def __add__(self, other):
        p0_new = self.p0 + other.p0
        p1_new = self.p1 + other.p1
        p2_new = self.p2 + other.p2
        p3_new = self.p3 + other.p3
        return LorentzVector(p0_new, p1_new, p2_new, p3_new)
    
    def __sub__(self, other):
        p0_new = self.p0 - other.p0
        p1_new = self.p1 - other.p1
        p2_new = self.p2 - other.p2
        p3_new = self.p3 - other.p3
        return LorentzVector(p0_new, p1_new, p2_new, p3_new)
    
    def __mul__(self, other):
        return np.dot(self.pmu*other.pmu, self.mt)
    
    def __rmul__(self, other):
        return np.dot(self.pmu*other.pmu, self.mt)
    
    def mass2(self):
        return np.dot(self.pmu**2, self.mt)
    
    def mass(self):
        return np.sqrt(np.dot(self.pmu**2, self.mt))
    
    def energy(self):
        return self.p0
    
    def cosine(self):
        return self.p3 / self.momentum()
    
    def phi(self):
        return arctan2(self.p2, self.p1)
    
    def theta(self):
        return arccos(self.p3 / self.momentum())
    
    def momentum(self):
        return self.momentum3.mag()
    
    def set_p4(self, p0, p1, p2, p3):
        self.p0 = p0
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.pmu = np.array([p0, p1, p2, p3])
        self.momentum3 = Vector3(self.p1, self.p2, self.p3)
    
    def get_3momentum(self):
        return Vector3(self.p1, self.p2, self.p3)
    
    def get_3velocity(self):
        return Vector3(self.p1/self.p0, self.p2/self.p0, self.p3/self.p0)



    
class Decay2Body:
    """
    Performs a weighted MC sampling of the differential 2-body decay width
    by choosing a final-state particle of interest (particle #1 by convention)
    and drawing random variates in the rest-frame of the parent; particle 1
    has a random angle on the 2-sphere and a random energy between some E_max and E_min.
    The weight is given by dGamma / dE_1 dOmega
    Finally, we boost to the lab frame with the appropriate Jacobian factor.
    """
    def __init__(self, p_parent: LorentzVector, m1, m2, n_samples=1000):
        self.mp = p_parent.mass()  # parent particle
        self.m1 = m1  # decay body 1
        self.m2 = m2  # decay body 2

        self.lv_p = p_parent

        self.n_samples = n_samples
        self.p1_cm_4vectors = []
        self.p1_lab_4vectors = []
        self.p2_cm_4vectors = []
        self.p2_lab_4vectors = []
        self.weights = np.array([])
    
    def set_new_decay(self, p_parent: LorentzVector, m1, m2):
        self.mp = p_parent.mass()  # parent particle
        self.m1 = m1  # decay body 1
        self.m2 = m2  # decay body 2

        self.lv_p = p_parent

        self.p1_cm_4vectors = []
        self.p1_lab_4vectors = []
        self.p2_cm_4vectors = []
        self.p2_lab_4vectors = []
        self.weights = np.array([])
    
    def decay(self):
        p_cm = power((self.mp**2 - (self.m2 - self.m1)**2)*(self.mp**2 - (self.m2 + self.m1)**2), 0.5)/(2*self.mp)
        e1_cm = sqrt(p_cm**2 + self.m1**2)
        e2_cm = sqrt(p_cm**2 + self.m2**2)

        # Draw random variates on the 2-sphere
        phi1_rnd = 2*pi*np.random.ranf(self.n_samples)
        theta1_rnd = arccos(1 - 2*np.random.ranf(self.n_samples))

        v_in = -self.lv_p.get_3velocity()

        self.p1_cm_4vectors = [LorentzVector(e1_cm,
                            p_cm*cos(phi1_rnd[i])*sin(theta1_rnd[i]),
                            p_cm*sin(phi1_rnd[i])*sin(theta1_rnd[i]),
                            p_cm*cos(theta1_rnd[i])) for i in range(self.n_samples)]
        self.p2_cm_4vectors = [LorentzVector(e2_cm,
                            -p_cm*cos(phi1_rnd[i])*sin(theta1_rnd[i]),
                            -p_cm*sin(phi1_rnd[i])*sin(theta1_rnd[i]),
                            -p_cm*cos(theta1_rnd[i])) for i in range(self.n_samples)]
        self.p1_lab_4vectors = [lorentz_boost(p1, v_in) for p1 in self.p1_cm_4vectors]
        self.p2_lab_4vectors = [lorentz_boost(p2, v_in) for p2 in self.p2_cm_4vectors]
        self.weights = np.ones(self.n_samples) / self.n_samples




def lorentz_boost(momentum: LorentzVector, v: Vector3):
    """
    Lorentz boost momentum to a new frame with velocity v
    :param momentum: four vector
    :param v: velocity of new frame, 3-dimention
    :return: boosted momentum
    """
    n = v.unit_vec().vec
    beta = v.mag()
    if beta == 0.0:
        return momentum
    gamma = 1/np.sqrt(1-beta**2)
    mat = np.array([[gamma, -gamma*beta*n[0], -gamma*beta*n[1], -gamma*beta*n[2]],
                    [-gamma*beta*n[0], 1+(gamma-1)*n[0]*n[0], (gamma-1)*n[0]*n[1], (gamma-1)*n[0]*n[2]],
                    [-gamma*beta*n[1], (gamma-1)*n[1]*n[0], 1+(gamma-1)*n[1]*n[1], (gamma-1)*n[1]*n[2]],
                    [-gamma*beta*n[2], (gamma-1)*n[2]*n[0], (gamma-1)*n[2]*n[1], 1+(gamma-1)*n[2]*n[2]]])
    boosted_p4 = mat @ momentum.pmu
    return LorentzVector(boosted_p4[0], boosted_p4[1], boosted_p4[2], boosted_p4[3])