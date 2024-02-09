# Runge-Kutta method for time evolving a charged particle through a magnetic field

from constants import *
from meson_physics import *

import numpy as np
from numpy import log, exp, pi, sqrt, power, sin, cos, arccos, heaviside
from scipy.stats import norm, expon


import matplotlib.pyplot as plt

# https://www.geeksforgeeks.org/runge-kutta-4th-order-method-solve-differential-equation/
# https://math.stackexchange.com/questions/14942/runge-kutta-method-for-newton-law
# https://primer-computational-mathematics.github.io/book/c_mathematics/numerical_methods/5_Runge_Kutta_method.html




def RungeKuttaCoupled(t, x, y, z, dt, dxdt, dydt, dzdt):
    k1 = dt*dydt(x, y, z)
    h1 = dt*dzdt(x, y, z)
    i1 = dt*dxdt(x, y, z)

    k2 = dt*dydt(x+i1/2., y+k1/2., z+h1/2.)
    h2 = dt*dzdt(x+i1/2., y+k1/2., z+h1/2.)
    i2 = dt*dxdt(x+i1/2., y+k1/2., z+h1/2.)

    k3 = dt*dydt(x+i2/2., y+k2/2., z+h2/2.)
    h3 = dt*dzdt(x+i2/2., y+k2/2., z+h2/2.)
    i3 = dt*dxdt(x+i2/2., y+k2/2., z+h2/2.)

    k4 = dt*dydt(x+i3, y+k3, z+h3)
    h4 = dt*dzdt(x+i3, y+k3, z+h3)
    i4 = dt*dxdt(x+i3, y+k3, z+h3)

    x = x + 1./6.*(i1 + 2*i2 + 2*i3 + i4)
    y = y + 1./6.*(k1 + 2*k2 + 2*k3 + k4)
    z = z + 1./6.*(h1 + 2*h2 + 2*h3 + h4)
    t = t + dt
    
    return t, x, y, z




class ChargedPionFluxMiniBooNE:
    def __init__(self, proton_energy=8000.0, meson_charge=1.0, solid_angle_cut=0.00924,
                 n_samples=10000, n_pot=18.75e20):
        self.n_samples = n_samples
        self.solid_angle_cut = solid_angle_cut
        self.charge = meson_charge
        self.n_pot = n_pot

        if meson_charge == 1.0:
            self.meson_type = "pi_plus"
        else:
            self.meson_type = "pi_minus"
        
        self.ep = proton_energy
        self.p_proton = np.sqrt(proton_energy**2 - M_P**2)

        self.sigma_times_n = 1e2 * self.sigmap(self.p_proton) * 1e-27 * 6.022e23 * 1.85 / 9.0  # sigma (mb->cm2) * number density cm^-3 -->meters

        self.x0 = np.array([])
        self.y0 = np.array([])
        self.z0 = np.array([])
        self.px0 = np.array([])
        self.py0 = np.array([])
        self.pz0 = np.array([])

        self.pip_p_post_horn = []
        self.pip_theta_post_horn = []
        self.pip_wgt_post_horn = []
        self.acceptance_wgt = []

        self.rksolver = RKHorn(horn_current=170.0, particle_mass=M_PI, charge=self.charge, step=0.1)
        self.meson_flux_pre_horn = None

    def sigmap(self, p):
        # returns total pion production cross section in mb
        A = 307.8
        B = 0.897
        C = -2.598
        D = -4.973
        n = 0.003
        return A + B*power(p/1.0e3,n) + C*log(p/1.0e3)*log(p/1.0e3) + D*log(p/1.0e3)

    def simulate_beam_spot(self):
        r1 = norm.rvs(size=self.n_samples)
        r2 = norm.rvs(size=self.n_samples)
        r3 = norm.rvs(size=self.n_samples)
        r4 = norm.rvs(size=self.n_samples)

        sigma_x = 1.51e-1  # cm
        sigma_y = 0.75e-1  # cm
        sigma_theta_x = 0.66e-3  # mrad
        sigma_theta_y = 0.40e-3  # mrad

        self.x0 = r1*sigma_x
        self.y0 = r2*sigma_y
        self.z0 = expon.rvs(scale=1/self.sigma_times_n, size=self.n_samples)

        self.px0 = sqrt(self.ep**2 - M_P**2)*r3*sigma_theta_x
        self.py0 = sqrt(self.ep**2 - M_P**2)*r4*sigma_theta_y
        self.pz0 = sqrt(self.ep**2 - M_P**2 - self.px0**2 - self.py0**2)

    def B(self, r):
        # B field in T for r in cm
        return heaviside(r - 2.2, 0.0) * (4*pi*1e-2) * 170 / (2*pi*r)
    
    def charged_meson_flux_mc(self, p_min, p_max, theta_min, theta_max):
        # Charged meson monte carlo flux simulation
        # Based on the Sanford-Wang and Feynman scaling parameterized proton prodution cross sections
        # momentum from [p_min, p_max] in GeV

        if self.meson_type not in ["pi_plus", "pi_minus", "k_plus", "K0S"]:
            raise Exception("meson_type not in list of available fluxes")
        
        meson_mass=M_PI
        meson_lifetime = PION_LIFETIME
        if self.meson_type == "k_plus" or self.meson_type == "K0S":
            meson_mass = M_K
            meson_lifetime = KAON_LIFETIME
        
        p_list = np.random.uniform(p_min, p_max, self.n_samples)
        theta_list = np.random.uniform(theta_min, theta_max, self.n_samples)

        xs_wgt = meson_production_d2SdpdOmega(p_list, theta_list, self.p_proton/1.0e3, meson_type=self.meson_type) * sin(theta_list)
        probability_decay = p_decay(p_list*1.0e3, meson_mass, meson_lifetime, 50)
        pi_plus_wgts = probability_decay * (2*pi*(theta_max-theta_min) * (p_max-p_min)) * self.n_pot * xs_wgt / self.n_samples / self.sigmap(self.p_proton)
        return np.array([p_list*1000.0, theta_list, pi_plus_wgts]).transpose()


    def focus_pions(self):
        self.pip_p_post_horn = []
        self.pip_theta_post_horn = []
        self.pip_wgt_post_horn = []
        self.acceptance_wgt = []
        phis = np.random.uniform(0.0, 2*pi, self.n_samples)
        for i, pflux in enumerate(self.meson_flux_pre_horn):
            p = pflux[0]
            theta = pflux[1]
            wgt = pflux[2]
            
            self.rksolver.set_new_particle(r0=[self.x0[i], self.y0[i], self.z0[i]],
                                           p0=[p*cos(phis[i])*sin(theta), p*sin(phis[i])*sin(theta), p*cos(theta)],
                                           charge=self.charge)
            self.rksolver.simulate(discard_history=True)

            p_final = sqrt(self.rksolver.px[-1]**2 + self.rksolver.py[-1]**2 + self.rksolver.pz[-1]**2)
            theta_final = arccos(self.rksolver.pz[-1] / p_final)

            self.pip_p_post_horn.append(p_final)
            self.pip_theta_post_horn.append(theta_final)
            self.pip_wgt_post_horn.append(wgt*(self.z0[i] <= 0.71))
            self.acceptance_wgt.append((self.z0[i] <= 0.71)*wgt*(theta_final < self.solid_angle_cut))
    
    def simulate(self):
        self.simulate_beam_spot()
        self.meson_flux_pre_horn = self.charged_meson_flux_mc(p_min=0.05, p_max=7.0,
                                        theta_min=0.0, theta_max=np.pi/4)
        print(self.meson_flux_pre_horn)

        self.focus_pions()


        


class RKHorn:
    """
    Takes in initial position r0 in meters
    initial momentum in MeV
    step size in nanoseconds
    Horn Current in kiloamps
    Elementary charge of particle
    """
    def __init__(self, r0=[0.0, 0.0, 0.0], p0 = [0.0, 0.0, 10.0], step=0.1, particle_mass=M_PI,
                    charge=1.0, horn_current=170.0):
        # initial conditions
        self.t = [0.0]
        self.x = [r0[0]]
        self.y = [r0[1]]
        self.z = [r0[2]]
        e0 = sqrt(p0[0]**2 + p0[1]**2 + p0[2]**2 + particle_mass**2)
        self.px = [p0[0]*C_LIGHT*1e-2 / e0]  # convert velocities to m / s
        self.py = [p0[1]*C_LIGHT*1e-2 / e0]
        self.pz = [p0[2]*C_LIGHT*1e-2 / e0]
        self.local_B = np.array([0.0, 0.0, 0.0])

        self.is_alive = True

        # solve params
        self.dt = step * 1e-9  # convert ns to s
        self.q = charge * CHARGE_COULOMBS  # convert to coulombs
        self.m = particle_mass / MEV_PER_KG  # convert to kg
        self.current = horn_current * 1000.0  # convert to Amps

        self.horn_start = 0.0  # meters
        self.horn_end = 1.9  # meters
        self.horn_radius_by_z = np.genfromtxt("data/mb_horn_radius_by_z.txt", delimiter=",")
        self.horn_max_radius = 0.3
    
    def inner_radius(self, z):
        return np.interp(z, 1e-2*self.horn_radius_by_z[:,0], 1e-2*self.horn_radius_by_z[:,1])

    def bfield(self, x, y, z):
        # takes in x, y, z in meters, returns in T
        r = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y,x)
        if r <= self.inner_radius(z) or r > self.horn_max_radius:  # minimum radius of horn
            return np.array([0.0, 0.0, 0.0])
        if z > self.horn_end:
            return np.array([0.0, 0.0, 0.0])
        if z < self.horn_start:
            return np.array([0.0, 0.0, 0.0])
        return (MU0 * self.current / (2 * pi * r)) * np.array([-sin(theta), cos(theta), 0.0])

    def dpxdtau(self, px, py, pz):
        return (self.q/self.m) * (self.local_B[2] * py - self.local_B[1] * pz)
    
    def dpydtau(self, px, py, pz):
        return (self.q/self.m) * (-self.local_B[2] * px + self.local_B[0] * pz)
    
    def dpzdtau(self, px, py, pz):
        return (self.q/self.m) * (self.local_B[1] * px - self.local_B[0] * py)

    def update(self):
        # Update B field
        self.local_B = self.bfield(self.x[-1], self.y[-1], self.z[-1])

        # Update velocity
        t_new, px_new, py_new, pz_new = RungeKuttaCoupled(self.t[-1], self.px[-1], self.py[-1], self.pz[-1],
                                                            self.dt, self.dpxdtau, self.dpydtau, self.dpzdtau)
        print("new px = ", px_new)

        self.px.append(px_new)
        self.py.append(py_new)
        self.pz.append(pz_new)
        self.t.append(t_new)
    
        # use velocity to update position
        # x(t + dt) = x(t) + dt * v(t+dt)
        self.x.append(self.x[-1] + self.dt * px_new)
        self.y.append(self.y[-1] + self.dt * py_new)
        self.z.append(self.z[-1] + self.dt * pz_new)

        #print("new x, y, z = ", self.x[-1], self.y[-1], self.z[-1])

    def set_new_particle(self, r0=[0.0, 0.0, 0.0], p0 =[0.0, 0.0, 10.0], particle_mass=M_PI, charge=1.0):
        # reset initial conditions
        self.t = [0.0]
        self.x = [r0[0]]
        self.y = [r0[1]]
        self.z = [r0[2]]
        e0 = sqrt(p0[0]**2 + p0[1]**2 + p0[2]**2 + particle_mass**2)
        self.px = [p0[0]*C_LIGHT*1e-2 / e0]  # convert velocities to m / s
        self.py = [p0[1]*C_LIGHT*1e-2 / e0]
        self.pz = [p0[2]*C_LIGHT*1e-2 / e0]
        self.local_B = np.array([0.0, 0.0, 0.0])

        self.q = charge * CHARGE_COULOMBS  # convert to coulombs
        self.m = particle_mass / MEV_PER_KG  # convert to kg

    def simulate(self, stop_z=2.0, discard_history=False):
        self.is_alive = True
        while self.is_alive:
            self.update()
            
            if discard_history and len(self.t) > 2:
                del self.t[0]
                del self.x[0]
                del self.y[0]
                del self.z[0]
                del self.px[0]
                del self.py[0]
                del self.pz[0]
            
            if len(self.t) > 10000:  # check time limit
                self.is_alive = False
            
            if abs(self.z[-1]) > stop_z:  # if we reach end of horn
                self.is_alive = False

            u = np.random.uniform(0.0,1.0)
            rnd_lifetime = -np.log(1-u) * pion_lifetime(self.px[-1]**2 + self.py[-1]**2 + self.pz[-1]**2)
            if rnd_lifetime < self.dt:
                self.is_alive = False

