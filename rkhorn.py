# Runge-Kutta method for time evolving a charged particle through a magnetic field

try:
    from .constants import *
    from .meson_physics import *
    from .bethe_bloch import *
    from .kinematics import *
except ImportError:
    from constants import *
    from meson_physics import *
    from bethe_bloch import *
    from kinematics import *

import numpy as np
from numpy import log, exp, pi, sqrt, power, sin, cos, arccos, heaviside
from scipy.stats import norm, expon

from tqdm import tqdm
import os

# All units converted internally to:
# seconds
# MeV
# meters


# Inspiration:
# https://www.geeksforgeeks.org/runge-kutta-4th-order-method-solve-differential-equation/
# https://math.stackexchange.com/questions/14942/runge-kutta-method-for-newton-law
# https://primer-computational-mathematics.github.io/book/c_mathematics/numerical_methods/5_Runge_Kutta_method.html




def runge_kutta_coupled(t, x, y, z, dt, dxdt, dydt, dzdt):
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
    def __init__(self, proton_energy=8938.0, meson_charge=1.0, solid_angle_cut=0.00924,
                 n_samples=10000, n_pot=18.75e20, horn_current=174.0, meson_species="pi_plus",
                 step_size_ns=0.1, verbose=False, absorption="reaction"):
        self.n_samples = n_samples
        self.kaon_sw_only = False
        self.solid_angle_cut = solid_angle_cut
        self.charge = meson_charge
        self.meson_type = meson_species
        self.n_pot = n_pot

        if meson_species == "pi_plus":
            self.charge = 1.0
            self.meson_mass = M_PI
            self.meson_lifetime = PION_LIFETIME
        if meson_species == "pi_minus":
            self.charge = -1.0
            self.meson_mass = M_PI
            self.meson_lifetime = PION_LIFETIME
        if meson_species == "k_plus":
            self.charge = 1.0
            self.meson_mass = M_K
            self.meson_lifetime = KAON_LIFETIME
        if meson_species == "k_minus":
            self.charge = -1.0
            self.meson_mass = M_K
            self.meson_lifetime = KAON_LIFETIME
        
        self.ep = proton_energy
        self.p_proton = np.sqrt(proton_energy**2 - M_P**2)

        # sigma (mb->cm2) * number density cm^-3 --> 1/meters
        self.sigma_reaction = proton_be_xs(self.p_proton/1.0e3, "prod")  # mb
        self.sigma_times_n = 1e2 * self.sigma_reaction * 1e-27 * AVOGADRO * DENSITY_BERYLLIUM / 9.012

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

        self.rksolver = RKHorn(horn_current=horn_current, particle_mass=self.meson_mass,
                               particle_charge=self.charge, step_size_ns=step_size_ns,
                               particle_lifetime=self.meson_lifetime, absorption=absorption)
        if self.meson_mass == M_K:
            self.rksolver.xs_scale = KAON_TO_PION_XS
        self.meson_flux_pre_horn = None

    def sigmap(self, p):
        a = 307.8
        b = 0.897
        c = -2.598
        d = -4.973
        n = 0.003
        return a + b*power(p/1.0e3,n) + c*log(p/1.0e3)*log(p/1.0e3) + d*log(p/1.0e3)

    def interaction_probability(self):
        return 1.0 - exp(-self.sigma_times_n * self.rksolver.target_length)

    def simulate_beam_spot(self):
        r1 = norm.rvs(size=self.n_samples)
        r2 = norm.rvs(size=self.n_samples)
        r3 = norm.rvs(size=self.n_samples)
        r4 = norm.rvs(size=self.n_samples)

        sigma_x = 1e-2 * 1.51e-1  # cm --> m
        sigma_y = 1e-2 * 0.75e-1  # cm --> m
        sigma_theta_x = 0.66e-3  # mrad
        sigma_theta_y = 0.40e-3  # mrad

        self.x0 = r1*sigma_x
        self.y0 = r2*sigma_y
        u = np.random.uniform(0.0, 1.0, self.n_samples)
        self.z0 = -log(1.0 - u*self.interaction_probability()) / self.sigma_times_n

        self.px0 = sqrt(self.ep**2 - M_P**2)*r3*sigma_theta_x
        self.py0 = sqrt(self.ep**2 - M_P**2)*r4*sigma_theta_y
        self.pz0 = sqrt(self.ep**2 - M_P**2 - self.px0**2 - self.py0**2)
    
    def charged_meson_flux_mc(self, p_min, p_max, theta_min, theta_max):
        # Charged meson monte carlo flux simulation
        # Based on the Sanford-Wang and Feynman scaling parameterized proton prodution cross sections
        # momentum from [p_min, p_max] in GeV

        if self.meson_type not in ["pi_plus", "pi_minus", "k_plus"]:
            raise Exception("meson_type not in list of available fluxes: pi_plus, pi_minus, k_plus")
        
        meson_mass=M_PI
        meson_lifetime = PION_LIFETIME
        if self.meson_type == "k_plus" or self.meson_type == "K0S":
            meson_mass = M_K
            meson_lifetime = KAON_LIFETIME
        
        p_list = np.random.uniform(p_min, p_max, self.n_samples)
        theta_list = np.random.uniform(theta_min, theta_max, self.n_samples)

        xs_wgt = meson_production_d2SdpdOmega(p_list, theta_list, self.p_proton/1.0e3, meson_type=self.meson_type, sw_only=self.kaon_sw_only) * sin(theta_list)
        #probability_decay = p_decay(p_list*1.0e3, meson_mass, meson_lifetime, 50)
        pi_plus_wgts = (2*pi*(theta_max-theta_min) * (p_max-p_min)) * self.n_pot * xs_wgt / self.n_samples \
                        * self.interaction_probability() / self.sigma_reaction
        return np.array([p_list*1000.0, theta_list, pi_plus_wgts]).transpose()

    def in_decay_volume(self, decay_pos, start_pos=None, direction=None):
        beam_stop, pipe_radius = 50.0, 0.9144
        coll_z1, coll_z2, coll_r1, coll_r2 = 2.59, 4.73, 0.30, 0.355
        z = decay_pos[2]
        if z > beam_stop or sqrt(decay_pos[0]**2 + decay_pos[1]**2) > pipe_radius:
            return 0.0
        if start_pos is not None and z > coll_z1 and direction[2] > 0.0:
            for z_check in [coll_z1, min(z, coll_z2)]:
                point = start_pos + direction * (z_check - start_pos[2]) / direction[2]
                aperture = coll_r1 + (coll_r2 - coll_r1) * (z_check - coll_z1) / (coll_z2 - coll_z1)
                if sqrt(point[0]**2 + point[1]**2) > aperture:
                    return 0.0
        return 1.0

    def focus_pions(self, discard_zeros=False, verbose=False):
        self.pip_p_post_horn = []
        self.pip_theta_post_horn = []
        self.pip_wgt_post_horn = []
        self.acceptance_wgt = []
        self.decay_positions = []
        self.decay_in_pipe_wgt = []
        self.pip_phi_post_horn = []
        phis = np.random.uniform(0.0, 2*pi, self.n_samples)
        for i, pflux in enumerate(tqdm(self.meson_flux_pre_horn)):
            p = pflux[0]
            theta = pflux[1]
            wgt = pflux[2]
            
            self.rksolver.set_new_particle(r0=[self.x0[i], self.y0[i], self.z0[i]],
                                           p0=[p*cos(phis[i])*sin(theta), p*sin(phis[i])*sin(theta), p*cos(theta)],
                                           charge=self.charge)
            self.rksolver.simulate(discard_history=True)

            p_final = sqrt(self.rksolver.px[-1]**2 + self.rksolver.py[-1]**2 + self.rksolver.pz[-1]**2)
            theta_final = arccos(self.rksolver.pz[-1] / p_final) if p_final > 0.0 else 0.0
            phi_final = np.arctan2(self.rksolver.py[-1], self.rksolver.px[-1])

            wgt = wgt * exp(-self.rksolver.optical_depth)

            # Handle the decays
            if self.rksolver.decay_pos is not None:  # if the decay already happened in the target or horn
                decay_in_vol = self.in_decay_volume(self.rksolver.decay_pos)

                if discard_zeros:
                    if decay_in_vol < 1.0:
                        continue
                
                self.decay_positions.append(self.rksolver.decay_pos)
                self.decay_in_pipe_wgt.append(decay_in_vol)
                self.pip_p_post_horn.append(p_final)
                self.pip_theta_post_horn.append(theta_final)
                self.pip_phi_post_horn.append(phi_final)
                self.pip_wgt_post_horn.append(wgt)
                self.acceptance_wgt.append(wgt*(theta_final < self.solid_angle_cut))

            else:
                # project out current vector to the MC decay position
                decay_length = expon.rvs(scale=(p_final*self.meson_lifetime*0.01*C_LIGHT/self.meson_mass))

                # get unit vector of current direction
                pvec = np.array([self.rksolver.px[-1],self.rksolver.py[-1],self.rksolver.pz[-1]])
                uvec = pvec / p_final
                rvec = np.array([self.rksolver.x[-1],self.rksolver.y[-1],self.rksolver.z[-1]])

                projected_decay_pos = rvec + decay_length * uvec
                
                decay_in_vol = self.in_decay_volume(projected_decay_pos, rvec, uvec)

                if discard_zeros:
                    if decay_in_vol < 1.0:
                        continue
                
                self.decay_positions.append(projected_decay_pos)
                self.decay_in_pipe_wgt.append(decay_in_vol)
                self.pip_p_post_horn.append(p_final)
                self.pip_theta_post_horn.append(theta_final)
                self.pip_phi_post_horn.append(phi_final)
                self.pip_wgt_post_horn.append(wgt)
                self.acceptance_wgt.append(wgt*(theta_final < self.solid_angle_cut))

    def simulate(self, discard_zeros=False, p_min=0.01, p_max=7.5, theta_min=0.0, theta_max=np.pi/2):
        self.simulate_beam_spot()
        self.meson_flux_pre_horn = self.charged_meson_flux_mc(p_min=p_min, p_max=p_max,
                                        theta_min=theta_min, theta_max=theta_max)

        self.focus_pions(discard_zeros)

    def decayed_flux_array(self):
        pos = np.array(self.decay_positions).reshape(-1, 3)
        return np.column_stack([self.pip_p_post_horn, self.pip_theta_post_horn,
                                pos[:,0], pos[:,1], pos[:,2], self.decay_in_pipe_wgt,
                                self.pip_wgt_post_horn, self.pip_phi_post_horn])
    
    def focus_pions_with_histories(self):
        self.pip_p_post_horn = []
        self.pip_theta_post_horn = []
        self.pip_wgt_post_horn = []
        self.acceptance_wgt = []
        phis = np.random.uniform(0.0, 2*pi, self.n_samples)

        x_histories = []
        y_histories = []
        z_histories = []

        for i, pflux in enumerate(self.meson_flux_pre_horn):
            p = pflux[0]
            theta = pflux[1]
            wgt = pflux[2]
            
            self.rksolver.set_new_particle(r0=[self.x0[i], self.y0[i], self.z0[i]],
                                           p0=[p*cos(phis[i])*sin(theta), p*sin(phis[i])*sin(theta), p*cos(theta)],
                                           charge=self.charge)
            self.rksolver.simulate(discard_history=False)

            p_final = sqrt(self.rksolver.px[-1]**2 + self.rksolver.py[-1]**2 + self.rksolver.pz[-1]**2)
            theta_final = arccos(self.rksolver.pz[-1] / p_final) if p_final > 0.0 else 0.0

            self.pip_p_post_horn.append(p_final)
            self.pip_theta_post_horn.append(theta_final)
            self.pip_wgt_post_horn.append(wgt)
            self.acceptance_wgt.append(wgt)

            x_histories.append(self.rksolver.x)
            y_histories.append(self.rksolver.y)
            z_histories.append(self.rksolver.z)
        
        return x_histories, y_histories, z_histories
    
    def simulate_with_histories(self):
        self.simulate_beam_spot()
        self.meson_flux_pre_horn = self.charged_meson_flux_mc(p_min=0.01, p_max=6.5,
                                        theta_min=0.0, theta_max=np.pi/2)
        
        x_hist, y_hist, z_hist = self.focus_pions_with_histories()

        return x_hist, y_hist, z_hist
    
    def del_zero_weights(self):
        pass



        


class RKHorn:
    """
    Takes in horn parameters
    step size in nanoseconds
    Horn Current in kiloamps
    Elementary charge of particle
    dat file of horn geometry (z in cm, r in cm)
    """
    def __init__(self, particle_mass=M_PI, particle_charge=1.0, horn_current=170.0,
                 horn_geometry_file="data/mb_horn_radius_by_z.txt", horn_outer_radius=0.3,
                 step_size_ns=0.1, target_radius=0.0048, target_length=0.711,
                 particle_lifetime=PION_LIFETIME, inner_conductor_thickness=0.003,
                 absorption="reaction"):
        self.target_radius = target_radius
        self.target_length = target_length

        self.conductor_thickness = inner_conductor_thickness

        self.absorption = absorption
        self.xs_scale = 1.0
        self.optical_depth = 0.0
        self.n_steps = 0
        self.proper_time_left = np.inf

        # initial conditions
        r0=[0.0, 0.0, 0.0]
        p0=[0.0, 0.0, 0.0]
        self.t = [0.0]
        self.x = [r0[0]]
        self.y = [r0[1]]
        self.z = [r0[2]]
        e0 = sqrt(p0[0]**2 + p0[1]**2 + p0[2]**2 + particle_mass**2)
        self.px = [p0[0]*C_LIGHT*1e-2 / e0]  # convert velocities to m / s
        self.py = [p0[1]*C_LIGHT*1e-2 / e0]
        self.pz = [p0[2]*C_LIGHT*1e-2 / e0]
        self.local_b = np.array([0.0, 0.0, 0.0])

        self.is_alive = True
        self.decay_pos = None

        # solve params
        self.dt = step_size_ns * 1e-9  # convert ns to s
        self.q = particle_charge * CHARGE_COULOMBS  # convert to coulombs
        self.m = particle_mass
        self.m_kg = particle_mass / MEV_PER_KG  # convert to kg
        self.tau = particle_lifetime
        self.current = horn_current * 1000.0  # convert to Amps

        fpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), horn_geometry_file)
        self.horn_radius_by_z = np.genfromtxt(fpath)
        self.horn_start = 0.0  # meters
        self.horn_end = max(1e-2*self.horn_radius_by_z[:,0])  # meters
        self.horn_max_radius = horn_outer_radius
    
    def inner_radius(self, z):
        return np.interp(z, 1e-2*self.horn_radius_by_z[:,0], 1e-2*self.horn_radius_by_z[:,1])

    def material_at(self, r, z):
        if (r < self.target_radius) and (0.0 <= z <= self.target_length):
            return BERYLLIUM
        if self.horn_start <= z <= self.horn_end:
            r_conductor = self.inner_radius(z)
            if r_conductor - self.conductor_thickness <= r <= r_conductor:
                return ALUMINUM
        return None

    def bfield(self, x, y, z):
        # takes in x, y, z in meters, returns in T
        r = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y,x)
        if (r <= self.inner_radius(z)) or (r > self.horn_max_radius):  # minimum radius of horn
            return np.array([0.0, 0.0, 0.0])
        if z > self.horn_end:
            return np.array([0.0, 0.0, 0.0])
        if z < self.horn_start:
            return np.array([0.0, 0.0, 0.0])
        return (MU0 * self.current / (2 * pi * r)) * np.array([-sin(theta), cos(theta), 0.0])

    # Momentum definitions, returning in units of MeV/s
    def dpxdtau(self, px, py, pz):
        return (MEV_PER_KG*self.q/sqrt(px**2 + py**2 + pz**2 + self.m**2)) * (self.local_b[2] * py - self.local_b[1] * pz)
    
    def dpydtau(self, px, py, pz):
        return (MEV_PER_KG*self.q/sqrt(px**2 + py**2 + pz**2 + self.m**2)) * (-self.local_b[2] * px + self.local_b[0] * pz)
    
    def dpzdtau(self, px, py, pz):
        return (MEV_PER_KG*self.q/sqrt(px**2 + py**2 + pz**2 + self.m**2)) * (self.local_b[1] * px - self.local_b[0] * py)

    def update(self):
        px, py, pz = self.px[-1], self.py[-1], self.pz[-1]
        energy = sqrt(px**2 + py**2 + pz**2 + self.m**2)
        half_step = 0.5 * self.dt * C_LIGHT*1e-2
        x_mid = self.x[-1] + half_step * px / energy
        y_mid = self.y[-1] + half_step * py / energy
        z_mid = self.z[-1] + half_step * pz / energy

        # Update B field
        self.local_b = self.bfield(x_mid, y_mid, z_mid)

        # Update momentum
        if (self.local_b[0] == 0.0) and (self.local_b[1] == 0.0) and (self.local_b[2] == 0.0):
            t_new, px_new, py_new, pz_new = self.t[-1] + self.dt, px, py, pz
        else:
            t_new, px_new, py_new, pz_new = runge_kutta_coupled(self.t[-1], px, py, pz,
                                                                self.dt, self.dpxdtau, self.dpydtau, self.dpzdtau)

        # Check if in material and apply Bethe-Bloch energy loss to the kinetic energy (KE)
        energy = sqrt(px_new**2 + py_new**2 + pz_new**2 + self.m**2)
        r = np.sqrt(x_mid**2 + y_mid**2)

        material = self.material_at(r, z_mid)
        if material is not None:
            p_new = sqrt(px_new**2 + py_new**2 + pz_new**2)
            step_cm = (p_new/energy)*C_LIGHT*self.dt  # (cm/s) * seconds

            if self.absorption is not None and p_new > 0.0:
                kind = "prod" if self.absorption == "reaction" else "ine"
                sigma_mb = self.xs_scale * pion_nucleus_xs(p_new/1.0e3, np.sign(self.q), material["name"], kind)
                self.optical_depth += (AVOGADRO * material["rho"] / material["A"]) * sigma_mb*1e-27 * step_cm

            # Bethe-Bloch energy loss
            delta_e = abs(dEdx(beta=p_new/energy, gamma=energy/self.m, M=self.m, Z=material["Z"], A=material["A"],
                              rho=material["rho"], Imin=material["I"])*step_cm)  # (MeV/cm) * cm

            p_after_loss = sqrt(power(np.clip(energy - delta_e, a_min=self.m, a_max=np.inf), 2) - self.m**2)

            px_new = (p_after_loss/p_new) * px_new
            py_new = (p_after_loss/p_new) * py_new
            pz_new = (p_after_loss/p_new) * pz_new
            energy = sqrt(p_after_loss**2 + self.m**2)

        self.px.append(px_new)
        self.py.append(py_new)
        self.pz.append(pz_new)
        self.t.append(t_new)

        self.x.append(x_mid + half_step * px_new / energy)
        self.y.append(y_mid + half_step * py_new / energy)
        self.z.append(z_mid + half_step * pz_new / energy)

    def set_new_particle(self, r0=[0.0, 0.0, 0.0], p0 =[0.0, 0.0, 10.0],
                         particle_mass=None, charge=1.0, particle_lifetime=None):
        # reset initial conditions
        if particle_mass is None:
            particle_mass = self.m
        if particle_lifetime is not None:
            self.tau = particle_lifetime
        self.is_alive = True
        self.decay_pos = None
        self.optical_depth = 0.0
        self.n_steps = 0
        self.proper_time_left = np.random.exponential(self.tau)
        self.t = [0.0]
        self.x = [r0[0]]
        self.y = [r0[1]]
        self.z = [r0[2]]
        e0 = sqrt(p0[0]**2 + p0[1]**2 + p0[2]**2 + particle_mass**2)
        self.px = [p0[0]]
        self.py = [p0[1]]
        self.pz = [p0[2]]
        self.local_b = np.array([0.0, 0.0, 0.0])

        self.q = charge * CHARGE_COULOMBS  # convert to coulombs
        self.m_kg = particle_mass / MEV_PER_KG  # convert to kg
        self.m = particle_mass

    def simulate(self, discard_history=False):
        self.is_alive = True
        while self.is_alive:
            self.update()
            
            self.n_steps += 1
            if discard_history and len(self.t) > 2:
                del self.t[0]
                del self.x[0]
                del self.y[0]
                del self.z[0]
                del self.px[0]
                del self.py[0]
                del self.pz[0]
            
            if self.n_steps > 10000:
                # check time limit
                self.is_alive = False
            
            if self.z[-1] > self.horn_end:
                # if we reach end of horn
                self.is_alive = False
            
            if self.z[-1] < 0.0:
                # if we go backwards
                self.is_alive = False
            
            if self.px[-1]**2 + self.py[-1]**2 + self.pz[-1]**2 <= 0.0:
                self.is_alive = False
                self.decay_pos = [self.x[-1], self.y[-1], self.z[-1]]

            p_now = sqrt(self.px[-1]**2 + self.py[-1]**2 + self.pz[-1]**2)
            self.proper_time_left -= self.dt * self.tau / lab_lifetime(p_now, self.m, self.tau)
            if self.proper_time_left < 0.0:
                self.is_alive = False
                self.decay_pos = [self.x[-1], self.y[-1], self.z[-1]]




def forced_detector_flux(meson_flux, m_meson=M_PI, m_lepton=M_MU, br=1.0, det_pos=(0.0, 0.0, 541.0)):
    """
    Neutrino flux at a point from two-body decays M -> lepton + nu, with no rejection sampling.

    Every decayed meson contributes exactly one weighted neutrino, the one that points at det_pos.
    In the meson rest frame the decay is isotropic, dN/dOmega* = 1/(4 pi), and E* = (M^2 - m_l^2)/(2M).
    For a meson with velocity beta, Lorentz factor gamma, at lab angle alpha to the line of sight:
        E_nu          = E* / (gamma (1 - beta cos(alpha)))
        dOmega*/dOmega = 1 / (gamma (1 - beta cos(alpha)))^2
    so the flux per unit area at distance L is  wgt * br / (4 pi L^2) * dOmega*/dOmega.

    meson_flux: array with the columns of ChargedPionFluxMiniBooNE.decayed_flux_array()
                (p [MeV], theta, x, y, z [m], decay wgt, wgt, phi). If the phi column is missing
                the azimuth is drawn at random.
    det_pos:    detector position in meters, same frame as the decay positions
    returns:    E_nu [MeV], flux weights [nu / cm^2 / (unit of wgt, e.g. POT)]
    """
    p, theta = meson_flux[:,0], meson_flux[:,1]
    r_decay = meson_flux[:,2:5]
    wgt = meson_flux[:,6]
    phi = meson_flux[:,7] if meson_flux.shape[1] > 7 else np.random.uniform(0.0, 2*pi, p.shape[0])

    energy = sqrt(p**2 + m_meson**2)
    gamma = energy / m_meson
    beta = p / energy

    line_of_sight = np.asarray(det_pos) - r_decay
    dist = sqrt(np.sum(line_of_sight**2, axis=1))  # meters
    u_meson = np.array([sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]).transpose()
    cos_alpha = np.sum(u_meson*line_of_sight, axis=1) / dist

    doppler = 1.0 / (gamma * (1.0 - beta*cos_alpha))
    e_nu_rest = (m_meson**2 - m_lepton**2) / (2*m_meson)

    return e_nu_rest*doppler, wgt * br * doppler**2 / (4*pi*power(100.0*dist, 2))




class ChargedMesonNeutrinoFlux:
    """
    Simulates the neutrino flux from charged mesons produced in the RKHorn.
    """
    def __init__(self, meson_flux, meson_type="pion", m_lepton=M_MU, energy_cut=140.0, n_samples=1):
        self.meson_flux = meson_flux
        if meson_type == "pion":
            self.mm = M_PI
            self.ckm = V_UD
            self.f_m = F_PI
            self.total_width = PION_WIDTH
        elif meson_type == "kaon":
            self.mm = M_K
            self.ckm = V_US
            self.f_m = F_K
            self.total_width = KAON_WIDTH
        else:
            raise Exception("Meson type not understood!", meson_type)
        self.m_lepton = m_lepton
        self.det_dist = 541
        self.dump_dist = 50
        self.det_length = 12
        self.det_radius = 6.106
        self.solid_angles = []
        self.energy_cut = energy_cut
        self.nsamples = n_samples
        self.energies = []
        self.cosines = []
        self.decay_pos = []
        self.weights = []
    
    def partial_width(self):
        return np.clip((G_F*self.f_m*self.m_lepton*self.ckm)**2 * self.mm * (1-(self.m_lepton/self.mm)**2)**2 / (8*pi), a_min=0.0, a_max=np.inf)
    
    def br(self):
        if self.m_lepton == M_MU:
            return BR_PI_MUNU if self.mm == M_PI else BR_K_MUNU
        return self.partial_width()/self.total_width

    def flux_at_detector(self):
        return forced_detector_flux(self.meson_flux, m_meson=self.mm, m_lepton=self.m_lepton, br=self.br(),
                                    det_pos=(0.0, 0.0, self.det_dist))

    def simulate_neutrinos(self):
        m_m = self.mm

        # Reproduce the neutrino flux abundance
        p_k = self.meson_flux[:,0]
        wgt_k = self.meson_flux[:,2]

        theta_rnd = arccos(1 - 2*np.random.ranf(self.nsamples))
        c_rnd = cos(theta_rnd)
        e_nu = (m_m**2 - M_MU**2)/(2*m_m)
        p_nu_z = e_nu * c_rnd

        beta = p_k / sqrt(p_k**2 + m_m**2)
        gamma = np.power(1-beta**2, -0.5)
        e_nu_lab = gamma*(e_nu + beta*p_nu_z)
        p_nu_z_lab = gamma*(p_nu_z + beta*e_nu)
        cos_lab = p_nu_z_lab / e_nu_lab

        # TODO: simulate decay position
        decay_lengths = METER_BY_MEV * p_k / self.partial_width() / self.mm
        u = []
        for length in decay_lengths:
            umax = exp(-2*self.dump_dist/length) * power(exp(self.dump_dist/length) - 1, 2) \
                if length > 1.0 else 1.0
            u.append(np.random.uniform(0.0, min(umax, 1.0)))
        
        x = decay_quantile(np.array(u), p_k, self.mm, self.partial_width())
        
        # Append decay positions and solid angle cosines for the geometric acceptance of each meson decay
        self.decay_pos.append(x)
        solid_angle_cosine = cos(arctan(self.det_length/(self.det_dist-x)/2))

        wgt_nu = self.br() * wgt_k * heaviside(cos_lab - solid_angle_cosine, 0.0)  * (2 * (max(e_nu_lab)-min(e_nu_lab)))
        return e_nu_lab, wgt_nu
    
    def simulate_neutrinos_from_decayed_flux(self):
        self.energies = []
        self.cosines = []
        self.decay_pos = []
        self.weights = []

        det_vector = np.array([0,0,self.det_dist])

        meson_p4 = LorentzVector(self.mm, 0.0, 0.0, 0.0)
        mc = Decay2Body(p_parent=LorentzVector(self.mm, 0.0, 0.0, 0.0), m1=0.0, m2=self.m_lepton, n_samples=self.nsamples)
        phis = np.random.uniform(0.0, 2*pi, self.meson_flux.shape[0]*self.nsamples)

        for i, m in tqdm(enumerate(self.meson_flux)):
            p_m = m[0]
            theta_m = m[1]
            x = m[2]
            y = m[3]
            z = m[4]

            flux_wgt = m[6]

            meson_p4.set_p4(sqrt(p_m**2 + self.mm**2), p_m*sin(theta_m)*cos(phis[i]),
                                                        p_m*sin(theta_m)*sin(phis[i]),
                                                        p_m*cos(theta_m))
            mc.set_new_decay(meson_p4, 0.0, self.m_lepton)
            mc.decay()

            nu_p4s = mc.p1_lab_4vectors

            for nu_p4 in nu_p4s:
                # determine the dot product between the neutrino production position vector and the detector vector
                lnu = det_vector - np.array([x,y,z])
                dist_to_det = np.sqrt(np.sum(lnu*lnu))
                pnu = np.array([nu_p4.p1, nu_p4.p2, nu_p4.p3])
                theta_nu_det = arccos(np.dot(lnu, pnu)/(nu_p4.momentum()*dist_to_det))

                self.energies.append(nu_p4.energy())
                self.cosines.append(nu_p4.cosine())
                self.weights.append(self.br()*flux_wgt*heaviside(self.det_radius/dist_to_det - theta_nu_det, 0.0)/self.nsamples)

            
        
        return np.array(self.energies), np.array(self.weights)