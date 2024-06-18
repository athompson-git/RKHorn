import numpy as np
from numpy import sin, pi
import matplotlib.pyplot as plt

from matplotlib import cm, ticker
from matplotlib.pylab import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

import sys
sys.path.append("../")
from alplib.charged_meson_3body import *
from rkhorn import *


# Generate K+ and Pi+ flux
nsamples = 294377
k_plus_flux = charged_meson_flux_mc("k_plus", 0.75, 6.5, 0.0, 0.4, n_samples=nsamples, n_pot=1.0)
pi_plus_flux = charged_meson_flux_mc("pi_plus", 0.0, 6.5, 0.03, 0.21, n_samples=nsamples, n_pot=0.3)

# import RKHorn-generated pion flux
rkhorn_pion_data = np.genfromtxt("data/bnb_focused_piplus_flux_rkhorn_1POT_v2.txt")
rkhorn_kaon_data = np.genfromtxt("data/bnb_focused_Kplus_flux_rkhorn_1POT.txt")

# columns: pi_momenta_post_horn, pi_thetas_post_horn, decay_x, decay_y, decay_z, decay_wgts, pi_weights_post_horn
# create minimal array for ChargedMesonNeutrinoFlux
pi_plus_flux_rkhorn = rkhorn_pion_data[rkhorn_pion_data[:,5] > 0]
K_plus_flux_rkhorn = rkhorn_kaon_data[rkhorn_kaon_data[:,5] > 0]


class ChargedMesonNeutrinoFlux:
    def __init__(self, meson_flux, meson_type="pion", m_lepton=M_MU, energy_cut=140.0, n_samples=1):
        self.meson_flux = meson_flux
        if meson_type == "pion":
            self.mm = M_PI
            self.ckm = V_UD
            self.fM = F_PI
            self.total_width = PION_WIDTH
        elif meson_type == "kaon":
            self.mm = M_K
            self.ckm = V_US
            self.fM = F_K
            self.total_width = KAON_WIDTH
        else:
            raise Exception("Meson type not understood!", meson_type)
        self.m_lepton = m_lepton
        self.det_dist = 541
        self.dump_dist = 50
        self.det_length = 12
        self.solid_angles = []
        self.energy_cut = energy_cut
        self.nsamples = n_samples
        self.energies = []
        self.cosines = []
        self.decay_pos = []
        self.weights = []
    
    def partial_width(self):
        return (G_F*self.fM*self.m_lepton*self.ckm)**2 * self.mm * (1-(self.m_lepton/self.mm)**2)**2 / (8*pi)
    
    def br(self):
        return self.partial_width()/self.total_width

    def simulate_neutrinos(self):
        M_M = self.mm

        # Reproduce the neutrino flux abundance
        p_k = self.meson_flux[:,0]
        wgt_k = self.meson_flux[:,2]

        phi_rnd = 2*np.pi*np.random.ranf(nsamples)
        theta_rnd = arccos(1 - 2*np.random.ranf(nsamples))
        c_rnd = cos(theta_rnd)
        e_nu = (M_M**2 - M_MU**2)/(2*M_M)
        p_nu_z = e_nu * c_rnd

        beta = p_k / sqrt(p_k**2 + M_M**2)
        gamma = np.power(1-beta**2, -0.5)
        e_nu_lab = gamma*(e_nu + beta*p_nu_z)
        p_nu_z_lab = gamma*(p_nu_z + beta*e_nu)
        cos_lab = p_nu_z_lab / e_nu_lab

        # TODO: simulate decay position
        decay_l = METER_BY_MEV * p_k / self.partial_width() / self.mm
        u = []
        for l in decay_l:
            umax = exp(-2*self.dump_dist/l) * power(exp(self.dump_dist/l) - 1, 2) \
                if l > 1.0 else 1.0
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

        for i, m in enumerate(self.meson_flux):
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
                self.weights.append(flux_wgt*heaviside(self.det_length/dist_to_det - theta_nu_det, 0.0)/self.nsamples)

            
        
        return np.array(self.energies), np.array(self.weights)





# Simulate K+ and pi+ sourced neutrino fluxes
kaon_nu_gen = ChargedMesonNeutrinoFlux(k_plus_flux, meson_type="kaon")
pion_nu_gen = ChargedMesonNeutrinoFlux(pi_plus_flux, meson_type="pion")
pion_nu_rkhorn = ChargedMesonNeutrinoFlux(pi_plus_flux_rkhorn, meson_type="pion", n_samples=5)
kaon_nu_rkhorn = ChargedMesonNeutrinoFlux(K_plus_flux_rkhorn, meson_type="kaon", n_samples=5)

nu_k_e, nu_k_wgts = kaon_nu_gen.simulate_neutrinos()
nu_pi_e, nu_pi_wgts = pion_nu_gen.simulate_neutrinos()
nu_pi_e_rk, nu_pi_wgts_rk = pion_nu_rkhorn.simulate_neutrinos_from_decayed_flux()
nu_K_e_rk, nu_K_wgts_rk = kaon_nu_rkhorn.simulate_neutrinos_from_decayed_flux()

print("Sum of weights RKHorn = {}".format(np.sum(nu_pi_wgts_rk)))
print("Sum of weights SW truncated = {}".format(np.sum(nu_pi_wgts)))

# Import the neutrino fluxes to compare
nu_pi_dat = np.genfromtxt("data/MiniBooNE_numu_from_pion_per_POT_per50MeV.txt")
nu_k_dat = np.genfromtxt("data/MiniBooNE_numu_from_Kaon_per_POT_per50MeV.txt")
nu_pi_dat[:,0] *= 1e3
nu_k_dat[:,0] *= 1e3

def flux_integral(data, bins):
    centers = (bins[1:] + bins[:-1])/2
    def dPhidE(enu):
        return 1.13e6 * np.interp(enu/1000, data[:,0]/1000, data[:,1])
    
    return np.array([quad(dPhidE, bins[i], bins[i+1])[0] for i in range(bins.shape[0]-1)])


energy_bins = np.arange(0, 3800, 50)
e_bin_centers = (energy_bins[1:] + energy_bins[:-1])/2


plt.hist(nu_k_e, weights=nu_k_wgts, bins=energy_bins, histtype='step', 
            label=r"FS $K^+ \to \mu^+ \nu_\mu$", color="red")
plt.hist(nu_pi_e, weights=nu_pi_wgts, bins=energy_bins, histtype='step',
            label=r"SW $\pi^+ \to \mu^+ \nu_\mu$", color="blue")
plt.hist(nu_pi_e_rk, weights=nu_pi_wgts_rk, bins=energy_bins, histtype='step',
            label=r"SW + RKHorn $\pi^+ \to \mu^+ \nu_\mu$", color="mediumpurple", ls='dashed')
plt.hist(nu_pi_dat[:,0], weights=nu_pi_dat[:,1], bins=energy_bins, histtype='step',
            label=r"MiniBooNE $\pi^+$ model", color='lightsteelblue')
plt.hist(nu_k_dat[:,0], weights=nu_k_dat[:,1], bins=energy_bins, histtype='step',
            label=r"MiniBooNE $K^+$ model", color='violet')
#plt.yscale('log')
plt.xlim((0.0,energy_bins[-1]))
#plt.ylim((1e-5,10))
plt.xlabel(r"$E_\nu$ [MeV]")
plt.ylabel(r"$N_\nu$ / POT")
plt.legend()
plt.show()
plt.close()


plt.hist(nu_k_e, weights=nu_k_wgts, bins=energy_bins, histtype='step', 
            label=r"FS $K^+ \to \mu^+ \nu_\mu$", color="red")
plt.hist(nu_pi_e, weights=nu_pi_wgts, bins=energy_bins, histtype='step',
            label=r"SW $\pi^+ \to \mu^+ \nu_\mu$", color="blue")
plt.hist(nu_pi_e_rk, weights=nu_pi_wgts_rk, bins=energy_bins, histtype='step',
            label=r"SW + RKHorn $\pi^+ \to \mu^+ \nu_\mu$", color="blue", ls='dashed')
plt.hist(nu_K_e_rk, weights=nu_K_wgts_rk, bins=energy_bins, histtype='step',
            label=r"SW + RKHorn $K^+ \to \mu^+ \nu_\mu$", color="red", ls='dashed')
plt.plot(nu_pi_dat[:,0], nu_pi_dat[:,1], marker='_', ls='none',
            label=r"MiniBooNE $\pi^+$ model", color='lightsteelblue')
plt.plot(nu_k_dat[:,0], nu_k_dat[:,1], marker='_', ls='none',
            label=r"MiniBooNE $K^+$ model", color='salmon')
plt.yscale('log')
plt.xlim((0.0,3800.0))
#plt.ylim((1e-5,10))
plt.xlabel(r"$E_\nu$ [MeV]")
plt.ylabel(r"$N_\nu$ / POT")
plt.legend()
plt.tight_layout()
plt.show()
plt.close()



# compare to rates / 50 MeV / POT with PhysRevD.84.114021

print("Sum of my pi+ -> nu weights = ", np.sum(nu_pi_wgts_rk))
print("Sum of my K+ -> nu weights = ", np.sum(nu_K_wgts_rk))
print("Sum of MB pi+ -> nu weights = ", np.sum(nu_pi_dat[:,1]))
print("Sum of MB K+ -> nu weights = ", np.sum(nu_k_dat[:,1]))

energy_bins_50 = np.arange(0.0, 3800, 50)
plt.hist(nu_pi_e_rk, weights=nu_pi_wgts_rk, bins=energy_bins_50, histtype='step',
            label=r"SW + RKHorn $\pi^+ \to \mu^+ \nu_\mu$", color="b")
plt.hist(nu_K_e_rk, weights=nu_K_wgts_rk, bins=energy_bins_50, histtype='step',
            label=r"SW + RKHorn $K^+ \to \mu^+ \nu_\mu$", color="r")
plt.plot(nu_pi_dat[:,0], nu_pi_dat[:,1], marker='_', ls='none',
            label=r"MiniBooNE $\pi^+$ model", color='lightsteelblue')
plt.plot(nu_k_dat[:,0], nu_k_dat[:,1], marker='_', ls='none',
            label=r"MiniBooNE $K^+$ model", color='salmon')
plt.yscale('log')
plt.xlim((0.0,3800.0))
plt.ylim((1e-9,2e-4))
plt.xlabel(r"$E_\nu$ [MeV]", fontsize=14)
plt.ylabel(r"$N_\nu$ / POT", fontsize=14)
plt.legend()
plt.tight_layout()
plt.show()
plt.close()
