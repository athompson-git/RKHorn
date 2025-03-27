from rkhorn import *
from rkhorn import ChargedPionFluxMiniBooNE
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.pylab import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

import numpy as np

#rkhornsim = ChargedPionFluxMiniBooNE(proton_energy=8.0e3 + M_P, meson_charge=1.0,
#                                     solid_angle_cut=0.00924, n_samples=1000000, n_pot=1,
#                                     horn_current=170.0)

# INPUT PARAMETERS
out_file = "data/bnb_focused_piminus_flux_rkhorn_1POT_BetheBloch.txt"
meson_type = "pi_minus"  # pi_minus, pi_plus, k_plus
if "minus" in meson_type:
    meson_charge = -1
else:
    meson_charge = 1
polarization = -1.0 # -1.0 for anti-nu mode, 1 for nu mode
n_samples = 1000000
p_proton = 8.89e3  # beam momentum
e_proton = sqrt(p_proton**2 + M_P**2)

rkhornsim = ChargedPionFluxMiniBooNE(proton_energy=e_proton, meson_charge=meson_charge,
                                     solid_angle_cut=0.00924, n_samples=n_samples, n_pot=1,
                                     horn_current=polarization*174.0, meson_species=meson_type)  # pi_plus, k_plus


rkhornsim.simulate(discard_zeros=True)


pi_thetas_pre_horn = rkhornsim.meson_flux_pre_horn[:,1]
pi_momenta_pre_horn = rkhornsim.meson_flux_pre_horn[:,0]
pi_weights_pre_horn = rkhornsim.meson_flux_pre_horn[:,2]

pi_thetas_post_horn = rkhornsim.pip_theta_post_horn
pi_momenta_post_horn = rkhornsim.pip_p_post_horn
pi_weights_post_horn = rkhornsim.pip_wgt_post_horn   # the number of pions that decayed pointing within detector solid angle

decay_positions = np.array(rkhornsim.decay_positions)
decay_wgts = np.array(rkhornsim.decay_in_pipe_wgt)
decay_x = decay_positions[:,0]
decay_y = decay_positions[:,1]
decay_z = decay_positions[:,2]


# PLOTTING THE DECAY HEAT MAP
z_values = np.linspace(0.0, 1.9, 50)
y_inner = rkhornsim.rksolver.inner_radius(z_values)
plt.hlines(y=[-0.3, 0.3], xmin=rkhornsim.rksolver.horn_start,
           xmax=rkhornsim.rksolver.horn_end, color="k", ls="dotted")
plt.vlines(x=[rkhornsim.rksolver.horn_start, rkhornsim.rksolver.horn_end],
           ymin=-0.3, ymax=0.3, color="k", ls="dotted")
plt.plot(z_values, y_inner, color="k")
plt.plot(z_values, -y_inner, color="k")

bins_z = np.linspace(0.0, 50.0, 100)
bins_y = np.linspace(-0.5, 0.5, 50)

# with weights
print(decay_wgts.shape, decay_z.shape, decay_y.shape)
plt.hist2d(decay_z, decay_y, weights=decay_wgts, bins=[bins_z, bins_y])
plt.ylabel(r"$y$ [m]", fontsize=12)
plt.xlabel(r"$z$ [m]", fontsize=12)
plt.colorbar()
plt.xlim((0.0, max(bins_z)))
plt.ylim((min(bins_y), max(bins_y)))
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.show()
plt.close()

# without weights
plt.hist2d(decay_z, decay_y, bins=[bins_z, bins_y])
plt.ylabel(r"$y$ [m]", fontsize=12)
plt.xlabel(r"$z$ [m]", fontsize=12)
plt.colorbar()
plt.xlim((0.0, max(bins_z)))
plt.ylim((min(bins_y), max(bins_y)))
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.show()
plt.close()

plt.hist(decay_z, bins=100)
plt.show()







# Plot energy/angle distributions
theta_bins = np.logspace(-4, np.log10(np.pi/4), 100)
momentum_bins = np.logspace(1, np.log10(7000.0), 100)

plt.hist(pi_thetas_post_horn, weights=pi_weights_post_horn, bins=theta_bins, histtype='step', color='b', density=False, label="Post-horn")
plt.hist(pi_thetas_post_horn, weights=pi_weights_post_horn*decay_wgts, bins=theta_bins, histtype='step', color='g', density=False, label="Post-horn (good decays)")
plt.hist(pi_thetas_pre_horn, weights=pi_weights_pre_horn, bins=theta_bins, density=False, color='k', histtype='step', label="Pre-horn")
plt.yscale('log')
plt.ylabel(r"$\pi^+$/POT", fontsize=14)
plt.xlabel(r"$\theta_\pi$ [rad]", fontsize=14)
plt.xscale('log')
plt.legend(loc="upper left")
plt.show()


plt.hist(pi_momenta_post_horn, weights=pi_weights_post_horn, bins=100, histtype='step', color='b', label="Post-horn")
plt.hist(pi_momenta_post_horn, weights=pi_weights_post_horn*decay_wgts, bins=100, histtype='step', color='g', label="Post-horn (good decays)")
plt.hist(pi_momenta_pre_horn, weights=pi_weights_pre_horn, bins=momentum_bins, color='k', histtype='step', label="Pre-horn")
plt.ylabel(r"$\pi^+$/POT", fontsize=14)
plt.xlabel(r"$p_\pi$ [MeV]", fontsize=14)
plt.legend(loc="upper left")
plt.show()

# dump to file
out_array = np.array([pi_momenta_post_horn, pi_thetas_post_horn, decay_x, decay_y, decay_z, decay_wgts, pi_weights_post_horn]).transpose()
np.savetxt(out_file, out_array)





# VISUALIZATION
rkhornsim2 = ChargedPionFluxMiniBooNE(proton_energy=8.0e3 + M_P, meson_charge=meson_charge,
                                     solid_angle_cut=0.00924, n_samples=20, n_pot=1,
                                     horn_current=polarization*300.0, meson_species=meson_type)

xs, ys, zs = rkhornsim2.simulate_with_histories()

horn_start = rkhornsim2.rksolver.horn_start
horn_end = rkhornsim2.rksolver.horn_end

# Plot the horn boundaries
z_values = np.linspace(0.0, 1.9, 50)
y_inner = rkhornsim2.rksolver.inner_radius(z_values)

fig = plt.figure(figsize=(20,10))
ax = fig.add_subplot(111)

plt.hlines(y=[-0.3, 0.3], xmin=horn_start, xmax=horn_end, color="k", ls="dotted")
plt.vlines(x=[horn_start, horn_end], ymin=-0.3, ymax=0.3, color="k", ls="dotted")
plt.plot(z_values, y_inner, color="k")
plt.plot(z_values, -y_inner, color="k")

for i in range(rkhornsim2.n_samples):
    plt.plot(zs[i], ys[i], color="k", markersize=1.0, alpha=0.5, marker='x', linewidth=1.0)


# plot the x component of the B field
z_vals_fine = np.linspace(0.0, 2.0, 200)
y_vals_fine = np.linspace(-0.5, 0.5, 200)
Z, Y = np.meshgrid(z_vals_fine, y_vals_fine)
BMAG = np.zeros_like(Z)
for i, z in enumerate(z_vals_fine):
    for j, y in enumerate(y_vals_fine):
        BMAG[i,j] = rkhornsim2.rksolver.bfield(0.0, Y[i,j], Z[i,j])[0]

plt.imshow(BMAG, extent=[z_vals_fine.min(), z_vals_fine.max(), y_vals_fine.min(), y_vals_fine.max()],
           origin='lower', aspect='auto', cmap='seismic')
cbar = fig.colorbar(label=r'$B_x$ [T]')
cbar.ax.yaxis.set_label_size(18)  # For vertical colorbar
cbar.ax.tick_params(labelsize=18)

ax.set_aspect('equal')

plt.ylabel(r"$y$ [m]", fontsize=18)
plt.xlabel(r"$z$ [m]", fontsize=18)
plt.ylim((-0.4,0.4))
plt.xlim((0.0, 2.0))
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.tight_layout()
plt.show()
plt.close()