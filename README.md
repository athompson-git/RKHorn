# RKHorn

[![DOI](https://zenodo.org/badge/754366396.svg)](https://doi.org/10.5281/zenodo.14219232)

A lightweight python library for simulating proton beam-induced charged pion fluxes and their radiative transport through a magnetic horn focusing system.

Contact the author: a.thompson@northwestern.edu

### How it works:
The RKHorn classes implement a very basic set of approximations to model the charged meson fluxes produced at the Booster Neutrino Beam (BNB) target and horn system, with a future possibility of modeling other target and horn systems. There are several main pieces to this model:
* We use the Sanford-Wang (SW) parameterization of the pion ($\pi^+ / \pi^-$) production cross section from 8 GeV protons on beryllium
* Similarly, we implement the Feynman-Scaling (FS) parameterizations for charged Kaon production
* The charged mesons are propagated through the target and horn according to simple equations of motion solved by Runge-Kutta finite timestep propagation
* When the mesons are in the target, Bethe-Bloch energy loss is applied at each timestep
* Each timestep update is governed by the Lorentz force law equation of motion, with an azimuthally symmetric magnetic field that is only non-zero between the inner and outer conductors of the horn

Let's begin with the modeling of the proton beam spot. We take the gaussian profile of the beam spot on the face of the Be target from [1](####references); we simulate $(x_0,y_0,z_0)$ positions and momenta $(p_x,p_y,p_z)$ of the colliding proton beam of total momentum $p = \sqrt{(T_p + m_p)^2 - m_p^2} \simeq 8.9$ GeV,
where $T_p = 8$ GeV kinetic energy, by drawing 

$$ x_0 \sim N(0,\sigma_{x}) $$
$$ y_0 \sim N(0, \sigma_{y}) $$
$$ z_0 \sim Expon(\lambda) $$
$$ p_x \sim p \times N(0, \sigma_{\theta_{x}}) $$
$$ p_y \sim p \times N(0, \sigma_{\theta_{y}}) $$
$$ p_z \sim \sqrt{p^2 - p_x^2 - p_y^2} $$

where $\sigma_{x} = 0.151$ cm, $\sigma_{y} = 0.75$ cm, $\sigma_{\theta_{x}} = 0.66$ mrad, and $\sigma_{\theta_{y}} = 0.40$ mrad, describing the spatial and angular
distribution of the beam spot. The $z_p$ depth into the target at which the proton scatters is determined from an exponential distribution defined by the
mean free path $\lambda$ of the proton described by its total cross section and the Be density, $\lambda = 1/(n \sigma)$.
In the above scheme, we have made some approximations, one for the angular distribution given that the angular spread is very small, and secondly,
we have constrained the total momentum to remain $p$.

Once we record these positions and momenta of the protons where they collide, we sample the distributions and weights of produced pions and kaons
using the differential cross sections under the Sanford-Wang and Feynman-scaling approximations.

Once the mesons are produced, we push their positions and momenta with a finite difference routine. This routine will depend on the meson's current position.
if it resides within the target material, then there is no magnetic field influence, but the meson does lose energy according to the Bethe-Bloch energy loss
per unit distance traveled, $dE/dx$. 

If the meson instead resides within the magnetic field region described by $\vec{B}(r,z)$, then we push its momenta and positions according to the Lorentz
force law,

$$ \frac{dp_m}{dt} = \frac{q}{E_m} (\vec{p_m} \times \vec{B}) $$


#### References:
* [1] "A measurement of hadron production cross-sections for the simulation of accelerator neutrino beams and a search for muon neutrino to electron neutrino oscillations in the $\Delta m^2$ ~ 1 eV$^2$ region" David W. Schmitz, Thesis. [https://inspirehep.net/literature/791454](https://inspirehep.net/literature/791454)
* [2] "Improved Parameterization of K+ Production in p-Be Collisions at Low Energy Using Feynman Scaling" Mariani, Cheng, Conrad, Shaevitz [Phys.Rev.D 84 (2011) 114021] [https://arxiv.org/abs/1110.0417v2](https://arxiv.org/abs/1110.0417v2)
* [3] "Neutrino flux prediction at MiniBooNE" (MiniBooNE Collaboration) [PRD 79, 072002 (2009)] [https://arxiv.org/abs/0806.1449](https://arxiv.org/abs/0806.1449) 

## Example Usage

There are two example files;
* Validation tests of the neutrino flux are performed in validate_neutrino_flux.ipynb
* See also example.ipynb for a closed notebook form of the following example


#### Simple usage of ```rkhorn``` classes

One calls the ChargedPionFluxMiniBooNE class to simulate the pi+ / pi- fluxes coming out of the BNB horn geometry. However, with this class the incomingg proton flux can be changed as well as the selected horn current in kilo-Amps;

```
import rkhorn
from rkhorn import ChargedPionFluxMiniBooNE

rkhornsim = ChargedPionFluxMiniBooNE(proton_energy=8.0e3 + M_P, meson_charge=1.0,
                                     solid_angle_cut=0.00924, n_samples=50000, n_pot=1,
                                     horn_current=300.0)
```

To simulate the beamspot, the production of charged pions of a chosen ```meson_charge=1``` or ```meson_charge=-1```, and the focusing, we call
```
rkhornsim.simulate()
```
which will simulate ```n_samples``` number of weighted Monte-Carlo samples of charged mesons through the horn system. To extract the final momenta, angles, and weights, use
```
pi_thetas_post_horn = rkhornsim.pip_theta_post_horn
pi_momenta_post_horn = rkhornsim.pip_p_post_horn
pi_weights_post_horn = rkhornsim.pip_wgt_post_horn
```
The weights will be normalized to the number of protons-on-target (POT), ```n_pot``` that the user specifies above.


### Plotting example (1)
We can visualize the tracks of the simulated mesons as they transport through the horn system. Let's start by simulating a small number of mesons, just ```n_samples=20```, but keep the histories of each pion by using ```simulate_with_histories()``` instead of ```simulate()```, which returns a list of lists for the x, y, and z positions over time:
```
rkhornsim2 = ChargedPionFluxMiniBooNE(proton_energy=8.0e3 + M_P, meson_charge=1.0,
                                     solid_angle_cut=0.00924, n_samples=20, n_pot=1,
                                     horn_current=300.0)

xs, ys, zs = rkhornsim2.simulate_with_histories()

horn_start = rkhornsim2.rksolver.horn_start
horn_end = rkhornsim2.rksolver.horn_end

# Plot the horn boundaries
z_values = np.linspace(0.0, 1.9, 50)
y_inner = rkhornsim2.rksolver.inner_radius(z_values)

fig = plt.figure()
ax = fig.add_subplot(111)

plt.hlines(y=[-0.3, 0.3], xmin=horn_start, xmax=horn_end, color="k", ls="dotted")
plt.vlines(x=[horn_start, horn_end], ymin=-0.3, ymax=0.3, color="k", ls="dotted")
plt.plot(z_values, y_inner, color="k")
plt.plot(z_values, -y_inner, color="k")

for i in range(rkhornsim2.n_samples):
    plt.plot(zs[i], ys[i], color="b", linewidth=1.0, alpha=0.5)

ax.set_aspect('equal')

plt.ylabel(r"$y$ [m]", fontsize=12)
plt.xlabel(r"$z$ [m]", fontsize=12)
plt.ylim((-0.4,0.4))
plt.xlim((-0.1, 3.0))
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.show()
```

Giving us the following nice plot:
![Meson transport through the BNB focusing horn](./plots/horn_transport_visualization.png)

Here we can see the pi+ tracks projected onto the y,z plane in blue. Their origins are always inside the BNB target inside the Horn boundary. We cut off the simulation either when their z values reach 3 meters, or when the time exceeds 10000 seconds or when the pion decays.

### Plotting example 2

We can visualize the statistics for many more pions by increasing the sample size to ```n_samples=50000```. This should take several seconds.

```
rkhornsim = ChargedPionFluxMiniBooNE(proton_energy=8.0e3 + M_P, meson_charge=1.0,
                                     solid_angle_cut=0.00924, n_samples=50000, n_pot=1,
                                     horn_current=300.0)


rkhornsim.simulate()
```

The distributions can then be plotted with
```

pi_thetas_pre_horn = rkhornsim.meson_flux_pre_horn[:,1]
pi_momenta_pre_horn = rkhornsim.meson_flux_pre_horn[:,0]
pi_weights_pre_horn = rkhornsim.meson_flux_pre_horn[:,2]

pi_thetas_post_horn = rkhornsim.pip_theta_post_horn
pi_momenta_post_horn = rkhornsim.pip_p_post_horn
pi_weights_post_horn = rkhornsim.pip_wgt_post_horn   # the number of pions that decayed pointing within detector solid angle

print(pi_weights_post_horn, pi_weights_pre_horn)

theta_bins = np.logspace(-4, np.log10(np.pi/4), 100)
momentum_bins = np.logspace(1, np.log10(7000.0), 100)

plt.hist(pi_thetas_post_horn, weights=pi_weights_post_horn, bins=theta_bins, histtype='step', color='b', density=False, label="Post-horn")
plt.hist(pi_thetas_pre_horn, weights=pi_weights_pre_horn, bins=theta_bins, density=False, color='k', histtype='step', label="Pre-horn")
plt.yscale('log')
plt.ylabel(r"$\pi^+$/POT", fontsize=14)
plt.xlabel(r"$\theta_\pi$ [rad]", fontsize=14)
plt.xscale('log')
plt.legend(loc="upper left")
plt.show()
```

![Meson transport through the BNB focusing horn](./plots/pre_vs_post_horn_angles_nu-mode_300kA.png)

This visualizes the effect of the focusing; the pions post-horn are now pointed to lower angles with respect to the beam line, and therefore their decay products (neutrinos) will have a higher flux within the detector solid angle downstream! Notice also that for angles below 1e-3 radians, these pions were created with angles that simply point down the target volume and never enter the magnetic field, so their initial and final angles are the same. We do not simulate the random walk as pions scatter through the target material (yet).

![Meson transport through the BNB focusing horn](./plots/pre_vs_post_horn_momenta_nu-mode_300kA.png)

Here we can also look at the pre-horn and post-horn momenta distribution. The post-horn momentta are softer, as expected.


## TODO
* adding pion scattering and absorption in material
* implementing custom horn geometries
* implementing alternative target materials
* simulating decay positions

## Implementing a custom horn geometry

This feature is not yet implemented in a user-friendly way, but you may hack the file ```data/mb_horn_radius_by_z.txt``` which encodes the shape of the inner conductor of the horn as a function of z position.