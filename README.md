# RKHorn

A lightweight python library for simulating proton beam-induced charged pion fluxes and their radiative transport through a magnetic horn focusing system.



## Example Usage

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
We can visualize the tracks of the simulated mesons as they transport through the horn system. Let's start by simulating a small number of mesons, just ```n_samples=20```:
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