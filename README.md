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