"""
Hand-checkable physics tests for RKHorn. Each test compares the code with a number that can be
worked out on paper or read from the references.

    python tests/test_physics.py        (or: pytest tests/test_physics.py)

References:
  [MB]  MiniBooNE, "Neutrino flux prediction at MiniBooNE", PRD 79, 072002 (arXiv:0806.1449)
  [FS]  Mariani et al., PRD 84, 114021 (arXiv:1110.0417)
"""

import contextlib
import io
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import numpy as np

from constants import *
from bethe_bloch import dEdx, BERYLLIUM
from meson_physics import proton_be_xs, pion_nucleus_xs, meson_production_d2SdpdOmega, lab_lifetime
from rkhorn import RKHorn, ChargedPionFluxMiniBooNE, ChargedMesonNeutrinoFlux, forced_detector_flux

P_BEAM = 8.89
E_BEAM = np.sqrt((1e3*P_BEAM)**2 + M_P**2)


def close(value, expected, rel):
    assert abs(value - expected) <= rel*abs(expected), "got {}, expected {} within {:.1%}".format(value, expected, rel)


def free_tracker(**kwargs):
    rk = RKHorn(target_length=-1.0, inner_conductor_thickness=0.0, absorption=None,
                particle_lifetime=1e9, **kwargs)
    return rk


def test_proton_beryllium_cross_sections():
    close(proton_be_xs(P_BEAM, "tot"), 285.4, 0.002)
    close(proton_be_xs(P_BEAM, "ine"), 212.4, 0.002)
    close(proton_be_xs(P_BEAM, "qel"), 34.9, 0.005)
    close(proton_be_xs(P_BEAM, "prod"), 177.5, 0.002)


def test_pion_cross_sections_match_figures():
    close(pion_nucleus_xs(2.0, +1, "Be", "ine"), 195.0, 0.06)
    close(pion_nucleus_xs(4.0, +1, "Be", "ine"), 180.0, 0.06)
    close(pion_nucleus_xs(4.0, -1, "Be", "ine"), 170.0, 0.06)
    close(pion_nucleus_xs(2.0, +1, "Al", "ine"), 400.0, 0.06)
    close(pion_nucleus_xs(4.0, +1, "Be", "qel"), 22.0, 0.10)
    close(pion_nucleus_xs(4.0, +1, "Al", "qel"), 52.0, 0.10)
    assert pion_nucleus_xs(1e-3, +1, "Be", "ine") >= 0.0 and pion_nucleus_xs(1e-3, +1, "Be", "prod") >= 0.0


def test_kaon_feynman_scaling_table():
    for p, theta, expected in [(1.52, 0.213, 9.37), (2.80, 0.106, 8.67), (4.30, 0.055, 4.73)]:
        close(meson_production_d2SdpdOmega(p, theta, P_BEAM, "k_plus"), expected, 0.01)


def test_cross_sections_never_negative():
    p = np.linspace(0.01, 8.8, 500)
    for species in ["pi_plus", "pi_minus", "k_plus"]:
        for sw_only in [False, True]:
            assert np.all(meson_production_d2SdpdOmega(p, 0.05, P_BEAM, species, sw_only=sw_only) >= 0.0)


def test_meson_multiplicities():
    rng = np.random.default_rng(1)
    n = 1000000
    p = rng.uniform(0.01, 8.5, n)
    theta = rng.uniform(0.0, np.pi/2, n)
    volume = 2*np.pi * (np.pi/2) * 8.49
    for species, expected in [("pi_plus", 0.8825), ("k_plus", 0.0689)]:
        d2s = meson_production_d2SdpdOmega(p, theta, P_BEAM, species) * np.sin(theta)
        sigma = np.sum(d2s[p*np.sin(theta) < 1.0]) * volume / n
        close(sigma / proton_be_xs(P_BEAM, "prod"), expected, 0.02)


def test_interaction_depth_and_probability():
    sim = ChargedPionFluxMiniBooNE(proton_energy=E_BEAM, n_samples=20000, n_pot=1)
    close(1.0/sim.sigma_times_n, 0.456, 0.005)
    close(sim.interaction_probability(), 0.790, 0.005)
    sim.simulate_beam_spot()
    assert sim.z0.min() >= 0.0 and sim.z0.max() <= sim.rksolver.target_length
    lam, length = 0.456, 0.711
    close(np.mean(sim.z0), lam - length*np.exp(-length/lam)/0.790, 0.03)


def test_field_strength_and_region():
    rk = free_tracker(horn_current=174.0)
    b = rk.bfield(0.0221, 0.0, 0.5)
    close(np.sqrt(np.sum(b**2)), 2e-7*174e3/0.0221, 1e-6)
    assert b[1] > 0.0 and b[0] == 0.0
    assert not np.any(rk.bfield(0.01, 0.0, 0.5))
    assert not np.any(rk.bfield(0.31, 0.0, 0.5))
    assert not np.any(rk.bfield(0.10, 0.0, 2.5))


def test_gyroradius_in_uniform_field():
    p, b_field, n_steps = 1000.0, 1.0, 200
    rk = free_tracker(step_size_ns=0.02)
    rk.bfield = lambda x, y, z: np.array([0.0, b_field, 0.0])
    rk.set_new_particle(r0=[0.0, 0.0, 0.0], p0=[0.0, 0.0, p], charge=1.0)
    for _ in range(n_steps):
        rk.update()
    beta = p / np.sqrt(p**2 + M_PI**2)
    path = beta * C_LIGHT*1e-2 * n_steps * rk.dt
    radius = (p/1e3) / (0.299792458 * b_field)
    turned = np.arctan2(abs(rk.px[-1]), rk.pz[-1])
    close(turned, path/radius, 1e-4)
    close(np.sqrt(rk.px[-1]**2 + rk.pz[-1]**2), p, 1e-9)


def test_horn_kick_and_focusing_sign():
    p, r0, length = 3000.0, 0.10, 0.30
    expected = 0.299792458 * (2e-7*174e3/r0) * length / (p/1e3)
    for charge, sign in [(1.0, -1.0), (-1.0, +1.0)]:
        rk = free_tracker(horn_current=174.0, step_size_ns=0.01, particle_charge=charge)
        rk.set_new_particle(r0=[r0, 0.0, 0.5], p0=[0.0, 0.0, p], charge=charge)
        while rk.z[-1] < 0.5 + length:
            rk.update()
        close(rk.px[-1]/rk.pz[-1], sign*expected, 0.03)


def test_step_size_convergence():
    angles = []
    for dt in [0.1, 0.005]:
        rk = free_tracker(horn_current=174.0, step_size_ns=dt)
        rk.set_new_particle(r0=[0.0, 0.0, 0.2], p0=[2000.0*np.sin(0.1), 0.0, 2000.0*np.cos(0.1)], charge=1.0)
        while 0.0 <= rk.z[-1] <= rk.horn_end:
            rk.update()
        angles.append(np.arctan2(rk.px[-1], rk.pz[-1]))
    assert abs(angles[0] - angles[1]) < 0.006


def test_bethe_bloch_minimum():
    beta_gamma = 3.5
    gamma = np.sqrt(1 + beta_gamma**2)
    value = dEdx(beta=beta_gamma/gamma, gamma=gamma, M=M_PI, Z=4, A=9.012)
    assert 2.95 <= value <= 3.15, value


def test_energy_loss_and_absorption_along_target():
    p0, length = 2000.0, 0.50
    rk = RKHorn(horn_current=174.0, step_size_ns=0.01, particle_lifetime=1e9,
                target_length=length, absorption="reaction")
    rk.set_new_particle(r0=[0.0, 0.0, 0.0], p0=[0.0, 0.0, p0], charge=1.0)
    while rk.z[-1] < length + 0.05:
        rk.update()
    e0 = np.sqrt(p0**2 + M_PI**2)
    loss = e0 - np.sqrt(rk.pz[-1]**2 + M_PI**2)
    close(loss, dEdx(beta=p0/e0, gamma=e0/M_PI, M=M_PI, Z=4, A=9.012) * 100*length, 0.02)
    n_be = AVOGADRO * BERYLLIUM["rho"] / BERYLLIUM["A"]
    close(rk.optical_depth, n_be * pion_nucleus_xs(p0/1e3, +1, "Be", "prod")*1e-27 * 100*length, 0.02)


def test_decay_in_flight_fraction():
    np.random.seed(3)
    n = 4000
    for mass, tau, p in [(M_PI, PION_LIFETIME, 300.0), (M_K, KAON_LIFETIME, 500.0)]:
        rk = RKHorn(particle_mass=mass, particle_lifetime=tau, target_length=-1.0,
                    inner_conductor_thickness=0.0, absorption=None)
        decayed = 0
        for _ in range(n):
            rk.set_new_particle(r0=[0.0, 0.0, 0.0], p0=[0.0, 0.0, p], charge=1.0)
            rk.simulate(discard_history=True)
            decayed += rk.decay_pos is not None
        expected = 1.0 - np.exp(-rk.horn_end / ((p/mass) * C_LIGHT*1e-2 * tau))
        close(decayed/n, expected, 0.08)
    close(lab_lifetime(300.0, M_PI, PION_LIFETIME), np.sqrt(300.0**2 + M_PI**2)/M_PI * PION_LIFETIME, 1e-12)


def test_decay_volume_geometry():
    sim = ChargedPionFluxMiniBooNE(proton_energy=E_BEAM, n_samples=1, n_pot=1)
    axis = np.array([0.0, 0.0, 1.0])
    start = np.array([0.0, 0.0, 2.0])
    assert sim.in_decay_volume(np.array([0.0, 0.0, 30.0]), start, axis) == 1.0
    assert sim.in_decay_volume(np.array([0.0, 0.0, 50.5]), start, axis) == 0.0
    assert sim.in_decay_volume(np.array([0.95, 0.0, 30.0]), start, axis) == 0.0
    off_axis = np.array([0.32, 0.0, 2.0])
    assert sim.in_decay_volume(np.array([0.32, 0.0, 30.0]), off_axis, axis) == 0.0
    assert sim.in_decay_volume(np.array([0.32, 0.0, 2.3]), off_axis, axis) == 1.0


def test_forced_detector_flux_on_axis():
    p, z = 2000.0, 41.0
    mesons = np.array([[p, 0.0, 0.0, 0.0, z, 1.0, 1.0, 0.0]])
    e_nu, flux = forced_detector_flux(mesons, m_meson=M_PI, br=1.0)
    energy = np.sqrt(p**2 + M_PI**2)
    boost = (energy + p) / M_PI
    close(e_nu[0], 29.79 * boost, 1e-3)
    close(flux[0], boost**2 / (4*np.pi*(100.0*(541.0 - z))**2), 1e-9)


def test_forced_detector_flux_against_rejection_sampling():
    np.random.seed(5)
    n_mesons, p = 40, 1000.0
    mesons = np.zeros((n_mesons, 8))
    mesons[:, 0] = p
    mesons[:, 1] = 0.004
    mesons[:, 4] = 20.0
    mesons[:, 5] = 1.0
    mesons[:, 6] = 1.0 / n_mesons
    mesons[:, 7] = np.random.uniform(0.0, 2*np.pi, n_mesons)
    _, flux = forced_detector_flux(mesons, m_meson=M_PI, br=BR_PI_MUNU)

    brute = ChargedMesonNeutrinoFlux(mesons, meson_type="pion", n_samples=2500)
    with contextlib.redirect_stderr(io.StringIO()):
        _, weights = brute.simulate_neutrinos_from_decayed_flux()
    area = np.pi * (100.0*brute.det_radius)**2
    close(np.sum(weights)/area, np.sum(flux), 0.06)


if __name__ == "__main__":
    tests = [(name, func) for name, func in sorted(globals().items()) if name.startswith("test_")]
    failed = 0
    for name, func in tests:
        try:
            func()
            print("PASS  " + name)
        except Exception as error:
            failed += 1
            print("FAIL  {}: {}".format(name, error))
    print("{} of {} tests passed".format(len(tests) - failed, len(tests)))
    sys.exit(1 if failed else 0)
