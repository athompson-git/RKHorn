#!/usr/bin/env python3
"""BSM particle fluxes at SBND from BNB charged mesons."""

import argparse
import os

import numpy as np
from scipy import integrate

from constants import (BR_PI_MUNU, F_K, F_PI, G_F, KAON_WIDTH, M_E, M_K,
                       M_MU, M_PI, METER_BY_MEV, PION_WIDTH, V_UD, V_US)
import flux_io
from rkhorn import forced_detector_flux

VEV = (np.sqrt(2.0)*G_F)**-0.5
KAPPA = 2.0/9.0
M_TOP = 172.7e3
VTD_VTS = 3.5e-4
SBND_Z = 110.0
FACE_SIDE = 4.0
EDGES = np.linspace(0.0, 7600.0, 153)
OUT_DIR = os.path.join(flux_io.DATA_DIR, "bsm_data")
SAMPLES = (("pi_plus", 174.0), ("k_plus", 174.0), ("pi_minus", -174.0))
PARENTS = {"pi": (M_PI, F_PI, V_UD, PION_WIDTH),
           "k": (M_K, F_K, V_US, KAON_WIDTH)}
ZP_MASSES = (10.0, 20.0, 30.0)
ZP_COUPLING = 1.0e-4
ZP_ALPHA_D = 0.5
ZP_NU_CHARGE = 1.0
S_MASSES = (10.0, 50.0, 100.0)
S_SIN_THETA = 1.0e-4
N_REST = 200000
N_THROW = 64

_PAULI = (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]),
          np.array([[1, 0], [0, -1]]))
_ZERO = np.zeros((2, 2))
_UNIT = np.eye(2)
GAMMA = np.array([np.block([[_UNIT, _ZERO], [_ZERO, -_UNIT]])]
                 + [np.block([[_ZERO, s], [-s, _ZERO]]) for s in _PAULI])
GAMMA5 = np.block([[_ZERO, _UNIT], [_UNIT, _ZERO]])
CHIRAL = np.eye(4) - GAMMA5
METRIC = np.array([1.0, -1.0, -1.0, -1.0])


def mdot(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Minkowski product of four-vectors of shape (n, 4)."""
    return np.sum(a*b*METRIC, axis=-1)


def slash(p: np.ndarray) -> np.ndarray:
    """Contracts four-vectors of shape (n, 4) with gamma."""
    return np.einsum("na,aij->nij", p*METRIC, GAMMA)


def spin_trace(lep: np.ndarray, m_l: float, nu: np.ndarray,
               v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
    """Spin sum of (u v1 v)(u v2 v)* for lepton and neutrino.

    Args:
        lep: Charged lepton four-momenta, shape (n, 4).
        m_l: Charged lepton mass in MeV.
        nu: Massless neutrino four-momenta, shape (n, 4).
        v1: First vertex matrices, shape (n, 4, 4).
        v2: Second vertex matrices, shape (n, 4, 4).

    Returns:
        Real part of the trace for each event.
    """
    bar = GAMMA[0] @ np.conj(np.swapaxes(v2, -1, -2)) @ GAMMA[0]
    chain = (slash(lep) + m_l*np.eye(4)) @ v1 @ slash(nu)
    return np.einsum("nij,nji->n", chain, bar).real


def kallen(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Kallen triangle function."""
    return a*a + b*b + c*c - 2.0*(a*b + a*c + b*c)


def isotropic(n: int, rng: np.random.Generator) -> np.ndarray:
    """Draws n unit vectors uniformly on the sphere."""
    cos_t = rng.uniform(-1.0, 1.0, n)
    sin_t = np.sqrt(1.0 - cos_t**2)
    phi = rng.uniform(0.0, 2.0*np.pi, n)
    return np.column_stack([sin_t*np.cos(phi), sin_t*np.sin(phi), cos_t])


def boost(p: np.ndarray, parent: np.ndarray) -> np.ndarray:
    """Boosts p from the parent rest frame to the frame of parent."""
    mass = np.sqrt(mdot(parent, parent))
    p_dot = np.sum(p[:, 1:]*parent[:, 1:], axis=1)
    energy = (p[:, 0]*parent[:, 0] + p_dot)/mass
    scale = p_dot/(mass*(parent[:, 0] + mass)) + p[:, 0]/mass
    return np.column_stack([energy, p[:, 1:] + scale[:, None]*parent[:, 1:]])


def two_body(parent: np.ndarray, m1: np.ndarray, m2: np.ndarray,
             rng: np.random.Generator) -> tuple[np.ndarray, ...]:
    """Decays parents isotropically into two particles.

    Args:
        parent: Parent four-momenta, shape (n, 4).
        m1: Mass of the first product in MeV.
        m2: Mass of the second product in MeV.
        rng: Random generator.

    Returns:
        Both four-momenta and the rest-frame momentum.
    """
    m0_sq = mdot(parent, parent)
    p_star = np.sqrt(np.maximum(kallen(m0_sq, m1**2, m2**2), 0.0)/m0_sq)/2.0
    vec = p_star[:, None]*isotropic(parent.shape[0], rng)
    p1 = np.column_stack([np.sqrt(p_star**2 + m1**2), vec])
    p2 = np.column_stack([np.sqrt(p_star**2 + m2**2), -vec])
    return boost(p1, parent), boost(p2, parent), p_star


def three_body(m0: float, m1: float, m2: float, m3: float, n: int,
               rng: np.random.Generator) -> tuple[np.ndarray, ...]:
    """Samples flat three-body phase space in the rest frame.

    Args:
        m0: Parent mass in MeV.
        m1: Mass of particle 1 in MeV.
        m2: Mass of particle 2 in MeV.
        m3: Mass of particle 3 in MeV.
        n: Number of events.
        rng: Random generator.

    Returns:
        Three four-momenta and the phase-space weights.
    """
    s_min, s_max = (m1 + m2)**2, (m0 - m3)**2
    s_pair = rng.uniform(s_min, s_max, n)
    rest = np.zeros((n, 4))
    rest[:, 0] = m0
    pair, p3, k_out = two_body(rest, np.sqrt(s_pair), np.full(n, m3), rng)
    p1, p2, k_in = two_body(pair, np.full(n, m1), np.full(n, m2), rng)
    weight = ((s_max - s_min)/(2.0*np.pi)*k_out/(4.0*np.pi*m0)
              * k_in/(4.0*np.pi*np.sqrt(s_pair)))
    return p1, p2, p3, weight


def scalar_msq(lep: np.ndarray, nu: np.ndarray, k: np.ndarray, m_l: float,
               legs: tuple[float, float, float] = (1.0, 1.0, 1.0)
               ) -> np.ndarray:
    """Spin-summed |M|^2 of M -> l nu S over (c sin(theta))^2.

    Args:
        lep: Charged lepton four-momenta.
        nu: Neutrino four-momenta.
        k: Scalar four-momenta.
        m_l: Charged lepton mass in MeV.
        legs: Switches for lepton, contact and meson emission.

    Returns:
        Squared amplitude in units of the meson constant c^2.
    """
    p = lep + nu + k
    m0_sq, m_s_sq = mdot(p, p), mdot(k, k)
    current = slash(p) @ CHIRAL
    pole_l = m_l/VEV/(2.0*mdot(lep, k) + m_s_sq)
    pole_m = (((1.0 + KAPPA)*m0_sq + 2.0*KAPPA*mdot(p, k))
              / (m_s_sq - 2.0*mdot(p, k))/VEV)
    vertex = (legs[0]*pole_l[:, None, None]
              * ((slash(lep + k) + m_l*np.eye(4)) @ current)
              + legs[1]*(2.0*KAPPA - 2.0)/VEV*current
              + legs[2]*pole_m[:, None, None]*(slash(p - k) @ CHIRAL))
    return spin_trace(lep, m_l, nu, vertex, vertex)


def vector_msq(lep: np.ndarray, nu: np.ndarray, k: np.ndarray, m_l: float,
               q_nu: float, chi: tuple[np.ndarray, np.ndarray] | None = None
               ) -> np.ndarray:
    """Spin-summed |M|^2 of M -> l nu Z' over (c g)^2.

    Args:
        lep: Charged lepton four-momenta.
        nu: Neutrino four-momenta.
        k: Vector four-momenta.
        m_l: Charged lepton mass in MeV.
        q_nu: Neutrino charge relative to the charged lepton.
        chi: Dark matter pair momenta for the decay tensor.

    Returns:
        Squared amplitude in units of the meson constant c^2.
    """
    m_v_sq = mdot(k, k)
    current = slash(lep + nu + k) @ CHIRAL
    left = ((slash(lep + k) + m_l*np.eye(4)) @ current
            / (2.0*mdot(lep, k) + m_v_sq)[:, None, None])
    right = (q_nu*current @ slash(nu + k)
             / (2.0*mdot(nu, k) + m_v_sq)[:, None, None])

    def trace(g1: np.ndarray, g2: np.ndarray) -> np.ndarray:
        """Trace with the vector index contracted with g1, g2."""
        return spin_trace(lep, m_l, nu, g1 @ left - right @ g1,
                          g2 @ left - right @ g2)

    t_metric = sum(METRIC[a]*trace(GAMMA[a], GAMMA[a]) for a in range(4))
    if chi is None:
        return -t_metric + trace(slash(k), slash(k))/m_v_sq
    m_chi_sq = mdot(chi[0], chi[0])
    norm = (4.0*m_v_sq + 8.0*m_chi_sq)/3.0
    t_chi = trace(slash(chi[0]), slash(chi[1]))
    return 4.0*(2.0*t_chi - 0.5*m_v_sq*t_metric)/norm


def meson_constant(parent: str) -> float:
    """Returns c^2/(2 M Gamma_tot) of the M -> l nu vertex."""
    mass, f_m, ckm, width = PARENTS[parent]
    return (G_F*f_m*ckm)**2/2.0/(2.0*mass*width)


def zprime_sample(parent: str, m_v: float, m_chi: float, g_v: float,
                  q_nu: float, n: int, rng: np.random.Generator,
                  m_l: float = M_MU, polarised: bool = True
                  ) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Rest-frame dark matter energies from M -> l nu Z'.

    Args:
        parent: Parent key, "pi" or "k".
        m_v: Vector mass in MeV.
        m_chi: Dark matter mass in MeV.
        g_v: Vector coupling to the charged lepton.
        q_nu: Neutrino charge relative to the charged lepton.
        n: Number of rest-frame events.
        rng: Random generator.
        m_l: Charged lepton mass in MeV.
        polarised: Keep the vector spin in its decay.

    Returns:
        Energies, probabilities, BR and its statistical error.
    """
    lep, nu, k, weight = three_body(PARENTS[parent][0], m_l, 0.0, m_v, n,
                                    rng)
    chi = two_body(k, np.full(n, m_chi), np.full(n, m_chi), rng)[:2]
    msq = vector_msq(lep, nu, k, m_l, q_nu, chi if polarised else None)
    rate = np.maximum(weight*msq*g_v**2*meson_constant(parent), 0.0)
    prob = np.tile(rate, 2)
    energy = np.concatenate([chi[0][:, 0], chi[1][:, 0]])
    return (energy, prob/prob.sum(), rate.mean(),
            rate.std()/np.sqrt(n))


def scalar_sample(parent: str, m_l: float, m_s: float, sin_theta: float,
                  n: int, rng: np.random.Generator
                  ) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Rest-frame scalar energies from M -> l nu S.

    Args:
        parent: Parent key, "pi" or "k".
        m_l: Charged lepton mass in MeV.
        m_s: Scalar mass in MeV.
        sin_theta: Sine of the Higgs mixing angle.
        n: Number of rest-frame events.
        rng: Random generator.

    Returns:
        Energies, rates per event, BR and its statistical error.
    """
    lep, nu, k, weight = three_body(PARENTS[parent][0], m_l, 0.0, m_s, n,
                                    rng)
    rate = (weight*scalar_msq(lep, nu, k, m_l)*sin_theta**2
            * meson_constant(parent))
    return k[:, 0], rate/n, rate.mean(), rate.std()/np.sqrt(n)


def scalar_ctau(m_s: float, sin_theta: float) -> float:
    """Proper decay length in m from S -> e+e- and mu+mu-."""
    width = sum(sin_theta**2*m_f**2*m_s/(8.0*np.pi*VEV**2)
                * (1.0 - 4.0*m_f**2/m_s**2)**1.5
                for m_f in (M_E, M_MU) if m_s > 2.0*m_f)
    return METER_BY_MEV/width if width > 0.0 else np.inf


def kaon_to_pi_scalar(m_s: float, sin_theta: float) -> tuple[float, float]:
    """Returns E_S in the kaon frame and BR(K+ -> pi+ S)."""
    amp = 3.0*VTD_VTS*M_TOP**2*M_K**2/(32.0*np.pi**2*VEV**3)
    lam = kallen(1.0, m_s**2/M_K**2, M_PI**2/M_K**2)
    width = sin_theta**2*amp**2*np.sqrt(lam)/(16.0*np.pi*M_K)
    return (M_K**2 + m_s**2 - M_PI**2)/(2.0*M_K), width/KAON_WIDTH


def forced(mes: np.ndarray, m0: float, m_x: float, e_star: np.ndarray,
           det: np.ndarray) -> tuple[np.ndarray, ...]:
    """Forces an isotropic decay product at a detector point.

    Args:
        mes: Meson table with the columns of flux_io.flux_array.
        m0: Meson mass in MeV.
        m_x: Product mass in MeV.
        e_star: Product energy in the meson frame, per meson.
        det: Detector point(s) in m, shape (3,) or (n, 3).

    Returns:
        Lab energy, dN/dOmega, distance [m], z cosine per root.
    """
    energy_m = np.sqrt(mes[:, 0]**2 + m0**2)
    gamma, beta = energy_m/m0, mes[:, 0]/energy_m
    sight = det - mes[:, 2:5]
    dist = np.sqrt(np.sum(sight**2, axis=1))
    axis = np.column_stack([np.sin(mes[:, 1])*np.cos(mes[:, 7]),
                            np.sin(mes[:, 1])*np.sin(mes[:, 7]),
                            np.cos(mes[:, 1])])
    b_cos = beta*np.sum(axis*sight, axis=1)/dist
    a, den = e_star/gamma, 1.0 - b_cos**2
    disc = a**2 - m_x**2*den
    root = np.sqrt(np.maximum(disc, 0.0))
    p_lab = np.array([a*b_cos + root, a*b_cos - root])/den
    e_lab = a + b_cos*p_lab
    slope = np.abs(p_lab - b_cos*e_lab)
    good = (disc > 0.0) & (p_lab > 0.0) & (e_lab > 0.0) & (slope > 0.0)
    p_star = np.sqrt(e_star**2 - m_x**2)
    with np.errstate(divide="ignore", invalid="ignore"):
        dn_domega = np.where(
            good, p_lab**2/(4.0*np.pi*p_star*gamma*slope), 0.0)
    return e_lab, dn_domega, dist, sight[:, 2]/dist


def detector_flux(mes: np.ndarray, m0: float, m_x: float,
                  e_star: np.ndarray, prob: np.ndarray, per_decay: float,
                  ctau: float, rng: np.random.Generator,
                  n_throw: int = N_THROW) -> dict[str, np.ndarray]:
    """Energy spectra at the SBND front face.

    Args:
        mes: Meson table with the columns of flux_io.flux_array.
        m0: Meson mass in MeV.
        m_x: Mass of the particle reaching SBND in MeV.
        e_star: Sample of energies in the meson frame.
        prob: Probabilities of the entries of e_star.
        per_decay: Particles per meson decay, BR included.
        ctau: Proper decay length in m, np.inf if stable.
        rng: Random generator.
        n_throw: Decays thrown per meson.

    Returns:
        Spectra, totals and their statistical errors.
    """
    n_mes, n_bin = mes.shape[0], EDGES.size - 1
    sums = {"centre": np.zeros(n_mes*n_bin), "face": np.zeros(n_mes*n_bin)}
    row = np.tile(np.arange(n_mes)*n_bin, 2)
    centre = np.array([0.0, 0.0, SBND_Z])
    for _ in range(n_throw):
        pick = e_star[rng.choice(e_star.size, size=n_mes, p=prob)]
        spot = np.column_stack(
            [rng.uniform(-0.5, 0.5, (n_mes, 2))*FACE_SIDE,
             np.full(n_mes, SBND_Z)])
        for name, det in (("centre", centre), ("face", spot)):
            e_lab, dn_domega, dist, cos_z = forced(mes, m0, m_x, pick, det)
            area = (1.0e-4 if name == "centre"
                    else FACE_SIDE**2*np.abs(cos_z))
            p_lab = np.sqrt(np.maximum(e_lab**2 - m_x**2, 1.0e-30))
            alive = np.exp(-dist*m_x/(p_lab*ctau))
            weight = (mes[:, 6]*per_decay/n_throw*area/dist**2
                      * dn_domega*alive)
            column = np.clip(np.digitize(e_lab, EDGES) - 1, 0, n_bin - 1)
            sums[name] += np.bincount((row + column.ravel()),
                                      weights=weight.ravel(),
                                      minlength=n_mes*n_bin)
    out = {}
    for name, label in (("centre", "flux_centre"), ("face", "face_rate")):
        table = sums[name].reshape(n_mes, n_bin)
        out[label] = table.sum(axis=0)
        out[f"{label}_err"] = np.sqrt(np.sum(table**2, axis=0))
        out[f"{label}_total"] = table.sum()
        out[f"{label}_total_err"] = np.sqrt(np.sum(table.sum(axis=1)**2))
    return out


def save_benchmark(name: str, flux: dict[str, np.ndarray], species: str,
                   current: float, info: dict) -> None:
    """Writes one benchmark flux to data/bsm_data.

    Args:
        name: File name without extension.
        flux: Output of detector_flux.
        species: Parent meson species.
        current: Horn current in kA.
        info: Model parameters to store with the flux.
    """
    os.makedirs(OUT_DIR, exist_ok=True)
    np.savez(
        os.path.join(OUT_DIR, f"{name}.npz"), energy_edges_mev=EDGES,
        flux_centre_per_pot_cm2=flux["flux_centre"],
        flux_centre_err=flux["flux_centre_err"],
        flux_centre_total_per_pot_cm2=flux["flux_centre_total"],
        flux_centre_total_err=flux["flux_centre_total_err"],
        face_rate_per_pot=flux["face_rate"],
        face_rate_err=flux["face_rate_err"],
        face_rate_total_per_pot=flux["face_rate_total"],
        face_rate_total_err=flux["face_rate_total_err"],
        parent_species=species, horn_current_ka=current,
        detector_distance_m=SBND_Z, detector_face_width_m=FACE_SIDE,
        detector_face_height_m=FACE_SIDE, **info)
    print(f"{name}: BR = {info['branching_ratio']:.3e}, centre = "
          f"{flux['flux_centre_total']:.3e} +- "
          f"{flux['flux_centre_total_err']:.1e} /POT/cm2, face = "
          f"{flux['face_rate_total']:.3e} +- "
          f"{flux['face_rate_total_err']:.1e} /POT")


def run_zprime(species: str, current: float, mes: np.ndarray,
               rng: np.random.Generator) -> None:
    """Computes and saves the dark matter fluxes of model A.

    Args:
        species: Parent meson species.
        current: Horn current in kA.
        mes: Meson table with the columns of flux_io.flux_array.
        rng: Random generator.
    """
    parent = species.split("_")[0]
    for m_v in ZP_MASSES:
        if m_v >= PARENTS[parent][0] - M_MU:
            continue
        m_chi = m_v/3.0
        ratio = m_chi**2/m_v**2
        width_chi = (ZP_ALPHA_D*m_v/3.0*(1.0 + 2.0*ratio)
                     * np.sqrt(1.0 - 4.0*ratio))
        width_nu = ZP_NU_CHARGE**2*ZP_COUPLING**2*m_v/(12.0*np.pi)
        br_chi = width_chi/(width_chi + width_nu)
        e_star, prob, br, br_err = zprime_sample(
            parent, m_v, m_chi, ZP_COUPLING, ZP_NU_CHARGE, N_REST, rng)
        flux = detector_flux(mes, PARENTS[parent][0], m_chi, e_star, prob,
                             2.0*br*br_chi, np.inf, rng)
        info = {"model": "zprime_dark_matter",
                "particle": "chi_plus_chibar",
                "channel_names": np.array([f"{parent} -> mu nu Zp"]),
                "channel_br": np.array([br]), "branching_ratio": br,
                "branching_ratio_err": br_err, "m_zprime_mev": m_v,
                "m_chi_mev": m_chi, "g_lepton": ZP_COUPLING,
                "nu_charge": ZP_NU_CHARGE, "alpha_dark": ZP_ALPHA_D,
                "br_zprime_to_chichi": br_chi,
                "zprime_width_mev": width_chi + width_nu,
                "particles_per_decay": 2.0}
        save_benchmark(f"zprime_dm_{species}_{current:+.0f}kA_"
                       f"mzp{m_v:.0f}MeV", flux, species, current, info)


def run_scalar(species: str, current: float, mes: np.ndarray,
               rng: np.random.Generator) -> None:
    """Computes and saves the scalar fluxes of model B.

    Args:
        species: Parent meson species.
        current: Horn current in kA.
        mes: Meson table with the columns of flux_io.flux_array.
        rng: Random generator.
    """
    parent = species.split("_")[0]
    for m_s in S_MASSES:
        ctau = scalar_ctau(m_s, S_SIN_THETA)
        if parent == "k":
            e_k, br = kaon_to_pi_scalar(m_s, S_SIN_THETA)
            names, brs, br_err = ["k -> pi S"], [br], 0.0
            e_star, prob = np.array([e_k]), np.array([1.0])
        else:
            leptons = [(tag, m_l) for tag, m_l in (("e", M_E), ("mu", M_MU))
                       if m_s < M_PI - m_l]
            parts = [scalar_sample(parent, m_l, m_s, S_SIN_THETA, N_REST,
                                   rng) for _, m_l in leptons]
            names = [f"pi -> {tag} nu S" for tag, _ in leptons]
            brs = [part[2] for part in parts]
            br_err = np.sqrt(sum(part[3]**2 for part in parts))
            e_star = np.concatenate([part[0] for part in parts])
            prob = np.concatenate([part[1] for part in parts])
            prob = np.maximum(prob, 0.0)/np.maximum(prob, 0.0).sum()
        flux = detector_flux(mes, PARENTS[parent][0], m_s, e_star, prob,
                             sum(brs), ctau, rng)
        info = {"model": "higgs_portal_scalar", "particle": "S",
                "channel_names": np.array(names),
                "channel_br": np.array(brs), "branching_ratio": sum(brs),
                "branching_ratio_err": br_err, "m_scalar_mev": m_s,
                "sin_theta": S_SIN_THETA, "ctau_m": ctau,
                "width_mev": METER_BY_MEV/ctau, "particles_per_decay": 1.0}
        save_benchmark(f"higgs_scalar_{species}_{current:+.0f}kA_"
                       f"ms{m_s:.0f}MeV", flux, species, current, info)


def load_mesons(species: str, current: float) -> np.ndarray:
    """Loads the decays in the beamline of one meson sample."""
    path = flux_io.pickle_path(species, current)
    return flux_io.flux_array(flux_io.load_meson_pickle(path))


def brute_force_face(mes: np.ndarray, m0: float, m_x: float,
                     e_star: np.ndarray, prob: np.ndarray, n_throw: int,
                     rng: np.random.Generator) -> tuple[float, float]:
    """Counts products hitting the SBND face without forcing.

    Args:
        mes: Meson table with the columns of flux_io.flux_array.
        m0: Meson mass in MeV.
        m_x: Product mass in MeV.
        e_star: Sample of energies in the meson frame.
        prob: Probabilities of the entries of e_star.
        n_throw: Decays thrown per meson.
        rng: Random generator.

    Returns:
        Products per meson-weight unit and statistical error.
    """
    n_mes = mes.shape[0]
    axis = np.column_stack([np.sin(mes[:, 1])*np.cos(mes[:, 7]),
                            np.sin(mes[:, 1])*np.sin(mes[:, 7]),
                            np.cos(mes[:, 1])])
    parent = np.column_stack([np.sqrt(mes[:, 0]**2 + m0**2),
                              mes[:, 0, None]*axis])
    total, total_sq = 0.0, 0.0
    for _ in range(n_throw):
        energy = e_star[rng.choice(e_star.size, size=n_mes, p=prob)]
        vec = (np.sqrt(energy**2 - m_x**2)[:, None]
               * isotropic(n_mes, rng))
        lab = boost(np.column_stack([energy, vec]), parent)
        reach = (SBND_Z - mes[:, 4])/np.where(lab[:, 3] > 0.0, lab[:, 3],
                                              np.nan)
        hit_x = mes[:, 2] + reach*lab[:, 1]
        hit_y = mes[:, 3] + reach*lab[:, 2]
        hit = ((np.abs(hit_x) < FACE_SIDE/2.0)
               & (np.abs(hit_y) < FACE_SIDE/2.0))
        total += mes[hit, 6].sum()/n_throw
        total_sq += np.sum(mes[hit, 6]**2)/n_throw**2
    return total, np.sqrt(total_sq)


def validate(rng: np.random.Generator) -> None:
    """Prints the kinematics, rate and flux checks."""
    lep, nu, k, weight = three_body(M_PI, M_MU, 0.0, 20.0, 400000, rng)

    def dalitz(s: float) -> float:
        """Width of the Dalitz plot at pair mass squared s."""
        return np.sqrt(kallen(s, M_MU**2, 0.0)
                       * kallen(M_PI**2, s, 400.0))/s

    exact = integrate.quad(dalitz, M_MU**2, (M_PI - 20.0)**2)[0]
    exact /= 128.0*np.pi**3*M_PI**2
    print(f"1 phase space pi->mu nu X(20): MC {weight.mean():.5e} +- "
          f"{weight.std()/weight.size**0.5:.1e}, Dalitz area {exact:.5e}")
    x = 400.0/M_PI**2
    closed = M_PI**2/(256.0*np.pi**3)*(1.0 - x**2 + 2.0*x*np.log(x))
    flat = three_body(M_PI, 0.0, 0.0, 20.0, 400000, rng)[3]
    print(f"1 phase space pi->0 0 X(20): MC {flat.mean():.5e} +- "
          f"{flat.std()/flat.size**0.5:.1e}, closed form {closed:.5e}")
    lab = np.tile([2000.0, 30.0, -40.0, 1900.0], (lep.shape[0], 1))
    lab[:, 0] = np.sqrt(M_PI**2 + np.sum(lab[:, 1:]**2, axis=1))
    total = sum(boost(part, lab) for part in (lep, nu, k))
    print(f"1 max |sum p - P| after boost: "
          f"{np.abs(total - lab).max():.2e} MeV, max mass shift "
          f"{np.abs(np.sqrt(mdot(k, k)) - 20.0).max():.2e} MeV")

    rest = np.zeros((4, 4))
    rest[:, 0] = M_PI
    mu, nu0, p_two = two_body(rest, np.full(4, M_MU), np.zeros(4), rng)
    current = slash(mu + nu0) @ CHIRAL
    width = (spin_trace(mu, M_MU, nu0, current, current)[0]
             * (G_F*F_PI*V_UD)**2/2.0*p_two[0]/(8.0*np.pi*M_PI**2))
    formula = ((G_F*F_PI*V_UD*M_MU)**2*M_PI*(1.0 - M_MU**2/M_PI**2)**2
               / (8.0*np.pi))
    print(f"2 SM pi->mu nu width: traces {width:.5e}, formula "
          f"{formula:.5e}, BR*width_tot {BR_PI_MUNU*PION_WIDTH:.5e} MeV")

    lep, nu, k, weight = three_body(M_K, M_MU, 0.0, 50.0, 20000, rng)
    m12, m23 = mdot(nu + k, nu + k), mdot(lep + k, lep + k)
    mk2, mm2, mv2 = M_K**2, M_MU**2, 2500.0
    d23 = m23 - mm2
    krnjaic_v = 4.0*mm2*(
        2.0 + (m12 + 2.0*mm2 - 2.0*mk2)/d23
        - (mk2 - mm2)*(mv2 + 2.0*mm2)/d23**2
        + 2.0*((mk2 - mm2)**2 + mv2*mm2)/(m12*d23)
        - mv2*(mk2 - mm2)/m12**2 + (m23 + mm2 - 2.0*mk2)/m12)
    ratio = vector_msq(lep, nu, k, M_MU, 1.0)/2.0/krnjaic_v
    print(f"2 vector |M|^2 / Krnjaic (A.7): min {ratio.min():.6f}, "
          f"max {ratio.max():.6f}")
    left = slash(lep + k) + M_MU*np.eye(4)
    left = left @ slash(lep + nu + k) @ CHIRAL
    left /= (2.0*mdot(lep, k) + mv2)[:, None, None]
    right = slash(lep + nu + k) @ CHIRAL @ slash(nu + k)
    right /= (2.0*mdot(nu, k) + mv2)[:, None, None]
    ward = slash(k) @ left - right @ slash(k)
    ward = spin_trace(lep, M_MU, nu, ward, ward)
    print(f"2 Ward identity k.M k.M* / |M|^2, max: "
          f"{np.abs(ward/vector_msq(lep, nu, k, M_MU, 1.0)).max():.1e}")
    krnjaic_s = 4.0*mm2/(2.0*mm2*d23**2)*(
        mk2*(m23 + mm2)**2 - m23*((m23 + mm2)**2 + m12*d23)
        + mv2*(m23**2 - mm2*mk2))
    ratio = (scalar_msq(lep, nu, k, M_MU, (1.0, 0.0, 0.0))/2.0
             / (krnjaic_s*(M_MU/VEV)**2))
    print(f"2 scalar lepton-leg |M|^2 / Krnjaic (A.8): min "
          f"{ratio.min():.6f}, max {ratio.max():.6f}")
    for m_v in (1.0, 10.0, 100.0, 300.0):
        br = zprime_sample("k", m_v, m_v/3.0, 1.0, 1.0, 200000, rng)[2:]
        print(f"2 BR(K->mu nu V)/alpha_V at {m_v:.0f} MeV: "
              f"{4.0*np.pi*br[0]:.4f} +- {4.0*np.pi*br[1]:.4f}")
    for m_v in ZP_MASSES:
        pol = zprime_sample("pi", m_v, m_v/3.0, 1.0, 1.0, 200000, rng)[2:]
        iso = zprime_sample("pi", m_v, m_v/3.0, 1.0, 1.0, 200000, rng,
                            polarised=False)[2:]
        br_e = [zprime_sample("pi", m_v, m_v/3.0, 1.0, q, 200000, rng,
                              m_l=M_E, polarised=False)[2]
                for q in (1.0, 0.0)]
        br_0 = zprime_sample("pi", m_v, m_v/3.0, 1.0, 0.0, 200000, rng,
                             polarised=False)[2]
        print(f"2 BR(pi->l nu V)/g^2 at {m_v:.0f} MeV: mu decay tensor "
              f"{pol[0]:.4e} +- {pol[1]:.1e}, mu pol sum {iso[0]:.4e} +- "
              f"{iso[1]:.1e}, mu q_nu=0 {br_0:.4e}, e q_nu=1 "
              f"{br_e[0]:.4e}, e q_nu=0 {br_e[1]:.4e}")
    for parent, m_l, tag in (("pi", M_MU, "pi mu"), ("pi", M_E, "pi e"),
                             ("k", M_MU, "K mu"), ("k", M_E, "K e")):
        br = zprime_sample(parent, 17.0, 5.0, 6.0e-5*0.30282, 0.0, 200000,
                           rng, m_l=m_l, polarised=False)[2]
        print(f"2 Dutta et al. Table II point, {tag}: BR = {br:.3e}")
    for e_soft in (1.0, 0.1, 0.01):
        rest = np.zeros((2000, 4))
        rest[:, 0] = M_PI
        m_soft = 0.5*e_soft
        s_pair = M_PI**2 + m_soft**2 - 2.0*M_PI*e_soft
        pair, k = two_body(rest, np.full(2000, np.sqrt(s_pair)),
                           np.full(2000, m_soft), rng)[:2]
        lep, nu = two_body(pair, np.full(2000, M_MU), np.zeros(2000),
                           rng)[:2]
        d_l = 2.0*mdot(lep, k) + m_soft**2
        d_n = 2.0*mdot(nu, k) + m_soft**2
        d_p = m_soft**2 - 2.0*M_PI*e_soft
        born = 8.0*M_MU**2*mdot(lep, nu)
        jay = 2.0*lep/d_l[:, None] - 2.0*nu/d_n[:, None]
        eik_v = born*(-mdot(jay, jay) + mdot(jay, k)**2/m_soft**2)
        eik_s = born*(2.0*M_MU**2/d_l + (1.0 + KAPPA)*M_PI**2/d_p
                      + KAPPA - 1.0)**2/VEV**2
        r_v = vector_msq(lep, nu, k, M_MU, 1.0)/eik_v
        r_s = scalar_msq(lep, nu, k, M_MU)/eik_s
        print(f"2 soft limit E_X = {e_soft} MeV: vector/eikonal "
              f"{r_v.mean():.4f} (rms {r_v.std():.4f}), scalar/pole "
              f"{r_s.mean():.4f} (rms {r_s.std():.4f})")
    for m_s in (1.0,) + S_MASSES:
        e_part = scalar_sample("pi", M_E, m_s, 1.0, 400000, rng)
        x = m_s**2/M_PI**2
        shape = (1.0 - 8.0*x + x**2)*(1.0 - x**2) - 12.0*x**2*np.log(x)
        dawson = (np.sqrt(2.0)*G_F*M_PI**4*(1.0 - KAPPA)**2*shape
                  / (96.0*np.pi**2*M_MU**2*(1.0 - M_MU**2/M_PI**2)**2))
        sm_width = ((G_F*F_PI*V_UD*M_MU)**2*M_PI
                    * (1.0 - M_MU**2/M_PI**2)**2/(8.0*np.pi))
        line = (f"2 BR(pi->e nu S)/sin^2 at {m_s:.0f} MeV: MC "
                f"{e_part[2]:.4e} +- {e_part[3]:.1e}, closed form "
                f"{dawson*sm_width/PION_WIDTH:.4e}")
        if m_s < M_PI - M_MU:
            parts = [scalar_sample("pi", M_MU, m_s, 1.0, 400000, rng)[2]]
            for legs in ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0),
                         (0.0, 0.0, 1.0)):
                lep, nu, k, weight = three_body(M_PI, M_MU, 0.0, m_s,
                                                400000, rng)
                parts.append((weight*scalar_msq(lep, nu, k, M_MU, legs)
                              ).mean()*meson_constant("pi"))
            line += (f"; mu mode all {parts[0]:.4e}, lepton only "
                     f"{parts[1]:.4e}, contact only {parts[2]:.4e}, "
                     f"meson only {parts[3]:.4e}")
        print(line)
        print(f"2 BR(K->pi S)/sin^2 at {m_s:.0f} MeV: "
              f"{kaon_to_pi_scalar(m_s, 1.0)[1]:.4e}, ctau*sin^2 = "
              f"{scalar_ctau(m_s, 1.0):.4e} m")

    mes = load_mesons("pi_plus", 174.0)
    e_nu = np.array([(M_PI**2 - M_MU**2)/(2.0*M_PI)])
    one = np.array([1.0])
    flux = detector_flux(mes, M_PI, 0.0, e_nu, one, BR_PI_MUNU, np.inf,
                         rng)
    ref_e, ref_w = forced_detector_flux(mes, m_meson=M_PI, br=BR_PI_MUNU,
                                        det_pos=(0.0, 0.0, SBND_Z))
    ref_hist = np.histogram(ref_e, bins=EDGES, weights=ref_w)[0]
    print(f"3 nu flux at (0,0,110): this code "
          f"{flux['flux_centre_total']:.6e}, forced_detector_flux "
          f"{ref_w.sum():.6e} /POT/cm2, max bin difference "
          f"{np.abs(flux['flux_centre'] - ref_hist).max():.1e}")
    brute = brute_force_face(mes, M_PI, 0.0, e_nu, one, 400, rng)
    print(f"3 nu through face: forced {flux['face_rate_total']:.4e}, "
          f"brute force {BR_PI_MUNU*brute[0]:.4e} +- "
          f"{BR_PI_MUNU*brute[1]:.1e} /POT, centre flux x area "
          f"{flux['flux_centre_total']*FACE_SIDE**2*1.0e4:.4e}")
    samples = (("S 100 MeV", 100.0, scalar_sample("pi", M_E, 100.0, 1.0,
                                                  N_REST, rng)),
               ("chi 20/3 MeV", 20.0/3.0, zprime_sample(
                   "pi", 20.0, 20.0/3.0, 1.0, 1.0, N_REST, rng)))
    for tag, m_x, (e_star, prob, _, _) in samples:
        prob = np.maximum(prob, 0.0)/np.maximum(prob, 0.0).sum()
        flux = detector_flux(mes, M_PI, m_x, e_star, prob, 1.0, np.inf,
                             rng)
        brute = brute_force_face(mes, M_PI, m_x, e_star, prob, 400, rng)
        print(f"3 {tag} through face per unit BR: forced "
              f"{flux['face_rate_total']:.4e}, brute force "
              f"{brute[0]:.4e} +- {brute[1]:.1e}")

    e_star, prob = samples[0][2][:2]
    prob = np.maximum(prob, 0.0)/np.maximum(prob, 0.0).sum()
    br_one = samples[0][2][2]
    for sin_theta in (1.0e-5, 1.0e-4, 1.0e-2, 3.0e-2, 0.1, 0.3, 1.0):
        ctau = scalar_ctau(100.0, sin_theta)
        flux = detector_flux(mes, M_PI, 100.0, e_star, prob,
                             br_one*sin_theta**2, ctau,
                             np.random.default_rng(7))
        print(f"4 S(100 MeV) sin_theta = {sin_theta:.0e}: ctau = "
              f"{ctau:.3e} m, face rate / sin^2 = "
              f"{flux['face_rate_total']/sin_theta**2:.4e} /POT")


def main() -> None:
    """Computes the benchmark fluxes or runs the checks."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validate", action="store_true",
                        help="print the validation numbers and exit")
    parser.add_argument("--seed", type=int, default=2024)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    if args.validate:
        validate(rng)
        return
    for species, current in SAMPLES:
        mes = load_mesons(species, current)
        run_zprime(species, current, mes, rng)
        run_scalar(species, current, mes, rng)


if __name__ == "__main__":
    main()
