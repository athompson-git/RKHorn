"""
Reproducible comparison of the RKHorn neutrino flux at MiniBooNE with the published flux
(MiniBooNE, PRD 79, 072002; digitized curves in data/public_digitized_fluxes).

    python tests/validate_flux.py                 # default statistics, fixed seed
    python tests/validate_flux.py --n 200000      # more statistics
    python tests/validate_flux.py --plot          # also draw the comparison
    python tests/validate_flux.py --n 500000 --theta-max 1.5708 --save   # regenerate the files in data/

Neutrino mode:      nu_mu    from pi+ (and, separately, K+), horn current +174 kA
Antineutrino mode:  nu_mubar from pi-,                       horn current -174 kA

Everything is in nu / POT / cm^2. The simulated flux is histogrammed in 50 MeV bins and the
published curve is integrated over the same bins, so the two are compared bin for bin.
"""

import argparse
import contextlib
import io
import os
import sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

PUBLISHED_TOTAL = {"numu": 5.19e-10, "numubar": 2.93e-10}

DATA_FILES = {"pi_plus": "bnb_focused_piplus_flux_rkhorn_1POT_BetheBloch.txt",
              "pi_minus": "bnb_focused_piminus_flux_rkhorn_1POT_BetheBloch.txt",
              "k_plus": "bnb_focused_Kplus_flux_rkhorn_1POT_BetheBloch.txt"}

BIN_EDGES = np.arange(0.0, 3050.0, 50.0)
BIN_WIDTH_GEV = 0.05
BANDS = [(0, 500), (500, 1000), (1000, 1500), (1500, 2000), (2000, 3000)]


def run_species(args):
    species, horn_current, n_samples, seed, step_ns, theta_max, absorption, save = args
    from rkhorn import ChargedPionFluxMiniBooNE, forced_detector_flux
    from constants import M_P, M_PI, M_K, BR_PI_MUNU, BR_K_MUNU

    np.random.seed(seed)
    e_proton = np.sqrt(8.89e3**2 + M_P**2)
    sim = ChargedPionFluxMiniBooNE(proton_energy=e_proton, n_samples=n_samples, n_pot=1,
                                   horn_current=horn_current, meson_species=species,
                                   step_size_ns=step_ns, absorption=absorption)
    with contextlib.redirect_stderr(io.StringIO()):
        sim.simulate(discard_zeros=True, theta_max=theta_max)

    mesons = sim.decayed_flux_array()
    if save:
        np.savetxt(os.path.join(ROOT, "data", DATA_FILES[species]), mesons, fmt="%.6e",
                   header="horn-focused {} per POT, horn current {} kA, {} mesons sampled; columns: ".format(
                              species, horn_current, n_samples)
                          + "p[MeV] theta[rad] decay_x[m] decay_y[m] decay_z[m] decay_in_volume weight[1/POT] phi[rad]")
    if species.startswith("k"):
        e_nu, w_nu = forced_detector_flux(mesons, m_meson=M_K, br=BR_K_MUNU)
    else:
        e_nu, w_nu = forced_detector_flux(mesons, m_meson=M_PI, br=BR_PI_MUNU)

    hist, _ = np.histogram(e_nu, weights=w_nu, bins=BIN_EDGES)
    hist_sq, _ = np.histogram(e_nu, weights=w_nu**2, bins=BIN_EDGES)
    produced = np.sum(sim.meson_flux_pre_horn[:,2])
    stopped = np.sum(mesons[mesons[:,0] <= 0.0, 6])
    return species, hist, hist_sq, produced, np.sum(mesons[:,6]), stopped


def published_per_bin(filename):
    data = np.genfromtxt(os.path.join(ROOT, "data", "public_digitized_fluxes", filename))
    data = data[np.argsort(data[:,0])]
    out = np.zeros(BIN_EDGES.shape[0] - 1)
    for i in range(out.shape[0]):
        e = np.linspace(BIN_EDGES[i], BIN_EDGES[i+1], 51)
        f = np.interp(e, data[:,0], data[:,1], left=0.0, right=0.0)
        out[i] = np.sum(0.5*(f[1:] + f[:-1])*np.diff(e)) * 1e-3
    return out


def report(label, sim, sim_sq, ref, published_total):
    err = np.sqrt(sim_sq)
    print("\n{}".format(label))
    print("  integrated 0-3 GeV:  RKHorn {:.3e} +- {:.1e}   digitized {:.3e}   paper total {:.3e}"
          .format(np.sum(sim), np.sqrt(np.sum(sim_sq)), np.sum(ref), published_total))
    print("  RKHorn / digitized = {:.3f}".format(np.sum(sim)/np.sum(ref)))
    i_sim, i_ref = np.argmax(sim), np.argmax(ref)
    print("  peak bin: RKHorn {:.0f}-{:.0f} MeV, {:.2e} /GeV   digitized {:.0f}-{:.0f} MeV, {:.2e} /GeV"
          .format(BIN_EDGES[i_sim], BIN_EDGES[i_sim+1], sim[i_sim]/BIN_WIDTH_GEV,
                  BIN_EDGES[i_ref], BIN_EDGES[i_ref+1], ref[i_ref]/BIN_WIDTH_GEV))
    print("  ratio by band [MeV]:")
    for lo, hi in BANDS:
        sel = (BIN_EDGES[:-1] >= lo) & (BIN_EDGES[1:] <= hi)
        ratio = np.sum(sim[sel])/np.sum(ref[sel])
        print("    {:5d}-{:5d}: {:.3f} +- {:.3f}".format(lo, hi, ratio, np.sqrt(np.sum(sim_sq[sel]))/np.sum(ref[sel])))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--n", type=int, default=40000, help="mesons sampled per species")
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--step-ns", type=float, default=0.1, help="transport time step in ns")
    parser.add_argument("--theta-max", type=float, default=0.6,
                        help="largest production angle sampled [rad]; wider angles do not reach the detector")
    parser.add_argument("--absorption", default="reaction", choices=["reaction", "inelastic", "none"],
                        help="hadronic absorption of mesons in the target and horn")
    parser.add_argument("--save", action="store_true", help="overwrite the meson flux files in data/")
    parser.add_argument("--plot", action="store_true")
    opts = parser.parse_args()

    absorption = None if opts.absorption == "none" else opts.absorption
    jobs = [("pi_plus", 174.0, opts.n, opts.seed, opts.step_ns, opts.theta_max, absorption, opts.save),
            ("k_plus", 174.0, opts.n, opts.seed + 1, opts.step_ns, opts.theta_max, absorption, opts.save),
            ("pi_minus", -174.0, opts.n, opts.seed + 2, opts.step_ns, opts.theta_max, absorption, opts.save)]
    with Pool(3) as pool:
        results = {r[0]: r[1:] for r in pool.map(run_species, jobs)}

    print("n = {} per species, seed = {}, step = {} ns, theta_max = {} rad, absorption = {}"
          .format(opts.n, opts.seed, opts.step_ns, opts.theta_max, opts.absorption))
    print("\nmesons per POT   produced(sampled range)   decayed in volume   of which stopped")
    for species in ["pi_plus", "k_plus", "pi_minus"]:
        _, _, produced, decayed, stopped = results[species]
        print("  {:9s}      {:.4f}                  {:.4f}              {:.4f}".format(species, produced, decayed, stopped))

    ref_numu = published_per_bin("numu_diff_flux_by_POT_GeV_cm2.txt")
    ref_numubar = published_per_bin("numubar_diff_flux_by_POT_GeV_cm2.txt")

    numu, numu_sq = results["pi_plus"][0], results["pi_plus"][1]
    report("nu_mu from pi+, neutrino mode", numu, numu_sq, ref_numu, 0.9672*PUBLISHED_TOTAL["numu"])
    report("nu_mubar from pi-, antineutrino mode", results["pi_minus"][0], results["pi_minus"][1],
           ref_numubar, 0.984*PUBLISHED_TOTAL["numubar"])

    kaon, kaon_sq = results["k_plus"][0], results["k_plus"][1]
    print("\nnu_mu from K+ -> mu nu, neutrino mode")
    print("  integrated 0-3 GeV: RKHorn {:.3e} +- {:.1e}   paper {:.3e} (2.65% of the nu_mu total, all energies)"
          .format(np.sum(kaon), np.sqrt(np.sum(kaon_sq)), 0.0265*PUBLISHED_TOTAL["numu"]))
    sel = (BIN_EDGES[:-1] >= 1000) & (BIN_EDGES[1:] <= 3000)
    print("  mean flux 1-3 GeV:  RKHorn {:.2e} /GeV   paper roughly 4e-12 /GeV (read off Fig. 29)"
          .format(np.sum(kaon[sel])/2.0))

    if opts.plot:
        import matplotlib.pyplot as plt
        centers = 0.5*(BIN_EDGES[1:] + BIN_EDGES[:-1])
        for sim, ref, color, name in [(numu, ref_numu, "b", r"$\nu_\mu$ from $\pi^+$"),
                                      (results["pi_minus"][0], ref_numubar, "r", r"$\bar{\nu}_\mu$ from $\pi^-$")]:
            plt.hist(centers, weights=sim/BIN_WIDTH_GEV, bins=BIN_EDGES, histtype="step", color=color, label="RKHorn " + name)
            plt.hist(centers, weights=ref/BIN_WIDTH_GEV, bins=BIN_EDGES, histtype="step", color=color, linestyle="dashed",
                     label="MiniBooNE " + name)
        plt.yscale("log")
        plt.ylim((1e-13, 1e-9))
        plt.xlabel(r"$E_\nu$ [MeV]")
        plt.ylabel(r"$\Phi$ [$\nu$ / POT / GeV / cm$^2$]")
        plt.legend()
        plt.show()


if __name__ == "__main__":
    main()
