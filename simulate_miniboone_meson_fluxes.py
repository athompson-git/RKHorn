#!/usr/bin/env python3
"""Simulates BNB meson samples and saves them in data/."""

import argparse
import contextlib
import io
from multiprocessing import Pool
import os

import numpy as np

from constants import M_P
import flux_io
from rkhorn import ChargedPionFluxMiniBooNE

P_PROTON = 8.89e3
NOMINAL_CURRENT = {"pi_plus": 174.0, "k_plus": 174.0, "pi_minus": -174.0}
TXT_FILES = {
    "pi_plus": "bnb_focused_piplus_flux_rkhorn_1POT_BetheBloch.txt",
    "pi_minus": "bnb_focused_piminus_flux_rkhorn_1POT_BetheBloch.txt",
    "k_plus": "bnb_focused_Kplus_flux_rkhorn_1POT_BetheBloch.txt",
}
TRACKS_FILE = os.path.join(flux_io.DATA_DIR, "tracks", "horn_3d_tracks.npz")

Job = tuple[str, float, int, int, float, bool]


def run_chunk(job: Job) -> dict:
    """Simulates one chunk of a meson sample.

    Args:
        job: Species, current [kA], mesons, seed, theta max, keep all.

    Returns:
        The meson record of the chunk.
    """
    species, horn_current, n_samples, seed, theta_max, keep_all = job
    np.random.seed(seed)
    sim = ChargedPionFluxMiniBooNE(
        proton_energy=np.sqrt(P_PROTON**2 + M_P**2), n_samples=n_samples,
        n_pot=1, horn_current=horn_current, meson_species=species)
    with contextlib.redirect_stderr(io.StringIO()):
        sim.simulate(discard_zeros=not keep_all, theta_max=theta_max)
    return sim.meson_record()


def simulate_samples(configs: list[tuple[str, float]], n_samples: int,
                     n_chunks: int, seed: int, theta_max: float,
                     keep_all: bool) -> list[dict]:
    """Simulates several meson samples in parallel.

    Args:
        configs: Species and horn current [kA] of each sample.
        n_samples: Mesons per sample.
        n_chunks: Parallel chunks per sample.
        seed: First random seed.
        theta_max: Largest production angle sampled, in rad.
        keep_all: Keep mesons that decay outside the beamline.

    Returns:
        One merged meson record per sample.
    """
    species_order = sorted(NOMINAL_CURRENT)
    jobs = [(species, current, n_samples//n_chunks,
             seed + 1000*species_order.index(species) + chunk, theta_max,
             keep_all)
            for species, current in configs for chunk in range(n_chunks)]
    with Pool(min(len(jobs), os.cpu_count())) as pool:
        chunks = pool.map(run_chunk, jobs)
    return [flux_io.merge_records(chunks[k*n_chunks:(k + 1)*n_chunks])
            for k in range(len(configs))]


def simulate_tracks(n_tracks: int, horn_current: float, seed: int) -> None:
    """Saves the trajectories of a few pions for the 3D figure.

    Args:
        n_tracks: Number of pions.
        horn_current: Horn current in kA.
        seed: Random seed.
    """
    np.random.seed(seed)
    sim = ChargedPionFluxMiniBooNE(
        proton_energy=np.sqrt(P_PROTON**2 + M_P**2), n_samples=n_tracks,
        n_pot=1, horn_current=horn_current, meson_species="pi_plus")
    sim.simulate_beam_spot()
    sim.meson_flux_pre_horn = sim.charged_meson_flux_mc(
        p_min=0.5, p_max=4.0, theta_min=0.03, theta_max=0.30)
    xs, ys, zs = sim.focus_pions_with_histories()
    os.makedirs(os.path.dirname(TRACKS_FILE), exist_ok=True)
    np.savez_compressed(
        TRACKS_FILE, x=np.concatenate(xs), y=np.concatenate(ys),
        z=np.concatenate(zs), lengths=[len(z) for z in zs],
        momenta=sim.meson_flux_pre_horn[:, 0], horn_current=horn_current)


def main() -> None:
    """Runs the requested simulations and writes the files."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", nargs="+", default=list(TXT_FILES),
                        choices=list(TXT_FILES))
    parser.add_argument("--horn-current", type=float, default=None,
                        help="kA; default +174 (pi+, K+), -174 (pi-)")
    parser.add_argument("--n", type=int, default=150000,
                        help="mesons sampled per species")
    parser.add_argument("--chunks", type=int, default=2,
                        help="parallel chunks per species")
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--theta-max", type=float, default=0.6,
                        help="largest production angle sampled [rad]")
    parser.add_argument("--format", default="pickle",
                        choices=["pickle", "txt", "both"])
    parser.add_argument("--in-volume-only", action="store_true",
                        help="drop mesons that decay outside the beamline")
    parser.add_argument("--tracks", type=int, default=0,
                        help="also save this many pion trajectories")
    opts = parser.parse_args()

    configs = [(s, NOMINAL_CURRENT[s] if opts.horn_current is None
                else opts.horn_current) for s in opts.species]
    records = simulate_samples(configs, opts.n, opts.chunks, opts.seed,
                               opts.theta_max, not opts.in_volume_only)
    for (species, current), record in zip(configs, records):
        if opts.format in ("pickle", "both"):
            path = flux_io.pickle_path(species, current)
            flux_io.save_meson_pickle(record, path)
            print(f"wrote {path}")
        if opts.format in ("txt", "both"):
            path = os.path.join(flux_io.DATA_DIR, TXT_FILES[species])
            flux_io.save_meson_txt(record, path)
            print(f"wrote {path}")
    if opts.tracks:
        simulate_tracks(opts.tracks, 174.0, opts.seed)
        print(f"wrote {TRACKS_FILE}")


if __name__ == "__main__":
    main()
