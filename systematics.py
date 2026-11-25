"""Flux envelopes from meson-production fit uncertainties."""

import argparse
import dataclasses
import os
from collections.abc import Iterator

import numpy as np

import flux_io
from constants import BR_K_MUNU, BR_PI_MUNU, M_K, M_P, M_PI
from meson_physics import meson_production_d2SdpdOmega
from rkhorn import ChargedPionFluxMiniBooNE, disc_averaged_flux

OUTPUT_DIR = os.path.join(flux_io.DATA_DIR, "systematics")

N_THROWS = 1000
SEED = 20261125
QUANTILES = (0.16, 0.5, 0.84)
N_VALIDATION = 2000000

P_EDGES = np.linspace(0.0, 7.5, 31)
THETA_EDGES = np.linspace(0.0, 0.6, 61)
NU_EDGES = np.arange(0.0, 3050.0, 50.0)

SW_FREE = (0, 1, 3, 4, 5, 6, 7)
SW_POSITIVE = (0, 1, 3, 5)
FS_FREE = (0, 1, 2, 3, 4, 5, 6)
FS_POSITIVE = (0, 3, 4, 5)

PI_PLUS_VALUES = (220.7, 1.080, 1.000, 1.978, 1.32, 5.572, 0.0868,
                  9.686, 1.0)
PI_PLUS_COV = (
    (1707.2, 1.146, -17.646, -15.968, -8.81, -0.7347, -60.816),
    (0.0396, -0.1072, -0.0993, 0.0325, 0.0007, -0.0777),
    (0.5945, 0.5049, 0.0655, 0.0025, 0.198),
    (0.4411, 0.0568, 0.0025, 0.2271),
    (0.2066, 0.0047, 0.1031),
    (0.0005, 0.0641),
    (16.0189,),
)

PI_MINUS_VALUES = (213.7, 0.9379, 5.454, 1.210, 1.284, 4.781, 0.07338,
                   8.329, 1.0)
PI_MINUS_COV = (
    (3688.9, 7.61, -15.666, -17.48, -11.329, -0.9925, -91.4),
    (0.0388, -0.0437, -0.0509, 0.0102, -0.0009, -0.1957),
    (0.0841, 0.0895, 0.0301, 0.0029, 0.2588),
    (0.0986, 0.0375, 0.0033, 0.3141),
    (0.1595, 0.0051, 0.1933),
    (0.0005, 0.064),
    (17.242,),
)

K_PLUS_VALUES = (11.70, 0.88, 4.77, 1.51, 2.21, 2.17, 1.51)
K_PLUS_COV = (
    (1.094, 0.0502, 2.99e-3, -0.0332, -0.0375, 0.125, 0.0743),
    (0.0161, 1.39e-3, -1.44e-3, -0.0126, 0.0322, 0.022),
    (7.47e-3, 2.06e-3, 1.93e-3, 0.0135, -3.34e-3),
    (3.46e-3, 2.03e-3, -4.11e-3, -6.28e-3),
    (0.0146, -0.0154, -0.0244),
    (0.182, 0.126),
    (0.159,),
)

K_PLUS_SB_VALUES = (11.29, 0.87, 4.75, 1.51, 2.21, 2.17, 1.51)
K_PLUS_SB_COV = (
    (0.84, 0.048, 0.0039, -0.032, -0.036, 0.12, 0.069),
    (0.016, 0.0014, -0.0015, -0.013, 0.032, 0.022),
    (0.0073, 0.0020, 0.0019, 0.014, -0.0029),
    (0.0034, 0.0020, -0.0039, -0.0060),
    (0.015, -0.015, -0.024),
    (0.18, 0.12),
    (0.15,),
)


def symmetric_matrix(upper_rows: tuple) -> np.ndarray:
    """Builds a symmetric matrix from its upper triangle.

    Args:
        upper_rows: Row i holds elements (i, i) to (i, n-1).

    Returns:
        Full symmetric matrix.
    """
    size = len(upper_rows)
    matrix = np.zeros((size, size))
    for i, row in enumerate(upper_rows):
        matrix[i, i:] = row
        matrix[i:, i] = row
    return matrix


@dataclasses.dataclass(frozen=True)
class FitModel:
    """Published production fit with its covariance."""

    meson_type: str
    code_values: tuple
    center: tuple
    covariance: np.ndarray
    free: tuple
    positive: tuple
    source: str


@dataclasses.dataclass(frozen=True)
class Sample:
    """Meson sample and the fit used to vary it."""

    name: str
    species: str
    horn_current: float
    model: str
    sigma_scale: float
    m_meson: float
    branching_ratio: float


MODELS = {
    "pi_plus": FitModel(
        "pi_plus", PI_PLUS_VALUES, PI_PLUS_VALUES,
        symmetric_matrix(PI_PLUS_COV), SW_FREE, SW_POSITIVE,
        "arXiv:0806.1449 Table V"),
    "pi_minus": FitModel(
        "pi_minus", PI_MINUS_VALUES, PI_MINUS_VALUES,
        symmetric_matrix(PI_MINUS_COV), SW_FREE, SW_POSITIVE,
        "arXiv:0806.1449 Table VI"),
    "k_plus": FitModel(
        "k_plus", K_PLUS_VALUES, K_PLUS_VALUES,
        symmetric_matrix(K_PLUS_COV), FS_FREE, FS_POSITIVE,
        "arXiv:0806.1449 Table VIII"),
    "k_plus_sciboone": FitModel(
        "k_plus", K_PLUS_VALUES, K_PLUS_SB_VALUES,
        symmetric_matrix(K_PLUS_SB_COV), FS_FREE, FS_POSITIVE,
        "arXiv:1110.0417 Tables IX and X"),
}

SAMPLES = (
    Sample("pi_plus", "pi_plus", 174.0, "pi_plus", 1.0,
           M_PI, BR_PI_MUNU),
    Sample("pi_minus", "pi_minus", -174.0, "pi_minus", 1.0,
           M_PI, BR_PI_MUNU),
    Sample("k_plus", "k_plus", 174.0, "k_plus", 1.0,
           M_K, BR_K_MUNU),
    Sample("pi_plus_horn_off", "pi_plus", 0.0, "pi_plus", 1.0,
           M_PI, BR_PI_MUNU),
    Sample("k_plus_sciboone", "k_plus", 174.0, "k_plus_sciboone", 1.0,
           M_K, BR_K_MUNU),
)


def covariance_factor(covariance: np.ndarray) -> tuple[np.ndarray, float]:
    """Factorises a covariance matrix for sampling.

    Args:
        covariance: Symmetric matrix, possibly rounded.

    Returns:
        Matrix A with A A^T = covariance, smallest eigenvalue.
    """
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    clipped = np.clip(eigenvalues, 0.0, None)
    return eigenvectors * np.sqrt(clipped), float(eigenvalues[0])


def is_physical(model: FitModel, throws: np.ndarray) -> np.ndarray:
    """Flags throws that keep the fall-off parameters positive.

    Args:
        model: Fit that defines the positive parameters.
        throws: Parameter vectors, one per row.

    Returns:
        Boolean mask over the rows.
    """
    return np.all(throws[:, list(model.positive)] > 0.0, axis=1)


def draw_throws(model: FitModel, n_throws: int,
                rng: np.random.Generator, sigma_scale: float = 1.0,
                reject: bool = True) -> tuple[np.ndarray, float]:
    """Draws parameter vectors from the fit covariance.

    Args:
        model: Fit to sample.
        n_throws: Number of accepted vectors.
        rng: Random generator.
        sigma_scale: Factor applied to all standard deviations.
        reject: Redraw throws that fail is_physical.

    Returns:
        Throws of shape (n_throws, n_params), rejected fraction.
    """
    factor, _ = covariance_factor(model.covariance)
    free = list(model.free)
    accepted = np.empty((0, len(model.center)))
    n_drawn = 0
    while accepted.shape[0] < n_throws:
        normal = rng.standard_normal((n_throws, len(free)))
        batch = np.tile(np.array(model.center), (n_throws, 1))
        batch[:, free] += sigma_scale * normal @ factor.T
        n_drawn += n_throws
        if reject:
            batch = batch[is_physical(model, batch)]
        accepted = np.vstack([accepted, batch])
    return accepted[:n_throws], 1.0 - accepted.shape[0] / n_drawn


def cross_section(model: FitModel, p: np.ndarray, theta: np.ndarray,
                  p_beam: float, values: np.ndarray) -> np.ndarray:
    """Evaluates the production cross section of a fit.

    Args:
        model: Fit that names the functional form.
        p: Meson momentum in GeV/c.
        theta: Meson angle to the beam in rad.
        p_beam: Proton momentum in GeV/c.
        values: Fit parameters.

    Returns:
        Cross section in mb/(GeV/c sr).
    """
    return meson_production_d2SdpdOmega(
        p, theta, p_beam, meson_type=model.meson_type,
        params=tuple(float(v) for v in values))


@dataclasses.dataclass(frozen=True)
class Histogram:
    """Binning of one weighted quantity of a meson sample."""

    name: str
    edges: np.ndarray
    rows: np.ndarray
    bins: np.ndarray
    weights: np.ndarray
    scale: np.ndarray

    def fill(self, ratio: np.ndarray) -> np.ndarray:
        """Histograms the reweighted sample.

        Args:
            ratio: Cross-section ratio of every meson.

        Returns:
            Bin contents times the bin scale.
        """
        n_bins = self.edges.shape[0] - 1
        counts = np.bincount(self.bins,
                             weights=self.weights * ratio[self.rows],
                             minlength=n_bins + 1)
        return counts[:n_bins] * self.scale


def make_histogram(name: str, edges: np.ndarray, values: np.ndarray,
                   weights: np.ndarray, rows: np.ndarray,
                   density: bool) -> Histogram:
    """Prepares the bin index of every entry.

    Args:
        name: Key prefix in the output file.
        edges: Bin edges.
        values: Quantity to histogram, one per entry.
        weights: Nominal weight of every entry.
        rows: Index of the meson behind every entry.
        density: Divide the contents by the bin width.

    Returns:
        Histogram ready to fill.
    """
    n_bins = edges.shape[0] - 1
    bins = np.searchsorted(edges, values, side="right") - 1
    bins[values == edges[-1]] = n_bins - 1
    bins[(bins < 0) | (bins >= n_bins)] = n_bins
    scale = 1.0 / np.diff(edges) if density else np.ones(n_bins)
    return Histogram(name, edges, rows, bins, weights, scale)


def build_histograms(record: dict, sample: Sample) -> list[Histogram]:
    """Sets up the five histograms of a meson sample.

    Args:
        record: Meson record from flux_io.load_meson_pickle.
        sample: Sample description.

    Returns:
        Pre-horn, post-horn and neutrino histograms.
    """
    data = {key: np.asarray(value, dtype=np.float64)
            for key, value in record.items() if key != "meta"}
    all_rows = np.arange(data["p0"].shape[0])
    decayed = np.flatnonzero(data["decay_in_volume"] > 0.0)
    e_nu, w_nu = disc_averaged_flux(
        flux_io.flux_array(record), m_meson=sample.m_meson,
        br=sample.branching_ratio)
    flux_rows = np.tile(decayed, e_nu.shape[0] // decayed.shape[0])
    return [
        make_histogram("pre_horn_p", P_EDGES, data["p0"] * 1.0e-3,
                       data["weight0"], all_rows, True),
        make_histogram("pre_horn_theta", THETA_EDGES, data["theta0"],
                       data["weight0"], all_rows, True),
        make_histogram("post_horn_p", P_EDGES, data["p"][decayed] * 1.0e-3,
                       data["weight"][decayed], decayed, True),
        make_histogram("post_horn_theta", THETA_EDGES,
                       data["theta"][decayed], data["weight"][decayed],
                       decayed, True),
        make_histogram("nu_flux", NU_EDGES, e_nu, w_nu, flux_rows, False),
    ]


def weight_ratios(record: dict, model: FitModel,
                  throws: np.ndarray) -> Iterator[np.ndarray]:
    """Computes the cross-section ratio of every meson.

    Args:
        record: Meson record.
        model: Fit to vary.
        throws: Parameter vectors, one per row.

    Yields:
        Varied over generated cross section, one array per throw.
    """
    p = np.asarray(record["p0"], dtype=np.float64) * 1.0e-3
    theta = np.asarray(record["theta0"], dtype=np.float64)
    p_beam = record["meta"]["p_proton_mev"] * 1.0e-3
    generated = cross_section(model, p, theta, p_beam, model.code_values)
    for values in throws:
        varied = cross_section(model, p, theta, p_beam, values)
        yield np.divide(varied, generated, out=np.zeros_like(varied),
                        where=generated > 0.0)


def summarise(histogram: Histogram, nominal: np.ndarray,
              thrown: np.ndarray) -> dict:
    """Reduces the thrown histograms to an envelope.

    Args:
        histogram: Binning that was filled.
        nominal: Contents for the central parameters.
        thrown: Contents per throw, shape (n_throws, n_bins).

    Returns:
        Output arrays keyed by name.
    """
    name = histogram.name
    lower, median, upper = np.quantile(thrown, QUANTILES, axis=0)
    width = 1.0 / histogram.scale
    total = np.sum(nominal * width)
    totals = np.sum(thrown * width, axis=1) / total
    total_lower, _, total_upper = np.quantile(totals, QUANTILES)
    return {
        f"{name}_edges": histogram.edges,
        f"{name}_nominal": nominal,
        f"{name}_lower": lower,
        f"{name}_upper": upper,
        f"{name}_median": median,
        f"{name}_std": np.std(thrown, axis=0, ddof=1),
        f"{name}_integral_nominal": total,
        f"{name}_integral_frac_std": np.std(totals, ddof=1),
        f"{name}_integral_frac_lower": total_lower - 1.0,
        f"{name}_integral_frac_upper": total_upper - 1.0,
    }


def production_normalisation(record: dict) -> float:
    """Returns mesons per POT per mb of integrated cross section.

    Args:
        record: Meson record with its "meta" dict.

    Returns:
        Interaction probability over reaction cross section.
    """
    meta = record["meta"]
    energy = np.sqrt(meta["p_proton_mev"]**2 + M_P**2)
    sim = ChargedPionFluxMiniBooNE(
        proton_energy=energy, n_samples=1, n_pot=meta["n_pot"],
        horn_current=meta["horn_current_ka"],
        meson_species=meta["species"])
    return (meta["n_pot"] * sim.interaction_probability()
            / sim.sigma_reaction)


def validate_reweighting(record: dict, model: FitModel,
                         values: np.ndarray, ratio: np.ndarray,
                         rng: np.random.Generator) -> dict:
    """Compares reweighted yield with a direct integral.

    Args:
        record: Meson record.
        model: Fit that was varied.
        values: Parameters of the throw to check.
        ratio: Cross-section ratio of every meson for it.
        rng: Random generator for the independent sample.

    Returns:
        Both yields per POT with their statistical errors.
    """
    meta = record["meta"]
    p_min, p_max = meta["p_range_gev"]
    theta_min, theta_max = meta["theta_range_rad"]
    p = rng.uniform(p_min, p_max, N_VALIDATION)
    theta = rng.uniform(theta_min, theta_max, N_VALIDATION)
    integrand = cross_section(model, p, theta,
                              meta["p_proton_mev"] * 1.0e-3,
                              values) * np.sin(theta)
    volume = 2.0 * np.pi * (theta_max - theta_min) * (p_max - p_min)
    factor = volume * production_normalisation(record)
    weights = np.asarray(record["weight0"], dtype=np.float64) * ratio
    n_mesons = weights.shape[0]
    return {
        "validation_reweighted": np.sum(weights),
        "validation_reweighted_error": np.sqrt(
            n_mesons * np.var(weights, ddof=1)),
        "validation_direct": factor * np.mean(integrand),
        "validation_direct_error": factor * np.std(integrand, ddof=1)
        / np.sqrt(N_VALIDATION),
        "validation_params": np.array(values),
    }


def run_sample(sample: Sample, n_throws: int, seed: int) -> dict | None:
    """Builds the envelope of one meson sample.

    Args:
        sample: Sample description.
        n_throws: Number of parameter throws.
        seed: Seed of the random generator.

    Returns:
        Output arrays keyed by name, or None without input.
    """
    path = flux_io.pickle_path(sample.species, sample.horn_current)
    if not os.path.exists(path):
        print(f"{sample.name}: {path} not found, skipped")
        return None
    record = flux_io.load_meson_pickle(path)
    model = MODELS[sample.model]
    rng = np.random.default_rng(seed)
    throws, rejected_fraction = draw_throws(model, n_throws, rng,
                                     sample.sigma_scale)
    histograms = build_histograms(record, sample)
    center = np.array(model.center)[np.newaxis, :]
    nominal_ratio = next(weight_ratios(record, model, center))
    thrown = [np.empty((n_throws, h.edges.shape[0] - 1))
              for h in histograms]
    first_ratio = None
    for i, ratio in enumerate(weight_ratios(record, model, throws)):
        if first_ratio is None:
            first_ratio = ratio
        for histogram, contents in zip(histograms, thrown):
            contents[i] = histogram.fill(ratio)

    _, min_eigenvalue = covariance_factor(model.covariance)
    flux = histograms[-1]
    n_cells = NU_EDGES.shape[0]
    per_meson = np.bincount(
        flux.rows * n_cells + flux.bins,
        weights=flux.weights * nominal_ratio[flux.rows],
        minlength=nominal_ratio.shape[0] * n_cells).reshape(-1, n_cells)
    flux_stat = np.sum(per_meson**2, axis=0)[:-1]
    output = {
        "sample": sample.name,
        "species": sample.species,
        "horn_current_ka": sample.horn_current,
        "source": model.source,
        "n_throws": n_throws,
        "rejected_fraction": rejected_fraction,
        "seed": seed,
        "n_mesons": record["meta"]["n_samples"],
        "quantiles": np.array(QUANTILES),
        "sigma_scale": sample.sigma_scale,
        "params_center": np.array(model.center),
        "params_generated": np.array(model.code_values),
        "params_free_index": np.array(model.free),
        "params_positive_index": np.array(model.positive),
        "params_covariance": model.covariance,
        "params_covariance_min_eigenvalue": min_eigenvalue,
        "params_throws": throws,
        "nu_flux_mc_stat_error": np.sqrt(flux_stat),
        "nu_flux_covariance": np.cov(thrown[-1], rowvar=False),
    }
    for histogram, contents in zip(histograms, thrown):
        output.update(summarise(histogram, histogram.fill(nominal_ratio),
                                contents))
    output.update(validate_reweighting(record, model, throws[0],
                                       first_ratio, rng))
    return output


def report(output: dict) -> None:
    """Prints the headline numbers of one envelope.

    Args:
        output: Result of run_sample.
    """
    print(f"\n{output['sample']} ({output['source']}, sigma x "
          f"{output['sigma_scale']:.0f}): {output['n_throws']} throws, "
          f"{100*output['rejected_fraction']:.1f}% rejected")
    for name in ("pre_horn_p", "post_horn_p", "nu_flux"):
        print(f"  {name:12s} integral {output[f'{name}_integral_nominal']:.4e}"
              f"  std {100*output[f'{name}_integral_frac_std']:5.1f}%"
              f"  16/84% {100*output[f'{name}_integral_frac_lower']:+6.1f}%"
              f" {100*output[f'{name}_integral_frac_upper']:+6.1f}%")
    direct = output["validation_direct"]
    reweighted = output["validation_reweighted"]
    print(f"  throw 0 yield: reweighted {reweighted:.5f} +- "
          f"{output['validation_reweighted_error']:.5f}, direct "
          f"{direct:.5f} +- {output['validation_direct_error']:.5f}, "
          f"ratio {reweighted/direct:.4f}")


def main() -> None:
    """Writes the envelope file of every available sample."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--throws", type=int, default=N_THROWS)
    parser.add_argument("--seed", type=int, default=SEED)
    opts = parser.parse_args()
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for offset, sample in enumerate(SAMPLES):
        output = run_sample(sample, opts.throws, opts.seed + offset)
        if output is None:
            continue
        report(output)
        np.savez(os.path.join(OUTPUT_DIR, f"envelope_{sample.name}.npz"),
                 **output)


if __name__ == "__main__":
    main()
