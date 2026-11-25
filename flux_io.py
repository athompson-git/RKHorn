"""Reads and writes simulated meson samples."""

import os
import pickle

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")
PICKLE_DIR = os.path.join(DATA_DIR, "pickle")

FLUX_COLUMNS = ("p", "theta", "decay_x", "decay_y", "decay_z",
                "decay_in_volume", "weight", "phi")
TXT_COLUMNS = ("p[MeV] theta[rad] decay_x[m] decay_y[m] decay_z[m] "
               "decay_in_volume weight[1/POT] phi[rad]")


def pickle_path(species: str, horn_current: float) -> str:
    """Returns the default pickle file of a meson sample.

    Args:
        species: Meson species, e.g. "pi_plus".
        horn_current: Horn current in kA.

    Returns:
        Path inside data/pickle.
    """
    name = f"mesons_{species}_{horn_current:+.0f}kA.pkl"
    return os.path.join(PICKLE_DIR, name)


def save_meson_pickle(record: dict, path: str) -> None:
    """Writes a meson record to a pickle file.

    Args:
        record: Output of ChargedPionFluxMiniBooNE.meson_record().
        path: File to write.
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "wb") as file:
        pickle.dump(record, file, protocol=pickle.HIGHEST_PROTOCOL)


def load_meson_pickle(path: str) -> dict:
    """Reads a meson record written by save_meson_pickle.

    Args:
        path: File to read.

    Returns:
        Arrays keyed by quantity, plus a "meta" dict.
    """
    with open(path, "rb") as file:
        return pickle.load(file)


def merge_records(records: list[dict]) -> dict:
    """Combines equal-size samples of the same configuration.

    Args:
        records: Meson records from independent runs.

    Returns:
        One record with the weights averaged over the runs.
    """
    merged = {key: np.concatenate([r[key] for r in records])
              for key in records[0] if key != "meta"}
    for key in ("weight", "weight0"):
        merged[key] = merged[key] / len(records)
    merged["meta"] = dict(records[0]["meta"])
    merged["meta"]["n_samples"] = sum(r["meta"]["n_samples"]
                                      for r in records)
    return merged


def flux_array(record: dict, in_volume_only: bool = True) -> np.ndarray:
    """Converts a meson record to the eight-column flux table.

    Args:
        record: Meson record.
        in_volume_only: Keep only decays in the open beamline.

    Returns:
        Rows with the columns named in FLUX_COLUMNS.
    """
    table = np.column_stack([record[key] for key in FLUX_COLUMNS])
    if in_volume_only:
        table = table[record["decay_in_volume"] > 0.0]
    return table.astype(np.float64)


def save_meson_txt(record: dict, path: str) -> None:
    """Writes the decays in the beamline as a text table.

    Args:
        record: Meson record.
        path: File to write.
    """
    meta = record["meta"]
    header = (f"horn-focused {meta['species']} per POT, horn current "
              f"{meta['horn_current_ka']} kA, {meta['n_samples']} mesons "
              f"sampled; columns: {TXT_COLUMNS}")
    np.savetxt(path, flux_array(record), fmt="%.6e", header=header)
