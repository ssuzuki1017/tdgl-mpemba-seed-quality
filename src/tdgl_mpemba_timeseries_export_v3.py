
"""
Time-series export and mean-curve plotting for the revised TDGL Mpemba study.

This script complements tdgl_mpemba_revised.py.

Why this exists
---------------
The revised simulation script saves sample-level summary CSV files. Those are
enough for transition-time statistics, survival curves, and scatter plots.
However, publication figures often benefit from mean time-series curves such as

    Q(t) = <phi^2>,
    P(t) = ordered-area fraction,
    C_max(t) = maximum ordered cluster size.

For phi6, this v3 script records additional initial seed quality measures:
    p_ordered_seed = fraction of sites with |phi(0)| > 0.5 phi_s,
    c_ordered_seed = max cluster size with |phi(0)| > 0.5 phi_s,
    seed_compactness = compactness of the largest barrier-crossing seed cluster,
    seed_perimeter = 4-neighbor perimeter of the largest barrier-crossing seed cluster,
    seed_perimeter_to_area = seed_perimeter / c_seed,
    ordered_seed_compactness = compactness of the largest ordered-like seed cluster.

This script runs the same revised simulations while storing measured time-series
arrays in compressed NPZ files, and then creates publication-oriented mean-curve
figures with optional standard-deviation or standard-error bands.

Typical workflow
----------------

1. Run a small phi4 time-series test:

    python tdgl_mpemba_timeseries_export.py run_phi4 \
      --out results_ts/phi4_timeseries.npz \
      --figdir figures_phi4_ts \
      --N 64 --nsamples 5 --tmax 10 --preeq_steps 500

2. Run a small phi6 time-series test:

    python tdgl_mpemba_timeseries_export.py run_phi6 \
      --out results_ts/phi6_timeseries.npz \
      --figdir figures_phi6_ts \
      --N 64 --nsamples 3 --tmax 50 --preeq_steps 500

3. Re-plot an existing NPZ:

    python tdgl_mpemba_timeseries_export.py plot_phi6 \
      --npz results_ts/phi6_timeseries.npz \
      --figdir figures_phi6_ts

Notes for production
--------------------
- Time-series arrays can become large. Use measure_every to reduce file size.
- For paper-quality final figures, increase N, nsamples, preeq_steps, and tmax.
- For first-order phi6 nucleation, mean curves are useful but should be shown
  alongside survival curves and sample-level statistics because nucleation is
  a stochastic rare-event process.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import matplotlib.pyplot as plt


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import tdgl_mpemba_revised as base  # noqa: E402


# -----------------------------------------------------------------------------
# Helper utilities
# -----------------------------------------------------------------------------

def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def save_current_figure(figdir: Path, stem: str) -> None:
    ensure_dir(figdir)
    plt.tight_layout()
    plt.savefig(figdir / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.savefig(figdir / f"{stem}.pdf", bbox_inches="tight")
    plt.close()


def array_stats(arr: np.ndarray, axis: int = 1, stderr: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """
    Return mean and spread along an axis.

    arr convention in this script:
        (n_labels, n_samples, n_times)

    For plotting, axis=1 averages over samples.
    """
    mean = np.nanmean(arr, axis=axis)
    std = np.nanstd(arr, axis=axis)
    if stderr:
        n = np.sum(np.isfinite(arr), axis=axis)
        spread = std / np.sqrt(np.maximum(n, 1))
    else:
        spread = std
    return mean, spread


def transition_time(times: np.ndarray, values: np.ndarray, threshold: float) -> float:
    return base.transition_time(times, values, threshold)


def persistent_crossing_time(
    times: np.ndarray,
    values: np.ndarray,
    threshold: float,
    min_consecutive: int,
) -> float:
    return base.persistent_crossing_time(times, values, threshold, min_consecutive=min_consecutive)


def survival_curve(samples: np.ndarray, times: np.ndarray) -> np.ndarray:
    return base.survival_curve(samples, times)


def largest_cluster_geometry_periodic(mask: np.ndarray) -> dict[str, float]:
    """
    Geometry of the largest connected cluster in a 2D periodic binary mask.

    Connectivity: 4-neighbor.
    Returned quantities:
      - area: number of sites in the largest cluster
      - bbox_area: area of the unwrapped bounding box
      - compactness: area / bbox_area
      - rg: radius of gyration of the largest cluster
      - perimeter: number of nearest-neighbor boundary edges of the cluster
      - perimeter_to_area: perimeter / area

    The perimeter is counted on the periodic square lattice. Each nearest-neighbor
    edge connecting a site inside the largest cluster to a site outside that cluster
    contributes one unit. This is a direct lattice proxy for interface length.

    Compactness is a simple operational proxy for whether a seed is spatially
    concentrated. A compact square-like seed has compactness close to 1, while a
    stringy or dispersed seed has a smaller value. For no cluster, area=0,
    compactness=0, perimeter=0, perimeter_to_area=nan, and rg=nan.
    """
    mask = np.asarray(mask, dtype=bool)
    N, M = mask.shape
    visited = np.zeros_like(mask, dtype=bool)

    best_area = 0
    best_coords = None
    best_sites = None

    for i0 in range(N):
        for j0 in range(M):
            if not mask[i0, j0] or visited[i0, j0]:
                continue

            stack = [(i0, j0, 0, 0)]
            visited[i0, j0] = True
            coords = []
            sites = []

            while stack:
                i, j, ui, uj = stack.pop()
                coords.append((ui, uj))
                sites.append((i, j))

                for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    ni = (i + di) % N
                    nj = (j + dj) % M
                    if mask[ni, nj] and not visited[ni, nj]:
                        visited[ni, nj] = True
                        stack.append((ni, nj, ui + di, uj + dj))

            area = len(coords)
            if area > best_area:
                best_area = area
                best_coords = coords
                best_sites = sites

    if best_area == 0 or best_coords is None or best_sites is None:
        return {
            "area": 0.0,
            "bbox_area": 0.0,
            "compactness": 0.0,
            "rg": np.nan,
            "perimeter": 0.0,
            "perimeter_to_area": np.nan,
        }

    arr = np.asarray(best_coords, dtype=float)
    min_xy = np.min(arr, axis=0)
    max_xy = np.max(arr, axis=0)
    bbox_lengths = max_xy - min_xy + 1.0
    bbox_area = float(bbox_lengths[0] * bbox_lengths[1])
    compactness = float(best_area / bbox_area) if bbox_area > 0 else 0.0

    center = np.mean(arr, axis=0)
    rg = float(np.sqrt(np.mean(np.sum((arr - center) ** 2, axis=1)))) if best_area > 0 else np.nan

    cluster_mask = np.zeros_like(mask, dtype=bool)
    for i, j in best_sites:
        cluster_mask[i, j] = True

    perimeter_edges = 0
    for i, j in best_sites:
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ni = (i + di) % N
            nj = (j + dj) % M
            if not cluster_mask[ni, nj]:
                perimeter_edges += 1

    perimeter = float(perimeter_edges)
    perimeter_to_area = float(perimeter / best_area) if best_area > 0 else np.nan

    return {
        "area": float(best_area),
        "bbox_area": bbox_area,
        "compactness": compactness,
        "rg": rg,
        "perimeter": perimeter,
        "perimeter_to_area": perimeter_to_area,
    }


def write_sample_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    ensure_dir(path.parent)
    keys = list(rows[0].keys())
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def metadata_to_jsonable(d: dict) -> str:
    def convert(v):
        if isinstance(v, np.ndarray):
            return v.tolist()
        if isinstance(v, (np.floating, np.integer)):
            return v.item()
        if isinstance(v, Path):
            return str(v)
        return v
    return json.dumps({k: convert(v) for k, v in d.items()}, indent=2)


def parse_labels(values: Optional[list[float]]) -> tuple[float, ...]:
    if values:
        return tuple(float(v) for v in values)
    return (1.05, 1.10, 1.20, 1.50, 2.00, 3.00)


# -----------------------------------------------------------------------------
# Time-series runners
# -----------------------------------------------------------------------------

@dataclass
class Phi4TSConfig:
    N: int = 128
    dx: float = 1.0
    dt: float = 0.05
    tmax: float = 20.0
    preeq_steps: int = 2000
    D_f: float = 0.0
    q_threshold: float = 0.5
    nsamples: int = 20
    seed: int = 1234
    init_mode: str = "preeq"
    remove_zero_mode_initial: bool = True
    D0: float = 0.05
    measure_every: int = 1


@dataclass
class Phi6TSConfig:
    N: int = 128
    dx: float = 1.0
    dt: float = 0.02
    tmax: float = 200.0
    preeq_steps: int = 3000
    a_f: float = 0.10
    b: float = 1.0
    c: float = 1.0
    D_f: float = 1.0e-4
    nsamples: int = 20
    seed: int = 5678
    init_mode: str = "preeq"
    remove_zero_mode_initial: bool = True
    D0: float = 0.02
    init_scheme: str = "noise_only"
    a_i_base: float = 0.30
    measure_every: int = 5
    q_fraction: float = 0.5
    ordered_fraction_threshold: float = 1.0e-3
    cluster_threshold: int = 20
    min_consecutive: int = 3


def run_phi4_timeseries(labels: Iterable[float], config: Phi4TSConfig, out_npz: Path) -> None:
    labels = tuple(float(x) for x in labels)
    params = base.make_phi4_initial_params(labels, D0=config.D0)

    rng = np.random.default_rng(config.seed)
    N, dx, dt = config.N, config.dx, config.dt
    k2 = base.make_k2(N, dx)

    nsteps = int(round(config.tmax / dt))
    measure_indices = np.arange(0, nsteps + 1, config.measure_every, dtype=int)
    times = measure_indices * dt

    Q = np.full((len(labels), config.nsamples, len(times)), np.nan, dtype=np.float32)
    t_tr = np.full((len(labels), config.nsamples), np.nan, dtype=float)

    sample_rows: list[dict] = []

    for ilabel, p in enumerate(params):
        for sample_id in range(config.nsamples):
            if config.init_mode == "preeq":
                if p.r_i is None:
                    raise ValueError("r_i required for phi4 pre-equilibration")
                phi = base.pre_equilibrate_phi4(
                    N=N,
                    dx=dx,
                    dt=dt,
                    nsteps=config.preeq_steps,
                    r_i=p.r_i,
                    D_i=p.D_i,
                    rng=rng,
                    remove_mean_at_end=config.remove_zero_mode_initial,
                )
            elif config.init_mode == "gaussian":
                if p.mu_i is None:
                    raise ValueError("mu_i required for Gaussian initialization")
                phi = base.gaussian_initial_condition(
                    N=N,
                    k2=k2,
                    mu_i=p.mu_i,
                    D_i=p.D_i,
                    rng=rng,
                    remove_zero=config.remove_zero_mode_initial,
                )
            else:
                raise ValueError("init_mode must be 'preeq' or 'gaussian'")

            meas_pos = 0
            for n in range(nsteps + 1):
                if n == measure_indices[meas_pos]:
                    Q[ilabel, sample_id, meas_pos] = np.mean(phi**2)
                    meas_pos += 1
                    if meas_pos >= len(measure_indices):
                        break

                nonlinear = -phi**3
                phi = base.semi_implicit_step(
                    phi=phi,
                    k2=k2,
                    dt=dt,
                    dx=dx,
                    linear_alpha=1.0,
                    nonlinear=nonlinear,
                    D=config.D_f,
                    rng=rng,
                )

            t_val = transition_time(times, Q[ilabel, sample_id], config.q_threshold)
            t_tr[ilabel, sample_id] = t_val

            sample_rows.append(
                {
                    "model": "phi4",
                    "init_mode": config.init_mode,
                    "N": N,
                    "dx": dx,
                    "dt": dt,
                    "tmax": config.tmax,
                    "preeq_steps": config.preeq_steps,
                    "D_f": config.D_f,
                    "label": p.label,
                    "r_i": p.r_i,
                    "D_i": p.D_i,
                    "sample_id": sample_id,
                    "q_threshold": config.q_threshold,
                    "t_tr": t_val,
                }
            )

    ensure_dir(out_npz.parent)
    metadata = {
        "model": "phi4",
        **asdict(config),
        "labels": list(labels),
    }

    np.savez_compressed(
        out_npz,
        model=np.array("phi4"),
        labels=np.array(labels, dtype=float),
        times=times.astype(float),
        Q=Q,
        t_tr=t_tr,
        metadata_json=np.array(metadata_to_jsonable(metadata)),
    )

    write_sample_csv(out_npz.with_suffix(".samples.csv"), sample_rows)
    print(f"Saved phi4 time-series data to {out_npz}")
    print(f"Saved sample summary to {out_npz.with_suffix('.samples.csv')}")


def run_phi6_timeseries(labels: Iterable[float], config: Phi6TSConfig, out_npz: Path) -> None:
    labels = tuple(float(x) for x in labels)
    params = base.make_phi6_initial_params(
        labels,
        D0=config.D0,
        scheme=config.init_scheme,
        a_i_base=config.a_i_base,
    )

    rng = np.random.default_rng(config.seed)
    N, dx, dt = config.N, config.dx, config.dt
    k2 = base.make_k2(N, dx)

    nsteps = int(round(config.tmax / dt))
    measure_indices = np.arange(0, nsteps + 1, config.measure_every, dtype=int)
    times = measure_indices * dt

    phi_b, phi_s = base.extrema_phi6(config.a_f, config.b, config.c)
    if math.isnan(phi_s):
        raise ValueError("No nonzero extrema. Need a_f < b^2/(4c).")

    ordered_threshold = 0.5 * phi_s
    q_threshold = config.q_fraction * phi_s**2

    shape = (len(labels), config.nsamples, len(times))
    Q = np.full(shape, np.nan, dtype=np.float32)
    P = np.full(shape, np.nan, dtype=np.float32)
    Cmax = np.full(shape, np.nan, dtype=np.float32)

    t_nuc_area = np.full((len(labels), config.nsamples), np.nan, dtype=float)
    t_nuc_cluster = np.full((len(labels), config.nsamples), np.nan, dtype=float)
    t_tr = np.full((len(labels), config.nsamples), np.nan, dtype=float)
    p_seed = np.full((len(labels), config.nsamples), np.nan, dtype=float)
    c_seed = np.full((len(labels), config.nsamples), np.nan, dtype=float)
    p_ordered_seed = np.full((len(labels), config.nsamples), np.nan, dtype=float)
    c_ordered_seed = np.full((len(labels), config.nsamples), np.nan, dtype=float)
    seed_compactness = np.full((len(labels), config.nsamples), np.nan, dtype=float)
    seed_bbox_area = np.full((len(labels), config.nsamples), np.nan, dtype=float)
    seed_rg = np.full((len(labels), config.nsamples), np.nan, dtype=float)
    seed_perimeter = np.full((len(labels), config.nsamples), np.nan, dtype=float)
    seed_perimeter_to_area = np.full((len(labels), config.nsamples), np.nan, dtype=float)
    ordered_seed_compactness = np.full((len(labels), config.nsamples), np.nan, dtype=float)
    ordered_seed_bbox_area = np.full((len(labels), config.nsamples), np.nan, dtype=float)
    ordered_seed_rg = np.full((len(labels), config.nsamples), np.nan, dtype=float)
    ordered_seed_perimeter = np.full((len(labels), config.nsamples), np.nan, dtype=float)
    ordered_seed_perimeter_to_area = np.full((len(labels), config.nsamples), np.nan, dtype=float)

    sample_rows: list[dict] = []

    for ilabel, p in enumerate(params):
        for sample_id in range(config.nsamples):
            if config.init_mode == "preeq":
                if p.a_i is None:
                    raise ValueError("a_i required for phi6 pre-equilibration")
                phi = base.pre_equilibrate_phi6(
                    N=N,
                    dx=dx,
                    dt=dt,
                    nsteps=config.preeq_steps,
                    a_i=p.a_i,
                    D_i=p.D_i,
                    rng=rng,
                    b=config.b,
                    c=config.c,
                    remove_mean_at_end=config.remove_zero_mode_initial,
                )
            elif config.init_mode == "gaussian":
                if p.mu_i is None:
                    raise ValueError("mu_i required for Gaussian initialization")
                phi = base.gaussian_initial_condition(
                    N=N,
                    k2=k2,
                    mu_i=p.mu_i,
                    D_i=p.D_i,
                    rng=rng,
                    remove_zero=config.remove_zero_mode_initial,
                )
            else:
                raise ValueError("init_mode must be 'preeq' or 'gaussian'")

            # Barrier-crossing seeds: local regions beyond the metastable barrier phi_b.
            seed_mask = np.abs(phi) > phi_b
            p_seed[ilabel, sample_id] = np.mean(seed_mask)
            seed_geom = largest_cluster_geometry_periodic(seed_mask)
            c_seed[ilabel, sample_id] = seed_geom["area"]
            seed_compactness[ilabel, sample_id] = seed_geom["compactness"]
            seed_bbox_area[ilabel, sample_id] = seed_geom["bbox_area"]
            seed_rg[ilabel, sample_id] = seed_geom["rg"]
            seed_perimeter[ilabel, sample_id] = seed_geom["perimeter"]
            seed_perimeter_to_area[ilabel, sample_id] = seed_geom["perimeter_to_area"]

            # Ordered-like seeds: initial regions already close to the ordered phase.
            # This uses the same threshold as the ordered-area fraction P(t),
            # namely |phi| > 0.5 * phi_s. These variables are useful because
            # barrier-crossing seeds can be numerous but still subcritical.
            ordered_seed_mask = np.abs(phi) > ordered_threshold
            p_ordered_seed[ilabel, sample_id] = np.mean(ordered_seed_mask)
            ordered_seed_geom = largest_cluster_geometry_periodic(ordered_seed_mask)
            c_ordered_seed[ilabel, sample_id] = ordered_seed_geom["area"]
            ordered_seed_compactness[ilabel, sample_id] = ordered_seed_geom["compactness"]
            ordered_seed_bbox_area[ilabel, sample_id] = ordered_seed_geom["bbox_area"]
            ordered_seed_rg[ilabel, sample_id] = ordered_seed_geom["rg"]
            ordered_seed_perimeter[ilabel, sample_id] = ordered_seed_geom["perimeter"]
            ordered_seed_perimeter_to_area[ilabel, sample_id] = ordered_seed_geom["perimeter_to_area"]

            meas_pos = 0
            for n in range(nsteps + 1):
                if n == measure_indices[meas_pos]:
                    ordered_mask = np.abs(phi) > ordered_threshold
                    Q[ilabel, sample_id, meas_pos] = np.mean(phi**2)
                    P[ilabel, sample_id, meas_pos] = np.mean(ordered_mask)
                    Cmax[ilabel, sample_id, meas_pos] = base.max_cluster_size_periodic(ordered_mask)
                    meas_pos += 1
                    if meas_pos >= len(measure_indices):
                        break

                nonlinear = config.b * phi**3 - config.c * phi**5
                phi = base.semi_implicit_step(
                    phi=phi,
                    k2=k2,
                    dt=dt,
                    dx=dx,
                    linear_alpha=-config.a_f,
                    nonlinear=nonlinear,
                    D=config.D_f,
                    rng=rng,
                )

            t_area = persistent_crossing_time(
                times,
                P[ilabel, sample_id],
                config.ordered_fraction_threshold,
                config.min_consecutive,
            )
            t_cluster = persistent_crossing_time(
                times,
                Cmax[ilabel, sample_id],
                config.cluster_threshold,
                config.min_consecutive,
            )
            t_transition = transition_time(times, Q[ilabel, sample_id], q_threshold)

            t_nuc_area[ilabel, sample_id] = t_area
            t_nuc_cluster[ilabel, sample_id] = t_cluster
            t_tr[ilabel, sample_id] = t_transition

            sample_rows.append(
                {
                    "model": "phi6",
                    "init_mode": config.init_mode,
                    "N": N,
                    "dx": dx,
                    "dt": dt,
                    "tmax": config.tmax,
                    "preeq_steps": config.preeq_steps,
                    "a_f": config.a_f,
                    "b": config.b,
                    "c": config.c,
                    "D_f": config.D_f,
                    "label": p.label,
                    "a_i": p.a_i,
                    "D_i": p.D_i,
                    "sample_id": sample_id,
                    "phi_b": phi_b,
                    "phi_s": phi_s,
                    "ordered_threshold": ordered_threshold,
                    "q_threshold": q_threshold,
                    "ordered_fraction_threshold": config.ordered_fraction_threshold,
                    "cluster_threshold": config.cluster_threshold,
                    "p_seed": p_seed[ilabel, sample_id],
                    "c_seed": c_seed[ilabel, sample_id],
                    "p_ordered_seed": p_ordered_seed[ilabel, sample_id],
                    "c_ordered_seed": c_ordered_seed[ilabel, sample_id],
                    "seed_compactness": seed_compactness[ilabel, sample_id],
                    "seed_bbox_area": seed_bbox_area[ilabel, sample_id],
                    "seed_rg": seed_rg[ilabel, sample_id],
                    "seed_perimeter": seed_perimeter[ilabel, sample_id],
                    "seed_perimeter_to_area": seed_perimeter_to_area[ilabel, sample_id],
                    "ordered_seed_compactness": ordered_seed_compactness[ilabel, sample_id],
                    "ordered_seed_bbox_area": ordered_seed_bbox_area[ilabel, sample_id],
                    "ordered_seed_rg": ordered_seed_rg[ilabel, sample_id],
                    "ordered_seed_perimeter": ordered_seed_perimeter[ilabel, sample_id],
                    "ordered_seed_perimeter_to_area": ordered_seed_perimeter_to_area[ilabel, sample_id],
                    "t_nuc_area": t_area,
                    "t_nuc_cluster": t_cluster,
                    "t_tr": t_transition,
                }
            )

    survival_cluster = np.full((len(labels), len(times)), np.nan, dtype=np.float32)
    for ilabel in range(len(labels)):
        survival_cluster[ilabel] = survival_curve(t_nuc_cluster[ilabel], times)

    ensure_dir(out_npz.parent)
    metadata = {
        "model": "phi6",
        **asdict(config),
        "labels": list(labels),
        "phi_b": phi_b,
        "phi_s": phi_s,
        "ordered_threshold": ordered_threshold,
        "q_threshold": q_threshold,
    }

    np.savez_compressed(
        out_npz,
        model=np.array("phi6"),
        labels=np.array(labels, dtype=float),
        times=times.astype(float),
        Q=Q,
        P=P,
        Cmax=Cmax,
        t_nuc_area=t_nuc_area,
        t_nuc_cluster=t_nuc_cluster,
        t_tr=t_tr,
        p_seed=p_seed,
        c_seed=c_seed,
        p_ordered_seed=p_ordered_seed,
        c_ordered_seed=c_ordered_seed,
        seed_compactness=seed_compactness,
        seed_bbox_area=seed_bbox_area,
        seed_rg=seed_rg,
        seed_perimeter=seed_perimeter,
        seed_perimeter_to_area=seed_perimeter_to_area,
        ordered_seed_compactness=ordered_seed_compactness,
        ordered_seed_bbox_area=ordered_seed_bbox_area,
        ordered_seed_rg=ordered_seed_rg,
        ordered_seed_perimeter=ordered_seed_perimeter,
        ordered_seed_perimeter_to_area=ordered_seed_perimeter_to_area,
        survival_cluster=survival_cluster,
        metadata_json=np.array(metadata_to_jsonable(metadata)),
    )

    write_sample_csv(out_npz.with_suffix(".samples.csv"), sample_rows)
    print(f"Saved phi6 time-series data to {out_npz}")
    print(f"Saved sample summary to {out_npz.with_suffix('.samples.csv')}")


# -----------------------------------------------------------------------------
# Plotting from NPZ
# -----------------------------------------------------------------------------

def load_npz(path: str | Path) -> dict:
    data = np.load(path, allow_pickle=False)
    out = {key: data[key] for key in data.files}
    if "metadata_json" in out:
        out["metadata"] = json.loads(str(out["metadata_json"].item()))
    else:
        out["metadata"] = {}
    return out


def plot_phi4_mean_curves(npz_path: Path, figdir: Path, stderr: bool = False) -> None:
    data = load_npz(npz_path)
    labels = data["labels"]
    times = data["times"]
    Q = data["Q"]
    metadata = data["metadata"]
    q_threshold = float(metadata.get("q_threshold", 0.5))

    Q_mean, Q_spread = array_stats(Q, axis=1, stderr=stderr)

    plt.figure(figsize=(6.5, 4.8))
    for i, label in enumerate(labels):
        plt.plot(times, Q_mean[i], label=f"label={label:g}")
        lo = Q_mean[i] - Q_spread[i]
        hi = Q_mean[i] + Q_spread[i]
        plt.fill_between(times, lo, hi, alpha=0.18)
    plt.axhline(q_threshold, linestyle="--", linewidth=1)
    plt.xlabel("time")
    plt.ylabel(r"$Q(t)=\langle\phi^2\rangle$")
    plt.legend()
    save_current_figure(figdir, "phi4_mean_Q_curves")

    t_tr = data["t_tr"]
    means = np.nanmean(t_tr, axis=1)
    stds = np.nanstd(t_tr, axis=1)
    if stderr:
        ns = np.sum(np.isfinite(t_tr), axis=1)
        stds = stds / np.sqrt(np.maximum(ns, 1))

    plt.figure(figsize=(6.2, 4.6))
    plt.errorbar(labels, means, yerr=stds, marker="o", capsize=3)
    plt.xlabel("initial-state label")
    plt.ylabel(r"transition time $t_{tr}$")
    save_current_figure(figdir, "phi4_timeseries_ttr_summary")


def plot_phi6_mean_curves(npz_path: Path, figdir: Path, stderr: bool = False) -> None:
    data = load_npz(npz_path)
    labels = data["labels"]
    times = data["times"]
    Q = data["Q"]
    P = data["P"]
    Cmax = data["Cmax"]
    metadata = data["metadata"]
    q_threshold = float(metadata.get("q_threshold", np.nan))
    ordered_fraction_threshold = float(metadata.get("ordered_fraction_threshold", 1e-3))
    cluster_threshold = float(metadata.get("cluster_threshold", 20))

    Q_mean, Q_spread = array_stats(Q, axis=1, stderr=stderr)
    P_mean, P_spread = array_stats(P, axis=1, stderr=stderr)
    C_mean, C_spread = array_stats(Cmax, axis=1, stderr=stderr)

    plt.figure(figsize=(6.5, 4.8))
    for i, label in enumerate(labels):
        plt.plot(times, Q_mean[i], label=f"label={label:g}")
        plt.fill_between(times, Q_mean[i] - Q_spread[i], Q_mean[i] + Q_spread[i], alpha=0.18)
    if np.isfinite(q_threshold):
        plt.axhline(q_threshold, linestyle="--", linewidth=1)
    plt.xlabel("time")
    plt.ylabel(r"$Q(t)=\langle\phi^2\rangle$")
    plt.legend()
    save_current_figure(figdir, "phi6_mean_Q_curves")

    plt.figure(figsize=(6.5, 4.8))
    for i, label in enumerate(labels):
        plt.plot(times, P_mean[i], label=f"label={label:g}")
        plt.fill_between(times, P_mean[i] - P_spread[i], P_mean[i] + P_spread[i], alpha=0.18)
    plt.axhline(ordered_fraction_threshold, linestyle="--", linewidth=1)
    plt.xlabel("time")
    plt.ylabel(r"ordered-area fraction $P(t)$")
    plt.legend()
    save_current_figure(figdir, "phi6_mean_P_curves")

    plt.figure(figsize=(6.5, 4.8))
    for i, label in enumerate(labels):
        plt.plot(times, C_mean[i], label=f"label={label:g}")
        plt.fill_between(times, C_mean[i] - C_spread[i], C_mean[i] + C_spread[i], alpha=0.18)
    plt.axhline(cluster_threshold, linestyle="--", linewidth=1)
    plt.xlabel("time")
    plt.ylabel(r"maximum ordered cluster size $C_{max}(t)$")
    plt.legend()
    save_current_figure(figdir, "phi6_mean_Cmax_curves")

    if "survival_cluster" in data:
        survival = data["survival_cluster"]
    else:
        survival = np.vstack([survival_curve(data["t_nuc_cluster"][i], times) for i in range(len(labels))])

    plt.figure(figsize=(6.5, 4.8))
    for i, label in enumerate(labels):
        plt.plot(times, survival[i], label=f"label={label:g}")
    plt.xlabel("time")
    plt.ylabel(r"survival probability $S_{surv}(t)$")
    plt.ylim(-0.02, 1.02)
    plt.legend()
    save_current_figure(figdir, "phi6_survival_curves_from_timeseries")

    t_nuc = data["t_nuc_cluster"]
    t_tr = data["t_tr"]
    tn_mean = np.nanmean(t_nuc, axis=1)
    tt_mean = np.nanmean(t_tr, axis=1)
    tn_spread = np.nanstd(t_nuc, axis=1)
    tt_spread = np.nanstd(t_tr, axis=1)
    if stderr:
        n1 = np.sum(np.isfinite(t_nuc), axis=1)
        n2 = np.sum(np.isfinite(t_tr), axis=1)
        tn_spread = tn_spread / np.sqrt(np.maximum(n1, 1))
        tt_spread = tt_spread / np.sqrt(np.maximum(n2, 1))

    plt.figure(figsize=(6.4, 4.8))
    plt.errorbar(labels, tn_mean, yerr=tn_spread, marker="o", capsize=3, label=r"$t_{nuc}$")
    plt.errorbar(labels, tt_mean, yerr=tt_spread, marker="s", capsize=3, label=r"$t_{tr}$")
    plt.xlabel("initial-state label")
    plt.ylabel("time")
    plt.legend()
    save_current_figure(figdir, "phi6_timeseries_tnuc_ttr_summary")


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Export TDGL time series and plot mean curves.")
    sub = parser.add_subparsers(dest="command", required=True)

    p4 = sub.add_parser("run_phi4", help="Run phi4 simulation and save time-series NPZ")
    p4.add_argument("--out", required=True)
    p4.add_argument("--figdir", default=None)
    p4.add_argument("--labels", nargs="*", type=float, default=None)
    p4.add_argument("--N", type=int, default=128)
    p4.add_argument("--dx", type=float, default=1.0)
    p4.add_argument("--dt", type=float, default=0.05)
    p4.add_argument("--tmax", type=float, default=20.0)
    p4.add_argument("--preeq_steps", type=int, default=2000)
    p4.add_argument("--Df", type=float, default=0.0)
    p4.add_argument("--D0", type=float, default=0.05)
    p4.add_argument("--nsamples", type=int, default=20)
    p4.add_argument("--seed", type=int, default=1234)
    p4.add_argument("--init_mode", choices=["preeq", "gaussian"], default="preeq")
    p4.add_argument("--measure_every", type=int, default=1)
    p4.add_argument("--q_threshold", type=float, default=0.5)
    p4.add_argument("--stderr", action="store_true")

    p6 = sub.add_parser("run_phi6", help="Run phi6 simulation and save time-series NPZ")
    p6.add_argument("--out", required=True)
    p6.add_argument("--figdir", default=None)
    p6.add_argument("--labels", nargs="*", type=float, default=None)
    p6.add_argument("--N", type=int, default=128)
    p6.add_argument("--dx", type=float, default=1.0)
    p6.add_argument("--dt", type=float, default=0.02)
    p6.add_argument("--tmax", type=float, default=200.0)
    p6.add_argument("--preeq_steps", type=int, default=3000)
    p6.add_argument("--af", type=float, default=0.10)
    p6.add_argument("--b", type=float, default=1.0)
    p6.add_argument("--c", type=float, default=1.0)
    p6.add_argument("--Df", type=float, default=1.0e-4)
    p6.add_argument("--D0", type=float, default=0.02)
    p6.add_argument("--nsamples", type=int, default=20)
    p6.add_argument("--seed", type=int, default=5678)
    p6.add_argument("--init_mode", choices=["preeq", "gaussian"], default="preeq")
    p6.add_argument("--init_scheme", choices=["noise_only", "mass_and_noise"], default="noise_only")
    p6.add_argument("--a_i_base", type=float, default=0.30)
    p6.add_argument("--measure_every", type=int, default=5)
    p6.add_argument("--q_fraction", type=float, default=0.5)
    p6.add_argument("--ordered_fraction_threshold", type=float, default=1.0e-3)
    p6.add_argument("--cluster_threshold", type=int, default=20)
    p6.add_argument("--min_consecutive", type=int, default=3)
    p6.add_argument("--stderr", action="store_true")

    pp4 = sub.add_parser("plot_phi4", help="Plot phi4 mean curves from NPZ")
    pp4.add_argument("--npz", required=True)
    pp4.add_argument("--figdir", required=True)
    pp4.add_argument("--stderr", action="store_true")

    pp6 = sub.add_parser("plot_phi6", help="Plot phi6 mean curves from NPZ")
    pp6.add_argument("--npz", required=True)
    pp6.add_argument("--figdir", required=True)
    pp6.add_argument("--stderr", action="store_true")

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "run_phi4":
        labels = parse_labels(args.labels)
        cfg = Phi4TSConfig(
            N=args.N,
            dx=args.dx,
            dt=args.dt,
            tmax=args.tmax,
            preeq_steps=args.preeq_steps,
            D_f=args.Df,
            q_threshold=args.q_threshold,
            nsamples=args.nsamples,
            seed=args.seed,
            init_mode=args.init_mode,
            D0=args.D0,
            measure_every=args.measure_every,
        )
        out = Path(args.out)
        run_phi4_timeseries(labels, cfg, out)
        if args.figdir:
            plot_phi4_mean_curves(out, Path(args.figdir), stderr=args.stderr)

    elif args.command == "run_phi6":
        labels = parse_labels(args.labels)
        cfg = Phi6TSConfig(
            N=args.N,
            dx=args.dx,
            dt=args.dt,
            tmax=args.tmax,
            preeq_steps=args.preeq_steps,
            a_f=args.af,
            b=args.b,
            c=args.c,
            D_f=args.Df,
            nsamples=args.nsamples,
            seed=args.seed,
            init_mode=args.init_mode,
            D0=args.D0,
            init_scheme=args.init_scheme,
            a_i_base=args.a_i_base,
            measure_every=args.measure_every,
            q_fraction=args.q_fraction,
            ordered_fraction_threshold=args.ordered_fraction_threshold,
            cluster_threshold=args.cluster_threshold,
            min_consecutive=args.min_consecutive,
        )
        out = Path(args.out)
        run_phi6_timeseries(labels, cfg, out)
        if args.figdir:
            plot_phi6_mean_curves(out, Path(args.figdir), stderr=args.stderr)

    elif args.command == "plot_phi4":
        plot_phi4_mean_curves(Path(args.npz), Path(args.figdir), stderr=args.stderr)

    elif args.command == "plot_phi6":
        plot_phi6_mean_curves(Path(args.npz), Path(args.figdir), stderr=args.stderr)

    else:
        raise ValueError(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()
