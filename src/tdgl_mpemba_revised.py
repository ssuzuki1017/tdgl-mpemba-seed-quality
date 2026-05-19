"""
Revised TDGL code for Mpemba-like transition-time studies.

Features
--------
1. Pre-equilibrated initial conditions for both continuous (phi^4) and
   first-order (phi^6) models.
2. Optional Gaussian initial ensemble for controlled comparisons.
3. Semi-implicit spectral time stepping.
4. Transition-time detection for Q(t) = <phi^2>.
5. Cluster-based nucleation detection for first-order transitions.
6. Survival curves and sample-level CSV logging.

The code is intentionally written as a research script rather than a polished
package. Start with small N/nsamples for exploration, then increase them for
production runs.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
import csv
import math

import numpy as np
import matplotlib.pyplot as plt


Array = np.ndarray


# -----------------------------------------------------------------------------
# Fourier utilities
# -----------------------------------------------------------------------------

def make_k2(N: int, dx: float) -> Array:
    """Return k^2 grid for a 2D periodic N x N lattice."""
    kx = 2.0 * np.pi * np.fft.fftfreq(N, d=dx)
    ky = 2.0 * np.pi * np.fft.fftfreq(N, d=dx)
    kx, ky = np.meshgrid(kx, ky, indexing="ij")
    return kx**2 + ky**2


def remove_zero_mode(phi: Array) -> Array:
    """Remove the spatial mean of the field."""
    return phi - np.mean(phi)


# -----------------------------------------------------------------------------
# Potential and extrema
# -----------------------------------------------------------------------------

def potential_phi4_final(phi: Array) -> Array:
    """V(phi) = -phi^2/2 + phi^4/4."""
    return -0.5 * phi**2 + 0.25 * phi**4


def potential_phi4_initial(phi: Array, r_i: float) -> Array:
    """V_i(phi) = r_i phi^2/2 + phi^4/4, with r_i > 0."""
    return 0.5 * r_i * phi**2 + 0.25 * phi**4


def potential_phi6(phi: Array, a: float, b: float = 1.0, c: float = 1.0) -> Array:
    """V(phi) = a phi^2/2 - b phi^4/4 + c phi^6/6."""
    return 0.5 * a * phi**2 - 0.25 * b * phi**4 + (c / 6.0) * phi**6


def extrema_phi6(a: float, b: float = 1.0, c: float = 1.0) -> Tuple[float, float]:
    """
    Return barrier position phi_b and stable ordered-phase position phi_s.

    For V(phi) = a phi^2/2 - b phi^4/4 + c phi^6/6.
    Nonzero extrema exist for a <= b^2/(4c).
    """
    disc = b**2 - 4.0 * a * c
    if disc <= 0:
        return math.nan, math.nan

    y_minus = (b - math.sqrt(disc)) / (2.0 * c)
    y_plus = (b + math.sqrt(disc)) / (2.0 * c)
    return math.sqrt(y_minus), math.sqrt(y_plus)


def energy_density(phi: Array, k2: Array, dx: float, model: str, **params: float) -> float:
    """
    Estimate spatial average of free-energy density.

    The gradient term is evaluated spectrally as 0.5 * <|grad phi|^2>.
    This is mainly useful for deterministic-noise checks.
    """
    phi_hat = np.fft.fft2(phi)
    grad_density = 0.5 * np.sum(k2 * np.abs(phi_hat) ** 2) / (phi.size**2)

    if model == "phi4_final":
        pot_density = np.mean(potential_phi4_final(phi))
    elif model == "phi4_initial":
        pot_density = np.mean(potential_phi4_initial(phi, params["r_i"]))
    elif model == "phi6":
        pot_density = np.mean(potential_phi6(phi, params["a"], params.get("b", 1.0), params.get("c", 1.0)))
    else:
        raise ValueError(f"Unknown model: {model}")

    return float(grad_density + pot_density)


# -----------------------------------------------------------------------------
# Time stepping
# -----------------------------------------------------------------------------

def semi_implicit_step(
    phi: Array,
    k2: Array,
    dt: float,
    dx: float,
    linear_alpha: float,
    nonlinear: Array,
    D: float,
    rng: np.random.Generator,
) -> Array:
    """
    One semi-implicit spectral step for

        dphi/dt = laplacian(phi) + linear_alpha * phi + nonlinear(phi) + noise.

    Fourier update:

        phi_hat^{n+1} = FFT[phi^n + dt*nonlinear + noise]
                         / [1 + dt*(k^2 - linear_alpha)]

    For stable initial phi4: linear_alpha = -r_i, nonlinear = -phi^3.
    For final phi4:         linear_alpha =  1,   nonlinear = -phi^3.
    For phi6:               linear_alpha = -a,   nonlinear = b phi^3 - c phi^5.
    """
    if D > 0.0:
        sigma = math.sqrt(2.0 * D * dt / (dx**2))
        noise = sigma * rng.normal(0.0, 1.0, size=phi.shape)
    else:
        noise = 0.0

    rhs = phi + dt * nonlinear + noise
    rhs_hat = np.fft.fft2(rhs)
    denom = 1.0 + dt * (k2 - linear_alpha)
    phi_next = np.fft.ifft2(rhs_hat / denom).real
    return phi_next


def pre_equilibrate_phi4(
    N: int,
    dx: float,
    dt: float,
    nsteps: int,
    r_i: float,
    D_i: float,
    rng: np.random.Generator,
    remove_mean_at_end: bool = True,
    initial_std: float = 1.0e-3,
) -> Array:
    """Pre-equilibrate a high-temperature phi^4 state with r_i > 0."""
    k2 = make_k2(N, dx)
    phi = initial_std * rng.normal(0.0, 1.0, size=(N, N))

    for _ in range(nsteps):
        nonlinear = -phi**3
        phi = semi_implicit_step(
            phi=phi,
            k2=k2,
            dt=dt,
            dx=dx,
            linear_alpha=-r_i,
            nonlinear=nonlinear,
            D=D_i,
            rng=rng,
        )

    if remove_mean_at_end:
        phi = remove_zero_mode(phi)
    return phi


def pre_equilibrate_phi6(
    N: int,
    dx: float,
    dt: float,
    nsteps: int,
    a_i: float,
    D_i: float,
    rng: np.random.Generator,
    b: float = 1.0,
    c: float = 1.0,
    remove_mean_at_end: bool = True,
    initial_std: float = 1.0e-3,
) -> Array:
    """Pre-equilibrate a high-temperature phi^6 state, usually with a_i > b^2/(4c)."""
    k2 = make_k2(N, dx)
    phi = initial_std * rng.normal(0.0, 1.0, size=(N, N))

    for _ in range(nsteps):
        nonlinear = b * phi**3 - c * phi**5
        phi = semi_implicit_step(
            phi=phi,
            k2=k2,
            dt=dt,
            dx=dx,
            linear_alpha=-a_i,
            nonlinear=nonlinear,
            D=D_i,
            rng=rng,
        )

    if remove_mean_at_end:
        phi = remove_zero_mode(phi)
    return phi


def gaussian_initial_condition(
    N: int,
    k2: Array,
    mu_i: float,
    D_i: float,
    rng: np.random.Generator,
    remove_zero: bool = True,
) -> Array:
    """
    Controlled Gaussian initial ensemble with S(k,0) ~ D_i/(mu_i+k^2).

    Use this mainly for comparison with linear theory, not as the only
    production initial condition.
    """
    white = rng.normal(0.0, 1.0, size=(N, N))
    white_hat = np.fft.fft2(white)
    amp = np.sqrt(D_i / (mu_i + k2))
    if remove_zero:
        amp[0, 0] = 0.0
    phi = np.fft.ifft2(amp * white_hat).real
    if remove_zero:
        phi = remove_zero_mode(phi)
    return phi


# -----------------------------------------------------------------------------
# Measurements
# -----------------------------------------------------------------------------

def transition_time(times: Array, values: Array, threshold: float) -> float:
    """First crossing time values >= threshold, with linear interpolation."""
    idx = np.where(values >= threshold)[0]
    if len(idx) == 0:
        return math.nan
    i = int(idx[0])
    if i == 0:
        return float(times[0])
    t0, t1 = times[i - 1], times[i]
    y0, y1 = values[i - 1], values[i]
    if y1 == y0:
        return float(t1)
    return float(t0 + (threshold - y0) * (t1 - t0) / (y1 - y0))


def persistent_crossing_time(
    times: Array,
    values: Array,
    threshold: float,
    min_consecutive: int = 3,
) -> float:
    """First time values stay above threshold for min_consecutive measurements."""
    above = values >= threshold
    count = 0
    for i, flag in enumerate(above):
        if flag:
            count += 1
            if count >= min_consecutive:
                return float(times[i - min_consecutive + 1])
        else:
            count = 0
    return math.nan


def max_cluster_size_periodic(mask: Array) -> int:
    """Maximum 4-connected cluster size on a 2D periodic lattice."""
    N, M = mask.shape
    visited = np.zeros_like(mask, dtype=bool)
    max_size = 0

    for i in range(N):
        for j in range(M):
            if not mask[i, j] or visited[i, j]:
                continue

            size = 0
            q = deque([(i, j)])
            visited[i, j] = True

            while q:
                x, y = q.popleft()
                size += 1

                neighbors = (
                    ((x + 1) % N, y),
                    ((x - 1) % N, y),
                    (x, (y + 1) % M),
                    (x, (y - 1) % M),
                )
                for nx, ny in neighbors:
                    if mask[nx, ny] and not visited[nx, ny]:
                        visited[nx, ny] = True
                        q.append((nx, ny))

            if size > max_size:
                max_size = size

    return int(max_size)


def survival_curve(t_nuc_samples: Array, times: Array) -> Array:
    """
    Survival probability S(t) = Pr(t_nuc > t).

    NaN nucleation times are treated as censored observations that survive
    through the entire observed time window.
    """
    out = np.empty_like(times, dtype=float)
    for i, t in enumerate(times):
        out[i] = np.mean(np.isnan(t_nuc_samples) | (t_nuc_samples > t))
    return out


# -----------------------------------------------------------------------------
# Configuration dataclasses
# -----------------------------------------------------------------------------

@dataclass
class InitialParam:
    label: float
    D_i: float
    r_i: Optional[float] = None
    a_i: Optional[float] = None
    mu_i: Optional[float] = None


@dataclass
class Phi4Config:
    N: int = 128
    dx: float = 1.0
    dt: float = 0.05
    tmax: float = 20.0
    preeq_steps: int = 2000
    D_f: float = 0.0
    q_threshold: float = 0.5
    nsamples: int = 20
    seed: int = 1234
    init_mode: str = "preeq"  # "preeq" or "gaussian"
    remove_zero_mode_initial: bool = True
    output_csv: str = "phi4_samples.csv"


@dataclass
class Phi6Config:
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
    init_mode: str = "preeq"  # "preeq" or "gaussian"
    remove_zero_mode_initial: bool = True
    measure_every: int = 5
    q_fraction: float = 0.5
    ordered_fraction_threshold: float = 1.0e-3
    cluster_threshold: int = 20
    min_consecutive: int = 3
    output_csv: str = "phi6_samples.csv"


# -----------------------------------------------------------------------------
# Parameter helpers
# -----------------------------------------------------------------------------

def make_phi4_initial_params(
    Ti_list: Iterable[float],
    D0: float = 0.05,
) -> List[InitialParam]:
    """
    Convenience mapping from labels Ti to r_i and D_i.

    r_i = 2(Ti-1) mimics distance from the continuous transition.
    D_i = D0 Ti controls initial noise strength.
    """
    params: List[InitialParam] = []
    for Ti in Ti_list:
        params.append(InitialParam(label=float(Ti), r_i=2.0 * (Ti - 1.0), D_i=D0 * Ti, mu_i=2.0 * (Ti - 1.0)))
    return params


def make_phi6_initial_params(
    Ti_list: Iterable[float],
    D0: float = 0.02,
    scheme: str = "noise_only",
    a_i_base: float = 0.30,
) -> List[InitialParam]:
    """
    Convenience initial parameters for first-order phi^6 runs.

    scheme="noise_only":
        a_i is fixed above the high-temperature spinodal, while D_i grows with Ti.
        This isolates the effect of stronger initial fluctuations.

    scheme="mass_and_noise":
        both a_i and D_i grow with Ti. This is a useful robustness check.
    """
    params: List[InitialParam] = []
    for Ti in Ti_list:
        if scheme == "noise_only":
            a_i = a_i_base
        elif scheme == "mass_and_noise":
            a_i = a_i_base + 2.0 * (Ti - 1.0)
        else:
            raise ValueError("scheme must be 'noise_only' or 'mass_and_noise'")
        params.append(InitialParam(label=float(Ti), a_i=float(a_i), D_i=D0 * Ti, mu_i=max(1.0e-8, a_i)))
    return params


# -----------------------------------------------------------------------------
# Main experiment runners
# -----------------------------------------------------------------------------

def run_phi4_experiment(config: Phi4Config, params: List[InitialParam]) -> Tuple[Array, Dict[float, dict]]:
    """Run revised continuous-transition phi^4 experiment."""
    rng = np.random.default_rng(config.seed)
    N, dx, dt = config.N, config.dx, config.dt
    k2 = make_k2(N, dx)
    nsteps = int(round(config.tmax / dt))
    times = np.arange(nsteps + 1) * dt

    rows = []
    results: Dict[float, dict] = {}

    for p in params:
        Q_all = []
        t_all = []

        for sample_id in range(config.nsamples):
            if config.init_mode == "preeq":
                if p.r_i is None:
                    raise ValueError("r_i is required for phi4 pre-equilibration")
                phi = pre_equilibrate_phi4(
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
                    raise ValueError("mu_i is required for Gaussian initialization")
                phi = gaussian_initial_condition(N, k2, p.mu_i, p.D_i, rng, config.remove_zero_mode_initial)
            else:
                raise ValueError("init_mode must be 'preeq' or 'gaussian'")

            Q = np.empty(nsteps + 1)
            Q[0] = np.mean(phi**2)

            for n in range(1, nsteps + 1):
                nonlinear = -phi**3
                phi = semi_implicit_step(
                    phi=phi,
                    k2=k2,
                    dt=dt,
                    dx=dx,
                    linear_alpha=1.0,
                    nonlinear=nonlinear,
                    D=config.D_f,
                    rng=rng,
                )
                Q[n] = np.mean(phi**2)

            t_tr = transition_time(times, Q, config.q_threshold)
            Q_all.append(Q)
            t_all.append(t_tr)

            rows.append({
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
                "t_tr": t_tr,
            })

        Q_all_arr = np.array(Q_all)
        t_all_arr = np.array(t_all)
        results[p.label] = {
            "Q_mean": np.nanmean(Q_all_arr, axis=0),
            "Q_std": np.nanstd(Q_all_arr, axis=0),
            "t_tr_samples": t_all_arr,
            "t_tr_mean": float(np.nanmean(t_all_arr)),
            "t_tr_std": float(np.nanstd(t_all_arr)),
        }

    write_csv(config.output_csv, rows)
    return times, results


def run_phi6_experiment(config: Phi6Config, params: List[InitialParam]) -> Tuple[Array, Dict[float, dict], dict]:
    """Run revised first-order phi^6 experiment with cluster nucleation detection."""
    rng = np.random.default_rng(config.seed)
    N, dx, dt = config.N, config.dx, config.dt
    k2 = make_k2(N, dx)
    nsteps = int(round(config.tmax / dt))
    measure_indices = np.arange(0, nsteps + 1, config.measure_every, dtype=int)
    measure_times = measure_indices * dt

    phi_b, phi_s = extrema_phi6(config.a_f, config.b, config.c)
    if math.isnan(phi_s):
        raise ValueError("No nonzero extrema. Need a_f < b^2/(4c).")

    ordered_threshold = 0.5 * phi_s
    q_threshold = config.q_fraction * phi_s**2

    rows = []
    results: Dict[float, dict] = {}

    for p in params:
        Q_all = []
        P_all = []
        Cmax_all = []
        t_nuc_cluster_all = []
        t_nuc_area_all = []
        t_tr_all = []
        p_seed_all = []
        c_seed_all = []

        for sample_id in range(config.nsamples):
            if config.init_mode == "preeq":
                if p.a_i is None:
                    raise ValueError("a_i is required for phi6 pre-equilibration")
                phi = pre_equilibrate_phi6(
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
                    raise ValueError("mu_i is required for Gaussian initialization")
                phi = gaussian_initial_condition(N, k2, p.mu_i, p.D_i, rng, config.remove_zero_mode_initial)
            else:
                raise ValueError("init_mode must be 'preeq' or 'gaussian'")

            seed_mask = np.abs(phi) > phi_b
            p_seed = float(np.mean(seed_mask))
            c_seed = max_cluster_size_periodic(seed_mask)

            Q_meas = np.empty(len(measure_indices))
            P_meas = np.empty(len(measure_indices))
            Cmax_meas = np.empty(len(measure_indices), dtype=float)

            meas_pos = 0
            for n in range(nsteps + 1):
                if n == measure_indices[meas_pos]:
                    ordered_mask = np.abs(phi) > ordered_threshold
                    Q_meas[meas_pos] = np.mean(phi**2)
                    P_meas[meas_pos] = np.mean(ordered_mask)
                    Cmax_meas[meas_pos] = max_cluster_size_periodic(ordered_mask)
                    meas_pos += 1
                    if meas_pos >= len(measure_indices):
                        break

                nonlinear = config.b * phi**3 - config.c * phi**5
                phi = semi_implicit_step(
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
                measure_times,
                P_meas,
                config.ordered_fraction_threshold,
                min_consecutive=config.min_consecutive,
            )
            t_cluster = persistent_crossing_time(
                measure_times,
                Cmax_meas,
                config.cluster_threshold,
                min_consecutive=config.min_consecutive,
            )
            t_tr = transition_time(measure_times, Q_meas, q_threshold)

            Q_all.append(Q_meas)
            P_all.append(P_meas)
            Cmax_all.append(Cmax_meas)
            t_nuc_area_all.append(t_area)
            t_nuc_cluster_all.append(t_cluster)
            t_tr_all.append(t_tr)
            p_seed_all.append(p_seed)
            c_seed_all.append(c_seed)

            rows.append({
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
                "p_seed": p_seed,
                "c_seed": c_seed,
                "t_nuc_area": t_area,
                "t_nuc_cluster": t_cluster,
                "t_tr": t_tr,
            })

        Q_all_arr = np.array(Q_all)
        P_all_arr = np.array(P_all)
        Cmax_all_arr = np.array(Cmax_all)
        t_cluster_arr = np.array(t_nuc_cluster_all)
        t_area_arr = np.array(t_nuc_area_all)
        t_tr_arr = np.array(t_tr_all)
        p_seed_arr = np.array(p_seed_all)
        c_seed_arr = np.array(c_seed_all)

        results[p.label] = {
            "Q_mean": np.nanmean(Q_all_arr, axis=0),
            "Q_std": np.nanstd(Q_all_arr, axis=0),
            "P_mean": np.nanmean(P_all_arr, axis=0),
            "P_std": np.nanstd(P_all_arr, axis=0),
            "Cmax_mean": np.nanmean(Cmax_all_arr, axis=0),
            "Cmax_std": np.nanstd(Cmax_all_arr, axis=0),
            "t_nuc_cluster_samples": t_cluster_arr,
            "t_nuc_area_samples": t_area_arr,
            "t_tr_samples": t_tr_arr,
            "t_nuc_cluster_mean": float(np.nanmean(t_cluster_arr)),
            "t_nuc_cluster_std": float(np.nanstd(t_cluster_arr)),
            "t_nuc_area_mean": float(np.nanmean(t_area_arr)),
            "t_nuc_area_std": float(np.nanstd(t_area_arr)),
            "t_tr_mean": float(np.nanmean(t_tr_arr)),
            "t_tr_std": float(np.nanstd(t_tr_arr)),
            "p_seed_mean": float(np.nanmean(p_seed_arr)),
            "p_seed_std": float(np.nanstd(p_seed_arr)),
            "c_seed_mean": float(np.nanmean(c_seed_arr)),
            "c_seed_std": float(np.nanstd(c_seed_arr)),
            "survival_cluster": survival_curve(t_cluster_arr, measure_times),
        }

    metadata = {
        **asdict(config),
        "phi_b": phi_b,
        "phi_s": phi_s,
        "ordered_threshold": ordered_threshold,
        "q_threshold": q_threshold,
        "measure_times": measure_times,
    }
    write_csv(config.output_csv, rows)
    return measure_times, results, metadata


# -----------------------------------------------------------------------------
# I/O and plotting helpers
# -----------------------------------------------------------------------------

def write_csv(path: str, rows: List[dict]) -> None:
    if not rows:
        return
    path_obj = Path(path)
    path_obj.parent.mkdir(parents=True, exist_ok=True)
    keys = list(rows[0].keys())
    with path_obj.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def print_phi4_summary(params: List[InitialParam], results: Dict[float, dict]) -> None:
    print("phi4 summary")
    print("label        t_tr_mean     t_tr_std")
    print("------------------------------------")
    for p in params:
        r = results[p.label]
        print(f"{p.label:8.3f}   {r['t_tr_mean']:10.4f}   {r['t_tr_std']:9.4f}")


def print_phi6_summary(params: List[InitialParam], results: Dict[float, dict], metadata: dict) -> None:
    print("phi6 summary")
    print(f"a_f={metadata['a_f']}, phi_b={metadata['phi_b']:.5f}, phi_s={metadata['phi_s']:.5f}")
    print("label      p_seed      c_seed    t_nuc_cl    t_nuc_sd      t_tr       t_tr_sd")
    print("----------------------------------------------------------------------------")
    for p in params:
        r = results[p.label]
        print(
            f"{p.label:8.3f}  "
            f"{r['p_seed_mean']:9.3e}  "
            f"{r['c_seed_mean']:8.2f}  "
            f"{r['t_nuc_cluster_mean']:10.3f}  "
            f"{r['t_nuc_cluster_std']:9.3f}  "
            f"{r['t_tr_mean']:9.3f}  "
            f"{r['t_tr_std']:9.3f}"
        )


def plot_phi4(times: Array, results: Dict[float, dict], q_threshold: float) -> None:
    plt.figure(figsize=(7, 5))
    for label, r in results.items():
        plt.plot(times, r["Q_mean"], label=f"label={label}")
    plt.axhline(q_threshold, linestyle="--", linewidth=1)
    plt.xlabel("time")
    plt.ylabel(r"$Q(t)=\langle\phi^2\rangle$")
    plt.legend()
    plt.tight_layout()
    plt.show()

    labels = list(results.keys())
    y = [results[x]["t_tr_mean"] for x in labels]
    err = [results[x]["t_tr_std"] for x in labels]
    plt.figure(figsize=(7, 5))
    plt.errorbar(labels, y, yerr=err, marker="o", capsize=3)
    plt.xlabel("initial-state label")
    plt.ylabel("transition time")
    plt.tight_layout()
    plt.show()


def plot_phi6(times: Array, results: Dict[float, dict], metadata: dict) -> None:
    plt.figure(figsize=(7, 5))
    for label, r in results.items():
        plt.plot(times, r["Q_mean"], label=f"label={label}")
    plt.axhline(metadata["q_threshold"], linestyle="--", linewidth=1)
    plt.xlabel("time")
    plt.ylabel(r"$Q(t)=\langle\phi^2\rangle$")
    plt.legend()
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(7, 5))
    for label, r in results.items():
        plt.plot(times, r["survival_cluster"], label=f"label={label}")
    plt.xlabel("time")
    plt.ylabel(r"$S_{surv}(t)=Pr(t_{nuc}>t)$")
    plt.legend()
    plt.tight_layout()
    plt.show()

    labels = list(results.keys())
    t_nuc = [results[x]["t_nuc_cluster_mean"] for x in labels]
    t_nuc_err = [results[x]["t_nuc_cluster_std"] for x in labels]
    t_tr = [results[x]["t_tr_mean"] for x in labels]
    t_tr_err = [results[x]["t_tr_std"] for x in labels]

    plt.figure(figsize=(7, 5))
    plt.errorbar(labels, t_nuc, yerr=t_nuc_err, marker="o", capsize=3, label="cluster nucleation time")
    plt.errorbar(labels, t_tr, yerr=t_tr_err, marker="s", capsize=3, label="transition time")
    plt.xlabel("initial-state label")
    plt.ylabel("time")
    plt.legend()
    plt.tight_layout()
    plt.show()

    p_seed = [results[x]["p_seed_mean"] for x in labels]
    p_seed_err = [results[x]["p_seed_std"] for x in labels]
    plt.figure(figsize=(7, 5))
    plt.errorbar(labels, p_seed, yerr=p_seed_err, marker="o", capsize=3)
    plt.xlabel("initial-state label")
    plt.ylabel("initial seed fraction")
    plt.tight_layout()
    plt.show()


# -----------------------------------------------------------------------------
# Example runs
# -----------------------------------------------------------------------------

def example_phi4() -> None:
    """Small exploratory phi4 run. Increase nsamples/N for production."""
    Ti_list = (1.05, 1.10, 1.20, 1.50, 2.00, 3.00)
    params = make_phi4_initial_params(Ti_list, D0=0.05)
    config = Phi4Config(
        N=64,
        dt=0.05,
        tmax=10.0,
        preeq_steps=500,
        D_f=0.0,
        nsamples=5,
        init_mode="preeq",
        output_csv="results/phi4_samples.csv",
    )
    times, results = run_phi4_experiment(config, params)
    print_phi4_summary(params, results)
    plot_phi4(times, results, config.q_threshold)


def example_phi6() -> None:
    """Small exploratory phi6 run. Increase nsamples/N/tmax for production."""
    Ti_list = (1.05, 1.10, 1.20, 1.50, 2.00, 3.00)
    params = make_phi6_initial_params(Ti_list, D0=0.02, scheme="noise_only", a_i_base=0.30)
    config = Phi6Config(
        N=64,
        dt=0.02,
        tmax=50.0,
        preeq_steps=500,
        a_f=0.10,
        D_f=1.0e-4,
        nsamples=3,
        init_mode="preeq",
        measure_every=10,
        cluster_threshold=10,
        output_csv="results/phi6_samples.csv",
    )
    times, results, metadata = run_phi6_experiment(config, params)
    print_phi6_summary(params, results, metadata)
    plot_phi6(times, results, metadata)


if __name__ == "__main__":
    # Choose one. These are intentionally small smoke-test settings.
    # For production, increase N, nsamples, preeq_steps, and tmax.
    # example_phi4()
    example_phi6()
