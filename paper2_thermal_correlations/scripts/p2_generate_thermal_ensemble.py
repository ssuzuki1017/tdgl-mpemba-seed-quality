#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Paper 2 smoke test:
Generate thermal initial ensembles using FDT-consistent stochastic Model-A TDGL.

Initial free energy:
    F_i[phi] = ∫ d^2x [ 1/2 |grad phi|^2 + r_i/2 phi^2 + u_i/4 phi^4 ]

Dynamics:
    d phi / dt = Gamma (laplacian phi - r_i phi - u_i phi^3) + eta

Noise:
    <eta(x,t) eta(x',t')> = 2 Gamma T_i delta(x-x') delta(t-t')

Discrete real-space noise increment:
    sqrt(2 Gamma T_i dt / dx^2) * N(0,1)

Outputs:
    - initial thermal fields
    - sample-level summary CSV
    - radial correlation function C(r)
    - radial structure factor S(k)
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import deque
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


def load_config(path: Path) -> dict:
    with path.open("r", encoding="utf-8-sig") as f:
        return json.load(f)


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def phi6_barrier_location(a: float, b: float, c: float) -> float:
    disc = b * b - 4.0 * a * c
    if disc <= 0:
        raise ValueError("phi6 barrier does not exist: b^2 - 4ac <= 0")
    phi_b_sq = (b - math.sqrt(disc)) / (2.0 * c)
    return math.sqrt(phi_b_sq)


def make_k2_grid(N: int, dx: float) -> np.ndarray:
    k = 2.0 * np.pi * np.fft.fftfreq(N, d=dx)
    kx, ky = np.meshgrid(k, k, indexing="ij")
    return kx * kx + ky * ky


def tdgl_model_a_step(
    phi: np.ndarray,
    r_i: float,
    u_i: float,
    T_i: float,
    gamma: float,
    dt: float,
    dx: float,
    k2: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Semi-implicit update for:
        dphi/dt = Gamma (laplacian phi - r_i phi - u_i phi^3) + eta

    Linear part is treated implicitly in Fourier space.
    Nonlinear part and noise are treated explicitly in real space.
    """
    d = 2
    noise_std = math.sqrt(2.0 * gamma * T_i * dt / (dx ** d))
    noise_increment = noise_std * rng.normal(size=phi.shape)

    nonlinear_increment = dt * (-gamma * u_i * phi ** 3)

    numerator_real = phi + nonlinear_increment + noise_increment
    numerator_k = np.fft.fft2(numerator_real)

    denom = 1.0 + dt * gamma * (k2 + r_i)
    phi_new = np.fft.ifft2(numerator_k / denom).real
    return phi_new


def precompute_realspace_radial_bins(N: int, dx: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    ii, jj = np.meshgrid(np.arange(N), np.arange(N), indexing="ij")
    di = np.minimum(ii, N - ii) * dx
    dj = np.minimum(jj, N - jj) * dx
    rr = np.sqrt(di * di + dj * dj)
    rbin = np.rint(rr / dx).astype(int)
    maxbin = int(rbin.max())
    counts = np.bincount(rbin.ravel(), minlength=maxbin + 1)
    r_values = np.arange(maxbin + 1) * dx
    return rbin, counts, r_values


def radial_autocorrelation(phi: np.ndarray, rbin: np.ndarray, counts: np.ndarray) -> np.ndarray:
    N = phi.shape[0]
    f = phi - float(np.mean(phi))
    fhat = np.fft.fft2(f)

    # ifft2(|F|^2) gives sum_x f(x) f(x+r); divide by N^2 for spatial average.
    ac2 = np.fft.ifft2(np.abs(fhat) ** 2).real / (N * N)

    weighted = np.bincount(rbin.ravel(), weights=ac2.ravel(), minlength=len(counts))
    out = np.full(len(counts), np.nan, dtype=float)
    valid = counts > 0
    out[valid] = weighted[valid] / counts[valid]
    return out


def precompute_kspace_radial_bins(N: int, dx: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    k = 2.0 * np.pi * np.fft.fftfreq(N, d=dx)
    kx, ky = np.meshgrid(k, k, indexing="ij")
    kk = np.sqrt(kx * kx + ky * ky)
    dk = 2.0 * np.pi / (N * dx)
    kbin = np.rint(kk / dk).astype(int)
    maxbin = int(kbin.max())
    counts = np.bincount(kbin.ravel(), minlength=maxbin + 1)
    k_values = np.arange(maxbin + 1) * dk
    return kbin, counts, k_values


def radial_structure_factor(phi: np.ndarray, kbin: np.ndarray, counts: np.ndarray) -> np.ndarray:
    N = phi.shape[0]
    f = phi - float(np.mean(phi))
    fhat = np.fft.fft2(f)

    # This normalization is sufficient for comparing ensembles in this smoke stage.
    s2 = (np.abs(fhat) ** 2) / (N * N)

    weighted = np.bincount(kbin.ravel(), weights=s2.ravel(), minlength=len(counts))
    out = np.full(len(counts), np.nan, dtype=float)
    valid = counts > 0
    out[valid] = weighted[valid] / counts[valid]
    return out


def estimate_xi_from_C(r_values: np.ndarray, C: np.ndarray, max_fraction: float = 0.25) -> float:
    """
    Crude exponential-fit correlation length:
        C(r)/C(0) ~ exp(-r/xi)

    This is only a smoke-test estimator.
    Later we can replace it with a more careful second-moment or fit-range-controlled estimator.
    """
    if len(C) < 5 or not np.isfinite(C[0]) or C[0] <= 0:
        return float("nan")

    y = C / C[0]
    rmax = r_values[-1] * max_fraction

    mask = (
        (r_values > 0)
        & (r_values <= rmax)
        & np.isfinite(y)
        & (y > 0.03)
        & (y < 0.95)
    )

    if int(np.sum(mask)) < 3:
        return float("nan")

    x = r_values[mask]
    logy = np.log(y[mask])

    slope, intercept = np.polyfit(x, logy, 1)
    if slope >= 0:
        return float("nan")

    return float(-1.0 / slope)


def cluster_metrics_periodic(mask: np.ndarray, dx: float) -> Dict[str, float]:
    """
    4-neighbor periodic connected components.

    Returns metrics for the whole mask and for the largest barrier-crossing cluster.
    """
    N = mask.shape[0]
    visited = np.zeros_like(mask, dtype=bool)

    total_sites = int(np.sum(mask))
    pseed = total_sites / float(N * N)

    if total_sites == 0:
        return {
            "pseed": pseed,
            "n_clusters": 0,
            "largest_cluster_size": 0,
            "largest_compactness": float("nan"),
            "largest_rg": float("nan"),
            "largest_perimeter": float("nan"),
            "largest_perimeter_area": float("nan"),
        }

    directions = [(1, 0), (-1, 0), (0, 1), (0, -1)]

    n_clusters = 0
    largest_size = 0
    largest_original_coords: List[Tuple[int, int]] = []
    largest_unwrapped_coords: List[Tuple[int, int]] = []

    for i0 in range(N):
        for j0 in range(N):
            if not mask[i0, j0] or visited[i0, j0]:
                continue

            n_clusters += 1

            q = deque()
            q.append((i0, j0, 0, 0))
            visited[i0, j0] = True

            original_coords: List[Tuple[int, int]] = []
            unwrapped_coords: List[Tuple[int, int]] = []

            while q:
                i, j, ui, uj = q.popleft()
                original_coords.append((i, j))
                unwrapped_coords.append((ui, uj))

                for di, dj in directions:
                    ni = (i + di) % N
                    nj = (j + dj) % N
                    if mask[ni, nj] and not visited[ni, nj]:
                        visited[ni, nj] = True
                        q.append((ni, nj, ui + di, uj + dj))

            size = len(original_coords)
            if size > largest_size:
                largest_size = size
                largest_original_coords = original_coords
                largest_unwrapped_coords = unwrapped_coords

    if largest_size == 0:
        compactness = float("nan")
        rg = float("nan")
        perimeter = float("nan")
        perimeter_area = float("nan")
    else:
        unwrap = np.array(largest_unwrapped_coords, dtype=float)
        umin, vmin = np.min(unwrap, axis=0)
        umax, vmax = np.max(unwrap, axis=0)
        bbox_area_sites = (umax - umin + 1.0) * (vmax - vmin + 1.0)
        compactness = float(largest_size / bbox_area_sites) if bbox_area_sites > 0 else float("nan")

        center = np.mean(unwrap, axis=0)
        rg = float(np.sqrt(np.mean(np.sum((unwrap - center) ** 2, axis=1))) * dx)

        comp_mask = np.zeros_like(mask, dtype=bool)
        for i, j in largest_original_coords:
            comp_mask[i, j] = True

        perimeter_edges = 0
        for i, j in largest_original_coords:
            for di, dj in directions:
                ni = (i + di) % N
                nj = (j + dj) % N
                if not comp_mask[ni, nj]:
                    perimeter_edges += 1

        perimeter = float(perimeter_edges * dx)
        area = float(largest_size * dx * dx)
        perimeter_area = float(perimeter / area) if area > 0 else float("nan")

    return {
        "pseed": pseed,
        "n_clusters": n_clusters,
        "largest_cluster_size": largest_size,
        "largest_compactness": compactness,
        "largest_rg": rg,
        "largest_perimeter": perimeter,
        "largest_perimeter_area": perimeter_area,
    }


def write_rows_csv(path: Path, rows: List[dict], fieldnames: List[str]) -> None:
    ensure_dir(path.parent)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def save_radial_mean_csv(
    path: Path,
    control_name: str,
    control_value: float,
    x_name: str,
    x_values: np.ndarray,
    curves: List[np.ndarray],
    y_name: str,
) -> None:
    ensure_dir(path.parent)

    arr = np.array(curves, dtype=float)
    mean = np.nanmean(arr, axis=0)
    if arr.shape[0] > 1:
        sem = np.nanstd(arr, axis=0, ddof=1) / math.sqrt(arr.shape[0])
    else:
        sem = np.full_like(mean, np.nan)

    rows = []
    for x, m, s in zip(x_values, mean, sem):
        rows.append({
            control_name: control_value,
            x_name: float(x),
            f"{y_name}_mean": float(m),
            f"{y_name}_sem": float(s),
            "n_samples": int(arr.shape[0]),
        })

    fieldnames = [control_name, x_name, f"{y_name}_mean", f"{y_name}_sem", "n_samples"]

    mode = "a" if path.exists() else "w"
    with path.open(mode, newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if mode == "w":
            writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="Path to JSON config.")
    args = parser.parse_args()

    config_path = Path(args.config)
    cfg = load_config(config_path)

    run_id = str(cfg["run_id"])
    seed = int(cfg["seed"])
    rng = np.random.default_rng(seed)

    N = int(cfg["N"])
    dx = float(cfg["dx"])
    dt = float(cfg["dt"])
    n_pre_steps = int(cfg["n_pre_steps"])
    n_samples = int(cfg["n_samples"])

    gamma = float(cfg["gamma"])
    u_i = float(cfg["u_i"])
    T_i = float(cfg["T_i"])
    r_i_values = [float(x) for x in cfg["r_i_values"]]
    initial_std = float(cfg["initial_std"])
    save_fields = bool(cfg.get("save_fields", True))

    final_phi6 = cfg["final_phi6"]
    phi_b = phi6_barrier_location(
        a=float(final_phi6["a_f"]),
        b=float(final_phi6["b"]),
        c=float(final_phi6["c"]),
    )

    out_root = Path(cfg["out_dir"]) / run_id
    processed_root = Path(cfg["processed_dir"]) / run_id
    ensure_dir(out_root)
    ensure_dir(processed_root)

    # Save effective config.
    with (processed_root / "config_effective.json").open("w", encoding="utf-8-sig") as f:
        json.dump(cfg, f, indent=2)

    k2 = make_k2_grid(N, dx)

    rbin_real, counts_real, r_values = precompute_realspace_radial_bins(N, dx)
    kbin, counts_k, k_values = precompute_kspace_radial_bins(N, dx)

    summary_rows: List[dict] = []

    C_mean_path = processed_root / "C_r_mean.csv"
    S_mean_path = processed_root / "S_k_mean.csv"
    if C_mean_path.exists():
        C_mean_path.unlink()
    if S_mean_path.exists():
        S_mean_path.unlink()

    print("=== Paper 2 thermal ensemble smoke test ===")
    print(f"run_id        : {run_id}")
    print(f"N             : {N}")
    print(f"n_samples     : {n_samples}")
    print(f"n_pre_steps   : {n_pre_steps}")
    print(f"T_i           : {T_i}")
    print(f"phi_b final   : {phi_b:.6f}")
    print("")

    for r_i in r_i_values:
        print(f"[r_i={r_i:.4f}] generating samples...")

        r_label = f"ri_{r_i:.4f}".replace(".", "p")
        field_dir = out_root / r_label
        ensure_dir(field_dir)

        C_curves: List[np.ndarray] = []
        S_curves: List[np.ndarray] = []

        for sample_idx in range(n_samples):
            phi = initial_std * rng.normal(size=(N, N))

            for _ in range(n_pre_steps):
                phi = tdgl_model_a_step(
                    phi=phi,
                    r_i=r_i,
                    u_i=u_i,
                    T_i=T_i,
                    gamma=gamma,
                    dt=dt,
                    dx=dx,
                    k2=k2,
                    rng=rng,
                )

            C = radial_autocorrelation(phi, rbin_real, counts_real)
            S = radial_structure_factor(phi, kbin, counts_k)
            xi_est = estimate_xi_from_C(r_values, C)

            C_curves.append(C)
            S_curves.append(S)

            mask = np.abs(phi) > phi_b
            metrics = cluster_metrics_periodic(mask, dx=dx)

            row = {
                "run_id": run_id,
                "r_i": r_i,
                "T_i": T_i,
                "sample": sample_idx,
                "N": N,
                "dx": dx,
                "dt": dt,
                "n_pre_steps": n_pre_steps,
                "phi_mean": float(np.mean(phi)),
                "phi_var": float(np.var(phi)),
                "phi_std": float(np.std(phi)),
                "phi_abs_mean": float(np.mean(np.abs(phi))),
                "phi_b": phi_b,
                "xi_expfit": xi_est,
                **metrics,
            }
            summary_rows.append(row)

            if save_fields:
                np.save(field_dir / f"field_{sample_idx:04d}.npy", phi.astype(np.float32))

        save_radial_mean_csv(
            path=C_mean_path,
            control_name="r_i",
            control_value=r_i,
            x_name="r",
            x_values=r_values,
            curves=C_curves,
            y_name="C",
        )

        save_radial_mean_csv(
            path=S_mean_path,
            control_name="r_i",
            control_value=r_i,
            x_name="k",
            x_values=k_values,
            curves=S_curves,
            y_name="S",
        )

        pseed_vals = [row["pseed"] for row in summary_rows if row["r_i"] == r_i]
        xi_vals = [row["xi_expfit"] for row in summary_rows if row["r_i"] == r_i and np.isfinite(row["xi_expfit"])]
        comp_vals = [row["largest_compactness"] for row in summary_rows if row["r_i"] == r_i and np.isfinite(row["largest_compactness"])]

        pseed_mean = float(np.mean(pseed_vals)) if pseed_vals else float("nan")
        xi_mean = float(np.mean(xi_vals)) if xi_vals else float("nan")
        comp_mean = float(np.mean(comp_vals)) if comp_vals else float("nan")

        print(
            f"  pseed_mean={pseed_mean:.4f}, "
            f"xi_expfit_mean={xi_mean:.4f}, "
            f"compactness_mean={comp_mean:.4f}"
        )

    summary_path = processed_root / "sample_summary.csv"
    fieldnames = [
        "run_id",
        "r_i",
        "T_i",
        "sample",
        "N",
        "dx",
        "dt",
        "n_pre_steps",
        "phi_mean",
        "phi_var",
        "phi_std",
        "phi_abs_mean",
        "phi_b",
        "xi_expfit",
        "pseed",
        "n_clusters",
        "largest_cluster_size",
        "largest_compactness",
        "largest_rg",
        "largest_perimeter",
        "largest_perimeter_area",
    ]
    write_rows_csv(summary_path, summary_rows, fieldnames)

    print("")
    print("Done.")
    print(f"Sample summary: {summary_path}")
    print(f"C(r) mean     : {C_mean_path}")
    print(f"S(k) mean     : {S_mean_path}")
    print(f"Fields        : {out_root}")


if __name__ == "__main__":
    main()
