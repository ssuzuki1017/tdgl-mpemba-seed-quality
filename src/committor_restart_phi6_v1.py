#!/usr/bin/env python
r"""
Committor-style restart test for the TDGL phi6 nucleation problem.

Purpose
-------
For each initial configuration, this script freezes the post-preparation field
phi(0), then restarts the post-quench dynamics many times with independent
post-quench noise. The estimated committor-like quantity is

    q_nuc(phi0) = fraction of restarts that nucleate by tmax.

This is not an exact transition-path committor, because it uses a finite
observation time and an operational cluster-nucleation criterion. However, it is
a much more direct test than label-averaged seed statistics.

Run from the repository root, after ensuring that these files are in src/:
  - tdgl_mpemba_revised.py
  - tdgl_mpemba_timeseries_export_v3.py

Smoke test:
  python .\src\committor_restart_phi6_v1.py --outdir .\analysis_committor_smoke --N 32 --tmax 20 --preeq_steps 20 --n-configs-per-label 1 --n-restarts 2 --labels 1.05 1.50

Production-lite:
  python .\src\committor_restart_phi6_v1.py --outdir .\analysis_committor_v1 --n-configs-per-label 8 --n-restarts 16

Heavier production:
  python .\src\committor_restart_phi6_v1.py --outdir .\analysis_committor_v2 --n-configs-per-label 12 --n-restarts 24
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

import matplotlib.pyplot as plt


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import tdgl_mpemba_revised as base  # noqa: E402


DEFAULT_LABELS = (1.05, 1.10, 1.20, 1.50, 2.00, 3.00)


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def largest_cluster_geometry_periodic(mask: np.ndarray) -> dict[str, float]:
    """Geometry of the largest 4-neighbor connected cluster in a periodic 2D mask."""
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

            if len(coords) > best_area:
                best_area = len(coords)
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
    rg = float(np.sqrt(np.mean(np.sum((arr - center) ** 2, axis=1))))

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
    return {
        "area": float(best_area),
        "bbox_area": bbox_area,
        "compactness": compactness,
        "rg": rg,
        "perimeter": perimeter,
        "perimeter_to_area": float(perimeter / best_area) if best_area > 0 else np.nan,
    }


def online_nucleation_time(
    phi0: np.ndarray,
    *,
    k2: np.ndarray,
    dt: float,
    tmax: float,
    measure_every: int,
    a_f: float,
    b: float,
    c: float,
    D_f: float,
    rng: np.random.Generator,
    ordered_threshold: float,
    cluster_threshold: int,
    min_consecutive: int,
    dx: float,
) -> float:
    """Run one post-quench restart and return first persistent cluster-nucleation time."""
    phi = np.array(phi0, copy=True)
    nsteps = int(round(tmax / dt))
    consecutive = 0

    for n in range(nsteps + 1):
        if n % measure_every == 0:
            ordered_mask = np.abs(phi) > ordered_threshold
            cmax = base.max_cluster_size_periodic(ordered_mask)
            if cmax >= cluster_threshold:
                consecutive += 1
                if consecutive >= min_consecutive:
                    first_n = n - (min_consecutive - 1) * measure_every
                    return float(first_n * dt)
            else:
                consecutive = 0

        if n >= nsteps:
            break

        nonlinear = b * phi**3 - c * phi**5
        phi = base.semi_implicit_step(
            phi=phi,
            k2=k2,
            dt=dt,
            dx=dx,
            linear_alpha=-a_f,
            nonlinear=nonlinear,
            D=D_f,
            rng=rng,
        )

    return np.nan


def safe_corr(x: np.ndarray, y: np.ndarray, method: str) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    if method == "spearman":
        x = pd.Series(x).rank(method="average").to_numpy(dtype=float)
        y = pd.Series(y).rank(method="average").to_numpy(dtype=float)
    return float(np.corrcoef(x, y)[0, 1])


def stderr(x: pd.Series) -> float:
    y = pd.to_numeric(x, errors="coerce").dropna()
    if len(y) <= 1:
        return np.nan
    return float(y.std(ddof=1) / math.sqrt(len(y)))


def save_figures(config_summary: pd.DataFrame, label_summary: pd.DataFrame, outdir: Path, dpi: int) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.9), constrained_layout=True)

    ax = axes[0, 0]
    ax.errorbar(
        label_summary["label"],
        label_summary["q_nuc_mean"],
        yerr=label_summary["q_nuc_se"],
        marker="o",
        capsize=3,
    )
    ax.set_xlabel("Initial-state label")
    ax.set_ylabel(r"Committor-like probability $q_{\rm nuc}$")
    ax.set_title("(a) label mean")
    ax.set_ylim(-0.05, 1.05)

    panels = [
        ("seed_compactness", r"Initial seed compactness $C_{\rm comp}$", "(b)"),
        ("seed_rg", r"Initial seed radius of gyration $R_g$", "(c)"),
        ("seed_perimeter_to_area", r"Initial seed perimeter / area", "(d)"),
    ]

    for ax, (xcol, xlabel, title) in zip(axes.flat[1:], panels):
        ax.scatter(config_summary[xcol], config_summary["q_nuc"], s=28)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"$q_{\rm nuc}$")
        ax.set_ylim(-0.05, 1.05)
        ax.set_title(title)

    for ext in ["png", "pdf"]:
        path = outdir / f"figS7_committor_seed_geometry.{ext}"
        fig.savefig(path, dpi=dpi if ext == "png" else None, bbox_inches="tight")
        print(f"[saved] {path}")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Committor-style restart test for phi6 TDGL nucleation.")
    parser.add_argument("--outdir", default="analysis_committor_v1")
    parser.add_argument("--labels", nargs="*", type=float, default=list(DEFAULT_LABELS))

    parser.add_argument("--N", type=int, default=64)
    parser.add_argument("--dx", type=float, default=1.0)
    parser.add_argument("--dt", type=float, default=0.02)
    parser.add_argument("--tmax", type=float, default=300.0)
    parser.add_argument("--preeq_steps", type=int, default=500)
    parser.add_argument("--a_f", "--af", dest="a_f", type=float, default=0.02)
    parser.add_argument("--b", type=float, default=1.0)
    parser.add_argument("--c", type=float, default=1.0)
    parser.add_argument("--D_f", "--Df", dest="D_f", type=float, default=0.009)
    parser.add_argument("--D0", type=float, default=0.02)
    parser.add_argument("--init_scheme", choices=["noise_only", "mass_and_noise"], default="noise_only")
    parser.add_argument("--a_i_base", type=float, default=0.30)
    parser.add_argument("--measure_every", type=int, default=5)
    parser.add_argument("--cluster_threshold", type=int, default=20)
    parser.add_argument("--min_consecutive", type=int, default=3)

    parser.add_argument("--n-configs-per-label", type=int, default=8)
    parser.add_argument("--n-restarts", type=int, default=16)
    parser.add_argument("--seed", type=int, default=424242)
    parser.add_argument("--dpi", type=int, default=300)

    args = parser.parse_args()

    outdir = Path(args.outdir).resolve()
    ensure_dir(outdir)

    labels = tuple(float(x) for x in args.labels)
    params = base.make_phi6_initial_params(
        labels,
        D0=args.D0,
        scheme=args.init_scheme,
        a_i_base=args.a_i_base,
    )

    rng_init = np.random.default_rng(args.seed)
    N, dx, dt = args.N, args.dx, args.dt
    k2 = base.make_k2(N, dx)

    phi_b, phi_s = base.extrema_phi6(args.a_f, args.b, args.c)
    if not np.isfinite(phi_s):
        raise ValueError("No stable nonzero phi6 extremum. Check a_f, b, c.")

    ordered_threshold = 0.5 * phi_s

    config_rows = []
    restart_rows = []

    total = len(params) * args.n_configs_per_label * args.n_restarts
    done = 0

    for ilabel, p in enumerate(params):
        for config_id in range(args.n_configs_per_label):
            phi0 = base.pre_equilibrate_phi6(
                N=N,
                dx=dx,
                dt=dt,
                nsteps=args.preeq_steps,
                a_i=p.a_i,
                D_i=p.D_i,
                rng=rng_init,
                b=args.b,
                c=args.c,
                remove_mean_at_end=True,
            )

            seed_mask = np.abs(phi0) > phi_b
            seed_geom = largest_cluster_geometry_periodic(seed_mask)
            p_seed = float(np.mean(seed_mask))

            ordered_seed_mask = np.abs(phi0) > ordered_threshold
            ordered_seed_geom = largest_cluster_geometry_periodic(ordered_seed_mask)
            p_ordered_seed = float(np.mean(ordered_seed_mask))

            restart_tnuc = []
            for restart_id in range(args.n_restarts):
                # Deterministic, independent seed for each restart.
                restart_seed = (
                    int(args.seed)
                    + 1000003 * (ilabel + 1)
                    + 1009 * (config_id + 1)
                    + 37 * (restart_id + 1)
                )
                rng_post = np.random.default_rng(restart_seed)

                t_nuc = online_nucleation_time(
                    phi0,
                    k2=k2,
                    dt=dt,
                    tmax=args.tmax,
                    measure_every=args.measure_every,
                    a_f=args.a_f,
                    b=args.b,
                    c=args.c,
                    D_f=args.D_f,
                    rng=rng_post,
                    ordered_threshold=ordered_threshold,
                    cluster_threshold=args.cluster_threshold,
                    min_consecutive=args.min_consecutive,
                    dx=dx,
                )
                restart_tnuc.append(t_nuc)
                restart_rows.append({
                    "label": p.label,
                    "config_id": config_id,
                    "restart_id": restart_id,
                    "restart_seed": restart_seed,
                    "t_nuc_cluster": t_nuc,
                    "nucleated": int(np.isfinite(t_nuc)),
                    "p_seed": p_seed,
                    "c_seed": seed_geom["area"],
                    "seed_compactness": seed_geom["compactness"],
                    "seed_rg": seed_geom["rg"],
                    "seed_perimeter": seed_geom["perimeter"],
                    "seed_perimeter_to_area": seed_geom["perimeter_to_area"],
                    "p_ordered_seed": p_ordered_seed,
                    "c_ordered_seed": ordered_seed_geom["area"],
                    "ordered_seed_compactness": ordered_seed_geom["compactness"],
                    "ordered_seed_rg": ordered_seed_geom["rg"],
                    "ordered_seed_perimeter": ordered_seed_geom["perimeter"],
                    "ordered_seed_perimeter_to_area": ordered_seed_geom["perimeter_to_area"],
                })

                done += 1
                if done % max(1, total // 20) == 0 or done == total:
                    print(f"[progress] {done}/{total} restarts")

            finite = np.isfinite(restart_tnuc)
            config_rows.append({
                "label": p.label,
                "config_id": config_id,
                "n_restarts": args.n_restarts,
                "n_nucleated": int(np.sum(finite)),
                "q_nuc": float(np.mean(finite)),
                "mean_t_nuc_observed": float(np.nanmean(restart_tnuc)) if np.any(finite) else np.nan,
                "p_seed": p_seed,
                "c_seed": seed_geom["area"],
                "seed_compactness": seed_geom["compactness"],
                "seed_bbox_area": seed_geom["bbox_area"],
                "seed_rg": seed_geom["rg"],
                "seed_perimeter": seed_geom["perimeter"],
                "seed_perimeter_to_area": seed_geom["perimeter_to_area"],
                "p_ordered_seed": p_ordered_seed,
                "c_ordered_seed": ordered_seed_geom["area"],
                "ordered_seed_compactness": ordered_seed_geom["compactness"],
                "ordered_seed_bbox_area": ordered_seed_geom["bbox_area"],
                "ordered_seed_rg": ordered_seed_geom["rg"],
                "ordered_seed_perimeter": ordered_seed_geom["perimeter"],
                "ordered_seed_perimeter_to_area": ordered_seed_geom["perimeter_to_area"],
            })

            # Incremental saves make long runs safer.
            pd.DataFrame(config_rows).to_csv(outdir / "committor_initial_configs.csv", index=False)
            pd.DataFrame(restart_rows).to_csv(outdir / "committor_restarts.csv", index=False)

    configs = pd.DataFrame(config_rows)
    restarts = pd.DataFrame(restart_rows)

    label_rows = []
    for label, g in configs.groupby("label"):
        label_rows.append({
            "label": float(label),
            "n_configs": int(len(g)),
            "n_restarts_total": int(g["n_restarts"].sum()),
            "q_nuc_mean": float(g["q_nuc"].mean()),
            "q_nuc_se": stderr(g["q_nuc"]),
            "q_nuc_min": float(g["q_nuc"].min()),
            "q_nuc_max": float(g["q_nuc"].max()),
            "seed_compactness_mean": float(g["seed_compactness"].mean()),
            "seed_rg_mean": float(g["seed_rg"].mean()),
            "seed_perimeter_to_area_mean": float(g["seed_perimeter_to_area"].mean()),
        })
    label_summary = pd.DataFrame(label_rows).sort_values("label")
    label_summary.to_csv(outdir / "committor_label_summary.csv", index=False)

    metrics = [
        "p_seed",
        "c_seed",
        "seed_compactness",
        "seed_rg",
        "seed_perimeter",
        "seed_perimeter_to_area",
        "p_ordered_seed",
        "c_ordered_seed",
    ]
    corr_rows = []
    for metric in metrics:
        corr_rows.append({
            "metric": metric,
            "pearson_with_q_nuc": safe_corr(configs[metric].to_numpy(), configs["q_nuc"].to_numpy(), "pearson"),
            "spearman_with_q_nuc": safe_corr(configs[metric].to_numpy(), configs["q_nuc"].to_numpy(), "spearman"),
            "n_configs": int(len(configs)),
        })
    pd.DataFrame(corr_rows).to_csv(outdir / "committor_geometry_correlations.csv", index=False)

    save_figures(configs, label_summary, outdir, dpi=args.dpi)

    readme = f"""# Committor-style restart test

This directory contains a finite-time committor-style restart analysis.

Parameters:

```text
N = {args.N}
dt = {args.dt}
tmax = {args.tmax}
preeq_steps = {args.preeq_steps}
a_f = {args.a_f}
D_f = {args.D_f}
D0 = {args.D0}
cluster_threshold = {args.cluster_threshold}
min_consecutive = {args.min_consecutive}
n_configs_per_label = {args.n_configs_per_label}
n_restarts = {args.n_restarts}
labels = {labels}
```

Definition:

```text
q_nuc(phi0) = fraction of independent post-quench restarts that nucleate by tmax
```

Outputs:

- `committor_initial_configs.csv`
- `committor_restarts.csv`
- `committor_label_summary.csv`
- `committor_geometry_correlations.csv`
- `figS7_committor_seed_geometry.png`
- `figS7_committor_seed_geometry.pdf`

Interpretation:

This is not an exact transition-path committor because the observation time is finite and the nucleation event is defined by an operational cluster threshold. It is nevertheless a direct restart test of whether initial seed geometry contains information about finite-time nucleation outcomes.
"""
    (outdir / "README_committor_restart_test.md").write_text(readme, encoding="utf-8")

    print(f"\nDone. Outputs written to: {outdir}")


if __name__ == "__main__":
    main()
