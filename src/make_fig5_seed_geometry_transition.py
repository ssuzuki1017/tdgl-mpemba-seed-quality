#!/usr/bin/env python
r"""
Generate Fig. 5 for the TDGL Mpemba seed-geometry manuscript.

Fig. 5 is intended to make the core mechanism visible in the main text:

(a) sample-level relation between seed amount and seed compactness;
(b) direct relation between cluster-nucleation time and global transition time.

Usage from the repository root:

  python .\src\make_fig5_seed_geometry_transition.py `
    --data-dir .\data `
    --outdir .\figures\main

Outputs:

  fig5_seed_geometry_transition.pdf
  fig5_seed_geometry_transition.png
  fig5_seed_geometry_transition_summary.csv
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

import matplotlib.pyplot as plt


def find_main_files(data_dir: Path) -> list[Path]:
    main_dir = data_dir / "main"
    files = []
    for seed in (5678, 9876):
        preferred = main_dir / f"phi6_main_Df9em3_seed{seed}_fullcols.samples.csv"
        if preferred.exists():
            files.append(preferred)
        else:
            candidates = sorted(main_dir.glob(f"*seed{seed}*fullcols*.samples.csv"))
            if candidates:
                files.append(candidates[0])
    if len(files) != 2:
        raise FileNotFoundError(
            "Could not find both full-column main files for seed 5678 and seed 9876 under data/main."
        )
    return files


def load_main_data(data_dir: Path) -> pd.DataFrame:
    frames = []
    for p in find_main_files(data_dir):
        df = pd.read_csv(p)
        df["source_file"] = p.name
        frames.append(df)
        print(f"[found] {p}")
    out = pd.concat(frames, ignore_index=True)

    required = ["label", "p_seed", "c_seed", "seed_compactness", "seed_rg", "t_nuc_cluster", "t_tr"]
    missing = [c for c in required if c not in out.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    for col in required:
        out[col] = pd.to_numeric(out[col], errors="coerce")

    out["nucleated"] = out["t_nuc_cluster"].notna()
    out["transitioned"] = out["t_tr"].notna()
    return out


def safe_corr(x: pd.Series, y: pd.Series, method: str = "pearson") -> float:
    a = pd.to_numeric(x, errors="coerce")
    b = pd.to_numeric(y, errors="coerce")
    ok = a.notna() & b.notna()
    a = a[ok].to_numpy(dtype=float)
    b = b[ok].to_numpy(dtype=float)
    if len(a) < 3 or np.std(a) == 0.0 or np.std(b) == 0.0:
        return np.nan
    if method == "spearman":
        a = pd.Series(a).rank(method="average").to_numpy(dtype=float)
        b = pd.Series(b).rank(method="average").to_numpy(dtype=float)
    return float(np.corrcoef(a, b)[0, 1])


def make_summary(df: pd.DataFrame) -> pd.DataFrame:
    both = df[df["t_nuc_cluster"].notna() & df["t_tr"].notna()].copy()
    rows = []

    rows.append({
        "quantity": "sample_count_total",
        "value": len(df),
    })
    rows.append({
        "quantity": "sample_count_nucleated",
        "value": int(df["t_nuc_cluster"].notna().sum()),
    })
    rows.append({
        "quantity": "sample_count_transitioned",
        "value": int(df["t_tr"].notna().sum()),
    })
    rows.append({
        "quantity": "sample_count_tnuc_and_ttr_observed",
        "value": len(both),
    })

    for xcol, ycol, name in [
        ("p_seed", "seed_compactness", "p_seed_vs_seed_compactness"),
        ("c_seed", "seed_compactness", "c_seed_vs_seed_compactness"),
        ("seed_rg", "seed_compactness", "seed_rg_vs_seed_compactness"),
    ]:
        rows.append({
            "quantity": f"{name}_pearson",
            "value": safe_corr(df[xcol], df[ycol], "pearson"),
        })
        rows.append({
            "quantity": f"{name}_spearman",
            "value": safe_corr(df[xcol], df[ycol], "spearman"),
        })

    rows.append({
        "quantity": "t_nuc_vs_t_tr_pearson",
        "value": safe_corr(both["t_nuc_cluster"], both["t_tr"], "pearson"),
    })
    rows.append({
        "quantity": "t_nuc_vs_t_tr_spearman",
        "value": safe_corr(both["t_nuc_cluster"], both["t_tr"], "spearman"),
    })

    if len(both) >= 2:
        delay = both["t_tr"] - both["t_nuc_cluster"]
        rows.append({
            "quantity": "transition_delay_mean_ttr_minus_tnuc",
            "value": float(delay.mean()),
        })
        rows.append({
            "quantity": "transition_delay_std_ttr_minus_tnuc",
            "value": float(delay.std(ddof=1)),
        })

    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate Fig. 5 seed amount/geometry and tnuc/ttr scatter.")
    parser.add_argument("--data-dir", default="data")
    parser.add_argument("--outdir", default="figures/main")
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    data_dir = Path(args.data_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_main_data(data_dir)
    both = df[df["t_nuc_cluster"].notna() & df["t_tr"].notna()].copy()

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.1), constrained_layout=True)

    ax = axes[0]
    for label, g in sorted(df.groupby("label"), key=lambda item: float(item[0])):
        ax.scatter(
            g["p_seed"],
            g["seed_compactness"],
            s=22,
            alpha=0.72,
            label=f"{float(label):.2f}",
        )
    ax.set_xlabel(r"Barrier-crossing seed fraction $p_{\rm seed}$")
    ax.set_ylabel(r"Seed compactness $C_{\rm comp}$")
    ax.set_title(r"(a) seed amount versus geometry")
    ax.legend(title="label", fontsize=7, title_fontsize=8, loc="best", ncol=2)

    ax = axes[1]
    ax.scatter(both["t_nuc_cluster"], both["t_tr"], s=24, alpha=0.72)
    if len(both) > 0:
        xmin = float(np.nanmin(both["t_nuc_cluster"]))
        xmax = float(np.nanmax(both["t_nuc_cluster"]))
        ymin = float(np.nanmin(both["t_tr"]))
        ymax = float(np.nanmax(both["t_tr"]))
        lo = min(xmin, ymin)
        hi = max(xmax, ymax)
        ax.plot([lo, hi], [lo, hi], linestyle="--", linewidth=1.0, label=r"$t_{\rm tr}=t_{\rm nuc}$")
        ax.legend(fontsize=8, loc="best")
    ax.set_xlabel(r"Cluster-nucleation time $t_{\rm nuc}$")
    ax.set_ylabel(r"Global transition time $t_{\rm tr}$")
    ax.set_title(r"(b) nucleation time versus transition time")

    for ext in ("png", "pdf"):
        path = outdir / f"fig5_seed_geometry_transition.{ext}"
        fig.savefig(path, dpi=args.dpi if ext == "png" else None, bbox_inches="tight")
        print(f"[saved] {path}")
    plt.close(fig)

    summary = make_summary(df)
    summary_path = outdir / "fig5_seed_geometry_transition_summary.csv"
    summary.to_csv(summary_path, index=False)
    print(f"[saved] {summary_path}")


if __name__ == "__main__":
    main()
