from __future__ import annotations

from pathlib import Path
import argparse
import math
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def find_candidates(data_dir: Path, seed: int) -> list[Path]:
    """Find strict main-condition candidates for a given seed.

    We intentionally require dt0p02 and ns50 so that the dt=0.01 robustness
    run is not accidentally used as the main data.
    """
    patterns = [
        f"data/main/*N64*dt0p02*tmax300*pre500*Df9em3*D00p02*ns50*seed{seed}*af0p02*b1*c1*cth20*mc3*labels_default.samples.csv",
        f"data/main/*N64*dt0p02*pre500*Df9em3*ns50*seed{seed}*cth20*.samples.csv",
        f"**/*N64*dt0p02*tmax300*pre500*Df9em3*D00p02*ns50*seed{seed}*af0p02*b1*c1*cth20*mc3*labels_default.samples.csv",
        f"**/*N64*dt0p02*pre500*Df9em3*ns50*seed{seed}*cth20*.samples.csv",
    ]

    matches: list[Path] = []
    for pat in patterns:
        matches.extend(data_dir.glob(pat))

    matches = sorted({p.resolve() for p in matches if p.exists() and p.suffix.lower() == ".csv"})
    return matches


def choose_main_file(data_dir: Path, seed: int) -> Path | None:
    matches = find_candidates(data_dir, seed)
    if not matches:
        return None

    # Priority:
    # 1. data/main
    # 2. non-run2 file, if multiple exist in data/main
    # 3. shortest filename
    def key(p: Path):
        parts = {part.lower() for part in p.parts}
        in_data_main = ("data" in parts and "main" in parts)
        has_run2 = "run2" in p.name
        return (
            0 if in_data_main else 1,
            1 if has_run2 else 0,
            len(p.name),
            p.name,
        )

    return sorted(matches, key=key)[0]


def show_candidate_help(data_dir: Path):
    print("\nCandidate CSV files found under data-dir:")
    candidates = sorted(data_dir.glob("**/*Df9em3*cth20*.samples.csv"))
    if not candidates:
        print("  none")
        return
    for p in candidates[:80]:
        rel = p.relative_to(data_dir) if p.is_relative_to(data_dir) else p
        print(" ", rel)


def safe_pearson(x: pd.Series, y: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce")
    y = pd.to_numeric(y, errors="coerce")
    if len(x) < 2 or x.nunique(dropna=True) < 2 or y.nunique(dropna=True) < 2:
        return float("nan")
    return float(x.corr(y, method="pearson"))


def safe_spearman_without_scipy(x: pd.Series, y: pd.Series) -> float:
    """Spearman correlation without scipy.

    pandas Series.corr(method='spearman') imports scipy in recent setups.
    This computes Spearman as Pearson correlation of average ranks.
    """
    x = pd.to_numeric(x, errors="coerce")
    y = pd.to_numeric(y, errors="coerce")
    if len(x) < 2 or x.nunique(dropna=True) < 2 or y.nunique(dropna=True) < 2:
        return float("nan")
    xr = x.rank(method="average")
    yr = y.rank(method="average")
    return float(xr.corr(yr, method="pearson"))


def main():
    parser = argparse.ArgumentParser(description="Create supplemental t_nuc vs t_tr scatter plot.")
    parser.add_argument("--data-dir", default=".", help="Repository/data root. Default: current directory.")
    parser.add_argument("--outdir", default="tnuc_ttr_supplement", help="Output directory.")
    parser.add_argument("--formats", default="png,pdf", help="Comma-separated output formats.")
    parser.add_argument("--dpi", type=int, default=300, help="PNG DPI.")
    parser.add_argument("--show-title", action="store_true", help="Add a plot title for debugging.")
    args = parser.parse_args()

    data_dir = Path(args.data_dir).resolve()
    outdir = Path(args.outdir)
    if not outdir.is_absolute():
        outdir = data_dir / outdir
    outdir.mkdir(parents=True, exist_ok=True)
    formats = [s.strip().lower() for s in args.formats.split(",") if s.strip()] or ["png", "pdf"]

    f5678 = choose_main_file(data_dir, 5678)
    f9876 = choose_main_file(data_dir, 9876)

    if f5678 is None or f9876 is None:
        print("Could not find both strict main CSV files.")
        print("Required condition: N64, dt0p02, pre500, Df9em3, ns50, cth20, seed5678/seed9876")
        print("seed5678:", f5678)
        print("seed9876:", f9876)
        show_candidate_help(data_dir)
        sys.exit(1)

    print("Using main files:")
    print(" ", f5678.relative_to(data_dir) if f5678.is_relative_to(data_dir) else f5678)
    print(" ", f9876.relative_to(data_dir) if f9876.is_relative_to(data_dir) else f9876)

    frames = []
    for p in [f5678, f9876]:
        df = pd.read_csv(p)
        df["source_file"] = p.name
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)

    required = ["label", "t_nuc_cluster", "t_tr"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    d = df.copy()
    d["t_nuc_cluster"] = pd.to_numeric(d["t_nuc_cluster"], errors="coerce")
    d["t_tr"] = pd.to_numeric(d["t_tr"], errors="coerce")

    observed = d.dropna(subset=["t_nuc_cluster", "t_tr"]).copy()
    if observed.empty:
        raise ValueError("No samples with both t_nuc_cluster and t_tr were found.")

    pearson = safe_pearson(observed["t_nuc_cluster"], observed["t_tr"])
    spearman = safe_spearman_without_scipy(observed["t_nuc_cluster"], observed["t_tr"])

    summary = pd.DataFrame({
        "n_total": [len(d)],
        "n_observed_both": [len(observed)],
        "n_nuc_observed": [int(d["t_nuc_cluster"].notna().sum())],
        "n_tr_observed": [int(d["t_tr"].notna().sum())],
        "pearson_corr": [pearson],
        "spearman_corr_no_scipy": [spearman],
        "mean_t_nuc": [observed["t_nuc_cluster"].mean()],
        "mean_t_tr": [observed["t_tr"].mean()],
        "mean_ttr_minus_tnuc": [(observed["t_tr"] - observed["t_nuc_cluster"]).mean()],
        "median_ttr_minus_tnuc": [(observed["t_tr"] - observed["t_nuc_cluster"]).median()],
    })
    summary.to_csv(outdir / "summary_tnuc_ttr_relation.csv", index=False)
    observed.to_csv(outdir / "samples_with_observed_tnuc_and_ttr.csv", index=False)

    fig, ax = plt.subplots(figsize=(5.8, 4.4))
    for label, g in sorted(observed.groupby("label"), key=lambda x: float(x[0])):
        ax.scatter(
            g["t_nuc_cluster"],
            g["t_tr"],
            s=30,
            alpha=0.75,
            label=fr"$T_i={float(label):.2f}$",
        )

    lo = float(min(observed["t_nuc_cluster"].min(), observed["t_tr"].min()))
    hi = float(max(observed["t_nuc_cluster"].max(), observed["t_tr"].max()))
    pad = 0.03 * (hi - lo) if hi > lo else 1.0
    lo -= pad
    hi += pad
    ax.plot([lo, hi], [lo, hi], linestyle="--", linewidth=1.0, label=r"$t_{\rm tr}=t_{\rm nuc}$")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel(r"Nucleation time $t_{\rm nuc}$")
    ax.set_ylabel(r"Transition time $t_{\rm tr}$")
    if args.show_title:
        ax.set_title(r"Relation between $t_{\rm nuc}$ and $t_{\rm tr}$")
    ax.legend(fontsize=7.5, loc="best")
    fig.tight_layout()

    for fmt in formats:
        out = outdir / f"figS7_tnuc_vs_ttr_scatter.{fmt}"
        fig.savefig(out, dpi=args.dpi if fmt in {"png", "jpg", "jpeg", "tif", "tiff"} else None, bbox_inches="tight")
        print("Saved:", out)
    plt.close(fig)

    print("Saved:", outdir / "summary_tnuc_ttr_relation.csv")
    print("Saved:", outdir / "samples_with_observed_tnuc_and_ttr.csv")
    print("\nSummary:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
