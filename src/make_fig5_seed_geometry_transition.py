#!/usr/bin/env python
r"""
Generate improved Fig. 5 for the TDGL Mpemba seed-geometry manuscript.

Fig. 5:
(a) sample-level relation between seed amount and compactness,
    with large label-mean markers connected by a trend line;
(b) direct relation between cluster-nucleation time and global transition time,
    annotated with correlation and mean growth delay.

Usage from the repository root:

  python .\src\make_fig5_seed_geometry_transition.py `
    --data-dir .\data `
    --outdir .\figures\main

Outputs:
  figures/main/fig5_seed_geometry_transition.pdf
  figures/main/fig5_seed_geometry_transition.png
  figures/main/fig5_seed_geometry_transition_summary.csv
  figures/main/fig5_seed_geometry_transition_numbers.tex
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


def fmt_float(x: float, ndigits: int = 3) -> str:
    if not np.isfinite(x):
        return "nan"
    return f"{x:.{ndigits}f}"


def make_summary(df: pd.DataFrame) -> pd.DataFrame:
    both = df[df["t_nuc_cluster"].notna() & df["t_tr"].notna()].copy()
    rows = []

    rows.append({"quantity": "sample_count_total", "value": len(df)})
    rows.append({"quantity": "sample_count_nucleated", "value": int(df["t_nuc_cluster"].notna().sum())})
    rows.append({"quantity": "sample_count_transitioned", "value": int(df["t_tr"].notna().sum())})
    rows.append({"quantity": "sample_count_tnuc_and_ttr_observed", "value": len(both)})

    for xcol, ycol, name in [
        ("p_seed", "seed_compactness", "p_seed_vs_seed_compactness"),
        ("c_seed", "seed_compactness", "c_seed_vs_seed_compactness"),
        ("seed_rg", "seed_compactness", "seed_rg_vs_seed_compactness"),
    ]:
        rows.append({"quantity": f"{name}_pearson", "value": safe_corr(df[xcol], df[ycol], "pearson")})
        rows.append({"quantity": f"{name}_spearman", "value": safe_corr(df[xcol], df[ycol], "spearman")})

    rows.append({"quantity": "t_nuc_vs_t_tr_pearson", "value": safe_corr(both["t_nuc_cluster"], both["t_tr"], "pearson")})
    rows.append({"quantity": "t_nuc_vs_t_tr_spearman", "value": safe_corr(both["t_nuc_cluster"], both["t_tr"], "spearman")})

    if len(both) >= 2:
        delay = both["t_tr"] - both["t_nuc_cluster"]
        rows.append({"quantity": "transition_delay_mean_ttr_minus_tnuc", "value": float(delay.mean())})
        rows.append({"quantity": "transition_delay_std_ttr_minus_tnuc", "value": float(delay.std(ddof=1))})
        rows.append({"quantity": "transition_delay_median_ttr_minus_tnuc", "value": float(delay.median())})

    label_means = (
        df.groupby("label", as_index=False)
        .agg(
            p_seed_mean=("p_seed", "mean"),
            p_seed_se=("p_seed", lambda x: x.std(ddof=1) / math.sqrt(len(x))),
            seed_compactness_mean=("seed_compactness", "mean"),
            seed_compactness_se=("seed_compactness", lambda x: x.std(ddof=1) / math.sqrt(len(x))),
            seed_rg_mean=("seed_rg", "mean"),
        )
        .sort_values("label")
    )
    for _, r in label_means.iterrows():
        rows.append({"quantity": f"label_{r['label']:.2f}_p_seed_mean", "value": float(r["p_seed_mean"])})
        rows.append({"quantity": f"label_{r['label']:.2f}_seed_compactness_mean", "value": float(r["seed_compactness_mean"])})

    return pd.DataFrame(rows)


def write_tex_macros(summary: pd.DataFrame, outpath: Path) -> None:
    vals = {row["quantity"]: float(row["value"]) for _, row in summary.iterrows() if pd.notna(row["value"])}

    pearson = vals.get("t_nuc_vs_t_tr_pearson", np.nan)
    spearman = vals.get("t_nuc_vs_t_tr_spearman", np.nan)
    delay_mean = vals.get("transition_delay_mean_ttr_minus_tnuc", np.nan)
    delay_std = vals.get("transition_delay_std_ttr_minus_tnuc", np.nan)
    n_both = int(vals.get("sample_count_tnuc_and_ttr_observed", 0))

    p_comp_spear = vals.get("p_seed_vs_seed_compactness_spearman", np.nan)
    p_comp_pear = vals.get("p_seed_vs_seed_compactness_pearson", np.nan)

    content = rf"""% Auto-generated by make_fig5_seed_geometry_transition.py
\newcommand{{\FigFiveTNucTTrPearson}}{{{fmt_float(pearson, 3)}}}
\newcommand{{\FigFiveTNucTTrSpearman}}{{{fmt_float(spearman, 3)}}}
\newcommand{{\FigFiveGrowDelayMean}}{{{fmt_float(delay_mean, 2)}}}
\newcommand{{\FigFiveGrowDelayStd}}{{{fmt_float(delay_std, 2)}}}
\newcommand{{\FigFiveNBoth}}{{{n_both}}}
\newcommand{{\FigFivePSeedCompactnessPearson}}{{{fmt_float(p_comp_pear, 3)}}}
\newcommand{{\FigFivePSeedCompactnessSpearman}}{{{fmt_float(p_comp_spear, 3)}}}
"""
    outpath.write_text(content, encoding="utf-8")
    print(f"[saved] {outpath}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate improved Fig. 5 seed amount/geometry and tnuc/ttr scatter.")
    parser.add_argument("--data-dir", default="data")
    parser.add_argument("--outdir", default="figures/main")
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    data_dir = Path(args.data_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_main_data(data_dir)
    both = df[df["t_nuc_cluster"].notna() & df["t_tr"].notna()].copy()

    label_means = (
        df.groupby("label", as_index=False)
        .agg(
            p_seed_mean=("p_seed", "mean"),
            p_seed_se=("p_seed", lambda x: x.std(ddof=1) / math.sqrt(len(x))),
            seed_compactness_mean=("seed_compactness", "mean"),
            seed_compactness_se=("seed_compactness", lambda x: x.std(ddof=1) / math.sqrt(len(x))),
        )
        .sort_values("label")
    )

    pearson = safe_corr(both["t_nuc_cluster"], both["t_tr"], "pearson")
    spearman = safe_corr(both["t_nuc_cluster"], both["t_tr"], "spearman")
    delay = both["t_tr"] - both["t_nuc_cluster"]
    delay_mean = float(delay.mean()) if len(delay) else np.nan
    delay_std = float(delay.std(ddof=1)) if len(delay) > 1 else np.nan

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.2), constrained_layout=True)

    ax = axes[0]
    for label, g in sorted(df.groupby("label"), key=lambda item: float(item[0])):
        ax.scatter(
            g["p_seed"],
            g["seed_compactness"],
            s=18,
            alpha=0.34,
            label=f"{float(label):.2f}",
            linewidths=0,
        )

    ax.errorbar(
        label_means["p_seed_mean"],
        label_means["seed_compactness_mean"],
        xerr=label_means["p_seed_se"],
        yerr=label_means["seed_compactness_se"],
        marker="o",
        markersize=7.5,
        linewidth=1.6,
        capsize=3,
        markeredgewidth=1.0,
        markeredgecolor="black",
        label="label mean",
        zorder=10,
    )

    # Direction arrow from lowest-label mean to highest-label mean.
    if len(label_means) >= 2:
        start = label_means.iloc[0]
        end = label_means.iloc[-1]
        ax.annotate(
            "",
            xy=(end["p_seed_mean"], end["seed_compactness_mean"]),
            xytext=(start["p_seed_mean"], start["seed_compactness_mean"]),
            arrowprops=dict(arrowstyle="->", lw=1.4),
            zorder=11,
        )
        ax.text(
            end["p_seed_mean"],
            end["seed_compactness_mean"],
            " larger label",
            fontsize=8,
            va="center",
        )

    ax.set_xlabel(r"Barrier-crossing seed fraction $p_{\rm seed}$")
    ax.set_ylabel(r"Seed compactness $C_{\rm comp}$")
    ax.set_title(r"(a) seed amount versus geometry")
    ax.legend(title="label", fontsize=6.5, title_fontsize=7.5, loc="best", ncol=2)

    ax = axes[1]
    ax.scatter(both["t_nuc_cluster"], both["t_tr"], s=22, alpha=0.56)
    if len(both) > 0:
        xmin = float(np.nanmin(both["t_nuc_cluster"]))
        xmax = float(np.nanmax(both["t_nuc_cluster"]))
        ymin = float(np.nanmin(both["t_tr"]))
        ymax = float(np.nanmax(both["t_tr"]))
        lo = min(xmin, ymin)
        hi = max(xmax, ymax)
        ax.plot([lo, hi], [lo, hi], linestyle="--", linewidth=1.0, label=r"$t_{\rm tr}=t_{\rm nuc}$")
        ax.legend(fontsize=8, loc="best")

    annotation = (
        rf"$r={pearson:.3f}$, $\rho={spearman:.3f}$" + "\n" +
        rf"$\langle\Delta t_{{\rm grow}}\rangle={delay_mean:.2f}\pm{delay_std:.2f}$"
    )
    ax.text(
        0.04,
        0.96,
        annotation,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8.5,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.75, linewidth=0.5),
    )

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

    write_tex_macros(summary, outdir / "fig5_seed_geometry_transition_numbers.tex")


if __name__ == "__main__":
    main()
