from __future__ import annotations

from pathlib import Path
import argparse
import math
import sys
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# TDGL Mpemba manuscript figures, v4 multipanel version
# ============================================================
# This script creates paper-ready multi-panel figures from the
# *.samples.csv files produced by tdgl_run_auto_v3.py.
#
# Main outputs:
#   fig2_nucleation_multipanel
#   fig3_seed_amount_multipanel
#   fig4_seed_quality_multipanel
#
# Supplementary outputs:
#   figS1_ordered_seed_multipanel
#   figS2_numerical_checks_multipanel
#   figS3_robustness_multipanel
#
# Notes:
#   - Survival curves use right-censored Kaplan-Meier treatment.
#   - No titles are added by default.
#   - PNG and PDF are exported by default.
# ============================================================


# -----------------------------
# File discovery
# -----------------------------

def find_matches(data_dir: Path, patterns: list[str]) -> list[Path]:
    matches: list[Path] = []
    for pat in patterns:
        matches.extend(data_dir.glob(pat))
        matches.extend(data_dir.glob(f"**/{pat}"))
    return sorted({p.resolve() for p in matches if p.exists() and p.suffix.lower() == ".csv"})


def choose_file(data_dir: Path, label: str, patterns: list[str]) -> Path | None:
    matches = find_matches(data_dir, patterns)
    if not matches:
        print(f"[not found] {label}")
        return None

    # Prefer run2 when available because the main analysis used run2 in this project.
    run2 = [p for p in matches if "run2" in p.name]
    chosen = sorted(run2 or matches, key=lambda p: (len(p.name), p.name))[0]
    print(f"[found] {label}: {chosen.name}")
    return chosen


def load_concat(paths: Iterable[Path | None]) -> pd.DataFrame:
    frames = []
    for p in paths:
        if p is not None and p.exists():
            df = pd.read_csv(p)
            df["source_file"] = p.name
            frames.append(df)
    if not frames:
        raise FileNotFoundError("No input files found for the requested dataset.")
    return pd.concat(frames, ignore_index=True)


# -----------------------------
# Statistics
# -----------------------------

def stderr(series: pd.Series) -> float:
    series = pd.to_numeric(series, errors="coerce").dropna()
    if len(series) <= 1:
        return np.nan
    return float(series.std(ddof=1) / math.sqrt(len(series)))


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    required = ["label", "t_nuc_cluster", "t_tr", "p_seed", "c_seed"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    rows = []
    for label, g in sorted(df.groupby("label"), key=lambda x: float(x[0])):
        tn = pd.to_numeric(g["t_nuc_cluster"], errors="coerce")
        tt = pd.to_numeric(g["t_tr"], errors="coerce")
        n = len(g)
        p_nuc = float(tn.notna().mean())
        p_tr = float(tt.notna().mean())

        def mean_col(name: str) -> float:
            return float(pd.to_numeric(g[name], errors="coerce").mean()) if name in g.columns else np.nan

        def se_col(name: str) -> float:
            return stderr(g[name]) if name in g.columns else np.nan

        rows.append({
            "label": float(label),
            "n": n,
            "p_seed_mean": mean_col("p_seed"),
            "p_seed_se": se_col("p_seed"),
            "c_seed_mean": mean_col("c_seed"),
            "c_seed_se": se_col("c_seed"),
            "seed_compactness_mean": mean_col("seed_compactness"),
            "seed_compactness_se": se_col("seed_compactness"),
            "seed_rg_mean": mean_col("seed_rg"),
            "seed_rg_se": se_col("seed_rg"),
            "p_ordered_seed_mean": mean_col("p_ordered_seed"),
            "p_ordered_seed_se": se_col("p_ordered_seed"),
            "c_ordered_seed_mean": mean_col("c_ordered_seed"),
            "c_ordered_seed_se": se_col("c_ordered_seed"),
            "p_nuc": p_nuc,
            "p_nuc_se": math.sqrt(p_nuc * (1.0 - p_nuc) / n) if n else np.nan,
            "p_tr": p_tr,
            "p_tr_se": math.sqrt(p_tr * (1.0 - p_tr) / n) if n else np.nan,
            "t_nuc_mean": float(tn.mean()) if tn.notna().any() else np.nan,
            "t_tr_mean": float(tt.mean()) if tt.notna().any() else np.nan,
        })
    return pd.DataFrame(rows)


# -----------------------------
# Plot helpers
# -----------------------------

def save_in_formats(fig: plt.Figure, basepath: Path, formats: list[str], dpi: int):
    for fmt in formats:
        out = basepath.with_suffix(f".{fmt}")
        fig.savefig(
            out,
            dpi=dpi if fmt.lower() in {"png", "jpg", "jpeg", "tif", "tiff"} else None,
            bbox_inches="tight",
        )
        print(f"[saved] {out}")


def set_label_ticks(ax: plt.Axes, labels: Iterable[float]):
    labels = list(labels)
    ax.set_xticks(labels)
    ax.set_xticklabels([f"{x:.2f}" for x in labels])


def add_panel_label(ax: plt.Axes, label: str):
    ax.text(
        0.02, 0.98, label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12,
        fontweight="bold",
    )


def plot_errorbar_panel(
    ax: plt.Axes,
    summary: pd.DataFrame,
    ycol: str,
    yerrcol: str,
    ylabel: str,
    panel_label: str,
    xlabel: str = "Initial-state label",
):
    if ycol not in summary.columns or summary[ycol].isna().all():
        ax.text(0.5, 0.5, f"Unavailable: {ycol}", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    yerr = summary[yerrcol] if yerrcol in summary.columns else None
    ax.errorbar(summary["label"], summary[ycol], yerr=yerr, marker="o", capsize=3)
    set_label_ticks(ax, summary["label"])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    add_panel_label(ax, panel_label)


def plot_multi_errorbar_panel(
    ax: plt.Axes,
    datasets: list[tuple[str, pd.DataFrame | None]],
    ycol: str,
    yerrcol: str,
    ylabel: str,
    panel_label: str,
    xlabel: str = "Initial-state label",
    legend_loc: str = "best",
):
    usable = [(name, s) for name, s in datasets if s is not None and ycol in s.columns and not s[ycol].isna().all()]
    if not usable:
        ax.text(0.5, 0.5, f"Unavailable: {ycol}", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    for name, summary in usable:
        yerr = summary[yerrcol] if yerrcol in summary.columns else None
        ax.errorbar(summary["label"], summary[ycol], yerr=yerr, marker="o", capsize=3, label=name)
    set_label_ticks(ax, usable[0][1]["label"])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(loc=legend_loc)
    add_panel_label(ax, panel_label)


def kaplan_meier_curve(event_times: np.ndarray, observed: np.ndarray) -> tuple[list[float], list[float]]:
    """Kaplan-Meier survival curve with right-censoring."""
    if len(event_times) == 0:
        return [0.0], [1.0]

    order = np.argsort(event_times)
    t = np.asarray(event_times)[order]
    e = np.asarray(observed, dtype=bool)[order]

    xs = [0.0]
    ys = [1.0]
    surv = 1.0
    n_risk = len(t)

    for current_t in np.unique(t):
        mask = (t == current_t)
        d_i = int(np.sum(e[mask]))
        c_i = int(np.sum(~e[mask]))

        xs.append(float(current_t))
        ys.append(float(surv))

        if d_i > 0 and n_risk > 0:
            surv *= (1.0 - d_i / n_risk)
            xs.append(float(current_t))
            ys.append(float(surv))

        n_risk -= (d_i + c_i)

    return xs, ys


def plot_survival_panel(
    ax: plt.Axes,
    df: pd.DataFrame,
    time_col: str,
    panel_label: str,
    tmax: float = 300.0,
    max_labels: int | None = None,
):
    if time_col not in df.columns:
        ax.text(0.5, 0.5, f"Unavailable: {time_col}", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    labels = sorted(df["label"].unique())
    if max_labels is not None and len(labels) > max_labels:
        # Keep low, middle-ish, high labels for readability.
        idxs = np.linspace(0, len(labels) - 1, max_labels).round().astype(int)
        labels = [labels[i] for i in idxs]

    for lab in labels:
        g = df[df["label"] == lab]
        raw = pd.to_numeric(g[time_col], errors="coerce")
        observed = raw.notna().to_numpy()
        event_times = raw.fillna(tmax).to_numpy(dtype=float)
        xs, ys = kaplan_meier_curve(event_times, observed)
        if xs[-1] < tmax:
            xs.append(float(tmax))
            ys.append(float(ys[-1]))
        ax.step(xs, ys, where="post", label=fr"$\lambda={float(lab):.2f}$")

    ax.set_xlabel(r"Time $t$")
    ax.set_ylabel(r"Survival probability $S(t)$")
    ax.set_xlim(left=0.0, right=tmax)
    ax.set_ylim(-0.02, 1.02)
    ax.legend(loc="best", fontsize=9)
    add_panel_label(ax, panel_label)


# -----------------------------
# Multi-panel figure makers
# -----------------------------

def make_fig2(df_main: pd.DataFrame, sum_main: pd.DataFrame, outdir: Path, formats: list[str], dpi: int, tmax: float, max_survival_labels: int | None):
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    plot_errorbar_panel(
        axes[0], sum_main, "p_nuc", "p_nuc_se",
        r"Nucleation probability $P_{\mathrm{nuc}}$", "(a)"
    )
    plot_survival_panel(axes[1], df_main, "t_nuc_cluster", "(b)", tmax=tmax, max_labels=max_survival_labels)
    fig.tight_layout()
    save_in_formats(fig, outdir / "fig2_nucleation_multipanel", formats, dpi)
    plt.close(fig)


def make_fig3(sum_main: pd.DataFrame, outdir: Path, formats: list[str], dpi: int):
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    plot_errorbar_panel(
        axes[0], sum_main, "p_seed_mean", "p_seed_se",
        r"Mean barrier-crossing seed fraction $p_{\mathrm{seed}}$", "(a)"
    )
    plot_errorbar_panel(
        axes[1], sum_main, "c_seed_mean", "c_seed_se",
        r"Mean largest barrier-crossing cluster size $c_{\mathrm{seed}}$", "(b)"
    )
    fig.tight_layout()
    save_in_formats(fig, outdir / "fig3_seed_amount_multipanel", formats, dpi)
    plt.close(fig)


def make_fig4(sum_main: pd.DataFrame, outdir: Path, formats: list[str], dpi: int):
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    plot_errorbar_panel(
        axes[0], sum_main, "seed_compactness_mean", "seed_compactness_se",
        "Mean seed compactness", "(a)"
    )
    plot_errorbar_panel(
        axes[1], sum_main, "seed_rg_mean", "seed_rg_se",
        r"Mean seed radius of gyration $R_g$", "(b)"
    )
    fig.tight_layout()
    save_in_formats(fig, outdir / "fig4_seed_quality_multipanel", formats, dpi)
    plt.close(fig)


def make_figS1(sum_main: pd.DataFrame, outdir: Path, formats: list[str], dpi: int):
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    plot_errorbar_panel(
        axes[0], sum_main, "p_ordered_seed_mean", "p_ordered_seed_se",
        r"Mean ordered-like seed fraction $p_{\mathrm{ordered}}$", "(a)"
    )
    plot_errorbar_panel(
        axes[1], sum_main, "c_ordered_seed_mean", "c_ordered_seed_se",
        r"Mean largest ordered-like cluster size $c_{\mathrm{ordered}}$", "(b)"
    )
    fig.tight_layout()
    save_in_formats(fig, outdir / "figS1_ordered_seed_multipanel", formats, dpi)
    plt.close(fig)


def make_figS2(sum_cth10, sum_main, sum_cth30, sum_dt001, outdir: Path, formats: list[str], dpi: int):
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    plot_multi_errorbar_panel(
        axes[0], [("Cth=10", sum_cth10), ("Cth=20", sum_main), ("Cth=30", sum_cth30)],
        "p_nuc", "p_nuc_se", r"Nucleation probability $P_{\mathrm{nuc}}$", "(a)", legend_loc="best"
    )
    plot_multi_errorbar_panel(
        axes[1], [("dt=0.02", sum_main), ("dt=0.01", sum_dt001)],
        "p_nuc", "p_nuc_se", r"Nucleation probability $P_{\mathrm{nuc}}$", "(b)", legend_loc="best"
    )
    fig.tight_layout()
    save_in_formats(fig, outdir / "figS2_numerical_checks_multipanel", formats, dpi)
    plt.close(fig)


def make_figS3(sum_main, sum_pre1000, sum_n128, outdir: Path, formats: list[str], dpi: int):
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    datasets = [("N64 pre500", sum_main), ("N64 pre1000", sum_pre1000), ("N128 pre500", sum_n128)]
    plot_multi_errorbar_panel(
        axes[0], datasets, "p_seed_mean", "p_seed_se",
        r"Mean barrier-crossing seed fraction $p_{\mathrm{seed}}$", "(a)", legend_loc="best"
    )
    plot_multi_errorbar_panel(
        axes[1], datasets, "seed_compactness_mean", "seed_compactness_se",
        "Mean seed compactness", "(b)", legend_loc="best"
    )
    fig.tight_layout()
    save_in_formats(fig, outdir / "figS3_robustness_multipanel", formats, dpi)
    plt.close(fig)


# -----------------------------
# README
# -----------------------------

def write_readme(outdir: Path, formats: list[str]):
    text = f"""# TDGL manuscript figures v4 multipanel

Generated by `tdgl_manuscript_figures_v4_multipanel.py`.

## Output formats

{', '.join(formats)}

## Main figures

- `fig2_nucleation_multipanel`: nucleation probability and right-censored survival curve
- `fig3_seed_amount_multipanel`: barrier-crossing seed fraction and largest cluster size
- `fig4_seed_quality_multipanel`: seed compactness and seed radius of gyration

## Supplementary figures

- `figS1_ordered_seed_multipanel`: ordered-like seed statistics
- `figS2_numerical_checks_multipanel`: cluster-threshold and time-step checks
- `figS3_robustness_multipanel`: pre-equilibration and system-size robustness

## Notes

- Survival curves use a Kaplan-Meier right-censoring treatment.
- The D_f=8e-3 compactness comparison is intentionally not included by default because D_f is a post-quench parameter and should not affect initial seed statistics.
- No titles are added; panel labels are included inside each subplot.
"""
    (outdir / "README_generated_figures_v4_multipanel.md").write_text(text, encoding="utf-8")


# -----------------------------
# Main
# -----------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate multi-panel manuscript figures from TDGL CSV outputs.")
    parser.add_argument("--data-dir", default=".", help="Directory containing *.samples.csv files. Default: current directory.")
    parser.add_argument("--outdir", default="manuscript_figures_v4", help="Output directory.")
    parser.add_argument("--formats", default="png,pdf", help="Comma-separated output formats. Example: png,pdf")
    parser.add_argument("--dpi", type=int, default=300, help="Raster DPI for PNG output.")
    parser.add_argument("--tmax", type=float, default=300.0, help="Right-censoring time for survival curves.")
    parser.add_argument("--max-survival-labels", type=int, default=0, help="Limit survival curves to this many labels. 0 means all labels.")
    args = parser.parse_args()

    data_dir = Path(args.data_dir).resolve()
    outdir = Path(args.outdir)
    if not outdir.is_absolute():
        outdir = data_dir / outdir
    outdir.mkdir(parents=True, exist_ok=True)

    formats = [x.strip().lower() for x in args.formats.split(",") if x.strip()]
    if not formats:
        formats = ["png", "pdf"]

    print(f"[data-dir] {data_dir}")
    print(f"[outdir]   {outdir}")

    # Main condition: Df=9e-3, N=64, pre500, cth20, seed 5678 + 9876.
    main_5678 = choose_file(data_dir, "main Df9e-3 seed5678 cth20", [
        "*N64*pre500*Df9em3*seed5678*cth20*_run2.samples.csv",
        "*N64*pre500*Df9em3*seed5678*cth20*.samples.csv",
    ])
    main_9876 = choose_file(data_dir, "main Df9e-3 seed9876 cth20", [
        "*N64*pre500*Df9em3*seed9876*cth20*_run2.samples.csv",
        "*N64*pre500*Df9em3*seed9876*cth20*.samples.csv",
    ])

    if main_5678 is None and main_9876 is None:
        print("\nNo main D_f=9e-3 CSV files found.")
        print("Run with, for example:")
        print('  python tdgl_manuscript_figures_v4_multipanel.py --data-dir "C:\\path\\to\\csv_folder"')
        print("\nVisible CSV files:")
        visible = sorted(data_dir.glob("*.csv"))
        if not visible:
            print("  none")
        else:
            for p in visible[:50]:
                print(" ", p.name)
        sys.exit(1)

    df_main = load_concat([main_5678, main_9876])
    sum_main = summarize(df_main)
    sum_main.to_csv(outdir / "summary_main_df9.csv", index=False)

    # Optional supplementary conditions.
    cth10 = choose_file(data_dir, "supplement Cth=10", ["*N64*pre500*Df9em3*cth10*.samples.csv"])
    cth30 = choose_file(data_dir, "supplement Cth=30", ["*N64*pre500*Df9em3*cth30*.samples.csv"])
    dt001 = choose_file(data_dir, "supplement dt=0.01", ["*N64*dt0p01*pre500*Df9em3*cth20*.samples.csv"])
    n128 = choose_file(data_dir, "supplement N=128", ["*N128*pre500*Df9em3*cth40*.samples.csv"])
    pre1000 = choose_file(data_dir, "supplement preeq=1000", ["*N64*pre1000*Df9em3*cth20*.samples.csv"])

    sum_cth10 = summarize(pd.read_csv(cth10)) if cth10 else None
    sum_cth30 = summarize(pd.read_csv(cth30)) if cth30 else None
    sum_dt001 = summarize(pd.read_csv(dt001)) if dt001 else None
    sum_n128 = summarize(pd.read_csv(n128)) if n128 else None
    sum_pre1000 = summarize(pd.read_csv(pre1000)) if pre1000 else None

    # Main figures.
    max_survival_labels = None if args.max_survival_labels <= 0 else args.max_survival_labels
    make_fig2(df_main, sum_main, outdir, formats, args.dpi, args.tmax, max_survival_labels)
    make_fig3(sum_main, outdir, formats, args.dpi)
    make_fig4(sum_main, outdir, formats, args.dpi)

    # Supplementary figures.
    make_figS1(sum_main, outdir, formats, args.dpi)
    if sum_cth10 is not None and sum_cth30 is not None and sum_dt001 is not None:
        make_figS2(sum_cth10, sum_main, sum_cth30, sum_dt001, outdir, formats, args.dpi)
    else:
        print("[skip] figS2_numerical_checks_multipanel: cth10/cth30/dt0.01 data not complete")

    if sum_pre1000 is not None and sum_n128 is not None:
        make_figS3(sum_main, sum_pre1000, sum_n128, outdir, formats, args.dpi)
    else:
        print("[skip] figS3_robustness_multipanel: pre1000/N128 data not complete")

    write_readme(outdir, formats)
    print(f"\nDone. Multi-panel figures written to: {outdir}")


if __name__ == "__main__":
    main()
