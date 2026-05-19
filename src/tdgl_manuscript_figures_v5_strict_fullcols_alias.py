from __future__ import annotations

from pathlib import Path
import argparse
import math
import sys
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


MAIN_REQUIRED_COLUMNS = [
    "label", "t_nuc_cluster", "t_tr",
    "p_seed", "c_seed",
    "seed_compactness", "seed_rg",
    "p_ordered_seed", "c_ordered_seed",
]

BASE_REQUIRED_COLUMNS = ["label", "t_nuc_cluster", "t_tr", "p_seed", "c_seed"]


def has_columns(path: Path, required: list[str]) -> bool:
    try:
        cols = pd.read_csv(path, nrows=0).columns.tolist()
    except Exception:
        return False
    return all(c in cols for c in required)


def list_candidates(data_dir: Path, pattern: str) -> list[Path]:
    return sorted({p.resolve() for p in data_dir.rglob(pattern) if p.exists() and p.suffix.lower() == ".csv"})


def choose_main_file(data_dir: Path, seed: int) -> Path:
    """Choose main phi6 Df=9e-3, dt=0.02, ns=50, cth20 file with full seed-quality columns.

    The selection is token-based rather than glob-pattern-based to avoid accidentally
    picking dt=0.01 robustness files. It refuses to continue unless a file containing
    seed_compactness and seed_rg is found.
    """
    include_tokens = [
        "N64", "dt0p02", "tmax300", "pre500", "Df9em3", "D00p02",
        "ns50", f"seed{seed}", "cth20", ".samples.csv",
    ]
    exclude_tokens = ["dt0p01", "cth10", "cth30", "pre1000", "N128"]

    candidates = []

    # First, allow short alias files. This avoids Windows MAX_PATH / OneDrive issues
    # caused by very long original filenames. Recommended aliases:
    #   data/main/phi6_main_Df9em3_seed5678_fullcols.samples.csv
    #   data/main/phi6_main_Df9em3_seed9876_fullcols.samples.csv
    alias_candidates = [
        data_dir / "data" / "main" / f"phi6_main_Df9em3_seed{seed}_fullcols.samples.csv",
        data_dir / "data" / "main" / f"phi6_main_seed{seed}_fullcols.samples.csv",
    ]
    for p in alias_candidates:
        if p.exists():
            candidates.append(p.resolve())

    # Then search the original long filenames.
    for p in data_dir.rglob("*.samples.csv"):
        name = p.name
        if all(tok in name for tok in include_tokens) and not any(tok in name for tok in exclude_tokens):
            candidates.append(p.resolve())

    candidates = sorted(set(candidates))

    if not candidates:
        print(f"[ERROR] No strict main candidate found for seed={seed}")
        print("Required filename tokens:", include_tokens)
        nearby = []
        for p in data_dir.rglob("*.samples.csv"):
            name = p.name
            if f"seed{seed}" in name and "Df9em3" in name and "cth20" in name:
                nearby.append(p)
        if nearby:
            print("Nearby candidates:")
            for p in sorted(nearby):
                print("  ", p.relative_to(data_dir))
        raise FileNotFoundError(f"Strict main file for seed={seed} was not found.")

    print(f"[candidates] strict main seed{seed}:")
    for p in candidates:
        full = has_columns(p, MAIN_REQUIRED_COLUMNS)
        base = has_columns(p, BASE_REQUIRED_COLUMNS)
        try:
            rel = p.relative_to(data_dir)
        except ValueError:
            rel = p
        print(f"  {rel}  base_columns={base}  full_seed_quality_columns={full}")

    full_candidates = [p for p in candidates if has_columns(p, MAIN_REQUIRED_COLUMNS)]
    if not full_candidates:
        print("\n[ERROR] Strict main candidates were found, but none contain seed-quality columns.")
        print("Needed columns:", MAIN_REQUIRED_COLUMNS)
        print("Use the *_run2.samples.csv files produced after adding seed_compactness/seed_rg,")
        print("or copy those run2 files into data/main and rerun.")
        raise ValueError(f"No full-column main file found for seed={seed}.")

    def score(p: Path) -> tuple[int, int, str]:
        # Higher is better: prefer run2, then files under data/robustness or data/main are both OK.
        run2 = 1 if "run2" in p.name else 0
        # Prefer robustness/run2 over old data/main if both exist, because old main files lack full columns.
        robustness = 1 if "data" in p.parts and "robustness" in p.parts else 0
        try:
            rel = str(p.relative_to(data_dir))
        except ValueError:
            rel = str(p)
        return (run2, robustness, rel)

    chosen = sorted(full_candidates, key=score, reverse=True)[0]
    print(f"[found] STRICT-FULL main seed{seed}: {chosen.relative_to(data_dir)}")

    if "dt0p01" in chosen.name:
        raise RuntimeError("Internal error: dt0p01 was selected as main data.")
    if not has_columns(chosen, MAIN_REQUIRED_COLUMNS):
        raise RuntimeError("Internal error: selected main file lacks seed-quality columns.")
    return chosen

def choose_optional_file(data_dir: Path, label: str, patterns: list[str], required_columns: list[str] | None = None) -> Path | None:
    candidates: list[Path] = []
    for pat in patterns:
        candidates.extend(list_candidates(data_dir, pat))
    candidates = sorted(set(candidates))
    if required_columns:
        candidates_with_cols = [p for p in candidates if has_columns(p, required_columns)]
        if candidates_with_cols:
            candidates = candidates_with_cols
    if not candidates:
        print(f"[not found] {label}")
        return None

    def score(p: Path) -> tuple[int, int, str]:
        run2 = 1 if "run2" in p.name else 0
        try:
            rel = str(p.relative_to(data_dir))
        except ValueError:
            rel = str(p)
        # Prefer files under data/robustness for checks, but not mandatory.
        robust = 1 if "data" in p.parts and "robustness" in p.parts else 0
        return (run2, robust, rel)

    chosen = sorted(candidates, key=score, reverse=True)[0]
    print(f"[found] {label}: {chosen.relative_to(data_dir)}")
    return chosen


def load_concat(paths: Iterable[Path]) -> pd.DataFrame:
    frames = []
    for p in paths:
        df = pd.read_csv(p)
        df["source_file"] = p.name
        frames.append(df)
    if not frames:
        raise FileNotFoundError("No input files supplied.")
    return pd.concat(frames, ignore_index=True)


def stderr(series: pd.Series) -> float:
    x = pd.to_numeric(series, errors="coerce").dropna()
    if len(x) <= 1:
        return np.nan
    return float(x.std(ddof=1) / math.sqrt(len(x)))


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    missing = [c for c in BASE_REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    rows = []
    for label, g in sorted(df.groupby("label"), key=lambda kv: float(kv[0])):
        tn = pd.to_numeric(g["t_nuc_cluster"], errors="coerce")
        tt = pd.to_numeric(g["t_tr"], errors="coerce")

        def mean_col(name: str) -> float:
            return float(pd.to_numeric(g[name], errors="coerce").mean()) if name in g.columns else np.nan

        def se_col(name: str) -> float:
            return stderr(g[name]) if name in g.columns else np.nan

        n = len(g)
        p_nuc = float(tn.notna().mean())
        p_tr = float(tt.notna().mean())
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
            "p_nuc_se": math.sqrt(p_nuc * (1 - p_nuc) / n) if n else np.nan,
            "p_tr": p_tr,
            "p_tr_se": math.sqrt(p_tr * (1 - p_tr) / n) if n else np.nan,
        })
    return pd.DataFrame(rows)


def choose_survival_labels(labels: list[float], max_labels: int | None) -> list[float]:
    labels = sorted(labels)
    if max_labels is None or max_labels <= 0 or max_labels >= len(labels):
        return labels
    if max_labels == 1:
        return [labels[len(labels)//2]]
    if max_labels == 2:
        return [labels[0], labels[-1]]
    # include low, middle, high as representative labels
    idx = np.linspace(0, len(labels)-1, max_labels).round().astype(int)
    return [labels[i] for i in sorted(set(idx))]


def set_label_ticks(ax, labels):
    labels = list(labels)
    ax.set_xticks(labels)
    ax.set_xticklabels([f"{x:.2f}" for x in labels], rotation=35, ha="right")
    ax.tick_params(axis="x", labelsize=8)


def save_figure(fig, out_base: Path, dpi: int):
    for ext in ["png", "pdf"]:
        path = out_base.with_suffix(f".{ext}")
        fig.savefig(path, dpi=dpi if ext == "png" else None, bbox_inches="tight")
        print(f"[saved] {path}")


def errorbar_panel(ax, summary: pd.DataFrame, ycol: str, yerrcol: str, ylabel: str, panel_label: str):
    ax.errorbar(summary["label"], summary[ycol], yerr=summary[yerrcol], marker="o", capsize=3)
    set_label_ticks(ax, summary["label"])
    ax.set_xlabel("Initial-state label")
    ax.set_ylabel(ylabel)
    ax.text(0.03, 0.95, panel_label, transform=ax.transAxes, va="top", ha="left", fontweight="bold")


def kaplan_meier_curve(times: np.ndarray, observed: np.ndarray) -> tuple[list[float], list[float]]:
    order = np.argsort(times)
    t = np.asarray(times)[order]
    e = np.asarray(observed, dtype=bool)[order]
    xs = [0.0]
    ys = [1.0]
    surv = 1.0
    n_risk = len(t)
    for current_t in np.unique(t):
        mask = t == current_t
        d_i = int(np.sum(e[mask]))
        c_i = int(np.sum(~e[mask]))
        xs.append(float(current_t))
        ys.append(float(surv))
        if d_i > 0 and n_risk > 0:
            surv *= 1.0 - d_i / n_risk
            xs.append(float(current_t))
            ys.append(float(surv))
        n_risk -= (d_i + c_i)
    return xs, ys


def plot_fig2(df_main: pd.DataFrame, sum_main: pd.DataFrame, outdir: Path, dpi: int, max_survival_labels: int | None, tmax: float):
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.7), constrained_layout=True)
    ax = axes[0]
    errorbar_panel(ax, sum_main, "p_nuc", "p_nuc_se", r"Nucleation probability $P_{\mathrm{nuc}}$", "(a)")

    ax = axes[1]
    labels = choose_survival_labels(sorted(df_main["label"].unique()), max_survival_labels)
    for label in labels:
        g = df_main[df_main["label"] == label]
        raw = pd.to_numeric(g["t_nuc_cluster"], errors="coerce")
        observed = raw.notna().to_numpy()
        times = raw.fillna(tmax).to_numpy(dtype=float)
        xs, ys = kaplan_meier_curve(times, observed)
        if xs[-1] < tmax:
            xs.append(tmax)
            ys.append(ys[-1])
        ax.step(xs, ys, where="post", label=fr"$T_i={float(label):.2f}$")
    ax.set_xlim(0, tmax)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(r"Time $t$")
    ax.set_ylabel(r"Survival probability $S(t)$")
    ax.text(0.03, 0.95, "(b)", transform=ax.transAxes, va="top", ha="left", fontweight="bold")
    ax.legend(fontsize=8, loc="best")
    save_figure(fig, outdir / "fig2_nucleation_multipanel", dpi)
    plt.close(fig)


def plot_fig3(sum_main: pd.DataFrame, outdir: Path, dpi: int):
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.7), constrained_layout=True)
    errorbar_panel(axes[0], sum_main, "p_seed_mean", "p_seed_se", r"Mean barrier-crossing seed fraction $p_{\mathrm{seed}}$", "(a)")
    errorbar_panel(axes[1], sum_main, "c_seed_mean", "c_seed_se", r"Mean largest barrier-crossing cluster size $c_{\mathrm{seed}}$", "(b)")
    save_figure(fig, outdir / "fig3_seed_amount_multipanel", dpi)
    plt.close(fig)


def plot_fig4(sum_main: pd.DataFrame, outdir: Path, dpi: int):
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.7), constrained_layout=True)
    errorbar_panel(axes[0], sum_main, "seed_compactness_mean", "seed_compactness_se", "Mean seed compactness", "(a)")
    errorbar_panel(axes[1], sum_main, "seed_rg_mean", "seed_rg_se", r"Mean seed radius of gyration $R_g$", "(b)")
    save_figure(fig, outdir / "fig4_seed_quality_multipanel", dpi)
    plt.close(fig)


def multi_errorbar(ax, datasets, ycol: str, yerrcol: str, ylabel: str, panel_label: str):
    first_summary = None
    for label, summary in datasets:
        if summary is None or ycol not in summary.columns or summary[ycol].isna().all():
            continue
        first_summary = summary if first_summary is None else first_summary
        ax.errorbar(summary["label"], summary[ycol], yerr=summary[yerrcol], marker="o", capsize=3, label=label)
    if first_summary is not None:
        set_label_ticks(ax, first_summary["label"])
    ax.set_xlabel("Initial-state label")
    ax.set_ylabel(ylabel)
    ax.text(0.03, 0.95, panel_label, transform=ax.transAxes, va="top", ha="left", fontweight="bold")
    ax.legend(fontsize=8, loc="best")


def plot_supplementals(sum_main, sum_cth10, sum_cth30, sum_dt001, sum_n128, sum_pre1000, outdir: Path, dpi: int):
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.7), constrained_layout=True)
    errorbar_panel(axes[0], sum_main, "p_ordered_seed_mean", "p_ordered_seed_se", r"Mean ordered-like seed fraction $p_{\mathrm{ordered}}$", "(a)")
    errorbar_panel(axes[1], sum_main, "c_ordered_seed_mean", "c_ordered_seed_se", r"Mean largest ordered-like cluster size $c_{\mathrm{ordered}}$", "(b)")
    save_figure(fig, outdir / "figS1_ordered_seed_multipanel", dpi)
    plt.close(fig)

    if sum_cth10 is not None and sum_cth30 is not None and sum_dt001 is not None:
        fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.7), constrained_layout=True)
        multi_errorbar(axes[0], [("Cth=10", sum_cth10), ("Cth=20", sum_main), ("Cth=30", sum_cth30)], "p_nuc", "p_nuc_se", r"Nucleation probability $P_{\mathrm{nuc}}$", "(a)")
        multi_errorbar(axes[1], [("dt=0.02", sum_main), ("dt=0.01", sum_dt001)], "p_nuc", "p_nuc_se", r"Nucleation probability $P_{\mathrm{nuc}}$", "(b)")
        save_figure(fig, outdir / "figS2_numerical_checks_multipanel", dpi)
        plt.close(fig)

    if sum_n128 is not None and sum_pre1000 is not None:
        fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.7), constrained_layout=True)
        multi_errorbar(axes[0], [("N64 pre500", sum_main), ("N64 pre1000", sum_pre1000), ("N128 pre500", sum_n128)], "p_seed_mean", "p_seed_se", r"Mean barrier-crossing seed fraction $p_{\mathrm{seed}}$", "(a)")
        multi_errorbar(axes[1], [("N64 pre500", sum_main), ("N64 pre1000", sum_pre1000), ("N128 pre500", sum_n128)], "seed_compactness_mean", "seed_compactness_se", "Mean seed compactness", "(b)")
        save_figure(fig, outdir / "figS3_robustness_multipanel", dpi)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Strict-main multi-panel manuscript figures for TDGL Mpemba project.")
    parser.add_argument("--data-dir", default=".")
    parser.add_argument("--outdir", default="manuscript_figures_v5")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--tmax", type=float, default=300.0)
    parser.add_argument("--max-survival-labels", type=int, default=0, help="0 means all labels; 3 means low/mid/high representatives.")
    args = parser.parse_args()

    data_dir = Path(args.data_dir).resolve()
    outdir = Path(args.outdir)
    if not outdir.is_absolute():
        outdir = data_dir / outdir
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[data-dir] {data_dir}")
    print(f"[outdir]   {outdir}")

    main_5678 = choose_main_file(data_dir, 5678)
    main_9876 = choose_main_file(data_dir, 9876)
    df_main = load_concat([main_5678, main_9876])
    sum_main = summarize(df_main)
    sum_main.to_csv(outdir / "summary_main_df9.csv", index=False)

    # Optional checks. Make patterns strict enough not to pick main files accidentally.
    cth10 = choose_optional_file(data_dir, "supplement Cth=10", ["*N64*dt0p02*tmax300*pre500*Df9em3*ns30*seed5678*cth10*.samples.csv", "*N64*dt0p02*tmax300*pre500*Df9em3*ns20*seed5678*cth10*.samples.csv"], BASE_REQUIRED_COLUMNS)
    cth30 = choose_optional_file(data_dir, "supplement Cth=30", ["*N64*dt0p02*tmax300*pre500*Df9em3*ns30*seed5678*cth30*.samples.csv"], BASE_REQUIRED_COLUMNS)
    dt001 = choose_optional_file(data_dir, "supplement dt=0.01", ["*N64*dt0p01*tmax300*pre500*Df9em3*ns30*seed5678*cth20*.samples.csv"], BASE_REQUIRED_COLUMNS)
    n128 = choose_optional_file(data_dir, "supplement N=128", ["*N128*dt0p02*tmax300*pre500*Df9em3*ns20*seed5678*cth40*.samples.csv"], BASE_REQUIRED_COLUMNS)
    pre1000 = choose_optional_file(data_dir, "supplement preeq=1000", ["*N64*dt0p02*tmax300*pre1000*Df9em3*ns30*seed5678*cth20*.samples.csv"], BASE_REQUIRED_COLUMNS)

    sum_cth10 = summarize(pd.read_csv(cth10)) if cth10 else None
    sum_cth30 = summarize(pd.read_csv(cth30)) if cth30 else None
    sum_dt001 = summarize(pd.read_csv(dt001)) if dt001 else None
    sum_n128 = summarize(pd.read_csv(n128)) if n128 else None
    sum_pre1000 = summarize(pd.read_csv(pre1000)) if pre1000 else None

    max_survival = None if args.max_survival_labels <= 0 else args.max_survival_labels
    plot_fig2(df_main, sum_main, outdir, args.dpi, max_survival, args.tmax)
    plot_fig3(sum_main, outdir, args.dpi)
    plot_fig4(sum_main, outdir, args.dpi)
    plot_supplementals(sum_main, sum_cth10, sum_cth30, sum_dt001, sum_n128, sum_pre1000, outdir, args.dpi)

    print(f"\nDone. Strict-full-column main multi-panel figures written to: {outdir}")


if __name__ == "__main__":
    main()
