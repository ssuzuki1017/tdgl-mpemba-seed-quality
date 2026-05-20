#!/usr/bin/env python
r"""
Extended seed-geometry checks for the TDGL Mpemba seed-quality project.

This script does two things:

1. Bootstrap/permutation trend checks for existing seed-geometry metrics:
   - seed_compactness
   - seed_rg
   - p_seed
   - c_seed

2. Perimeter-to-area analysis if perimeter information is already present in the CSV files.
   The script auto-detects common perimeter column names such as:
   - seed_perimeter
   - seed_boundary
   - seed_interface_length
   - seed_perimeter_to_area
   - seed_perimeter_area_ratio

If no perimeter column is present, the script still completes and writes a short note explaining
that perimeter-to-area cannot be reconstructed from the current aggregate CSV files alone.

Example:
  python .\src\analyze_seed_geometry_trends_v1.py --data-dir .\data --outdir .\analysis_seed_geometry --bootstrap 2000 --permutations 20000
"""

from __future__ import annotations

from pathlib import Path
import argparse
import math
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

import matplotlib.pyplot as plt


REQUIRED_COLUMNS = [
    "label",
    "t_nuc_cluster",
    "p_seed",
    "c_seed",
    "seed_compactness",
    "seed_rg",
]

METRICS = [
    "p_seed",
    "c_seed",
    "seed_compactness",
    "seed_rg",
]

PERIMETER_COLUMNS = [
    "seed_perimeter",
    "seed_boundary",
    "seed_interface_length",
    "seed_boundary_length",
]

PERIMETER_RATIO_COLUMNS = [
    "seed_perimeter_to_area",
    "seed_perimeter_area_ratio",
    "seed_interface_to_area",
    "seed_boundary_to_area",
]


def find_main_files(data_dir: Path) -> list[Path]:
    candidates = []
    main_dir = data_dir / "main"
    for seed in [5678, 9876]:
        exact = main_dir / f"phi6_main_Df9em3_seed{seed}_fullcols.samples.csv"
        if exact.exists():
            candidates.append(exact)
            continue
        fallback = sorted(main_dir.glob(f"*seed{seed}*fullcols*.samples.csv"))
        if fallback:
            candidates.append(fallback[0])
            continue
    if len(candidates) < 2:
        raise FileNotFoundError(
            "Could not find both full-column main files under data/main. "
            "Expected phi6_main_Df9em3_seed5678_fullcols.samples.csv and "
            "phi6_main_Df9em3_seed9876_fullcols.samples.csv."
        )
    return candidates


def load_main_data(paths: Iterable[Path]) -> pd.DataFrame:
    frames = []
    for p in paths:
        df = pd.read_csv(p)
        df["source_file"] = p.name
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)

    missing = [c for c in REQUIRED_COLUMNS if c not in out.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    out["label"] = pd.to_numeric(out["label"], errors="coerce")
    out["nucleated"] = pd.to_numeric(out["t_nuc_cluster"], errors="coerce").notna().astype(int)

    for col in METRICS:
        out[col] = pd.to_numeric(out[col], errors="coerce")

    return out


def stderr(x: pd.Series) -> float:
    y = pd.to_numeric(x, errors="coerce").dropna()
    if len(y) <= 1:
        return np.nan
    return float(y.std(ddof=1) / math.sqrt(len(y)))


def label_summary(df: pd.DataFrame, metrics: list[str]) -> pd.DataFrame:
    rows = []
    for label, g in sorted(df.groupby("label"), key=lambda kv: float(kv[0])):
        row = {
            "label": float(label),
            "n": int(len(g)),
            "nucleated": int(g["nucleated"].sum()),
            "p_nuc": float(g["nucleated"].mean()),
        }
        row["p_nuc_se"] = math.sqrt(row["p_nuc"] * (1.0 - row["p_nuc"]) / row["n"]) if row["n"] else np.nan
        for m in metrics:
            row[f"{m}_mean"] = float(pd.to_numeric(g[m], errors="coerce").mean())
            row[f"{m}_se"] = stderr(g[m])
        rows.append(row)
    return pd.DataFrame(rows)


def spearman_rho(x: np.ndarray, y: np.ndarray) -> float:
    x = pd.Series(x).rank(method="average").to_numpy(dtype=float)
    y = pd.Series(y).rank(method="average").to_numpy(dtype=float)
    if np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def bootstrap_trend(df: pd.DataFrame, metric: str, rng: np.random.Generator, n_boot: int) -> np.ndarray:
    labels = sorted(df["label"].dropna().unique())
    values = []
    groups = {lab: df[df["label"] == lab] for lab in labels}

    for _ in range(n_boot):
        means = []
        ok = True
        for lab in labels:
            g = groups[lab]
            if len(g) == 0:
                ok = False
                break
            sample = g.sample(n=len(g), replace=True, random_state=int(rng.integers(0, 2**32 - 1)))
            means.append(pd.to_numeric(sample[metric], errors="coerce").mean())
        if ok:
            values.append(spearman_rho(np.asarray(labels, dtype=float), np.asarray(means, dtype=float)))
    return np.asarray(values, dtype=float)


def permutation_p_value(summary: pd.DataFrame, metric_mean_col: str, rng: np.random.Generator, n_perm: int) -> tuple[float, float]:
    labels = summary["label"].to_numpy(dtype=float)
    y = summary[metric_mean_col].to_numpy(dtype=float)
    obs = spearman_rho(labels, y)
    null = []
    for _ in range(n_perm):
        null.append(spearman_rho(labels, rng.permutation(y)))
    null = np.asarray(null, dtype=float)
    p = float((np.sum(np.abs(null) >= abs(obs)) + 1) / (len(null) + 1))
    return obs, p


def trend_tests(df: pd.DataFrame, summary: pd.DataFrame, metrics: list[str], rng: np.random.Generator, n_boot: int, n_perm: int) -> pd.DataFrame:
    rows = []
    labels = sorted(df["label"].dropna().unique())
    low_label = labels[0]
    high_label = labels[-1]

    for metric in metrics:
        mean_col = f"{metric}_mean"
        obs, perm_p = permutation_p_value(summary, mean_col, rng, n_perm=n_perm)
        boot = bootstrap_trend(df, metric, rng, n_boot=n_boot)
        low_mean = float(summary.loc[summary["label"] == low_label, mean_col].iloc[0])
        high_mean = float(summary.loc[summary["label"] == high_label, mean_col].iloc[0])
        rows.append({
            "metric": metric,
            "spearman_rho_label_means": obs,
            "permutation_p_two_sided": perm_p,
            "bootstrap_rho_mean": float(np.nanmean(boot)) if len(boot) else np.nan,
            "bootstrap_rho_ci_low": float(np.nanpercentile(boot, 2.5)) if len(boot) else np.nan,
            "bootstrap_rho_ci_high": float(np.nanpercentile(boot, 97.5)) if len(boot) else np.nan,
            "low_label": float(low_label),
            "high_label": float(high_label),
            "mean_low_label": low_mean,
            "mean_high_label": high_mean,
            "high_minus_low": high_mean - low_mean,
        })
    return pd.DataFrame(rows)


def detect_perimeter_metric(df: pd.DataFrame) -> tuple[pd.DataFrame, str | None, str]:
    df = df.copy()

    for col in PERIMETER_RATIO_COLUMNS:
        if col in df.columns:
            df["seed_perimeter_to_area_metric"] = pd.to_numeric(df[col], errors="coerce")
            return df, "seed_perimeter_to_area_metric", f"Using existing ratio column: {col}"

    for col in PERIMETER_COLUMNS:
        if col in df.columns:
            perimeter = pd.to_numeric(df[col], errors="coerce")
            if "seed_area" in df.columns:
                area = pd.to_numeric(df["seed_area"], errors="coerce")
                area_name = "seed_area"
            else:
                area = pd.to_numeric(df["c_seed"], errors="coerce")
                area_name = "c_seed"
            ratio = perimeter / area.replace(0, np.nan)
            df["seed_perimeter_to_area_metric"] = ratio
            return df, "seed_perimeter_to_area_metric", f"Computed ratio from perimeter column {col} and area proxy {area_name}."

    return df, None, "No perimeter or perimeter-to-area column found."


def save_trend_figure(summary: pd.DataFrame, metrics: list[str], out_base: Path, dpi: int) -> None:
    """Save a compact 2-column/3x2 grid version of the seed-geometry trend figure."""
    n = len(metrics)
    ncols = 2 if n <= 4 else 3
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.55 * ncols, 2.75 * nrows), constrained_layout=True)
    axes_flat = np.ravel(axes)

    pretty = {
        "p_seed": r"$p_{\rm seed}$",
        "c_seed": r"$c_{\rm seed}$",
        "seed_compactness": r"$C_{\rm comp}$",
        "seed_rg": r"$R_g$",
        "seed_perimeter_to_area_metric": r"$P_{\rm seed}/A_{\rm seed}$",
    }

    panel_letters = "abcdefghijklmnopqrstuvwxyz"
    for idx, (ax, metric) in enumerate(zip(axes_flat, metrics)):
        ax.errorbar(
            summary["label"],
            summary[f"{metric}_mean"],
            yerr=summary[f"{metric}_se"],
            marker="o",
            capsize=3,
            linewidth=1.2,
        )
        ax.set_xlabel("Initial-state label")
        ax.set_ylabel(pretty.get(metric, metric))
        ax.set_xticks(summary["label"])
        ax.set_xticklabels([f"{x:.2f}" for x in summary["label"]], rotation=30, ha="right")
        ax.text(
            0.03,
            0.94,
            f"({panel_letters[idx]})",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontweight="bold",
        )

    for ax in axes_flat[n:]:
        ax.axis("off")

    for ext in ["png", "pdf"]:
        path = out_base.with_suffix(f".{ext}")
        fig.savefig(path, dpi=dpi if ext == "png" else None, bbox_inches="tight")
        print(f"[saved] {path}")
    plt.close(fig)

def main() -> None:
    parser = argparse.ArgumentParser(description="Extended seed-geometry trend tests.")
    parser.add_argument("--data-dir", default="data")
    parser.add_argument("--outdir", default="analysis_seed_geometry")
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--permutations", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    data_dir = Path(args.data_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)

    main_files = find_main_files(data_dir)
    print("[main files]")
    for p in main_files:
        print(f"  {p}")

    df = load_main_data(main_files)
    df, perimeter_metric, perimeter_message = detect_perimeter_metric(df)

    metrics = METRICS.copy()
    if perimeter_metric is not None:
        metrics.append(perimeter_metric)
    else:
        note = outdir / "PERIMETER_TO_AREA_NOT_AVAILABLE.md"
        note.write_text(
            "# Perimeter-to-area ratio was not computed\n\n"
            f"{perimeter_message}\n\n"
            "The current full-column aggregate CSV files contain seed compactness and seed radius of gyration, "
            "but they do not contain the boundary length or perimeter of the largest barrier-crossing seed cluster.\n\n"
            "A true perimeter-to-area ratio cannot be reconstructed from these aggregate columns alone. "
            "To add this metric, modify the simulation/analysis code at the point where the largest seed cluster mask is available, "
            "compute its boundary length, and write a new column such as `seed_perimeter_to_area` or `seed_perimeter`.\n",
            encoding="utf-8",
        )
        print(f"[note] {note}")

    summary = label_summary(df, metrics)
    summary.to_csv(outdir / "seed_geometry_label_summary.csv", index=False)

    tests = trend_tests(df, summary, metrics, rng, n_boot=args.bootstrap, n_perm=args.permutations)
    tests.to_csv(outdir / "seed_geometry_trend_tests.csv", index=False)

    keep_cols = ["label", "nucleated"] + metrics
    df[keep_cols + ["source_file"]].to_csv(outdir / "seed_geometry_sample_data_used.csv", index=False)

    save_trend_figure(summary, metrics, outdir / "figS6_seed_geometry_trend_tests", args.dpi)

    readme = f"""# Extended seed-geometry trend checks

Input files:

{chr(10).join(f"- {p}" for p in main_files)}

Rows: {len(df)}
Observed nucleation events: {int(df["nucleated"].sum())}

Perimeter status:

{perimeter_message}

Outputs:

- `seed_geometry_label_summary.csv`
- `seed_geometry_trend_tests.csv`
- `seed_geometry_sample_data_used.csv`
- `figS6_seed_geometry_trend_tests.png`
- `figS6_seed_geometry_trend_tests.pdf`

Interpretation guide:

- A robust compactness trend should have a negative or positive Spearman trend consistent with the plotted label means and a bootstrap CI that does not straddle zero.
- A permutation p value is a label-level trend check, not a full mechanistic proof.
- If perimeter-to-area is unavailable, do not claim that interface length was tested.
- If perimeter-to-area becomes available later, rerun this script and it will include the perimeter metric automatically.
"""
    (outdir / "README_seed_geometry_trend_checks.md").write_text(readme, encoding="utf-8")

    print(f"\nDone. Outputs written to: {outdir}")


if __name__ == "__main__":
    main()
