#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
Sample-level seed-quality analysis for the TDGL Mpemba seed-quality project.

Purpose
-------
This script adds a reviewer-facing analysis that tests whether seed geometry
(seed compactness / radius of gyration) provides information about nucleation
beyond seed amount (p_seed / c_seed).

It uses only numpy, pandas, and matplotlib, matching the repository dependencies.

Typical use from repository root
--------------------------------
python .\src\analyze_seed_quality_sample_level.py --data-dir .\data --outdir .\analysis_seed_quality --bootstrap 500

Outputs
-------
analysis_seed_quality/
  fig5_sample_level_seed_quality.png
  fig5_sample_level_seed_quality.pdf
  seed_quality_model_comparison.csv
  seed_quality_logistic_coefficients.csv
  seed_quality_binned_probabilities.csv
  README_seed_quality_analysis.md
"""

from __future__ import annotations

from pathlib import Path
import argparse
import math
import sys
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# APS/PRE-friendly font embedding for PDF output.
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42


REQUIRED_COLUMNS = [
    "label",
    "t_nuc_cluster",
    "t_tr",
    "p_seed",
    "c_seed",
    "seed_compactness",
    "seed_rg",
]


def find_main_files(data_dir: Path) -> list[Path]:
    """Find the two full-column main CSV files."""
    candidates = []
    for seed in [5678, 9876]:
        patterns = [
            f"main/phi6_main_Df9em3_seed{seed}_fullcols.samples.csv",
            f"data/main/phi6_main_Df9em3_seed{seed}_fullcols.samples.csv",
        ]
        found = None
        for pat in patterns:
            p = data_dir / pat
            if p.exists():
                found = p
                break
        if found is None:
            # Fallback: recursive search.
            matches = sorted(data_dir.rglob(f"phi6_main_Df9em3_seed{seed}_fullcols.samples.csv"))
            if matches:
                found = matches[0]
        if found is None:
            raise FileNotFoundError(
                f"Could not find full-column main file for seed={seed}. "
                f"Expected data/main/phi6_main_Df9em3_seed{seed}_fullcols.samples.csv"
            )
        candidates.append(found.resolve())
    return candidates


def load_main(files: Iterable[Path]) -> pd.DataFrame:
    frames = []
    for p in files:
        df = pd.read_csv(p)
        df["source_file"] = p.name
        missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
        if missing:
            raise ValueError(f"{p} is missing required columns: {missing}")
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    out["nucleated"] = pd.to_numeric(out["t_nuc_cluster"], errors="coerce").notna().astype(int)
    out["transitioned"] = pd.to_numeric(out["t_tr"], errors="coerce").notna().astype(int)
    return out


def as_numeric_frame(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    x = pd.DataFrame(index=df.index)
    for c in columns:
        x[c] = pd.to_numeric(df[c], errors="coerce")
    return x


def standardize(x: pd.DataFrame) -> tuple[np.ndarray, dict[str, tuple[float, float]]]:
    stats: dict[str, tuple[float, float]] = {}

    # Intercept-only models pass an empty predictor frame.  Returning an
    # (n, 0) design block lets the caller still build a valid intercept-only
    # matrix with np.column_stack([ones, Xs]).
    if len(x.columns) == 0:
        return np.empty((len(x), 0), dtype=float), stats

    cols = []
    for c in x.columns:
        arr = x[c].to_numpy(dtype=float)
        mu = float(np.nanmean(arr))
        sd = float(np.nanstd(arr, ddof=0))
        if not np.isfinite(sd) or sd <= 0:
            sd = 1.0
        stats[c] = (mu, sd)
        cols.append((arr - mu) / sd)
    return np.column_stack(cols), stats


def sigmoid(z: np.ndarray) -> np.ndarray:
    z = np.clip(z, -35.0, 35.0)
    return 1.0 / (1.0 + np.exp(-z))


def fit_logistic(X: np.ndarray, y: np.ndarray, ridge: float = 1e-8, max_iter: int = 100) -> dict:
    """Logistic regression by Newton/IRLS. X must already include intercept."""
    beta = np.zeros(X.shape[1], dtype=float)
    penalty = np.eye(X.shape[1]) * ridge
    penalty[0, 0] = 0.0

    for _ in range(max_iter):
        eta = X @ beta
        mu = sigmoid(eta)
        W = np.clip(mu * (1.0 - mu), 1e-9, None)
        z = eta + (y - mu) / W
        XtW = X.T * W
        H = XtW @ X + penalty
        rhs = XtW @ z
        try:
            beta_new = np.linalg.solve(H, rhs)
        except np.linalg.LinAlgError:
            beta_new = np.linalg.lstsq(H, rhs, rcond=None)[0]
        if np.max(np.abs(beta_new - beta)) < 1e-8:
            beta = beta_new
            break
        beta = beta_new

    eta = X @ beta
    mu = np.clip(sigmoid(eta), 1e-12, 1 - 1e-12)
    ll = float(np.sum(y * np.log(mu) + (1 - y) * np.log(1 - mu)))
    k = int(X.shape[1])
    aic = 2 * k - 2 * ll

    W = np.clip(mu * (1.0 - mu), 1e-9, None)
    H = (X.T * W) @ X + penalty
    try:
        cov = np.linalg.inv(H)
        se = np.sqrt(np.maximum(np.diag(cov), 0.0))
    except np.linalg.LinAlgError:
        se = np.full_like(beta, np.nan)

    return {"beta": beta, "se": se, "loglik": ll, "aic": aic, "pred": mu}


def normal_two_sided_p(z: float) -> float:
    if not np.isfinite(z):
        return np.nan
    return math.erfc(abs(z) / math.sqrt(2.0))


def fit_model(df: pd.DataFrame, predictors: list[str], outcome: str = "nucleated") -> tuple[dict, pd.DataFrame]:
    # Build the analysis frame explicitly so that an empty predictor list
    # works for the intercept-only model.
    use = as_numeric_frame(df, predictors + [outcome]).dropna()
    if len(use) == 0:
        raise ValueError(f"No valid rows remain for predictors={predictors} and outcome={outcome}.")
    y = use[outcome].to_numpy(dtype=float)
    x_raw = use[predictors] if predictors else pd.DataFrame(index=use.index)
    Xs, stats = standardize(x_raw)
    X = np.column_stack([np.ones(len(use)), Xs])
    fit = fit_logistic(X, y)
    names = ["intercept"] + predictors
    rows = []
    for i, name in enumerate(names):
        beta = float(fit["beta"][i])
        se = float(fit["se"][i]) if i < len(fit["se"]) else np.nan
        z = beta / se if se and np.isfinite(se) and se > 0 else np.nan
        rows.append({
            "term": name,
            "coef_standardized": beta,
            "std_error": se,
            "z_approx": z,
            "p_approx": normal_two_sided_p(z),
            "odds_ratio_per_1sd": float(np.exp(np.clip(beta, -30.0, 30.0))) if np.isfinite(beta) else np.nan,
        })
    coef_df = pd.DataFrame(rows)
    fit["n"] = int(len(use))
    fit["events"] = int(y.sum())
    fit["predictors"] = predictors
    fit["stats"] = stats
    return fit, coef_df


def binomial_se(p: float, n: int) -> float:
    if n <= 0:
        return np.nan
    return math.sqrt(p * (1 - p) / n)


def quantile_bins(s: pd.Series, q: int) -> pd.Series:
    # duplicates="drop" prevents failure when quantile cutpoints are repeated.
    return pd.qcut(pd.to_numeric(s, errors="coerce"), q=q, duplicates="drop")


def binned_probability(df: pd.DataFrame, variable: str, q: int = 4) -> pd.DataFrame:
    tmp = df[[variable, "nucleated"]].copy()
    tmp = tmp.dropna()
    tmp["bin"] = quantile_bins(tmp[variable], q)
    rows = []
    for b, g in tmp.groupby("bin", observed=True):
        n = len(g)
        p = float(g["nucleated"].mean())
        rows.append({
            "analysis": f"{variable}_quantile",
            "bin": str(b),
            "x_mean": float(g[variable].mean()),
            "n": n,
            "events": int(g["nucleated"].sum()),
            "p_nuc": p,
            "p_nuc_se": binomial_se(p, n),
        })
    return pd.DataFrame(rows)


def compactness_within_pseed_bins(df: pd.DataFrame, q: int = 3) -> pd.DataFrame:
    tmp = df[["p_seed", "seed_compactness", "nucleated"]].copy().dropna()
    tmp["pseed_bin"] = quantile_bins(tmp["p_seed"], q)
    rows = []
    for pbin, g in tmp.groupby("pseed_bin", observed=True):
        med = float(g["seed_compactness"].median())
        for group_name, gg in [
            ("low_compactness", g[g["seed_compactness"] <= med]),
            ("high_compactness", g[g["seed_compactness"] > med]),
        ]:
            n = len(gg)
            if n == 0:
                continue
            p = float(gg["nucleated"].mean())
            rows.append({
                "analysis": "compactness_median_split_within_pseed_tercile",
                "bin": str(pbin),
                "group": group_name,
                "x_mean": float(gg["seed_compactness"].mean()),
                "p_seed_mean": float(gg["p_seed"].mean()),
                "n": n,
                "events": int(gg["nucleated"].sum()),
                "p_nuc": p,
                "p_nuc_se": binomial_se(p, n),
            })
    return pd.DataFrame(rows)


def plot_results(df: pd.DataFrame, binned: pd.DataFrame, within: pd.DataFrame, outdir: Path, dpi: int):
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.0), constrained_layout=True)

    # Panel A: P_nuc vs p_seed quantiles
    ax = axes[0, 0]
    a = binned[binned["analysis"] == "p_seed_quantile"]
    ax.errorbar(a["x_mean"], a["p_nuc"], yerr=a["p_nuc_se"], marker="o", capsize=3)
    ax.set_xlabel(r"Mean $p_{\mathrm{seed}}$ in quantile bin")
    ax.set_ylabel(r"Nucleation probability $P_{\mathrm{nuc}}$")
    ax.text(0.03, 0.95, "(a)", transform=ax.transAxes, va="top", ha="left", fontweight="bold")

    # Panel B: P_nuc vs compactness quantiles
    ax = axes[0, 1]
    b = binned[binned["analysis"] == "seed_compactness_quantile"]
    ax.errorbar(b["x_mean"], b["p_nuc"], yerr=b["p_nuc_se"], marker="o", capsize=3)
    ax.set_xlabel("Mean seed compactness in quantile bin")
    ax.set_ylabel(r"Nucleation probability $P_{\mathrm{nuc}}$")
    ax.text(0.03, 0.95, "(b)", transform=ax.transAxes, va="top", ha="left", fontweight="bold")

    # Panel C: high vs low compactness within p_seed terciles
    ax = axes[1, 0]
    for i, (bin_name, g) in enumerate(within.groupby("bin", observed=True)):
        x_low = i - 0.12
        x_high = i + 0.12
        low = g[g["group"] == "low_compactness"]
        high = g[g["group"] == "high_compactness"]
        if len(low):
            ax.errorbar([x_low], low["p_nuc"], yerr=low["p_nuc_se"], marker="o", capsize=3, label="low compactness" if i == 0 else None)
        if len(high):
            ax.errorbar([x_high], high["p_nuc"], yerr=high["p_nuc_se"], marker="s", capsize=3, label="high compactness" if i == 0 else None)
    ax.set_xticks(range(len(within["bin"].drop_duplicates())))
    ax.set_xticklabels(list(within["bin"].drop_duplicates()), rotation=25, ha="right")
    ax.set_xlabel(r"$p_{\mathrm{seed}}$ tercile")
    ax.set_ylabel(r"Nucleation probability $P_{\mathrm{nuc}}$")
    ax.legend(fontsize=8)
    ax.text(0.03, 0.95, "(c)", transform=ax.transAxes, va="top", ha="left", fontweight="bold")

    # Panel D: raw scatter for compactness and p_seed, with nucleation marker encoded by y jitter
    ax = axes[1, 1]
    tmp = df[["p_seed", "seed_compactness", "nucleated"]].dropna().copy()
    rng = np.random.default_rng(12345)
    y = tmp["nucleated"].to_numpy(dtype=float) + rng.normal(0, 0.03, len(tmp))
    ax.scatter(tmp["seed_compactness"], y, s=14, alpha=0.55)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["not nucleated", "nucleated"])
    ax.set_xlabel("Seed compactness")
    ax.set_ylabel("Outcome")
    ax.text(0.03, 0.95, "(d)", transform=ax.transAxes, va="top", ha="left", fontweight="bold")

    for ext in ["png", "pdf"]:
        path = outdir / f"fig5_sample_level_seed_quality.{ext}"
        fig.savefig(path, dpi=dpi if ext == "png" else None, bbox_inches="tight")
        print(f"[saved] {path}")
    plt.close(fig)


def write_readme(outdir: Path, model_df: pd.DataFrame):
    best = model_df.sort_values("aic").iloc[0]
    text = f"""# Sample-level seed-quality analysis

This analysis tests whether seed geometry adds predictive information about nucleation
beyond seed amount.

## Outcome

The binary outcome is:

```text
nucleated = observed t_nuc_cluster by the final observation time
```

## Models

All logistic-regression predictors are standardized before fitting. The reported odds
ratios are therefore per one standard deviation of the predictor.

## Main model-comparison result

Best AIC model:

```text
{best['model']}
```

AIC table:

```text
{model_df.to_string(index=False)}
```

## Interpretation rule

For the manuscript, use the result conservatively:

- If `amount_plus_quality` or `all_predictors` improves AIC over `seed_amount`,
  it supports the claim that seed geometry contains predictive information beyond
  seed amount.
- If not, weaken the claim to: the label-averaged seed geometry changes
  systematically, but the present sample-level statistics do not independently
  establish compactness as a predictor of nucleation.

"""
    (outdir / "README_seed_quality_analysis.md").write_text(text, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default="data", help="Repository data directory or repository root.")
    parser.add_argument("--outdir", default="analysis_seed_quality")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--bootstrap", type=int, default=0, help="Reserved for future use; current script uses asymptotic SEs.")
    args = parser.parse_args()

    data_dir = Path(args.data_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    files = find_main_files(data_dir)
    print("[main files]")
    for p in files:
        print(" ", p)

    df = load_main(files)
    df.to_csv(outdir / "sample_level_main_data_used.csv", index=False)

    print(f"[data] n={len(df)}, nucleated={int(df['nucleated'].sum())}, transitioned={int(df['transitioned'].sum())}")

    models = {
        "intercept_only": [],
        "label_only": ["label"],
        "seed_amount": ["p_seed", "c_seed"],
        "seed_quality": ["seed_compactness", "seed_rg"],
        "amount_plus_quality": ["p_seed", "c_seed", "seed_compactness", "seed_rg"],
        "all_predictors": ["label", "p_seed", "c_seed", "seed_compactness", "seed_rg"],
    }

    model_rows = []
    coef_frames = []
    best_aic = None

    for name, preds in models.items():
        fit, coef = fit_model(df, preds)
        model_rows.append({
            "model": name,
            "predictors": ", ".join(preds) if preds else "(intercept only)",
            "n": fit["n"],
            "events": fit["events"],
            "k": len(preds) + 1,
            "loglik": fit["loglik"],
            "aic": fit["aic"],
        })
        coef.insert(0, "model", name)
        coef_frames.append(coef)

    model_df = pd.DataFrame(model_rows).sort_values("aic")
    model_df["delta_aic"] = model_df["aic"] - model_df["aic"].min()
    model_df.to_csv(outdir / "seed_quality_model_comparison.csv", index=False)
    pd.concat(coef_frames, ignore_index=True).to_csv(outdir / "seed_quality_logistic_coefficients.csv", index=False)

    binned = pd.concat([
        binned_probability(df, "p_seed", q=4),
        binned_probability(df, "c_seed", q=4),
        binned_probability(df, "seed_compactness", q=4),
        binned_probability(df, "seed_rg", q=4),
    ], ignore_index=True)
    within = compactness_within_pseed_bins(df, q=3)
    binned.to_csv(outdir / "seed_quality_binned_probabilities.csv", index=False)
    within.to_csv(outdir / "seed_quality_within_pseed_bins.csv", index=False)

    plot_results(df, binned, within, outdir, args.dpi)
    write_readme(outdir, model_df)

    print("[saved]", outdir / "seed_quality_model_comparison.csv")
    print("[saved]", outdir / "seed_quality_logistic_coefficients.csv")
    print("[saved]", outdir / "seed_quality_binned_probabilities.csv")
    print("[saved]", outdir / "seed_quality_within_pseed_bins.csv")
    print("[saved]", outdir / "README_seed_quality_analysis.md")
    print("\nDone.")


if __name__ == "__main__":
    main()
