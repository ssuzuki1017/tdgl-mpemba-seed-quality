"""
Auto-naming launcher for TDGL Mpemba simulations.

Place this file in the same directory as:
  - tdgl_mpemba_timeseries_export_v3.py
  - tdgl_mpemba_revised.py

This wrapper automatically creates output filenames and figure directories from
the run parameters, so you do not have to type file names manually.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


DEFAULT_LABELS = (1.05, 1.10, 1.20, 1.50, 2.00, 3.00)


def safe_float(value: float | str) -> str:
    """Convert floats into filename-safe strings."""
    if isinstance(value, str):
        try:
            value = float(value)
        except ValueError:
            return value.replace(".", "p").replace("-", "m").replace("+", "").replace(" ", "")

    v = float(value)
    if (abs(v) > 0 and abs(v) < 1e-2) or abs(v) >= 1e3:
        s = f"{v:.0e}"  # e.g. 1e-04
        mant, exp = s.split("e")
        exp_int = int(exp)
        mant = mant.replace("-", "m").replace(".", "p")
        return f"{mant}em{abs(exp_int)}" if exp_int < 0 else f"{mant}ep{exp_int}"

    return f"{v:.6g}".replace("-", "m").replace(".", "p")


def label_tag(labels: list[float]) -> str:
    if tuple(labels) == DEFAULT_LABELS:
        return "labels_default"
    return "labels_" + "_".join(safe_float(x) for x in labels)


def build_stem(args: argparse.Namespace) -> str:
    parts = [
        args.model,
        f"N{args.N}",
        f"dt{safe_float(args.dt)}",
        f"tmax{safe_float(args.tmax)}",
        f"pre{args.preeq_steps}",
        f"Df{safe_float(args.Df)}",
        f"D0{safe_float(args.D0)}",
        f"ns{args.nsamples}",
        f"init{args.init_mode}",
        f"me{args.measure_every}",
        f"seed{args.seed}",
    ]

    if args.model == "phi6":
        parts.extend([
            f"af{safe_float(args.af)}",
            f"b{safe_float(args.b)}",
            f"c{safe_float(args.c)}",
            f"scheme{args.init_scheme}",
            f"aibase{safe_float(args.a_i_base)}",
            f"qfrac{safe_float(args.q_fraction)}",
            f"pth{safe_float(args.ordered_fraction_threshold)}",
            f"cth{args.cluster_threshold}",
            f"mc{args.min_consecutive}",
        ])
    else:
        parts.append(f"qth{safe_float(args.q_threshold)}")

    if args.labels:
        parts.append(label_tag(args.labels))

    if args.tag:
        clean_tag = args.tag.replace(" ", "_").replace("/", "_").replace("\\", "_").replace(":", "_")
        parts.append(clean_tag)

    return "_".join(parts)


def unique_path(path: Path, overwrite: bool) -> Path:
    if overwrite or not path.exists():
        return path

    i = 2
    while True:
        candidate = path.with_name(f"{path.stem}_run{i}{path.suffix}")
        if not candidate.exists():
            return candidate
        i += 1


def build_command(args: argparse.Namespace, out_npz: Path, figdir: Path | None) -> list[str]:
    script = Path(__file__).resolve().parent / "tdgl_mpemba_timeseries_export_v3.py"

    if args.model == "phi4":
        cmd = [
            sys.executable, str(script), "run_phi4",
            "--out", str(out_npz),
            "--N", str(args.N),
            "--dx", str(args.dx),
            "--dt", str(args.dt),
            "--tmax", str(args.tmax),
            "--preeq_steps", str(args.preeq_steps),
            "--Df", str(args.Df),
            "--D0", str(args.D0),
            "--nsamples", str(args.nsamples),
            "--seed", str(args.seed),
            "--init_mode", args.init_mode,
            "--measure_every", str(args.measure_every),
            "--q_threshold", str(args.q_threshold),
        ]
    else:
        cmd = [
            sys.executable, str(script), "run_phi6",
            "--out", str(out_npz),
            "--N", str(args.N),
            "--dx", str(args.dx),
            "--dt", str(args.dt),
            "--tmax", str(args.tmax),
            "--preeq_steps", str(args.preeq_steps),
            "--af", str(args.af),
            "--b", str(args.b),
            "--c", str(args.c),
            "--Df", str(args.Df),
            "--D0", str(args.D0),
            "--nsamples", str(args.nsamples),
            "--seed", str(args.seed),
            "--init_mode", args.init_mode,
            "--init_scheme", args.init_scheme,
            "--a_i_base", str(args.a_i_base),
            "--measure_every", str(args.measure_every),
            "--q_fraction", str(args.q_fraction),
            "--ordered_fraction_threshold", str(args.ordered_fraction_threshold),
            "--cluster_threshold", str(args.cluster_threshold),
            "--min_consecutive", str(args.min_consecutive),
        ]

    if figdir is not None:
        cmd.extend(["--figdir", str(figdir)])

    if args.labels:
        cmd.append("--labels")
        cmd.extend(str(x) for x in args.labels)

    if args.stderr:
        cmd.append("--stderr")

    return cmd


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run TDGL simulations with automatic output names.")
    p.add_argument("model", choices=["phi4", "phi6"])

    p.add_argument("--root", default="runs", help="Root directory for all auto-named outputs")
    p.add_argument("--results_dir", default=None, help="Override NPZ output directory")
    p.add_argument("--figures_dir", default=None, help="Override figure output directory")
    p.add_argument("--no_figs", action="store_true")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--dry_run", action="store_true")
    p.add_argument("--tag", default=None, help="Optional tag added to filename")
    p.add_argument("--labels", nargs="*", type=float, default=None)

    p.add_argument("--N", type=int, default=64)
    p.add_argument("--dx", type=float, default=1.0)
    p.add_argument("--dt", type=float, default=None)
    p.add_argument("--tmax", type=float, default=None)
    p.add_argument("--preeq_steps", type=int, default=500)
    p.add_argument("--Df", type=float, default=None)
    p.add_argument("--D0", type=float, default=None)
    p.add_argument("--nsamples", type=int, default=5)
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--init_mode", choices=["preeq", "gaussian"], default="preeq")
    p.add_argument("--measure_every", type=int, default=None)
    p.add_argument("--stderr", action="store_true")

    # phi4
    p.add_argument("--q_threshold", type=float, default=0.5)

    # phi6
    p.add_argument("--af", type=float, default=0.10)
    p.add_argument("--b", type=float, default=1.0)
    p.add_argument("--c", type=float, default=1.0)
    p.add_argument("--init_scheme", choices=["noise_only", "mass_and_noise"], default="noise_only")
    p.add_argument("--a_i_base", type=float, default=0.30)
    p.add_argument("--q_fraction", type=float, default=0.5)
    p.add_argument("--ordered_fraction_threshold", type=float, default=1.0e-3)
    p.add_argument("--cluster_threshold", type=int, default=10)
    p.add_argument("--min_consecutive", type=int, default=3)

    return p


def fill_defaults(args: argparse.Namespace) -> argparse.Namespace:
    args.labels = args.labels if args.labels else list(DEFAULT_LABELS)

    if args.model == "phi4":
        args.dt = 0.05 if args.dt is None else args.dt
        args.tmax = 10.0 if args.tmax is None else args.tmax
        args.Df = 0.0 if args.Df is None else args.Df
        args.D0 = 0.05 if args.D0 is None else args.D0
        args.seed = 1234 if args.seed is None else args.seed
        args.measure_every = 1 if args.measure_every is None else args.measure_every

    if args.model == "phi6":
        args.dt = 0.02 if args.dt is None else args.dt
        args.tmax = 100.0 if args.tmax is None else args.tmax
        args.Df = 1.0e-4 if args.Df is None else args.Df
        args.D0 = 0.02 if args.D0 is None else args.D0
        args.seed = 5678 if args.seed is None else args.seed
        args.measure_every = 5 if args.measure_every is None else args.measure_every

    return args


def main() -> None:
    args = fill_defaults(build_parser().parse_args())

    root = Path(args.root)
    results_dir = Path(args.results_dir) if args.results_dir else root / "results_ts"
    figures_root = Path(args.figures_dir) if args.figures_dir else root / "figures"

    stem = build_stem(args)
    out_npz = unique_path(results_dir / f"{stem}.npz", overwrite=args.overwrite)
    figdir = None if args.no_figs else figures_root / out_npz.stem

    out_npz.parent.mkdir(parents=True, exist_ok=True)
    if figdir is not None:
        figdir.mkdir(parents=True, exist_ok=True)

    cmd = build_command(args, out_npz, figdir)

    print("Output NPZ:")
    print(f"  {out_npz}")
    if figdir is not None:
        print("Figure directory:")
        print(f"  {figdir}")
    print("\nCommand:")
    print("  " + " ".join(f'"{x}"' if " " in x else x for x in cmd))
    print()

    if args.dry_run:
        print("Dry run only. Nothing executed.")
        return

    completed = subprocess.run(cmd)
    if completed.returncode != 0:
        raise SystemExit(completed.returncode)

    print("\nFinished.")
    print(f"NPZ: {out_npz}")
    print(f"Sample CSV: {out_npz.with_suffix('.samples.csv')}")
    if figdir is not None:
        print(f"Figures: {figdir}")


if __name__ == "__main__":
    main()
