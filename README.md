# TDGL Mpemba seed-quality study

This repository contains simulation code, analysis scripts, data tables, and manuscript source files for the study:

**Seed quality and nucleation-controlled Mpemba-like transition times in time-dependent Ginzburg--Landau models**

## Summary

The study investigates Mpemba-like transition times in spatially extended TDGL models. The main conclusion is that in a first-order `phi6` model, hotter initial-state labels generate more barrier-crossing seeds, but these seeds become less compact and more spatially extended. Nucleation probability is therefore not determined by seed amount alone.

## Repository contents

- `src/`: simulation and figure-generation scripts
- `data/`: CSV outputs used for the manuscript figures
- `figures/`: final figure PDFs
- `manuscript/`: REVTeX manuscript and supplemental material sources
- `checks/`: parameter, reference, and submission-readiness checks
- `README_REPRODUCIBILITY.md`: commands and data-to-figure mapping

## Quick start

```bash
pip install -r requirements.txt
```

Generate figures from CSV files:

```bash
python src/tdgl_phi4_baseline_figure.py --data-dir data --outdir figures/main
python src/tdgl_manuscript_figures_v4_multipanel.py --data-dir data --outdir figures
```

See `README_REPRODUCIBILITY.md` for full details.

## License

A license has not yet been selected. See `LICENSE_TODO.md`.
