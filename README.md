# TDGL Mpemba seed-quality study

This repository contains simulation code, data tables, figure-generation scripts, figures, and manuscript sources for the study:

**Seed quality and nucleation-controlled Mpemba-like transition times in time-dependent Ginzburg--Landau models**

## Overview

This project studies Mpemba-like transition times in spatially extended time-dependent Ginzburg--Landau models.

We compare:

1. a continuous-transition `phi4` TDGL baseline model, and
2. a first-order `phi6` TDGL model in which the post-quench state is metastable and transitions proceed through stochastic nucleation.

The main result is that, in the first-order model, larger initial-state labels generate more barrier-crossing seeds, but those seeds become less compact and more spatially extended. Therefore, the nucleation probability is not determined solely by the amount of barrier-crossing seed; seed quality, especially compactness and spatial extent, is also important.

## Repository structure

```text
tdgl-mpemba-seed-quality/
  README.md
  README_REPRODUCIBILITY.md
  CITATION.cff
  LICENSE
  requirements.txt

  manuscript/
    manuscript_revtex.tex
    supplemental_material.tex
    references_verified.bib

  src/
    tdgl_run_auto_v3.py
    tdgl_phi4_baseline_figure.py
    tdgl_manuscript_figures_v5_strict_fullcols_alias.py
    make_tnuc_ttr_supplement_v3.py

  data/
    main/
    processed/
    robustness/

  figures/
    main/
    supplemental/
```

## Requirements

The scripts use the standard scientific Python stack:

```text
numpy
pandas
matplotlib
```

Install dependencies with:

```bash
pip install -r requirements.txt
```

## Reproducibility

For numerical parameters, expected input files, figure-generation commands, and manuscript compilation instructions, see:

```text
README_REPRODUCIBILITY.md
```

If the CSV files are already included in `data/`, re-running the full simulations is not required to reproduce the manuscript figures.

## Main outputs

Main manuscript figures:

```text
figures/main/phi4_baseline_transition_time.pdf
figures/main/fig2_nucleation_multipanel.pdf
figures/main/fig3_seed_amount_multipanel.pdf
figures/main/fig4_seed_quality_multipanel.pdf
```

Supplemental figures:

```text
figures/supplemental/figS1_ordered_seed_multipanel.pdf
figures/supplemental/figS2_numerical_checks_multipanel.pdf
figures/supplemental/figS3_robustness_multipanel.pdf
```

## Figure generation smoke test

From the repository root, the manuscript figures can be regenerated from the included CSV files with:

```powershell
$root = (Get-Location).Path
$smoke = Join-Path $root "_smoke_figures"

python .\src\tdgl_phi4_baseline_figure.py `
  --data-dir .\data\main `
  --outdir "$smoke\main"

python .\src\tdgl_manuscript_figures_v5_strict_fullcols_alias.py `
  --data-dir .\data `
  --outdir "$smoke\v5" `
  --max-survival-labels 3
```

The second command writes both PNG and PDF versions of the multi-panel manuscript figures.

## Manuscript

The manuscript and supplemental material are written in REVTeX format.

Main manuscript:

```text
manuscript/manuscript_revtex.tex
```

Supplemental material:

```text
manuscript/supplemental_material.tex
```

Bibliography:

```text
manuscript/references_verified.bib
```

## Data availability

The simulation data and analysis scripts used to generate the figures are available in this repository.

For an archival release, create a GitHub release and archive it with Zenodo. After obtaining a DOI, update the manuscript data availability statement accordingly.

## License

Code in this repository is licensed under the MIT License.

Unless otherwise noted, manuscript text, figures, and data are made available under the Creative Commons Attribution 4.0 International License.
