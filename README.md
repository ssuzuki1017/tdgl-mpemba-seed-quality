# TDGL Mpemba seed-quality study

This repository contains simulation code, data tables, figure-generation scripts, figures, and manuscript sources for the study:

**Seed quality and nucleation-controlled transition times in TDGL models motivated by phase-transition Mpemba effects**

## Overview

This project studies transition times in spatially extended time-dependent Ginzburg--Landau models motivated by phase-transition Mpemba-like relaxation.

We compare:

1. a continuous-transition `phi4` TDGL baseline model, and
2. a first-order `phi6` TDGL model in which the post-quench state is metastable and transitions proceed through stochastic nucleation.

The main result is that, in the first-order model, larger initial-state labels generate more barrier-crossing seed candidates, but those seeds become less compact and more spatially extended. Therefore, the amount of barrier-crossing seed is not by itself sufficient to characterize nucleation-controlled transition times. Seed geometry, especially compactness and spatial extent, should be analyzed separately.

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
    manuscript_revtex.pdf
    supplemental_material.pdf

  src/
    tdgl_run_auto_v3.py
    tdgl_phi4_baseline_figure.py
    tdgl_manuscript_figures_v5_strict_fullcols_alias.py
    analyze_seed_quality_sample_level.py
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

For expected input files, figure-generation commands, smoke-test instructions, and manuscript compilation notes, see:

```text
README_REPRODUCIBILITY.md
```

If the CSV files are already included in `data/`, re-running the full stochastic simulations is not required to reproduce the manuscript figures.

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
figures/supplemental/figS5_sample_level_seed_quality.pdf
```

The supplemental material also discusses the relation between nucleation and global transition times.

## Figure-generation smoke test

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

The sample-level seed-quality analysis used for Supplemental Fig. S5 can be regenerated with:

```powershell
python .\src\analyze_seed_quality_sample_level.py `
  --data-dir .\data `
  --outdir .\analysis_seed_quality
```

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

Compiled PDFs:

```text
manuscript/manuscript_revtex.pdf
manuscript/supplemental_material.pdf
```

## Data availability

The simulation data and analysis scripts used to generate the figures are available in this repository.

For an archival release, create a GitHub release and archive it with Zenodo. After obtaining a DOI, update the manuscript data availability statement accordingly.

## License

Code in this repository is licensed under the MIT License.

Unless otherwise noted, manuscript text, figures, and data are made available under the Creative Commons Attribution 4.0 International License.
