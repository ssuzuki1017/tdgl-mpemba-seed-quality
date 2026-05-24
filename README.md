# TDGL seed-geometry study

This repository contains simulation code, processed data tables, figure-generation scripts, manuscript figures, and LaTeX sources for the manuscript:

**Seed geometry in nucleation-controlled transition times in time-dependent Ginzburg--Landau models**

## Overview

This project studies transition times in spatially extended time-dependent Ginzburg--Landau (TDGL) models motivated by phase-transition Mpemba-like protocols.

The main comparison is between:

1. a continuous-transition `phi4` TDGL baseline model, and
2. a first-order `phi6` TDGL model in which the post-quench state is metastable and the transition proceeds through stochastic nucleation.

The central result is conservative: seed amount alone is not a sufficient reduced description of nucleation-controlled transition times. In the first-order model, larger initial-state labels generate more barrier-crossing seed material, but the largest seed clusters become less compact and more spatially extended. Seed geometry is therefore an independent morphology diagnostic that should be tracked separately from seed amount. The simple scalar geometry measures used here are not claimed to be complete committor-level predictors.

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
    tdgl_phi4_baseline_figure.py
    tdgl_manuscript_figures_v5_strict_fullcols_alias.py
    make_fig5_seed_geometry_transition.py
    analyze_seed_quality_sample_level.py
    analyze_seed_geometry_trends_v1.py
    committor_restart_phi6_v1.py
    tdgl_run_auto_v3.py
    tdgl_mpemba_timeseries_export_v3.py
    tdgl_mpemba_revised.py

  data/
    main/
    robustness/
    processed/
      committor_v1/
      committor_v2/

  figures/
    main/
    supplemental/

  checks/
```

## Main manuscript figures

| Figure | File | Description |
|---|---|---|
| Fig. 1 | `figures/main/phi4_baseline_transition_time.pdf` | Continuous-transition `phi4` baseline transition time |
| Fig. 2 | `figures/main/fig2_nucleation_multipanel.pdf` | First-order `phi6` nucleation probability and survival curves |
| Fig. 3 | `figures/main/fig3_seed_amount_multipanel.pdf` | Barrier-crossing seed amount: `p_seed` and `c_seed` |
| Fig. 4 | `figures/main/fig4_seed_quality_multipanel.pdf` | Seed compactness and radius of gyration |
| Fig. 5 | `figures/main/fig5_seed_geometry_transition.pdf` | Sample-level seed amount--geometry relation and `t_nuc` versus `t_tr` |

## Supplemental figures

| Figure | File | Description |
|---|---|---|
| Fig. S1 | `figures/supplemental/figS1_ordered_seed_multipanel.pdf` | Ordered-like seed statistics |
| Fig. S2 | `figures/supplemental/figS2_numerical_checks_multipanel.pdf` | Numerical and operational checks for nucleation probability |
| Fig. S3 | `figures/supplemental/figS3_robustness_multipanel.pdf` | System-size and pre-equilibration robustness checks |
| Fig. S4 | `figures/supplemental/figS7_tnuc_vs_ttr_scatter.pdf` | Cluster-nucleation time versus global transition time |
| Fig. S5 | `figures/supplemental/figS5_sample_level_seed_quality.pdf` | Sample-level seed-quality diagnostics |
| Fig. S6 | `figures/supplemental/figS6_seed_geometry_trend_tests.pdf` | Bootstrap/permutation trend checks including perimeter-to-area |
| Fig. S7 | `figures/supplemental/figS7_committor_seed_geometry.pdf` | Finite-time committor-style restart analysis |

## Reproducibility

See [`README_REPRODUCIBILITY.md`](README_REPRODUCIBILITY.md) for numerical parameters, expected input files, figure-generation commands, analysis commands, and manuscript compilation instructions.

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

## Manuscript compilation

The manuscript source is written in REVTeX format.

Main manuscript:

```bash
cd manuscript
pdflatex manuscript_revtex.tex
bibtex manuscript_revtex
pdflatex manuscript_revtex.tex
pdflatex manuscript_revtex.tex
```

Supplemental material:

```bash
pdflatex supplemental_material.tex
pdflatex supplemental_material.tex
```

## Data and code availability

The simulation data, processed tables, random seeds, figure-generation scripts, manuscript figures, and manuscript sources are archived at Zenodo: [10.5281/zenodo.20357154](https://doi.org/10.5281/zenodo.20357154). They are also available in this GitHub repository.

## Citation

If you use this repository, cite the archived Zenodo release, DOI [10.5281/zenodo.20357154](https://doi.org/10.5281/zenodo.20357154), and the associated manuscript. Repository citation metadata are provided in [`CITATION.cff`](CITATION.cff).

## License

Code in this repository is licensed under the MIT License. Unless otherwise noted, manuscript text, figures, and data are made available under the Creative Commons Attribution 4.0 International License (CC BY 4.0).
