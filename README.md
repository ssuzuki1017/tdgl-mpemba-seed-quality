# TDGL Mpemba seed-quality study

This repository contains simulation code, analysis scripts, data tables, figure files, and manuscript sources for the study:

**Seed quality and nucleation-controlled Mpemba-like transition times in time-dependent Ginzburg--Landau models**

## Overview

This project studies Mpemba-like transition times in spatially extended time-dependent Ginzburg--Landau (TDGL) models.

The main comparison is between:

1. a continuous-transition `phi4` TDGL baseline model, and
2. a first-order `phi6` TDGL model in which the post-quench state is metastable and the transition proceeds through stochastic nucleation.

The central result is that, in the first-order model, larger temperature-like initial-state labels generate more barrier-crossing seeds, but those seeds become less compact and more spatially extended. Therefore, the nucleation probability is not determined solely by the amount of barrier-crossing seed; seed quality, especially compactness and spatial extent, is also important.

## Repository structure

```text
tdgl-mpemba-seed-quality/
  README.md
  README_REPRODUCIBILITY.md
  CITATION.cff
  requirements.txt
  LICENSE_TODO.md

  manuscript/
    manuscript_revtex.tex
    supplemental_material.tex
    references_verified.bib
    manuscript_revtex.pdf                 # optional compiled output
    supplemental_material.pdf              # optional compiled output

  src/
    tdgl_run_auto_v3.py
    tdgl_manuscript_figures_v4_multipanel.py
    tdgl_phi4_baseline_figure.py
    tdgl_manuscript_figures_v3.py          # optional legacy figure script
    tdgl_mpemba_revised.py                 # optional exploratory script

  data/
    phi4_samples.csv
    phi6_*.samples.csv
    summary_*.csv                          # optional analysis summaries

  figures/
    phi4_baseline_transition_time.pdf
    fig2_nucleation_multipanel.pdf
    fig3_seed_amount_multipanel.pdf
    fig4_seed_quality_multipanel.pdf
    figS1_p_ordered_seed.pdf
    figS2_c_ordered_seed.pdf
    figS3_cluster_threshold_nucleation_probability.pdf
    figS4_dt_dependence_nucleation_probability.pdf
    figS5_p_seed_robustness.pdf
    figS6_seed_compactness_robustness.pdf

  checks/
    TABLE_I_PARAMETER_CHECK.md
    REFERENCES_FINAL_CHECK.md
    PRE_SUBMISSION_CHECKLIST.md
    GITHUB_PREPARATION_GUIDE.md
```

File names may differ slightly depending on the final local export. The essential requirement is that the data files, scripts, manuscript sources, and figure files needed to reproduce the manuscript figures are all included.

## Main manuscript figures

| Figure | File | Description |
|---|---|---|
| Fig. 1 | `figures/phi4_baseline_transition_time.pdf` | Continuous-transition `phi4` baseline transition time |
| Fig. 2 | `figures/fig2_nucleation_multipanel.pdf` | First-order `phi6` nucleation probability and survival curves |
| Fig. 3 | `figures/fig3_seed_amount_multipanel.pdf` | Barrier-crossing seed amount: `p_seed` and `c_seed` |
| Fig. 4 | `figures/fig4_seed_quality_multipanel.pdf` | Seed quality: compactness and radius of gyration |

## Supplemental figures

| Figure | File | Description |
|---|---|---|
| Fig. S1 | `figures/figS1_p_ordered_seed.pdf` | Ordered-like seed fraction |
| Fig. S2 | `figures/figS2_c_ordered_seed.pdf` | Largest ordered-like seed cluster size |
| Fig. S3 | `figures/figS3_cluster_threshold_nucleation_probability.pdf` | Cluster-threshold robustness |
| Fig. S4 | `figures/figS4_dt_dependence_nucleation_probability.pdf` | Time-step robustness |
| Fig. S5 | `figures/figS5_p_seed_robustness.pdf` | Robustness of `p_seed` to system size and pre-equilibration |
| Fig. S6 | `figures/figS6_seed_compactness_robustness.pdf` | Robustness of seed compactness to system size and pre-equilibration |

## Reproducibility

See [`README_REPRODUCIBILITY.md`](README_REPRODUCIBILITY.md) for details on the numerical parameters, expected input files, figure-generation commands, and manuscript compilation.

## Requirements

The Python scripts were developed with the standard scientific Python stack:

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

To compile in Overleaf:

1. Upload the repository contents or the `manuscript/` and `figures/` folders.
2. Set `manuscript_revtex.tex` as the main document for the main manuscript.
3. Set `supplemental_material.tex` as the main document for the supplemental material.
4. Use pdfLaTeX.

The bibliography file is:

```text
manuscript/references_verified.bib
```

## Data availability

If this repository is public, the data availability statement in the manuscript can be written as:

```text
The simulation data and analysis scripts used to generate the figures are available at
https://github.com/ssuzuki1017/tdgl-mpemba-seed-quality.
```

If the repository remains private during submission, use:

```text
The simulation data and analysis scripts used to generate the figures are available from the author upon reasonable request.
```

For a final archival version, creating a GitHub release and linking it to Zenodo is recommended.

## Citation

If you use this repository, please cite the manuscript and/or the repository metadata in `CITATION.cff`.

## License

A license has not yet been finalized. See `LICENSE_TODO.md`.

Recommended choices:

- Code: MIT License or BSD-3-Clause License
- Manuscript, figures, and data: CC BY 4.0

Before making the repository public, replace `LICENSE_TODO.md` with the selected license file and clarify the reuse terms in this README.

## License

Code in this repository is licensed under the MIT License.

Unless otherwise noted, manuscript text, figures, and data are made available under the Creative Commons Attribution 4.0 International License (CC BY 4.0).
