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
  LICENSE
  requirements.txt

  manuscript/
    manuscript_revtex.tex
    supplemental_material.tex
    references_verified.bib

  src/
    tdgl_run_auto_v3.py
    tdgl_manuscript_figures_v4_multipanel.py
    tdgl_phi4_baseline_figure.py
    legacy/

  data/
    main/
      phi4_samples.csv
      phi6_N64_dt0p02_tmax300_pre500_Df9em3_D00p02_ns50_*seed5678*cth20*.samples.csv
      phi6_N64_dt0p02_tmax300_pre500_Df9em3_D00p02_ns50_*seed9876*cth20*.samples.csv
    robustness/
      additional CSV files for threshold, time-step, system-size, pre-equilibration, and exploratory checks

  figures/
    main/
      phi4_baseline_transition_time.pdf
      fig2_nucleation_multipanel.pdf
      fig3_seed_amount_multipanel.pdf
      fig4_seed_quality_multipanel.pdf
    supplemental/
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
    GITHUB_DIRECTORY_CONTENTS.md
```

## Main manuscript figures

| Figure | File | Description |
|---|---|---|
| Fig. 1 | `figures/main/phi4_baseline_transition_time.pdf` | Continuous-transition `phi4` baseline transition time |
| Fig. 2 | `figures/main/fig2_nucleation_multipanel.pdf` | First-order `phi6` nucleation probability and survival curves |
| Fig. 3 | `figures/main/fig3_seed_amount_multipanel.pdf` | Barrier-crossing seed amount: `p_seed` and `c_seed` |
| Fig. 4 | `figures/main/fig4_seed_quality_multipanel.pdf` | Seed quality: compactness and radius of gyration |

## Supplemental figures

| Figure | File | Description |
|---|---|---|
| Fig. S1 | `figures/supplemental/figS1_p_ordered_seed.pdf` | Ordered-like seed fraction |
| Fig. S2 | `figures/supplemental/figS2_c_ordered_seed.pdf` | Largest ordered-like seed cluster size |
| Fig. S3 | `figures/supplemental/figS3_cluster_threshold_nucleation_probability.pdf` | Cluster-threshold robustness |
| Fig. S4 | `figures/supplemental/figS4_dt_dependence_nucleation_probability.pdf` | Time-step robustness |
| Fig. S5 | `figures/supplemental/figS5_p_seed_robustness.pdf` | Robustness of `p_seed` to system size and pre-equilibration |
| Fig. S6 | `figures/supplemental/figS6_seed_compactness_robustness.pdf` | Robustness of seed compactness to system size and pre-equilibration |

## Reproducibility

See [`README_REPRODUCIBILITY.md`](README_REPRODUCIBILITY.md) for numerical parameters, expected input files, figure-generation commands, and manuscript compilation.

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

1. Upload the repository contents, including the `manuscript/` and `figures/` folders.
2. Set `manuscript/manuscript_revtex.tex` as the main document for the main manuscript.
3. Set `manuscript/supplemental_material.tex` as the main document for the supplemental material.
4. Use pdfLaTeX.

The bibliography file is:

```text
manuscript/references_verified.bib
```

## Data availability

The simulation data and analysis scripts used to generate the figures are available in this repository.

For a final archival version, creating a GitHub release and linking it to Zenodo is recommended. After obtaining a Zenodo DOI, the data availability statement in the manuscript should be updated to cite the DOI.

## Citation

If you use this repository, please cite the manuscript and/or the repository metadata in `CITATION.cff`.

## License

Code in this repository is licensed under the MIT License.

Unless otherwise noted, manuscript text, figures, and data are made available under the Creative Commons Attribution 4.0 International License (CC BY 4.0).
