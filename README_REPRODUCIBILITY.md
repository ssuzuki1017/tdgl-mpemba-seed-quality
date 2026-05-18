# README_REPRODUCIBILITY

Project: Seed quality and nucleation-controlled Mpemba-like transition times in TDGL models  
Date: 2026-05-18

## Purpose

This file records how to reproduce the numerical data summaries and figures used in the manuscript.

## Python environment

Recommended:

```bash
python >= 3.10
pip install numpy pandas matplotlib
```

The scripts use only standard-library modules plus:

- `numpy`
- `pandas`
- `matplotlib`

## Main simulation scripts

- `tdgl_run_auto_v3.py`  
  Main automated simulation runner used for `phi6` runs.
- `tdgl_mpemba_revised.py`  
  Revised TDGL simulation code.
- `tdgl_phi4_baseline_figure.py`  
  Generates the `phi4` baseline figure from `phi4_samples.csv`.
- `tdgl_manuscript_figures_v4_multipanel.py`  
  Generates the main multi-panel figures and supplemental figures.

## Main `phi6` data files

The main `phi6` results use two independent random seeds:

```text
phi6_N64_dt0p02_tmax300_pre500_Df9em3_D00p02_ns50_initpreeq_me5_seed5678_af0p02_b1_c1_schemenoise_only_aibase0p3_qfrac0p5_pth1em3_cth20_mc3_labels_default_run2.samples.csv
phi6_N64_dt0p02_tmax300_pre500_Df9em3_D00p02_ns50_initpreeq_me5_seed9876_af0p02_b1_c1_schemenoise_only_aibase0p3_qfrac0p5_pth1em3_cth20_mc3_labels_default_run2.samples.csv
```

Together these give 100 stochastic realizations per initial-state label.

## Main `phi6` command pattern

The two main runs correspond to:

```powershell
python .	dgl_run_auto_v3.py phi6 --N 64 --nsamples 50 --tmax 300 --preeq_steps 500 --af 0.02 --Df 9e-3 --D0 0.02 --cluster_threshold 20 --seed 5678
python .	dgl_run_auto_v3.py phi6 --N 64 --nsamples 50 --tmax 300 --preeq_steps 500 --af 0.02 --Df 9e-3 --D0 0.02 --cluster_threshold 20 --seed 9876
```

## Robustness check command patterns

```powershell
# cluster_threshold = 10
python .	dgl_run_auto_v3.py phi6 --N 64 --nsamples 30 --tmax 300 --preeq_steps 500 --af 0.02 --Df 9e-3 --D0 0.02 --cluster_threshold 10 --seed 5678

# cluster_threshold = 30
python .	dgl_run_auto_v3.py phi6 --N 64 --nsamples 30 --tmax 300 --preeq_steps 500 --af 0.02 --Df 9e-3 --D0 0.02 --cluster_threshold 30 --seed 5678

# dt = 0.01
python .	dgl_run_auto_v3.py phi6 --N 64 --nsamples 30 --tmax 300 --dt 0.01 --preeq_steps 500 --af 0.02 --Df 9e-3 --D0 0.02 --cluster_threshold 20 --seed 5678

# N = 128 check
python .	dgl_run_auto_v3.py phi6 --N 128 --nsamples 20 --tmax 300 --preeq_steps 500 --af 0.02 --Df 9e-3 --D0 0.02 --cluster_threshold 40 --seed 5678

# pre-equilibration = 1000 check
python .	dgl_run_auto_v3.py phi6 --N 64 --nsamples 30 --tmax 300 --preeq_steps 1000 --af 0.02 --Df 9e-3 --D0 0.02 --cluster_threshold 20 --seed 5678
```

## Figure generation

Run from the directory containing the `.samples.csv` files:

```powershell
python .	dgl_phi4_baseline_figure.py --formats png,pdf
python .	dgl_manuscript_figures_v4_multipanel.py --formats png,pdf
```

The final manuscript package uses the review-revised figures with simplified x-axis tick labels.

## Manuscript compilation

Main manuscript:

```text
manuscript_revtex_review_revised.tex
references_verified.bib
figures/*.pdf
```

Supplemental Material:

```text
supplemental_material_review_revised.tex
figures/*.pdf
```

Compile in Overleaf or locally with a typical sequence:

```bash
pdflatex manuscript_revtex_review_revised.tex
bibtex manuscript_revtex_review_revised
pdflatex manuscript_revtex_review_revised.tex
pdflatex manuscript_revtex_review_revised.tex

pdflatex supplemental_material_review_revised.tex
pdflatex supplemental_material_review_revised.tex
```

## Figure mapping

| manuscript figure | file |
|---|---|
| Fig. 1 | `phi4_baseline_transition_time.pdf` |
| Fig. 2 | `fig2_nucleation_multipanel.pdf` |
| Fig. 3 | `fig3_seed_amount_multipanel.pdf` |
| Fig. 4 | `fig4_seed_quality_multipanel.pdf` |
| Fig. S1 | `figS1_p_ordered_seed.pdf` |
| Fig. S2 | `figS2_c_ordered_seed.pdf` |
| Fig. S3 | `figS3_cluster_threshold_nucleation_probability.pdf` |
| Fig. S4 | `figS4_dt_dependence_nucleation_probability.pdf` |
| Fig. S5 | `figS5_p_seed_robustness.pdf` |
| Fig. S6 | `figS6_seed_compactness_robustness.pdf` |

## Data availability wording

Current manuscript wording:

```text
The simulation data and analysis scripts used to generate the figures are available from the author upon reasonable request.
```

If a GitHub or Zenodo repository is made public, replace this with the repository URL or DOI.
