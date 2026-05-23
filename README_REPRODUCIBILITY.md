# Reproducibility notes

This document describes how to reproduce the figures and analyses for:

**Seed geometry in nucleation-controlled transition times in time-dependent Ginzburg--Landau models**

All commands below are intended to be run from the repository root unless otherwise noted.

## 1. Environment

Recommended environment:

```text
Python >= 3.10
numpy
pandas
matplotlib
```

Install dependencies:

```bash
pip install -r requirements.txt
```

The scripts do not require GPU computation.

## 2. Main numerical parameters

The main first-order `phi6` calculations use:

| Parameter | Value |
|---|---|
| Model | first-order `phi6` TDGL |
| System size | `N = 64` |
| Grid spacing | `dx = 1` |
| Time step | `dt = 0.02` |
| Final time | `tmax = 300` |
| Pre-equilibration steps | `preeq_steps = 500` |
| Final potential parameter | `a_f = 0.02` |
| Landau coefficients | `b = c = 1` |
| Initial fluctuation scale | `D_i = D_0 T_i`, `D_0 = 0.02` |
| Post-quench noise strength | `D_f = 9e-3` |
| Initial-state labels | `1.05, 1.10, 1.20, 1.50, 2.00, 3.00` |
| Barrier-crossing seed threshold | `abs(phi) > phi_b` |
| Ordered-like seed threshold | `abs(phi) > 0.5 phi_s` |
| Nucleation cluster threshold | `C_th = 20` |
| Persistence condition | `n_cons = 3` |

The main `phi6` results combine two independent random seeds:

```text
seed = 5678, nsamples = 50 per label
seed = 9876, nsamples = 50 per label
combined = 100 realizations per label
total = 600 realizations
```

The continuous-transition `phi4` baseline uses:

```text
N = 64
dt = 0.05
tmax = 10
preeq_steps = 500
5 realizations per label
```

Robustness checks include:

```text
cluster_threshold = 10 and 30
dt = 0.01
N = 128
preeq_steps = 1000
```

## 3. Data files

The `data/` directory contains the CSV files used to generate the manuscript figures.

Expected main data files:

```text
data/main/phi4_samples.csv
data/main/phi6_main_Df9em3_seed5678_fullcols.samples.csv
data/main/phi6_main_Df9em3_seed9876_fullcols.samples.csv
```

Expected robustness data are under:

```text
data/robustness/
```

Processed restart-analysis tables are under:

```text
data/processed/committor_v1/
data/processed/committor_v2/
```

Important columns for `phi6` sample files:

```text
label
t_nuc_cluster
t_tr
p_seed
c_seed
p_ordered_seed
c_ordered_seed
seed_compactness
seed_rg
seed_perimeter
seed_perimeter_to_area
```

Important columns for `phi4_samples.csv`:

```text
label
t_tr
```

## 4. Reproducing the main manuscript figures

### Fig. 1: continuous-transition `phi4` baseline

```bash
python src/tdgl_phi4_baseline_figure.py --data-dir data/main --outdir figures/main
```

Expected outputs:

```text
figures/main/phi4_baseline_transition_time.pdf
figures/main/phi4_baseline_transition_time.png
figures/main/phi4_baseline_summary.csv
```

### Figs. 2--4 and supplemental multipanel figures

```bash
python src/tdgl_manuscript_figures_v5_strict_fullcols_alias.py --data-dir data --outdir figures --max-survival-labels 3
```

Expected main outputs:

```text
figures/main/fig2_nucleation_multipanel.pdf
figures/main/fig2_nucleation_multipanel.png
figures/main/fig3_seed_amount_multipanel.pdf
figures/main/fig3_seed_amount_multipanel.png
figures/main/fig4_seed_quality_multipanel.pdf
figures/main/fig4_seed_quality_multipanel.png
```

Expected supplemental outputs include:

```text
figures/supplemental/figS1_ordered_seed_multipanel.pdf
figures/supplemental/figS2_numerical_checks_multipanel.pdf
figures/supplemental/figS3_robustness_multipanel.pdf
```

### Fig. 5: sample-level seed geometry and transition timing

```bash
python src/make_fig5_seed_geometry_transition.py --data-dir data --outdir figures/main
```

Expected outputs:

```text
figures/main/fig5_seed_geometry_transition.pdf
figures/main/fig5_seed_geometry_transition.png
figures/main/fig5_seed_geometry_transition_summary.csv
figures/main/fig5_seed_geometry_transition_numbers.tex
```

## 5. Supplemental analyses

### Sample-level seed-quality analysis

```bash
python src/analyze_seed_quality_sample_level.py --data-dir data --outdir analysis_seed_quality --bootstrap 500
```

Expected outputs include:

```text
analysis_seed_quality/fig5_sample_level_seed_quality.pdf
analysis_seed_quality/seed_quality_model_comparison.csv
analysis_seed_quality/seed_quality_logistic_coefficients.csv
```

The figure used in the supplement is archived as:

```text
figures/supplemental/figS5_sample_level_seed_quality.pdf
figures/supplemental/figS5_sample_level_seed_quality.png
```

### Seed-geometry bootstrap/permutation analysis

```bash
python src/analyze_seed_geometry_trends_v1.py --data-dir data --outdir analysis_seed_geometry_perimeter_grid --bootstrap 2000 --permutations 20000
```

Expected outputs include:

```text
analysis_seed_geometry_perimeter_grid/figS6_seed_geometry_trend_tests.pdf
analysis_seed_geometry_perimeter_grid/seed_geometry_trend_tests.csv
```

The figure used in the supplement is archived as:

```text
figures/supplemental/figS6_seed_geometry_trend_tests.pdf
figures/supplemental/figS6_seed_geometry_trend_tests.png
```

### Finite-time committor-style restart analysis

The larger restart analysis used in the current supplement is reproduced with:

```bash
python src/committor_restart_phi6_v1.py --outdir analysis_committor_v2 --n-configs-per-label 12 --n-restarts 24
```

Expected scale:

```text
72 frozen initial configurations
1728 restart trajectories
```

Expected outputs:

```text
analysis_committor_v2/committor_initial_configs.csv
analysis_committor_v2/committor_restarts.csv
analysis_committor_v2/committor_label_summary.csv
analysis_committor_v2/committor_geometry_correlations.csv
analysis_committor_v2/figS7_committor_seed_geometry.pdf
analysis_committor_v2/figS7_committor_seed_geometry.png
```

Archived processed tables are stored under:

```text
data/processed/committor_v2/
```

## 6. Re-running simulations

If the processed CSV files are already included in `data/`, re-running simulations is not required to reproduce the manuscript figures.

Example main `phi6` simulation commands:

```bash
python src/tdgl_run_auto_v3.py phi6 --N 64 --nsamples 50 --tmax 300 --preeq_steps 500 --af 0.02 --Df 0.009 --D0 0.02 --cluster_threshold 20 --seed 5678

python src/tdgl_run_auto_v3.py phi6 --N 64 --nsamples 50 --tmax 300 --preeq_steps 500 --af 0.02 --Df 0.009 --D0 0.02 --cluster_threshold 20 --seed 9876
```

If command-line options differ, inspect:

```bash
python src/tdgl_run_auto_v3.py --help
```

## 7. Manuscript compilation

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

The bibliography file is:

```text
manuscript/references_verified.bib
```

## 8. Figure mapping

| Manuscript figure | Source file |
|---|---|
| Fig. 1 | `figures/main/phi4_baseline_transition_time.pdf` |
| Fig. 2 | `figures/main/fig2_nucleation_multipanel.pdf` |
| Fig. 3 | `figures/main/fig3_seed_amount_multipanel.pdf` |
| Fig. 4 | `figures/main/fig4_seed_quality_multipanel.pdf` |
| Fig. 5 | `figures/main/fig5_seed_geometry_transition.pdf` |
| Fig. S1 | `figures/supplemental/figS1_ordered_seed_multipanel.pdf` |
| Fig. S2 | `figures/supplemental/figS2_numerical_checks_multipanel.pdf` |
| Fig. S3 | `figures/supplemental/figS3_robustness_multipanel.pdf` |
| Fig. S4 | `figures/supplemental/figS7_tnuc_vs_ttr_scatter.pdf` |
| Fig. S5 | `figures/supplemental/figS5_sample_level_seed_quality.pdf` |
| Fig. S6 | `figures/supplemental/figS6_seed_geometry_trend_tests.pdf` |
| Fig. S7 | `figures/supplemental/figS7_committor_seed_geometry.pdf` |

## 9. Notes

The initial-state labels are protocol labels rather than calibrated thermodynamic temperatures.

The finite-time committor-style restart analysis is an operational restart diagnostic, not an exact transition-path committor.

The main conclusion is not that compactness alone predicts nucleation, but that seed-amount-only descriptions are incomplete for the present nucleation-controlled transition-time protocol.
