# Reproducibility notes

This document describes the files and commands needed to reproduce the manuscript figures for:

**Seed quality and nucleation-controlled Mpemba-like transition times in time-dependent Ginzburg--Landau models**

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

The main first-order `phi6` calculations use the following parameters.

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
| Barrier-crossing seed threshold | `|phi| > phi_b` |
| Ordered-like seed threshold | `|phi| > 0.5 phi_s` |
| Nucleation cluster threshold | `C_th = 20` |
| Persistence condition | `n_cons = 3` |

Main `phi6` results combine two independent random seeds:

```text
seed = 5678, nsamples = 50 per label
seed = 9876, nsamples = 50 per label
combined = 100 realizations per label
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

The `data/` directory should contain the CSV files used to generate the figures.

Minimum expected files:

```text
data/
  phi4_samples.csv

  phi6_N64_dt0p02_tmax300_pre500_Df9em3_D00p02_ns50_*seed5678*cth20*.samples.csv
  phi6_N64_dt0p02_tmax300_pre500_Df9em3_D00p02_ns50_*seed9876*cth20*.samples.csv

  phi6_N64_dt0p02_tmax300_pre500_Df9em3_*cth10*.samples.csv
  phi6_N64_dt0p02_tmax300_pre500_Df9em3_*cth30*.samples.csv
  phi6_N64_dt0p01_tmax300_pre500_Df9em3_*cth20*.samples.csv
  phi6_N128_dt0p02_tmax300_pre500_Df9em3_*cth40*.samples.csv
  phi6_N64_dt0p02_tmax300_pre1000_Df9em3_*cth20*.samples.csv
```

The exact filenames may be longer because they encode the numerical parameters. The figure scripts search for the relevant patterns automatically.

Important columns for `phi6` sample files include:

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
```

Important columns for `phi4_samples.csv` include:

```text
label
t_tr
```

## 4. Reproducing figures

Run commands from the repository root.

### 4.1 phi4 baseline figure

```bash
python src/tdgl_phi4_baseline_figure.py --data-dir data --outdir figures
```

Expected output:

```text
figures/phi4_baseline_transition_time.pdf
figures/phi4_baseline_transition_time.png
figures/phi4_baseline_summary.csv
```

### 4.2 phi6 main and supplemental figures

```bash
python src/tdgl_manuscript_figures_v4_multipanel.py --data-dir data --outdir figures
```

Expected main outputs:

```text
figures/fig2_nucleation_multipanel.pdf
figures/fig3_seed_amount_multipanel.pdf
figures/fig4_seed_quality_multipanel.pdf
```

Expected supplemental outputs:

```text
figures/figS1_p_ordered_seed.pdf
figures/figS2_c_ordered_seed.pdf
figures/figS3_cluster_threshold_nucleation_probability.pdf
figures/figS4_dt_dependence_nucleation_probability.pdf
figures/figS5_p_seed_robustness.pdf
figures/figS6_seed_compactness_robustness.pdf
```

If the survival-curve panel in Fig. 2 is visually crowded, use:

```bash
python src/tdgl_manuscript_figures_v4_multipanel.py --data-dir data --outdir figures --max-survival-labels 3
```

## 5. Re-running simulations

If the CSV files are already included in `data/`, re-running simulations is not required to reproduce the figures.

If simulations need to be re-run, use the main run script, for example:

```bash
python src/tdgl_run_auto_v3.py phi6 --N 64 --nsamples 50 --tmax 300 --preeq_steps 500 --af 0.02 --Df 9e-3 --D0 0.02 --cluster_threshold 20 --seed 5678
```

Second independent run:

```bash
python src/tdgl_run_auto_v3.py phi6 --N 64 --nsamples 50 --tmax 300 --preeq_steps 500 --af 0.02 --Df 9e-3 --D0 0.02 --cluster_threshold 20 --seed 9876
```

Robustness examples:

```bash
python src/tdgl_run_auto_v3.py phi6 --N 64 --nsamples 30 --tmax 300 --preeq_steps 500 --af 0.02 --Df 9e-3 --D0 0.02 --cluster_threshold 10 --seed 5678

python src/tdgl_run_auto_v3.py phi6 --N 64 --nsamples 30 --tmax 300 --preeq_steps 500 --af 0.02 --Df 9e-3 --D0 0.02 --cluster_threshold 30 --seed 5678

python src/tdgl_run_auto_v3.py phi6 --N 64 --nsamples 30 --tmax 300 --preeq_steps 500 --af 0.02 --Df 9e-3 --D0 0.02 --cluster_threshold 20 --dt 0.01 --seed 5678
```

The exact supported command-line options depend on the final version of `tdgl_run_auto_v3.py`. If a command fails, run:

```bash
python src/tdgl_run_auto_v3.py --help
```

## 6. Manuscript compilation

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

The manuscript was compiled using REVTeX 4.2 with pdfLaTeX.

In Overleaf:

1. Upload the repository contents.
2. Set `manuscript/manuscript_revtex.tex` as the main document.
3. Compile with pdfLaTeX.
4. For the supplement, set `manuscript/supplemental_material.tex` as the main document and compile separately.

## 7. Figure mapping

| Manuscript figure | Source file |
|---|---|
| Fig. 1 | `figures/phi4_baseline_transition_time.pdf` |
| Fig. 2 | `figures/fig2_nucleation_multipanel.pdf` |
| Fig. 3 | `figures/fig3_seed_amount_multipanel.pdf` |
| Fig. 4 | `figures/fig4_seed_quality_multipanel.pdf` |
| Fig. S1 | `figures/figS1_p_ordered_seed.pdf` |
| Fig. S2 | `figures/figS2_c_ordered_seed.pdf` |
| Fig. S3 | `figures/figS3_cluster_threshold_nucleation_probability.pdf` |
| Fig. S4 | `figures/figS4_dt_dependence_nucleation_probability.pdf` |
| Fig. S5 | `figures/figS5_p_seed_robustness.pdf` |
| Fig. S6 | `figures/figS6_seed_compactness_robustness.pdf` |

## 8. Data availability statement

For a private repository:

```text
The simulation data and analysis scripts used to generate the figures are available from the author upon reasonable request.
```

For a public GitHub repository:

```text
The simulation data and analysis scripts used to generate the figures are available at
https://github.com/ssuzuki1017/tdgl-mpemba-seed-quality.
```

For an archival release, create a GitHub release and connect it to Zenodo, then replace the GitHub URL with the Zenodo DOI.

## 9. Checklist before public release

- [ ] Remove or replace `LICENSE_TODO.md`.
- [ ] Add a final `LICENSE` file.
- [ ] Confirm that no private paths or personal local directories remain in scripts or README files.
- [ ] Confirm that all figure files can be regenerated from files in `data/`.
- [ ] Confirm that `manuscript_revtex.tex` and `supplemental_material.tex` compile in Overleaf.
- [ ] Confirm that `CITATION.cff` contains the final title and author information.
- [ ] Decide whether to keep the repository private until publication or make it public before submission.
