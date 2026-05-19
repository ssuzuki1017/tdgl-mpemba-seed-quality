# Reproducibility guide

This document describes how to reproduce the smoke-test figure outputs for the TDGL Mpemba seed-quality study.

The intended workflow is:

1. install the Python dependencies,
2. verify that the required CSV files are present,
3. regenerate the phi4 baseline figure,
4. regenerate the phi6 main and supplemental multi-panel figures,
5. regenerate the sample-level seed-quality diagnostic used for Supplemental Fig. S5,
6. confirm that the expected files are produced.

The full stochastic simulations do not need to be rerun if the CSV files are already included under `data/`.

## 1. Environment

From the repository root:

```powershell
python --version
python -m pip install -r requirements.txt
```

The required Python packages are listed in `requirements.txt`:

```text
numpy
pandas
matplotlib
```

## 2. Required data files

The phi4 baseline figure expects:

```text
data/main/phi4_samples.csv
```

The phi6 v5 multi-panel figure script expects full-column main files containing the seed-quality columns:

```text
data/main/phi6_main_Df9em3_seed5678_fullcols.samples.csv
data/main/phi6_main_Df9em3_seed9876_fullcols.samples.csv
```

These files must include at least the following columns:

```text
label
t_nuc_cluster
t_tr
p_seed
c_seed
seed_compactness
seed_rg
p_ordered_seed
c_ordered_seed
```

The repository may also contain older long-name main CSV files without `seed_compactness` and `seed_rg`. Those files are useful historical outputs, but they are not sufficient for regenerating the v5 seed-quality figures.

Optional robustness files are read from `data/robustness/` when present.

To list available CSV files:

```powershell
Get-ChildItem .\data -Recurse -Filter "*.csv" | Select-Object FullName
```

## 3. Smoke-test output directory

Create a temporary output directory so the committed `figures/` directory is not modified:

```powershell
$root = (Get-Location).Path
$smoke = Join-Path $root "_smoke_figures"

Remove-Item $smoke -Recurse -Force -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Path $smoke | Out-Null
```

## 4. Regenerate the phi4 baseline figure

```powershell
python .\src\tdgl_phi4_baseline_figure.py `
  --data-dir .\data\main `
  --outdir "$smoke\main"
```

Expected outputs include:

```text
_smoke_figures/main/phi4_baseline_summary.csv
_smoke_figures/main/phi4_baseline_transition_time.png
_smoke_figures/main/phi4_baseline_transition_time.pdf
_smoke_figures/main/README_phi4_baseline_figure.md
```

## 5. Regenerate the phi6 v5 multi-panel figures

```powershell
python .\src\tdgl_manuscript_figures_v5_strict_fullcols_alias.py `
  --data-dir .\data `
  --outdir "$smoke\v5" `
  --max-survival-labels 3
```

Important notes:

- Do not pass `--formats`; this script always writes both `.png` and `.pdf`.
- `--data-dir .\data` is the recommended form.
- The script should select the full-column alias files in `data/main/` for seed `5678` and seed `9876`.

Expected outputs include:

```text
_smoke_figures/v5/fig2_nucleation_multipanel.png
_smoke_figures/v5/fig2_nucleation_multipanel.pdf
_smoke_figures/v5/fig3_seed_amount_multipanel.png
_smoke_figures/v5/fig3_seed_amount_multipanel.pdf
_smoke_figures/v5/fig4_seed_quality_multipanel.png
_smoke_figures/v5/fig4_seed_quality_multipanel.pdf
_smoke_figures/v5/figS1_ordered_seed_multipanel.png
_smoke_figures/v5/figS1_ordered_seed_multipanel.pdf
_smoke_figures/v5/figS2_numerical_checks_multipanel.png
_smoke_figures/v5/figS2_numerical_checks_multipanel.pdf
_smoke_figures/v5/figS3_robustness_multipanel.png
_smoke_figures/v5/figS3_robustness_multipanel.pdf
_smoke_figures/v5/summary_main_df9.csv
```

Check the output list with:

```powershell
Get-ChildItem "$smoke\v5" -File | Select-Object Name
```

## 6. Regenerate the sample-level seed-quality analysis

The sample-level seed-quality diagnostic used for Supplemental Fig. S5 can be regenerated with:

```powershell
python .\src\analyze_seed_quality_sample_level.py `
  --data-dir .\data `
  --outdir .\analysis_seed_quality
```

Expected outputs include:

```text
analysis_seed_quality/fig5_sample_level_seed_quality.png
analysis_seed_quality/fig5_sample_level_seed_quality.pdf
analysis_seed_quality/seed_quality_model_comparison.csv
analysis_seed_quality/seed_quality_logistic_coefficients.csv
analysis_seed_quality/seed_quality_binned_probabilities.csv
analysis_seed_quality/seed_quality_within_pseed_bins.csv
analysis_seed_quality/sample_level_main_data_used.csv
analysis_seed_quality/README_seed_quality_analysis.md
```

For the supplemental material, the figure files are copied as:

```text
figures/supplemental/figS5_sample_level_seed_quality.png
figures/supplemental/figS5_sample_level_seed_quality.pdf
```

## 7. Compile the manuscript

The manuscript uses REVTeX. A typical local build is:

```powershell
cd .\manuscript

pdflatex manuscript_revtex.tex
bibtex manuscript_revtex
pdflatex manuscript_revtex.tex
pdflatex manuscript_revtex.tex

pdflatex supplemental_material.tex
pdflatex supplemental_material.tex
```

Expected PDF outputs:

```text
manuscript/manuscript_revtex.pdf
manuscript/supplemental_material.pdf
```

If `pdflatex` or `bibtex` is not recognized, install a LaTeX distribution such as MiKTeX or TeX Live and reopen the terminal so that the PATH is updated.

## 8. Git hygiene

The `_smoke_figures/` and `analysis_seed_quality/` directories are temporary output directories and do not need to be committed unless you intentionally want to archive intermediate analysis tables.

LaTeX auxiliary files should normally not be committed:

```text
*.aux
*.bbl
*.blg
*.log
*.out
*Notes.bib
*.synctex.gz
```

After the test:

```powershell
git status
```

If temporary outputs appear as untracked files, remove them:

```powershell
Remove-Item .\_smoke_figures -Recurse -Force -ErrorAction SilentlyContinue
Remove-Item .\analysis_seed_quality -Recurse -Force -ErrorAction SilentlyContinue
```

## 9. Pass criteria

The smoke test passes if:

1. dependencies install without error,
2. `data/main/phi4_samples.csv` is detected,
3. the phi4 baseline PNG/PDF files are generated,
4. the two full-column phi6 main alias files are detected,
5. the six phi6 multi-panel PNG/PDF figure pairs are generated,
6. `summary_main_df9.csv` is generated,
7. the sample-level seed-quality diagnostic figure and CSV tables are generated,
8. the main manuscript and supplemental PDFs compile without undefined citations or references,
9. temporary output directories and LaTeX auxiliary files are not committed.
