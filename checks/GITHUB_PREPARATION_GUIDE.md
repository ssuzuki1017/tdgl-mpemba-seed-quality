# GitHub preparation guide

Recommended repository name:

```text
tdgl-mpemba-seed-quality
```

## Recommended visibility

Start as **private** while cleaning files. Switch to public when:

- the manuscript has been submitted or posted to arXiv,
- the data files included are the intended final public dataset,
- the README and license decision are finalized.

## Repository structure

```text
tdgl-mpemba-seed-quality/
  README.md
  README_REPRODUCIBILITY.md
  CITATION.cff
  requirements.txt
  LICENSE_TODO.md
  manuscript/
    manuscript_revtex_review_revised.tex
    supplemental_material_review_revised.tex
    references_verified.bib
  src/
    tdgl_run_auto_v3.py
    tdgl_mpemba_revised.py
    tdgl_phi4_baseline_figure.py
    tdgl_manuscript_figures_v4_multipanel.py
  data/
    phi4_samples.csv
    phi6_main_seed5678.samples.csv
    phi6_main_seed9876.samples.csv
    robustness/
      ...
  figures/
    main/
      ...
    supplemental/
      ...
  checks/
    TABLE_I_PARAMETER_CHECK.md
    REFERENCES_FINAL_CHECK.md
    PRE_SUBMISSION_CHECKLIST.md
```

## Git commands

From the repository root:

```bash
git init
git add .
git commit -m "Initial reproducibility package for TDGL Mpemba seed-quality study"
git branch -M main
git remote add origin https://github.com/<USER>/<REPO>.git
git push -u origin main
```

## License decision

Do not publish without deciding a license.

Common options:

- Code: MIT License or BSD-3-Clause License
- Manuscript/figures/data: CC BY 4.0 is common, but confirm whether you want this
- No license: others can view but do not have clear permission to reuse

A placeholder `LICENSE_TODO.md` is included instead of selecting a license on your behalf.

## Before public release

- Remove intermediate/failed-run files that are not needed for reproduction.
- Confirm that all included CSVs are intended to be public.
- Confirm that paths in README are relative, not local Windows paths.
- Confirm that the manuscript's Data Availability statement matches the repository status.
- Consider creating a Zenodo DOI after the repository is public and stable.
