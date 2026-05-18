# Table I parameter check

Date: 2026-05-18

## Result

Table I in `manuscript_revtex_review_revised.tex` is consistent with the simulation filenames and CSV sample counts available in the project directory.

## Main `phi6` data set

Input CSVs:

- `phi6_N64_dt0p02_tmax300_pre500_Df9em3_D00p02_ns50_initpreeq_me5_seed5678_af0p02_b1_c1_schemenoise_only_aibase0p3_qfrac0p5_pth1em3_cth20_mc3_labels_default_run2.samples.csv`
- `phi6_N64_dt0p02_tmax300_pre500_Df9em3_D00p02_ns50_initpreeq_me5_seed9876_af0p02_b1_c1_schemenoise_only_aibase0p3_qfrac0p5_pth1em3_cth20_mc3_labels_default_run2.samples.csv`

Combined sample count:

| label | samples |
|---:|---:|
| 1.05 | 100 |
| 1.10 | 100 |
| 1.20 | 100 |
| 1.50 | 100 |
| 2.00 | 100 |
| 3.00 | 100 |


This confirms: **100 realizations per label** for the main `phi6` figures.

## `phi4` baseline

Input CSV: `phi4_samples.csv`

| label | samples |
|---:|---:|
| 1.05 | 5 |
| 1.10 | 5 |
| 1.20 | 5 |
| 1.50 | 5 |
| 2.00 | 5 |
| 3.00 | 5 |

This confirms: **5 realizations per label** for the `phi4` baseline.

## Robustness checks

| check | CSV | rows | samples per label |
|---|---|---:|---:|
| cluster threshold check Cth=10 | `phi6_N64_dt0p02_tmax300_pre500_Df9em3_D00p02_ns30_initpreeq_me5_seed5678_af0p02_b1_c1_schemenoise_only_aibase0p3_qfrac0p5_pth1em3_cth10_mc3_labels_default.samples.csv` | 180 | 30 |
| cluster threshold check Cth=30 | `phi6_N64_dt0p02_tmax300_pre500_Df9em3_D00p02_ns30_initpreeq_me5_seed5678_af0p02_b1_c1_schemenoise_only_aibase0p3_qfrac0p5_pth1em3_cth30_mc3_labels_default.samples.csv` | 180 | 30 |
| time step check dt=0.01 | `phi6_N64_dt0p01_tmax300_pre500_Df9em3_D00p02_ns30_initpreeq_me5_seed5678_af0p02_b1_c1_schemenoise_only_aibase0p3_qfrac0p5_pth1em3_cth20_mc3_labels_default.samples.csv` | 180 | 30 |
| N=128 check | `phi6_N128_dt0p02_tmax300_pre500_Df9em3_D00p02_ns20_initpreeq_me5_seed5678_af0p02_b1_c1_schemenoise_only_aibase0p3_qfrac0p5_pth1em3_cth40_mc3_labels_default.samples.csv` | 120 | 20 |
| preeq=1000 check | `phi6_N64_dt0p02_tmax300_pre1000_Df9em3_D00p02_ns30_initpreeq_me5_seed5678_af0p02_b1_c1_schemenoise_only_aibase0p3_qfrac0p5_pth1em3_cth20_mc3_labels_default.samples.csv` | 180 | 30 |


## Table I values checked

| item | value | status |
|---|---:|---|
| Main system size | `N=64` | OK |
| Grid spacing | `dx=1` | OK |
| Main time step | `dt=0.02` | OK |
| Main final time | `tmax=300` | OK |
| Main pre-equilibration | `preeq_steps=500` | OK |
| Final potential | `a_f=0.02`, `b=c=1` | OK |
| Initial noise scale | `D_i=D0 T_i`, `D0=0.02` | OK |
| Main post-quench noise | `D_f=9e-3` | OK |
| Initial-state labels | `1.05, 1.10, 1.20, 1.50, 2.00, 3.00` | OK |
| Nucleation threshold | `C_th=20` | OK |
| Persistence condition | `n_cons=3` | OK |
| Main `phi6` sample count | `100 realizations per label` | OK |
| `phi4` sample count | `5 realizations per label` | OK |
| Threshold checks | `30 realizations per label` | OK |
| Time-step check | `30 realizations per label` | OK |
| `N=128` check | `20 realizations per label` | OK |

## Note

The `preeq=1000` check also has 30 realizations per label. It is discussed in the robustness text and used in the supplemental robustness figure.
