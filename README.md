# Chatterjee's xi coefficient with psychological data: simulations and real-data analyses

Code and results behind the manuscript *When Does Chatterjee's xi Work with Psychological Data, and
When Does It Not? A Simulation and Real-Data Guide*. Everything reported in the manuscript can be
regenerated from this repository.

## What is here

| Folder | Content |
|---|---|
| `simulation/R/` | data-generating mechanisms, calibration of the dependence strength, helpers of each study |
| `simulation/scripts/` | one script per analysis (see below) |
| `simulation/results/` | raw results (CSV) of every analysis, one file per shard, with R/package versions and seeds |
| `paper_brm/` | the tables (`tables_en/`) and figures (`figures/*_en.pdf`) of the manuscript, generated from the results |
| `.github/workflows/` | the GitHub Actions workflows that ran the simulations and the real-data analyses |

## Scripts and manuscript studies

| Study | Scripts | Results |
|---|---|---|
| Study 1 (level and power of the xi_n test with discrete data) | `20_study1.R` | `results/study1/` |
| Study 1, tie-breaking noise (Table 4; earlier modular simulations) | `03_h2b.R`, `05_h2d.R`, `06_h2d_kvarying.R` | `results/03_h2b.csv`, `05_h2d.csv`, `06_h2d_kvarying.csv` |
| Study 2 (xi_n against Spearman, Pearson, dCor and tau*) | `21_study2.R` | `results/study2/` |
| Study 3 (every variable pair of `bfi`, `sat.act` and NHANES 2017-2018) | `30_realdata_screen.R`, `31_nhanes_screen.R`, `32_realdata_summary.R` | `results/realdata/`, `results/nhanes/`, `results/realdata_summary/` |
| Tables and figures | `40_make_tables.R` | `paper_brm/tables_en/`, `paper_brm/figures/` |

Scripts `00` to `12` are earlier, one-factor-at-a-time simulations of the same project; only `03`,
`05` and `06` are used in the manuscript. `simulation/README.md` documents every script, the seeds, the
software versions, how to run each block and how the Actions workflows are organized.

## Data

- `psych::bfi` and `psych::sat.act`: in the CRAN package `psych` (`data(bfi, package = "psych")`).
- NHANES 2017-2018 (US CDC/NCHS, public-use files): `DEMO_J`, `DPQ_J`, `SLQ_J`, `BMX_J`, downloaded
  without registration from `https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/`
  by `31_nhanes_screen.R` and `40_make_tables.R`; the size and md5 sum of each file are recorded in
  `results/nhanes/*_env.txt`.

## Reproducing

R (4.5 or later) with the packages `XICOR`, `energy`, `TauStar`, `psych` and `foreign` (the `energy`
package needs the system library `libgsl`). Run each script from the `simulation/` directory, e.g.
`Rscript scripts/20_study1.R typeI_asym 100`; each invocation is one shard of a block and writes
its CSV to `results/`. All scripts use fixed seeds, and two independent runs of the main blocks of
Study 2 produced byte-identical files.
