# Reproducibility package — simulation

R code and raw results for the Monte Carlo study of Chatterjee's xi
coefficient. Every script has a fixed seed, writes its results as CSV under
`results/`, and is run from this directory (`cd simulation`).

## Software

- Scripts 00-12: R 4.5.2; packages `XICOR` 0.4.1, `energy` 1.7.12,
  `TauStar` 1.1.9, `psych` 2.6.5 (versions currently used).
- Study 1: produced on GitHub Actions (`ubuntu-latest`, R 4.6.1, `XICOR` 0.4.1);
  the exact versions and the seed of each shard are in its `*_env.txt`.
- `dcor()` from `energy` builds an n x n distance matrix (O(n^2) memory), so the
  permutation-based dCor tests are run with small R on small machines.

## Layout

- `R/dgp.R` data-generating mechanisms (five latent families, five dependence
  structures); `R/calibration.R` calibration of the noise level `sigma` to a
  target population xi (monotonicity is checked, not assumed);
  `R/discretize.R` discretization and composite scores for the H2 modules;
  `R/statistics.R` loads `XICOR`, `energy`, `TauStar`; `R/study1_lib.R`,
  `R/study2_lib.R` helpers of Studies 1 and 2.
- `scripts/` one script per module; `results/` one CSV per module.

## Scripts and manuscript tables

Table numbers refer to `paper/paper.tex` as originally submitted.

| Script | Content | Table | Results file |
|---|---|---|---|
| `00_validation.R` | exact xi for the equicorrelated normal vs. Monte Carlo | 3 | `00_validation_normal.csv` |
| `01_h1.R`, `01b_h1_highregime.R` | finite-sample bias and MSE vs. n, moderate and high dependence | 4 | `01_h1.csv`, `01b_h1_highregime.csv` |
| `02_h2a.R` | ties in the response | 5 | `02_h2a.csv` |
| `03_h2b.R` | ties in the predictor, nested tie-breaking variance decomposition | 6 | `03_h2b.csv` |
| `04_h2c.R` | joint discretization | 7 | `04_h2c.csv` |
| `05_h2d.R`, `06_h2d_kvarying.R` | composite score of k items | 8, 9 | `05_h2d.csv`, `06_h2d_kvarying.csv` |
| `07_h3.R`, `07b_h3_typeI_power.R` | marginal families: bias, variance, type I error, power | 10 | `07_h3.csv`, `07b_h3_typeI_power.csv` |
| `08_h3_invariance.R` | invariance under monotone transformations | 11 | `08_h3_invariance.csv` |
| `09_h4.R` | power by dependence geometry vs. dCor and tau* | 12 | `09_h4_type1.csv`, `09_h4_power.csv` |
| `10_h5a.R` | power at fixed dependence | 13 | `10_h5a.csv` |
| `11_h5b.R` | local alternatives | 14 | `11_h5b.csv` |
| `12_realdata_bfi.R` | real-data illustration, `psych::bfi` | 15 | `12_realdata_bfi.csv` |

`07b_h3_typeI_power.R` takes one family per call (memory), e.g.
`for f in normal uniform exponential t mixture; do Rscript scripts/07b_h3_typeI_power.R $f; done`.

## Seeds

Every script calls `set.seed(20260906)` once at the top, except
`12_realdata_bfi.R`, whose random tie-breaks are not seeded (its xi values
reproduce only up to Monte Carlo noise in the third decimal). Study 1 uses one
seed per block and shard: `20261100 + 10000 * block_index + shard_index`
(`block_index` follows the order typeI_asym, typeI_perm, power_asym,
power_perm; `shard_index` is n, or 1/2/3 for lineal/cuadratica/senoidal).
Population cut points are drawn from a separate fixed stream (`set.seed(1)`).

## Study 1: type I error and power with discrete margins

`scripts/20_study1.R <block> <shard>` runs one block and shard and writes
`results/study1/<block>_<shard>.csv` (plus an `_env.txt` with R and package
versions and the seed).

| Block | Shard | Computes |
|---|---|---|
| `typeI_asym` | n in 50, 100, 250, 500, 1000 | type I error of the tie-aware asymptotic test (`XICOR`, `ties = TRUE`) and of the naive test that assumes a continuous response, under exact independence; 3 X types x 8 Y types; R = 5000 |
| `typeI_perm` | n in 100, 250, 500 | same datasets plus a permutation test (B = 499); X in {continuous, binary 50/50} x Y in {continuous, binary 90/10, binary 95/5, 7-category with floor effect}; R = 1000 |
| `power_asym` | structure (lineal, cuadratica, senoidal) | asymptotic power when dependence is calibrated on the latent variable to xi = 0.15; 3 X types x 8 Y types x 5 n; R = 2000 |
| `power_perm` | n in 100, 250 | asymptotic and permutation power (B = 199); 3 structures x 2 X types x 4 Y types; R = 500 |

Margin types. X: continuous; binary 50/50; 5-category skewed
(.05/.10/.15/.30/.40). Y: continuous; binary 50/50, 80/20, 90/10, 95/5;
5-category skewed; 7-category with floor effect
(.35/.25/.15/.10/.07/.05/.03); composite of 5 six-category items. Categories
are obtained by cutting the latent variable at population cut points, so
category counts are random from sample to sample, as in real data.

Quick check of the pipeline: `STUDY1_SMOKE=1 Rscript scripts/20_study1.R typeI_asym 100`
(R = 20, writes to a temporary directory). `STUDY1_R=<R>` overrides the number of
replications and `STUDY1_OUT=<dir>` the output directory.

### Running on GitHub Actions

`.github/workflows/study1.yml` runs the 13 shards as parallel jobs. A push to
branch `brm-revision` runs a smoke version (R = 20); pushing a tag named
`study1-v*` runs the full version. Each shard uploads its CSV as the artifact
`study1-<mode>-<block>-<shard>`; download them with `gh run download <run-id>`
and store the files of a full run in `results/study1/`. The files currently in
`results/study1/` come from run 40297517916 (tag `study1-v1`, commit `91f012d`).

## Study 2: xi_n against competing measures

`scripts/21_study2.R <block> <shard>` (workflow `.github/workflows/study2.yml`,
tag `study2-v*` for the full run) compares the tie-aware asymptotic test of xi_n
with Spearman and Pearson (two-sided t tests), distance correlation
(`energy::dcor.test`, permutation test, B = 199) and tau* (`TauStar::tauStarTest`,
mode chosen automatically from the data), always on the same simulated datasets.
X is continuous; Y is continuous, binary 50/50, binary 90/10, 5-category skewed or
a composite of 5 six-category items (same definitions as Study 1).

**Neutral calibration.** Dependence is not matched on any coefficient's scale but on
the signal share `s = Var(noiseless structure) / Var(latent Y)`, i.e. the R^2 of the
true generating function: `sigma^2 = Var(f) (1 - s) / s`, with `Var(f)` from a large
noiseless sample (`R/study2_lib.R`). The same `s` therefore means the same amount of
signal for the linear, quadratic, sinusoidal and circular structures.

| Block | Shard | Computes |
|---|---|---|
| `typeI` | n in 50, 100, 200, 400 | type I error of all five tests under independence, 5 Y types; R = 1000 |
| `power_cheap` | structure (lineal, cuadratica, senoidal, circular) | power of xi, Spearman, Pearson; s in .05, .13, .26, .50; n in 50, 100, 200, 400, 800; 5 Y types; R = 2000 |
| `power_full` | structure | power of all five tests; s in .13, .26 (Cohen's medium and large R^2); n in 50, 100, 200, 400; 4 Y types; R = 500 |
| `power_circ_hi` | circular | as `power_full` for the circular structure only, at s in .70, .85, .95: for a circular relationship no test has power at s <= .50, so that range is uninformative |

`warn_rate_tau` in the output is the share of replications in which
`tauStarTest` raised a numerical warning (counted, not hidden); `na_rate_tau` is
the share in which it failed. Seeds: `20262100 + 10000 * block_index + shard_index`
(+ `round(1000 * s)` inside the power blocks). Environment variables:
`STUDY2_SMOKE=1`, `STUDY2_R`, `STUDY2_OUT`, as in Study 1.

The files in `results/study2/` come from run 40325393994 (tag `study2-v2`, commit
`9d72d30`). An earlier full run of the first twelve shards (tag `study2-v1`) produced
byte-identical CSV files, which confirms that the fixed seeds make the Actions runs
exactly reproducible.

## Real-data screening

`scripts/30_realdata_screen.R <shard> <n_shards>` (workflow
`.github/workflows/realdata.yml`, tag `realdata-v*` for the full run, 12 shards)
analyses **every** pair of variables of two public datasets distributed with the
CRAN package `psych` (`data(bfi, package = "psych")`, `data(sat.act, package = "psych")`),
so that no pair is chosen by its result: the 25 `bfi` items, the 5 Big Five composite
scores (reverse-keyed item means), gender, education and age (503 pairs after removing
item-own-composite pairs) and the 6 variables of `sat.act` (15 pairs). Each pair uses its
complete cases and its own seed (`20263000 + pair index`), so results do not depend on how
the pairs are split into shards.

For every pair: Spearman's rho with a Fisher-type CI and test, Pearson's r with CI and
test, xi(X,Y) and xi(Y,X) averaged over 20 random tie-breaks with their SD, the asymptotic
tie-aware p-value of each direction (median over tie-breaks) and the share of tie-breaks
with p < .05, distance correlation with a permutation test (B = 499), and tau* with its
test (`tauStarTest`). No confidence interval is given for xi: the bootstrap fails for xi_n
and no meaningful interval exists in general (Dalitz et al., 2024). The decision rule used in
the summary (Benjamini-Hochberg-adjusted p < .05 within each measure, and unadjusted) is
fixed in the header of the script. Environment variables: `REALDATA_SMOKE=1`,
`REALDATA_OUT`, `REALDATA_MAXPAIRS`.

