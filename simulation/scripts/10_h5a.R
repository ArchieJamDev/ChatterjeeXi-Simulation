# Module: H5a -- desempeno practico a dependencia fija (manuscript
# Table "tab:h5a")
#
# Reference cell: forma = normal, tipo de variable = continua,
# estructura = lineal. Crosses n x xi_0 (fixed), including an exact
# independence condition (xi_0 = 0, generated directly -- NOT via a
# large sigma, per Section 3.2's precaution). Uses xi_n's own fast
# asymptotic test throughout (no permutation), since this table is
# about xi_n alone, not a cross-statistic comparison (that's H4).
#
# NOTE: run from the simulation/ directory, e.g.
#   cd simulation && Rscript scripts/10_h5a.R
source("R/dgp.R")
source("R/statistics.R")
source("R/calibration.R")

set.seed(20260906)

n_values <- c(20, 50, 100, 500, 2000)
xi_targets <- c(0.1, 0.3, 0.5, 0.8)
R <- 1000
alpha <- 0.05

## Calibrate sigma for each xi_0 target (structure=lineal, family=normal)
mono <- check_monotonicity("lineal", "normal", seed = 1)
sigmas <- sapply(xi_targets, function(xi0) {
  calibrate_sigma("lineal", "normal", xi_target = xi0, mono_check = mono, seed = 1)$sigma
})
names(sigmas) <- as.character(xi_targets)
cat("Calibrated sigmas:\n"); print(round(sigmas, 4))

## Power grid: rows = n, columns = xi_0 (including exact independence)
col_names <- c("xi0_0_exact", paste0("xi0_", xi_targets))
power_grid <- matrix(NA_real_, nrow = length(n_values), ncol = length(col_names),
                      dimnames = list(paste0("n=", n_values), col_names))

for (i in seq_along(n_values)) {
  n <- n_values[i]

  # xi_0 = 0, exact independence (X and Y generated independently, not via large sigma)
  rej <- vapply(seq_len(R), function(r) {
    X <- rnorm(n); Y <- rnorm(n)
    xicor(X, Y, pvalue = TRUE)$pval < alpha
  }, logical(1))
  power_grid[i, 1] <- mean(rej)

  # xi_0 > 0, calibrated sigma
  for (j in seq_along(xi_targets)) {
    sigma <- sigmas[j]
    rej <- vapply(seq_len(R), function(r) {
      xy <- generate_xy(n, "lineal", "normal", sigma)
      xicor(xy$X, xy$Y, pvalue = TRUE)$pval < alpha
    }, logical(1))
    power_grid[i, j + 1] <- mean(rej)
  }
  cat("n =", n, "done.\n")
}

cat("\n=== H5a: power grid (rows = n, cols = xi_0) ===\n")
print(round(power_grid, 3))

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(power_grid, file.path(out_dir, "10_h5a.csv"))
