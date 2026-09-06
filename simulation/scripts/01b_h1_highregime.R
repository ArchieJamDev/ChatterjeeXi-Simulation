# H1 follow-up: adds a HIGH xi_0 regime, requested by external review,
# to actually test H1's claim of "reduccion particularmente marcada
# cuando xi sea alto" -- the original 01_h1.R only tested xi~0.493
# (moderate), never a high regime, so that part of H1 was unverified.
# Also adds Monte Carlo SE for the bias estimate at each n (point 1.5
# of the review: point estimates alone don't establish "indistinguible
# de cero").
#
# NOTE: run from the simulation/ directory, e.g.
#   cd simulation && Rscript scripts/01b_h1_highregime.R
source("R/dgp.R")
source("R/statistics.R")
source("R/calibration.R")

set.seed(20260906)

n_values <- c(20, 50, 100, 500, 2000, 2800, 3500)
R <- 1000
xi_target_high <- 0.85

mono <- check_monotonicity("lineal", "normal", seed = 1)
cal <- calibrate_sigma("lineal", "normal", xi_target = xi_target_high, mono_check = mono, seed = 1)
sigma <- cal$sigma
cat("Calibrated sigma (high regime) =", round(sigma, 4), "-> xi_hat (n=2e4) =", round(cal$xi_hat, 4), "\n")

xi_true_reps <- replicate(5, {
  xy <- generate_xy(2e5, "lineal", "normal", sigma)
  xicor(xy$X, xy$Y)
})
xi_true <- mean(xi_true_reps)
cat("High-precision population xi (n=2e5 x5):", round(xi_true, 5),
    "(sd across reps:", round(sd(xi_true_reps), 5), ")\n")

results <- data.frame(
  n = n_values,
  cota_teorica = (n_values - 2) / (n_values + 1),
  media_xi_hat = NA_real_,
  sesgo = NA_real_,
  se_mc_sesgo = NA_real_,
  mse = NA_real_
)

for (i in seq_along(n_values)) {
  n <- n_values[i]
  xi_hats <- vapply(seq_len(R), function(r) {
    xy <- generate_xy(n, "lineal", "normal", sigma)
    xicor(xy$X, xy$Y)
  }, numeric(1))
  results$media_xi_hat[i] <- mean(xi_hats)
  results$sesgo[i]        <- mean(xi_hats) - xi_true
  results$se_mc_sesgo[i]  <- sd(xi_hats) / sqrt(R) # MC standard error of the mean bias
  results$mse[i]          <- mean((xi_hats - xi_true)^2)
  cat("n =", n, "done.\n")
}

cat("\n=== H1 high regime: sesgo de muestra finita de xi_n vs. n (xi_true =", round(xi_true, 4), ") ===\n")
print(round(results, 5))

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(results, file.path(out_dir, "01b_h1_highregime.csv"), row.names = FALSE)
