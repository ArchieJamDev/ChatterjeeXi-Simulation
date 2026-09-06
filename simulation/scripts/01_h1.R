# Module: H1 -- efecto del tamano de muestra (manuscript Table "tab:h1")
#
# Reference cell: forma = normal, tipo de variable = continua,
# estructura = lineal, regimen = moderado (xi_0 ~ 0.5).
#
# NOTE: run from the simulation/ directory, e.g.
#   cd simulation && Rscript scripts/01_h1.R
source("R/dgp.R")
source("R/statistics.R")
source("R/calibration.R")

set.seed(20260906)

## 1. Calibrate sigma for structure=lineal, family=normal, xi_0 ~ 0.5
mono <- check_monotonicity("lineal", "normal", seed = 1)
cal <- calibrate_sigma("lineal", "normal", xi_target = 0.5, mono_check = mono, seed = 1)
sigma <- cal$sigma
cat("Calibrated sigma =", round(sigma, 4), "-> xi_hat (n=2e4) =", round(cal$xi_hat, 4), "\n")

## 2. Population xi under this exact sigma, high-precision (n=1e6, averaged over 5 reps)
xi_true_reps <- replicate(5, {
  xy <- generate_xy(2e5, "lineal", "normal", sigma)
  xicor(xy$X, xy$Y)
})
xi_true <- mean(xi_true_reps)
cat("High-precision population xi (n=2e5 x5):", round(xi_true, 5),
    "(sd across reps:", round(sd(xi_true_reps), 5), ")\n")

## 3. Monte Carlo over n, R replications each
n_values <- c(20, 50, 100, 500, 2000)
R <- 1000

results <- data.frame(
  n = n_values,
  cota_teorica = (n_values - 2) / (n_values + 1),
  media_xi_hat = NA_real_,
  sesgo = NA_real_,
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
  results$mse[i]           <- mean((xi_hats - xi_true)^2)
  cat("n =", n, "done.\n")
}

cat("\n=== H1: sesgo de muestra finita de xi_n vs. n (xi_true =", round(xi_true, 4), ") ===\n")
print(round(results, 4))

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(results, file.path(out_dir, "01_h1.csv"), row.names = FALSE)
saveRDS(list(sigma = sigma, xi_true = xi_true), file.path(out_dir, "01_h1_calibration.rds"))
