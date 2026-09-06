# Module: H3 -- estabilidad ante la forma de la distribucion (manuscript
# Table "tab:h3")
#
# Reference cell: n de referencia = 100, estructura = lineal, regimen =
# moderado. Each marginal family is calibrated SEPARATELY (sigma
# depends on family, Section 3.2) to hit xi_target ~ 0.5 -- xi_n's own
# population reference. dCor's population value under that SAME sigma
# is whatever it turns out to be (Section 3.2: "sin referencia
# poblacional compartida entre estadisticos"), not forced to match.
#
# NOTE: run from the simulation/ directory, e.g.
#   cd simulation && Rscript scripts/07_h3.R
source("R/dgp.R")
source("R/statistics.R")
source("R/calibration.R")

set.seed(20260906)

n_ref <- 100
R <- 1000
xi_target <- 0.5
families <- c("normal", "uniform", "exponential", "t", "mixture")

results <- data.frame(
  familia = families,
  sigma = NA_real_,
  xi_true_calibrated = NA_real_,
  xi_n_media = NA_real_,
  dcor_media = NA_real_
)

for (i in seq_along(families)) {
  fam <- families[i]
  mono <- check_monotonicity("lineal", fam, seed = 1)
  cal <- calibrate_sigma("lineal", fam, xi_target = xi_target, mono_check = mono, seed = 1)
  sigma <- cal$sigma

  xi_hats <- dcor_hats <- numeric(R)
  for (r in seq_len(R)) {
    xy <- generate_xy(n_ref, "lineal", fam, sigma)
    xi_hats[r] <- xicor(xy$X, xy$Y)
    dcor_hats[r] <- dcor(xy$X, xy$Y)
  }

  results$sigma[i] <- sigma
  results$xi_true_calibrated[i] <- cal$xi_hat
  results$xi_n_media[i] <- mean(xi_hats)
  results$dcor_media[i] <- mean(dcor_hats)
  cat(fam, "done. sigma =", round(sigma, 4), "\n")
}

cat("\n=== H3: xi_n vs dCor across marginal families (n =", n_ref, ", target xi = 0.5) ===\n")
print(round(results[, -1], 4))
cat("\nRange of xi_n means across families:", round(diff(range(results$xi_n_media)), 4), "\n")
cat("Range of dCor means across families:", round(diff(range(results$dcor_media)), 4), "\n")

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(results, file.path(out_dir, "07_h3.csv"), row.names = FALSE)
