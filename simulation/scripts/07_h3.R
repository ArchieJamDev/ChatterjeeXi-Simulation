# Module: H3 -- estabilidad ante la forma de la distribucion (manuscript
# Table "tab:h3"): bias and variance across marginal families.
#
# Reference cell: n de referencia = 100, estructura = lineal, regimen =
# moderado. Each marginal family is calibrated SEPARATELY (sigma
# depends on family, Section 3.2) to hit xi_target ~ 0.5 -- xi_n's own
# population reference. dCor's population value under that SAME sigma
# is whatever it turns out to be (Section 3.2), not forced to match.
#
# Type I error / power (per external review, added separately in
# 07b_h3_typeI_power.R -- one family per lightweight invocation, to
# keep memory use low on this machine) are NOT computed here.
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
  dcor_true_calibrated = NA_real_,
  xi_n_media = NA_real_, xi_n_sesgo = NA_real_, xi_n_var = NA_real_,
  dcor_media = NA_real_, dcor_sesgo = NA_real_, dcor_var = NA_real_
)

for (i in seq_along(families)) {
  fam <- families[i]
  mono <- check_monotonicity("lineal", fam, seed = 1)
  cal <- calibrate_sigma("lineal", fam, xi_target = xi_target, mono_check = mono, seed = 1)
  sigma <- cal$sigma
  results$sigma[i] <- sigma
  results$xi_true_calibrated[i] <- cal$xi_hat

  ## NOTE: dcor() materializes an n x n distance matrix (O(n^2) memory);
  ## n=2000 is plenty for a stable population estimate (~32MB) and safe
  ## on this machine's limited RAM -- n=2e4 here previously caused an
  ## OOM kill (~3.2GB for that one matrix alone).
  xy_pop <- generate_xy(2000, "lineal", fam, sigma)
  results$dcor_true_calibrated[i] <- dcor(xy_pop$X, xy_pop$Y)
  rm(xy_pop); gc(verbose = FALSE)

  xi_hats <- dcor_hats <- numeric(R)
  for (r in seq_len(R)) {
    xy <- generate_xy(n_ref, "lineal", fam, sigma)
    xi_hats[r] <- xicor(xy$X, xy$Y)
    dcor_hats[r] <- dcor(xy$X, xy$Y)
  }
  results$xi_n_media[i] <- mean(xi_hats)
  results$xi_n_sesgo[i] <- mean(xi_hats) - results$xi_true_calibrated[i]
  results$xi_n_var[i]   <- var(xi_hats)
  results$dcor_media[i] <- mean(dcor_hats)
  results$dcor_sesgo[i] <- mean(dcor_hats) - results$dcor_true_calibrated[i]
  results$dcor_var[i]   <- var(dcor_hats)

  rm(xi_hats, dcor_hats); gc(verbose = FALSE)
  cat(fam, "done. sigma =", round(sigma, 4), "\n")
}

cat("\n=== H3: bias and variance across marginal families ===\n")
print(round(results[, -1], 5))

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(results, file.path(out_dir, "07_h3.csv"), row.names = FALSE)
