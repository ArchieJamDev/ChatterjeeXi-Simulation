# Module: H4 -- potencia segun geometria de la dependencia (manuscript
# Table "tab:h4")
#
# Reference cell: n de referencia = 100, forma = normal, tipo de
# variable = continua. Each structure is calibrated SEPARATELY to a
# shared xi_target = 0.10, chosen because it is the largest value
# achievable by ALL FIVE structures -- heterocedastica's own ceiling
# is only ~0.124 (at lambda ~ 0; recall lambda, not sigma, is its
# true dependence-strength knob, Sec 3.2/dgp.R), and circular's is
# ~0.23, both far below the ~0.97 ceiling of lineal/cuadratica/senoidal.
#
# Type I error is computed ONCE per statistic under exact independence
# (X _|_ Y, structure-agnostic by definition -- there is no
# "structure" without dependence) and reported identically across all
# five columns, per Section 3.3's requirement to verify comparable
# Type I error before comparing power.
#
# NOTE: run from the simulation/ directory, e.g.
#   cd simulation && Rscript scripts/09_h4.R
source("R/dgp.R")
source("R/statistics.R")
source("R/calibration.R")

set.seed(20260906)

n_ref <- 100
R <- 300
alpha <- 0.05
xi_target <- 0.10
structures <- c("lineal", "cuadratica", "senoidal", "circular", "heterocedastica")

## 1. Type I error (structure-agnostic, exact independence)
rej_xi <- rej_dcor <- rej_taustar <- logical(R)
for (r in seq_len(R)) {
  X <- rnorm(n_ref); Y <- rnorm(n_ref) # exact independence
  rej_xi[r]      <- xicor(X, Y, pvalue = TRUE)$pval < alpha
  rej_dcor[r]    <- dcor.test(X, Y, R = 99)$p.value < alpha
  rej_taustar[r] <- as.numeric(tauStarTest(X, Y, resamples = 99)$pVal) < alpha
}
type1 <- c(xi_n = mean(rej_xi), dcor = mean(rej_dcor), taustar = mean(rej_taustar))
cat("Type I error (shared across structures):\n"); print(round(type1, 4))

## 2. Power per structure, each calibrated to xi_target = 0.10
sigma_ranges <- list(
  lineal = c(1e-3, 20), cuadratica = c(1e-3, 20), senoidal = c(1e-3, 20),
  circular = c(1e-3, 20), heterocedastica = c(1e-4, 0.999) # lambda in [0,1)
)

power_results <- data.frame(estructura = structures, sigma_cal = NA_real_,
                             xi_hat_cal = NA_real_, xi_n = NA_real_, dcor = NA_real_, taustar = NA_real_)

for (i in seq_along(structures)) {
  st <- structures[i]
  mono <- check_monotonicity(st, "normal", sigma_grid = exp(seq(log(sigma_ranges[[st]][1]),
                                                                 log(sigma_ranges[[st]][2]), length.out = 15)), seed = 1)
  cal <- calibrate_sigma(st, "normal", xi_target = xi_target, mono_check = mono,
                          sigma_range = sigma_ranges[[st]], seed = 1)
  sigma <- cal$sigma
  power_results$sigma_cal[i] <- sigma
  power_results$xi_hat_cal[i] <- cal$xi_hat

  rej_xi <- rej_dcor <- rej_taustar <- logical(R)
  for (r in seq_len(R)) {
    xy <- generate_xy(n_ref, st, "normal", sigma)
    rej_xi[r]      <- xicor(xy$X, xy$Y, pvalue = TRUE)$pval < alpha
    rej_dcor[r]    <- dcor.test(xy$X, xy$Y, R = 99)$p.value < alpha
    rej_taustar[r] <- as.numeric(tauStarTest(xy$X, xy$Y, resamples = 99)$pVal) < alpha
  }
  power_results$xi_n[i] <- mean(rej_xi)
  power_results$dcor[i] <- mean(rej_dcor)
  power_results$taustar[i] <- mean(rej_taustar)
  cat(st, "done. sigma_cal =", round(sigma, 4), " xi_hat_cal =", round(cal$xi_hat, 4), "\n")
}

cat("\n=== H4: power by structure (n =", n_ref, ", xi_target =", xi_target, ") ===\n")
print(round(power_results[, -1], 4))

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(data.frame(t(type1)), file.path(out_dir, "09_h4_type1.csv"), row.names = FALSE)
write.csv(power_results, file.path(out_dir, "09_h4_power.csv"), row.names = FALSE)
