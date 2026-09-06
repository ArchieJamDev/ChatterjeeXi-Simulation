# Module: H5b -- alternativas locales (manuscript Table "tab:h5b")
#
# Following Auddy, Deb & Nandy (2024): the population xi is set to
# shrink with n at a specified rate, xi(n) = K * (n/n0)^(-gamma), with
# gamma = 1/4 (their xi_n detection boundary for testing independence)
# or gamma = 1/2 (the faster-shrinking, parametrically-optimal rate
# that xi_n is NOT guaranteed to detect). Both sequences are anchored
# to the SAME starting value xi(n0=100) = 0.15, so they are directly
# comparable -- only the decay rate as n grows differs.
#
# Expected qualitative pattern: power under gamma=1/4 should stay
# roughly stable/non-trivial as n grows (this is xi_n's own detection
# boundary); power under gamma=1/2 should visibly decay toward the
# nominal alpha as n grows (signal shrinks faster than xi_n's test can
# keep up with).
#
# NOTE: run from the simulation/ directory, e.g.
#   cd simulation && Rscript scripts/11_h5b.R
source("R/dgp.R")
source("R/statistics.R")
source("R/calibration.R")

set.seed(20260906)

n0 <- 100
xi0 <- 0.15
n_values <- c(100, 500, 2000)
gammas <- c("1/4" = 1/4, "1/2" = 1/2)
structures <- c("lineal", "senoidal")
R <- 1000
alpha <- 0.05

xi_target_fun <- function(n, gamma) xi0 * (n / n0)^(-gamma)

results <- expand.grid(estructura = structures, n = n_values, tasa = names(gammas),
                        stringsAsFactors = FALSE)
results$xi_target <- mapply(function(n, g) xi_target_fun(n, gammas[g]), results$n, results$tasa)
results$sigma_cal <- NA_real_
results$potencia <- NA_real_

mono_cache <- list()

for (i in seq_len(nrow(results))) {
  st <- results$estructura[i]; n <- results$n[i]; xi_t <- results$xi_target[i]
  if (is.null(mono_cache[[st]])) mono_cache[[st]] <- check_monotonicity(st, "normal", seed = 1)
  cal <- calibrate_sigma(st, "normal", xi_target = xi_t, mono_check = mono_cache[[st]], seed = 1)
  results$sigma_cal[i] <- cal$sigma

  rej <- vapply(seq_len(R), function(r) {
    xy <- generate_xy(n, st, "normal", cal$sigma)
    xicor(xy$X, xy$Y, pvalue = TRUE)$pval < alpha
  }, logical(1))
  results$potencia[i] <- mean(rej)
  cat(st, "n=", n, "tasa=", results$tasa[i], "xi_target=", round(xi_t, 4),
      "sigma=", round(cal$sigma, 4), "potencia=", round(results$potencia[i], 3), "\n")
}

cat("\n=== H5b: power under local alternatives (anchored xi(n0=100)=0.15) ===\n")
print(results, row.names = FALSE)

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(results, file.path(out_dir, "11_h5b.csv"), row.names = FALSE)
