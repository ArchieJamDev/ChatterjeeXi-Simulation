# Module: Validacion (manuscript Section 4.1 / Table "Validacion")
#
# Gate: reproduce the closed-form xi(X,Y) for the bivariate normal
# case from Ansari & Fuchs (2026), Eq. (33) / Example 5.5, before
# proceeding to interpret H3's non-elliptical cells. If this doesn't
# hold, there is a bug in dgp.R / statistics.R, not a substantive
# finding about xi.
#
# Ansari & Fuchs (2026), Eq. 33 (their (Y,X) notation matches our
# xi_n(X,Y) = "how well Y is determined by X" once X,Y are relabeled
# consistently -- verified: at rho=1, formula gives xi=1; at rho=0,
# xi=0, as required):
#   xi(X,Y) = (3/pi) * asin((1+rho^2)/2) - 1/2
# for standard bivariate normal (X,Y) with correlation rho.

# NOTE: run from the simulation/ directory, e.g.
#   cd simulation && Rscript scripts/00_validation.R
source("R/dgp.R")
source("R/statistics.R")
source("R/calibration.R")

xi_closed_form_normal <- function(rho) (3 / pi) * asin((1 + rho^2) / 2) - 1 / 2

set.seed(20260906)

n_pop <- 1e5
sigmas <- c(3, 2, 1.5, 1, 0.7, 0.5, 0.3, 0.15, 0.05)
rhos <- 1 / sqrt(1 + sigmas^2) # Corr(X,Y) for Y = X + sigma*E, X,E ~ N(0,1) indep.

results <- data.frame(
  sigma = sigmas,
  rho = rhos,
  xi_closed_form = xi_closed_form_normal(rhos),
  xi_mc_estimate = NA_real_
)

for (i in seq_along(sigmas)) {
  xy <- generate_xy(n_pop, "lineal", "normal", sigmas[i])
  results$xi_mc_estimate[i] <- xicor(xy$X, xy$Y)
}

results$abs_diff <- abs(results$xi_closed_form - results$xi_mc_estimate)

cat("=== Validation: xi_n vs. Ansari & Fuchs (2026) closed form, bivariate normal ===\n")
print(round(results, 4))
cat("\nMax abs. difference (closed form vs. n=", n_pop, " Monte Carlo):",
    round(max(results$abs_diff), 4), "\n")
cat("Gate", ifelse(max(results$abs_diff) < 0.01, "PASSED", "FAILED -- do not proceed to H3"), "\n")

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(results, file.path(out_dir, "00_validation_normal.csv"), row.names = FALSE)
