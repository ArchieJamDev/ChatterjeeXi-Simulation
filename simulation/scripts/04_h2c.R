# Module: H2c -- discretizacion conjunta (manuscript Table "tab:h2c")
#
# Reference cell: n de referencia = 100, forma = normal, estructura =
# lineal, regimen = moderado (reuses 01_h1's calibrated sigma/xi_true).
# Both X and Y are discretized to the same number of categories.
#
# Since X has ties here too, the same nested (dataset x B tie-breaks)
# design as H2b is used internally for correctness -- but H2c's table
# reports a single pooled total variance (not the between/within
# split, which is H2b's own result), matching the simpler two-row
# table design.
#
# NOTE: run from the simulation/ directory, e.g.
#   cd simulation && Rscript scripts/04_h2c.R
source("R/dgp.R")
source("R/discretize.R")
source("R/statistics.R")

set.seed(20260906)

cal <- readRDS("results/01_h1_calibration.rds")
sigma <- cal$sigma
n_ref <- 100
R_outer <- 300
B <- 20
alpha <- 0.05

conditions <- list(
  Continua = NA,
  Binaria  = 2,
  `Likert-3` = 3,
  `Likert-5` = 5,
  `Likert-7` = 7
)

results <- data.frame(
  condicion = names(conditions),
  varianza_empirica = NA_real_,
  potencia = NA_real_
)

for (i in seq_along(conditions)) {
  c_cat <- conditions[[i]]
  xi_mat <- matrix(NA_real_, nrow = R_outer, ncol = B)
  rej_mat <- matrix(NA, nrow = R_outer, ncol = B)

  for (r in seq_len(R_outer)) {
    xy <- generate_xy(n_ref, "lineal", "normal", sigma)
    Xobs <- if (is.na(c_cat)) xy$X else discretize(xy$X, c_cat)
    Yobs <- if (is.na(c_cat)) xy$Y else discretize(xy$Y, c_cat)
    for (b in seq_len(B)) {
      fit <- xicor(Xobs, Yobs, pvalue = TRUE)
      xi_mat[r, b] <- fit$xi
      rej_mat[r, b] <- fit$pval < alpha
    }
  }

  results$varianza_empirica[i] <- var(as.vector(xi_mat))
  results$potencia[i] <- mean(rej_mat)
  cat(names(conditions)[i], "done.\n")
}

cat("\n=== H2c: discretizacion conjunta (n =", n_ref, ", sigma =", round(sigma, 4),
    ", xi_true =", round(cal$xi_true, 4), ") ===\n")
print(results, row.names = FALSE)

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(results, file.path(out_dir, "04_h2c.csv"), row.names = FALSE)
