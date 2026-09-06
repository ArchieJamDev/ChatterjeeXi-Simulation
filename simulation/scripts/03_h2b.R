# Module: H2b -- empates en el predictor, diseno anidado (manuscript
# Table "tab:h2b" and Section 3.2 "Discretizacion, desempate anidado...")
#
# Reference cell: n de referencia = 100, forma = normal, estructura =
# lineal, regimen = moderado (reuses 01_h1's calibrated sigma/xi_true).
#
# Nested design: (i) generate a latent dataset, (ii) discretize X
# (Y stays continuous), (iii) compute xi_n under B independent random
# tie-breaks of that SAME discretized dataset, (iv) decompose total
# variance into between-dataset (ordinary sampling variability) and
# within-dataset (variability induced by tie-breaking alone) -- this
# is the decomposition planned in Section 3.2/3.3 as its own result.
# Empirically confirmed beforehand: xicor() re-randomizes the tie-
# break on every call given a fixed tied dataset (var > 0 across
# repeated calls), which is what makes this nested design meaningful.
#
# NOTE: run from the simulation/ directory, e.g.
#   cd simulation && Rscript scripts/03_h2b.R
source("R/dgp.R")
source("R/discretize.R")
source("R/statistics.R")

set.seed(20260906)

cal <- readRDS("results/01_h1_calibration.rds")
sigma <- cal$sigma
n_ref <- 100
R_outer <- 300 # independent datasets
B <- 20        # tie-break replicates per dataset (R_outer*B = 1000 calls/condition)
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
  var_entre = NA_real_,
  var_dentro = NA_real_,
  potencia = NA_real_
)

for (i in seq_along(conditions)) {
  c_cat <- conditions[[i]]
  xi_mat <- matrix(NA_real_, nrow = R_outer, ncol = B)
  rej_mat <- matrix(NA, nrow = R_outer, ncol = B)

  for (r in seq_len(R_outer)) {
    xy <- generate_xy(n_ref, "lineal", "normal", sigma)
    Xobs <- if (is.na(c_cat)) xy$X else discretize(xy$X, c_cat)
    for (b in seq_len(B)) {
      fit <- xicor(Xobs, xy$Y, pvalue = TRUE)
      xi_mat[r, b] <- fit$xi
      rej_mat[r, b] <- fit$pval < alpha
    }
  }

  dataset_means <- rowMeans(xi_mat)
  dataset_vars  <- apply(xi_mat, 1, var) # 0 for "Continua" (no ties to re-break)

  results$var_entre[i]  <- var(dataset_means)
  results$var_dentro[i] <- mean(dataset_vars)
  results$potencia[i]   <- mean(rej_mat)
  cat(names(conditions)[i], "done.\n")
}

cat("\n=== H2b: empates en el predictor, descomposicion de varianza (n =", n_ref,
    ", sigma =", round(sigma, 4), ", xi_true =", round(cal$xi_true, 4), ") ===\n")
print(results, row.names = FALSE)

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(results, file.path(out_dir, "03_h2b.csv"), row.names = FALSE)
