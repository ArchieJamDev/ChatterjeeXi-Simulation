# Module: H3 sub-experimento aislado -- invariancia ante transformaciones
# monotonas (manuscript Table "tab:h3invar", Section 3.2 "Invariancia
# ante transformaciones monotonas")
#
# Unlike the main H3 module (which compares DIFFERENT generative
# families, each its own additive-noise model), this isolates the
# invariance property itself: ONE base joint distribution (X,Y), with
# strictly monotone (nonlinear) transforms applied to X only, keeping
# the underlying dependence fixed. xi_n, Spearman, and Kendall are
# rank-based and should be numerically IDENTICAL (up to Monte Carlo /
# tie-break noise) across transforms, by construction -- this is not a
# hypothesis to test, it follows directly from how those statistics are
# defined (sorting/ranking is unaffected by a monotone transform of the
# variable being ranked). dCor is expected to drift because it uses
# raw distances.
#
# NOTE: run from the simulation/ directory, e.g.
#   cd simulation && Rscript scripts/08_h3_invariance.R
source("R/dgp.R")
source("R/statistics.R")

set.seed(20260906)

cal <- readRDS("results/01_h1_calibration.rds")
sigma <- cal$sigma
n_ref <- 200
R <- 300

transforms <- list(
  `Sin transformar` = function(x) x,
  `Transf. monótona 1 (cubo)` = function(x) x^3,
  `Transf. monótona 2 (exp)` = function(x) exp(x)
)

results <- data.frame(
  transformacion = names(transforms),
  xi_n = NA_real_, spearman = NA_real_, kendall = NA_real_, dcor = NA_real_
)

xi_mat <- spearman_mat <- kendall_mat <- dcor_mat <- matrix(NA_real_, nrow = R, ncol = length(transforms))

for (r in seq_len(R)) {
  xy <- generate_xy(n_ref, "lineal", "normal", sigma)
  for (j in seq_along(transforms)) {
    Xt <- transforms[[j]](xy$X)
    xi_mat[r, j] <- xicor(Xt, xy$Y)
    spearman_mat[r, j] <- cor(Xt, xy$Y, method = "spearman")
    kendall_mat[r, j] <- cor(Xt, xy$Y, method = "kendall")
    dcor_mat[r, j] <- dcor(Xt, xy$Y)
  }
}

results$xi_n <- colMeans(xi_mat)
results$spearman <- colMeans(spearman_mat)
results$kendall <- colMeans(kendall_mat)
results$dcor <- colMeans(dcor_mat)

cat("=== H3 sub-experiment: invariance to monotone transforms of X (n =", n_ref, ") ===\n")
print(t(round(results[, -1], 4)))
cat("\nRange across transforms -- xi_n:", round(diff(range(results$xi_n)), 5),
    " Spearman:", round(diff(range(results$spearman)), 5),
    " Kendall:", round(diff(range(results$kendall)), 5),
    " dCor:", round(diff(range(results$dcor)), 5), "\n")

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(results, file.path(out_dir, "08_h3_invariance.csv"), row.names = FALSE)
