# Module: H2d -- puntaje compuesto (manuscript Table "tab:h2d")
#
# Reference cell: n de referencia = 100, forma = normal, estructura =
# lineal, regimen = moderado (reuses 01_h1's calibrated sigma/xi_true).
# theta = Y (the response); k=10 Likert-6 items measure Y with
# independent item-level measurement error each; compares xi_n(X, .)
# for a single raw item vs. the k-item composite average.
#
# NOTE: run from the simulation/ directory, e.g.
#   cd simulation && Rscript scripts/05_h2d.R
source("R/dgp.R")
source("R/discretize.R")
source("R/statistics.R")

set.seed(20260906)

cal <- readRDS("results/01_h1_calibration.rds")
sigma <- cal$sigma
n_ref <- 100
R <- 1000
alpha <- 0.05
k_items <- 10
likert_cat <- 6

xi_single <- xi_composite <- numeric(R)
rej_single <- rej_composite <- logical(R)
distinct_single <- distinct_composite <- integer(R)

for (r in seq_len(R)) {
  xy <- generate_xy(n_ref, "lineal", "normal", sigma)
  cs <- composite_score(xy$Y, k = k_items, categories = likert_cat, item_noise_sd = 1)

  fit_s <- xicor(xy$X, cs$single_item, pvalue = TRUE)
  fit_c <- xicor(xy$X, cs$composite, pvalue = TRUE)

  xi_single[r]  <- fit_s$xi
  xi_composite[r] <- fit_c$xi
  rej_single[r] <- fit_s$pval < alpha
  rej_composite[r] <- fit_c$pval < alpha
  distinct_single[r] <- length(unique(cs$single_item))
  distinct_composite[r] <- length(unique(cs$composite))
}

results <- data.frame(
  condicion = c("Item unico (Likert-6)", "Puntaje compuesto (k=10)"),
  varianza_empirica = c(var(xi_single), var(xi_composite)),
  potencia = c(mean(rej_single), mean(rej_composite)),
  num_valores_distintos = c(round(mean(distinct_single), 1), round(mean(distinct_composite), 1))
)

cat("=== H2d: item unico vs. puntaje compuesto (n =", n_ref, ", sigma =", round(sigma, 4),
    ", xi_true =", round(cal$xi_true, 4), ") ===\n")
print(results, row.names = FALSE)

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(results, file.path(out_dir, "05_h2d.csv"), row.names = FALSE)
