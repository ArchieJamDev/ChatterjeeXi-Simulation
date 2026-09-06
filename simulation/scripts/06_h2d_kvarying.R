# Follow-up to H2d: does the number of composite items k matter, and
# does the benefit keep growing for large k (300-400 item inventories
# exist in real personality/interest assessment)? (manuscript Table
# "tab:h2d-kvar")
#
# Reference cell: n de referencia = 100, forma = normal, estructura =
# lineal, regimen = moderado (reuses 01_h1's calibrated sigma/xi_true).
# n here is always the number of SUBJECTS (n=100 people); k is the
# number of items averaged per subject into one composite score --
# two independent dimensions, never conflated.
#
# NOTE: run from the simulation/ directory, e.g.
#   cd simulation && Rscript scripts/06_h2d_kvarying.R
source("R/dgp.R")
source("R/discretize.R")
source("R/statistics.R")

set.seed(20260906)

cal <- readRDS("results/01_h1_calibration.rds")
sigma <- cal$sigma
n_ref <- 100
R <- 500
alpha <- 0.05
k_values <- c(1, 2, 3, 5, 10, 20, 40, 100, 200, 400)

out <- data.frame(k = k_values, varianza_empirica = NA_real_, potencia = NA_real_,
                   valores_distintos = NA_real_)

for (i in seq_along(k_values)) {
  k <- k_values[i]
  xi_vals <- numeric(R); rej <- logical(R); ndist <- integer(R)
  for (r in seq_len(R)) {
    xy <- generate_xy(n_ref, "lineal", "normal", sigma)
    cs <- composite_score(xy$Y, k = k, categories = 6, item_noise_sd = 1)
    fit <- xicor(xy$X, cs$composite, pvalue = TRUE)
    xi_vals[r] <- fit$xi
    rej[r] <- fit$pval < alpha
    ndist[r] <- length(unique(cs$composite))
  }
  out$varianza_empirica[i] <- var(xi_vals)
  out$potencia[i] <- mean(rej)
  out$valores_distintos[i] <- mean(ndist)
  cat("k =", k, "done.\n")
}

cat("\n=== H2d follow-up: effect of k (n =", n_ref, "subjects, sigma =", round(sigma, 4), ") ===\n")
print(out, row.names = FALSE)
cat("\nNote: 'valores_distintos' saturates near n =", n_ref,
    "(the number of subjects), not near the composite's own theoretical range 5k+1 --",
    "this is why variance stops shrinking once virtually every subject has a unique score.\n")

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(out, file.path(out_dir, "06_h2d_kvarying.csv"), row.names = FALSE)
