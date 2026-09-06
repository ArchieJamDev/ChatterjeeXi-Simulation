# Module: H2a -- empates en la respuesta (manuscript Table "tab:h2a")
#
# Reference cell: n de referencia = 100 (typical psychology sample
# size), forma = normal, estructura = lineal, regimen = moderado
# (reuses the sigma/xi_true calibrated in 01_h1.R for consistency
# across modules).
#
# X stays continuous; Y is discretized into c categories (H2a design:
# ties only in the response, not the covariate -- xi_n's tie-breaking
# randomization from ties in X, Sec 2.1, does not apply here).
#
# NOTE: run from the simulation/ directory, e.g.
#   cd simulation && Rscript scripts/02_h2a.R
source("R/dgp.R")
source("R/discretize.R")
source("R/statistics.R")

set.seed(20260906)

cal <- readRDS("results/01_h1_calibration.rds")
sigma <- cal$sigma
n_ref <- 100
R <- 1000
alpha <- 0.05

conditions <- list(
  Continua = NA, # NA = no discretization
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
  xi_hats <- numeric(R)
  rejections <- logical(R)
  for (r in seq_len(R)) {
    xy <- generate_xy(n_ref, "lineal", "normal", sigma)
    Yobs <- if (is.na(c_cat)) xy$Y else discretize(xy$Y, c_cat)
    fit <- xicor(xy$X, Yobs, pvalue = TRUE)
    xi_hats[r] <- fit$xi
    rejections[r] <- fit$pval < alpha
  }
  results$varianza_empirica[i] <- var(xi_hats)
  results$potencia[i] <- mean(rejections)
  cat(names(conditions)[i], "done.\n")
}

cat("\n=== H2a: empates en la respuesta (n =", n_ref, ", sigma =", round(sigma,4),
    ", xi_true =", round(cal$xi_true, 4), ") ===\n")
print(results, row.names = FALSE)

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(results, file.path(out_dir, "02_h2a.csv"), row.names = FALSE)
