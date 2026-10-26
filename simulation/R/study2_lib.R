# Study 2 helpers (BRM revision): xi_n against competing measures of association
# (Spearman, Pearson, distance correlation, tau*) in type I error and power, with
# discrete (binary / ordinal / composite) responses and several sample sizes.
#
# NEUTRAL CALIBRATION. Dependence strength is NOT matched on any coefficient's own
# scale. It is set through the "signal share"
#     s = Var(noiseless structure) / Var(latent Y)  (= sigma-free explained-variance
#                                                    share of the true function),
# i.e. sigma^2 = Var(f) * (1 - s) / s, with Var(f) taken from a large noiseless
# sample. s is the R^2 of the true generating function, so the same s means the same
# amount of "signal" for the linear, quadratic, sinusoidal and circular structures.
# (Pearson's r^2 equals s only for the linear structure, by design.)
#
# Requires study1_lib.R (make_margin, min_count) and dgp.R.

STRUCTS2 <- c("lineal", "cuadratica", "senoidal", "circular")
S_ALL <- c(0.05, 0.13, 0.26, 0.50)       # 0.13 and 0.26 = Cohen's medium and large R^2
Y_TYPES2 <- c("cont", "bin50", "bin10", "ord5s", "comp5")
Y_FULL <- c("cont", "bin50", "ord5s", "comp5")
STATS <- c("xi", "spearman", "pearson", "dcor", "tau")

#' Noise level for a target signal share s on the LATENT scale.
sigma_for_share <- function(structure, s, n_big = 2e5) {
  set.seed(1)
  v <- var(generate_xy(n_big, structure, "normal", 0)$Y)   # sigma = 0: Y = f(X)
  sqrt(v * (1 - s) / s)
}

#' p-values of the five tests for one sample. xi: tie-aware asymptotic test (XICOR
#' default); Spearman and Pearson: two-sided t tests; dcor: permutation test (B
#' permutations, energy::dcor.test); tau: TauStar::tauStarTest (mode "auto").
#' dcor and tau are only computed when full = TRUE (they are much more expensive).
all_tests <- function(x, y, B, full) {
  out <- setNames(rep(NA_real_, length(STATS) + 1), c(STATS, "tauwarn"))
  if (length(unique(x)) < 2 || length(unique(y)) < 2) return(out)
  out["xi"] <- XICOR::xicor(x, y, pvalue = TRUE)$pval
  out["spearman"] <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)$p.value)
  out["pearson"] <- suppressWarnings(cor.test(x, y)$p.value)
  if (full) {
    out["dcor"] <- tryCatch(energy::dcor.test(x, y, R = B)$p.value, error = function(e) NA_real_)
    nw <- 0
    out["tau"] <- tryCatch(
      withCallingHandlers(as.numeric(TauStar::tauStarTest(x, y)$pVal),
                          warning = function(w) { nw <<- nw + 1; invokeRestart("muffleWarning") }),
      error = function(e) NA_real_)
    out["tauwarn"] <- as.numeric(nw > 0)   # numerical warnings of tauStarTest, counted not hidden
  }
  out
}

#' R replications of one cell; gen() returns list(x=, y=) of observed margins.
#' Returns one row with the rejection rate of every test at alpha = .05 and .01
#' (NA when a test was not computed; NA p-values count as non-rejections).
run_cell2 <- function(gen, R, full, B = 199) {
  P <- matrix(NA_real_, R, length(STATS) + 1, dimnames = list(NULL, c(STATS, "tauwarn")))
  for (r in seq_len(R)) {
    d <- gen()
    P[r, ] <- all_tests(d$x, d$y, B, full)
  }
  rej <- function(a) sapply(STATS, function(s) {
    p <- P[, s]; if (all(is.na(p))) NA_real_ else sum(p < a, na.rm = TRUE) / R })
  r05 <- rej(0.05); r01 <- rej(0.01)
  names(r05) <- paste0("rej05_", STATS); names(r01) <- paste0("rej01_", STATS)
  data.frame(R = R, B = if (full) B else NA_integer_, na_rate_xi = mean(is.na(P[, "xi"])),
             na_rate_tau = if (full) mean(is.na(P[, "tau"])) else NA_real_,
             warn_rate_tau = if (full) mean(P[, "tauwarn"], na.rm = TRUE) else NA_real_,
             t(r05), t(r01))
}
