# Study 1 helpers (BRM revision): type I error and power of xi_n's tests when X
# and/or Y are discrete, with REALISTIC (random) category counts.
#
# Observed margins are built by cutting the LATENT variable at POPULATION cut
# points (fixed once per margin from a large latent sample), so category counts
# vary from sample to sample as in real data. (Cutting at each sample's own
# quantiles would fix the counts deterministically and distorts the type I error
# of the asymptotic test under imbalance: see the pilot analysis in the repo log.)
#
# Dependence is calibrated on the LATENT continuous variable, before
# discretizing (same rule as in the main design, Section 3.2).

X_TYPES <- c("cont", "bin50", "ord5s")
Y_TYPES <- c("cont", "bin50", "bin20", "bin10", "bin05", "ord5s", "ord7f", "comp5")

# Category probabilities (first category = lowest latent values). "s" = skewed,
# "f" = floor effect; both mimic typical item distributions.
MARGIN_PROBS <- list(
  bin50 = c(.50, .50), bin20 = c(.80, .20), bin10 = c(.90, .10), bin05 = c(.95, .05),
  ord5s = c(.05, .10, .15, .30, .40),
  ord7f = c(.35, .25, .15, .10, .07, .05, .03)
)
COMP_K <- 5      # items per composite score
COMP_CATS <- 6   # response categories per item (1-6, as in bfi)

#' Build a function latent -> observed for a margin type. `pop_latent` is a large
#' sample from the latent margin, used only to fix the population cut points.
make_margin <- function(type, pop_latent) {
  if (type == "cont") return(function(z) z)
  if (type == "comp5") {
    mu <- mean(pop_latent); s <- sd(pop_latent)
    pitem <- (pop_latent - mu) / s + rnorm(length(pop_latent))
    cuts <- quantile(pitem, seq_len(COMP_CATS - 1) / COMP_CATS, names = FALSE)
    return(function(z) {
      th <- (z - mu) / s
      items <- vapply(seq_len(COMP_K),
                      function(j) findInterval(th + rnorm(length(th)), cuts) + 1L,
                      integer(length(th)))
      rowMeans(items)
    })
  }
  pr <- MARGIN_PROBS[[type]]
  cuts <- quantile(pop_latent, cumsum(pr)[-length(pr)], names = FALSE)
  function(z) findInterval(z, cuts) + 1L
}

#' Smallest non-empty category count (minority count k) for integer-coded margins.
min_count <- function(y) {
  if (!is.integer(y)) return(NA_real_)
  tb <- tabulate(y)
  min(tb[tb > 0])
}

#' xi_n tie-aware asymptotic p-value (XICOR default, ties = TRUE) and the NAIVE
#' p-value that assumes continuous Y (sd = sqrt(2/(5n))), from the SAME xi.
#' Returns c(p_tie, p_naive, xi); NAs if X or Y is constant.
asym_test <- function(x, y) {
  n <- length(x)
  if (length(unique(y)) < 2 || length(unique(x)) < 2) return(c(NA_real_, NA_real_, NA_real_))
  r <- XICOR::xicor(x, y, pvalue = TRUE)
  c(r$pval, 1 - pnorm(sqrt(n) * r$xi / sqrt(2 / 5)), r$xi)
}

#' Permutation p-value for xi_n: p = (1 + #{xi* >= xi_obs}) / (B + 1), shuffling Y
#' (X ties are re-broken at random in every permutation). Uses XICOR's internal
#' calculateXI(), checked to exist so a package change fails loudly.
perm_test <- function(x, y, B, xi_obs) {
  stopifnot(exists("calculateXI", envir = asNamespace("XICOR")))
  f <- get("calculateXI", envir = asNamespace("XICOR"))
  rp <- vapply(seq_len(B), function(b) f(x, sample(y)), numeric(1))
  (1 + sum(rp >= xi_obs)) / (B + 1)
}

#' Run R replications of one cell. `gen()` returns list(x=, y=) of OBSERVED margins.
#' B > 0 adds the permutation test on the same datasets. Returns a one-row data.frame.
run_cell <- function(gen, R, B = 0) {
  p <- matrix(NA_real_, R, 3); pp <- rep(NA_real_, R); k <- rep(NA_real_, R)
  for (r in seq_len(R)) {
    d <- gen()
    p[r, ] <- asym_test(d$x, d$y)
    k[r] <- min_count(d$y)
    if (B > 0 && !is.na(p[r, 1])) pp[r] <- perm_test(d$x, d$y, B, p[r, 3])
  }
  ok <- !is.na(p[, 1])
  data.frame(
    R = R, B = B, na_rate = mean(!ok),
    rej05 = sum(p[ok, 1] < 0.05) / R, rej01 = sum(p[ok, 1] < 0.01) / R,
    rej05_naive = sum(p[ok, 2] < 0.05) / R,
    rej05_perm = if (B > 0) sum(pp[ok] < 0.05) / R else NA_real_,
    rej01_perm = if (B > 0) sum(pp[ok] < 0.01) / R else NA_real_,
    mean_xi = mean(p[ok, 3]), sd_xi = sd(p[ok, 3]),
    mean_min_count = mean(k, na.rm = TRUE)
  )
}

#' Calibrate the noise level sigma of a structure on the LATENT scale to a target xi.
calibrate_latent <- function(structure, xi_target) {
  mono <- check_monotonicity(structure, "normal", seed = 1)
  cal <- calibrate_sigma(structure, "normal", xi_target = xi_target, mono_check = mono,
                         tol = 0.005, seed = 1)
  list(sigma = cal$sigma, xi_hat = cal$xi_hat, method = cal$method)
}
