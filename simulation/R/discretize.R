# Discretization for H2a/b/c (manuscript Section 3.2)
#
# Cuts a continuous variable into c equal-probability ordered
# categories (1..c). Uses empirical quantiles of the observed sample
# rather than a theoretical CDF: for X (whose family is known exactly)
# the two coincide asymptotically; for Y = f(X) + sigma*E, no simple
# closed-form CDF exists in general, so the empirical-quantile version
# is the only one that is always well defined and is what the
# implementation actually uses for both X and Y, for consistency.

#' Discretize a numeric vector into `categories` equal-probability
#' ordered bins, coded 1..categories.
discretize <- function(x, categories) {
  probs <- seq(0, 1, length.out = categories + 1)
  breaks <- quantile(x, probs = probs, type = 7, names = FALSE)
  breaks[1] <- -Inf
  breaks[length(breaks)] <- Inf
  as.integer(cut(x, breaks = breaks, include.lowest = TRUE, labels = FALSE))
}

#' Randomly re-break ties within a discretized vector's rank
#' computation is handled by xi_n itself (Sec 2.1); this helper just
#' re-applies `discretize` with a fresh RNG state, used to generate the
#' B nested tie-break replicates for H2b/H2c (manuscript Section 3.2).
redraw_ties_note <- "Tie-breaking for xi_n happens inside xicor(); no separate helper needed here."

#' H2d: composite score from k Likert items measuring the same latent
#' trait theta, each with independent measurement error (manuscript
#' Section 3.2). Returns the item-average, plus the raw single-item
#' version (first item) for the H2d comparison table.
composite_score <- function(theta, k = 10, categories = 6, item_noise_sd = 1) {
  n <- length(theta)
  items <- sapply(seq_len(k), function(j) {
    discretize(theta + rnorm(n, 0, item_noise_sd), categories)
  })
  list(single_item = items[, 1], composite = rowMeans(items))
}
