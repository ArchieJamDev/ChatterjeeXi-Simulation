# Calibration of the dependence regime xi_0 via the noise level sigma
# (manuscript Section 3.2, "Calibracion del regimen de dependencia").
#
# Three precautions from the design, implemented literally:
#  1. Monotonicity of sigma -> xi(X,Y) is NOT assumed a priori; it is
#     checked empirically on a grid before choosing bisection vs. a
#     global grid search.
#  2. For H2 conditions, calibration happens on the latent continuous
#     variable, BEFORE discretizing.
#  3. Each statistic (xi_n, dCor, tau*) gets its OWN population
#     reference under the SAME (structure, family, sigma) cell -- this
#     file only calibrates sigma to hit a target xi_n population value;
#     the resulting population dCor/tau* under that same sigma are
#     whatever they are, and are estimated separately (population_stats()
#     below), never forced to match xi_0.

#' Monte Carlo estimate of the population value of a statistic under
#' one (structure, family, sigma) cell, using a single very large
#' sample (n_pop) as a cheap plug-in estimator of the population
#' functional. stat_fun(X,Y) -> scalar (e.g. function(x,y) xicor(x,y)).
population_stat <- function(structure, family, sigma, stat_fun,
                             n_pop = 1e5, lambda = 0.5, df = 3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  xy <- generate_xy(n_pop, structure, family, sigma, lambda = lambda, df = df)
  stat_fun(xy$X, xy$Y)
}

#' Check monotonicity of sigma -> xi(X,Y) on a grid, for one
#' (structure, family) combination. Returns the grid and a logical
#' flag; used to decide bisection vs. grid search before calibrating.
check_monotonicity <- function(structure, family, sigma_grid = 10^seq(-2, 1, length.out = 15),
                                n_pop = 2e4, lambda = 0.5, df = 3, seed = 1) {
  xi_vals <- vapply(sigma_grid, function(s) {
    population_stat(structure, family, s, function(x, y) xicor(x, y),
                     n_pop = n_pop, lambda = lambda, df = df, seed = seed)
  }, numeric(1))
  list(sigma = sigma_grid, xi = xi_vals,
       monotone = all(diff(xi_vals) <= 1e-6)) # non-increasing in sigma
}

#' Calibrate sigma to hit a target xi_0, for one (structure, family)
#' combination. Uses bisection if the sigma->xi curve is verified
#' monotone (non-increasing); otherwise falls back to picking the grid
#' point closest to xi_0 (a global search), as specified in the design.
calibrate_sigma <- function(structure, family, xi_target, lambda = 0.5, df = 3,
                             sigma_range = c(1e-3, 20), tol = 0.01,
                             n_pop = 2e4, max_iter = 40, seed = 1,
                             mono_check = NULL) {
  if (is.null(mono_check)) mono_check <- check_monotonicity(structure, family, lambda = lambda, df = df, seed = seed)

  xi_of_sigma <- function(s) population_stat(structure, family, s, function(x, y) xicor(x, y),
                                              n_pop = n_pop, lambda = lambda, df = df, seed = seed)

  if (isTRUE(mono_check$monotone)) {
    lo <- sigma_range[1]; hi <- sigma_range[2]
    for (i in seq_len(max_iter)) {
      mid <- sqrt(lo * hi) # geometric bisection: sigma spans orders of magnitude
      xi_mid <- xi_of_sigma(mid)
      if (abs(xi_mid - xi_target) < tol) return(list(sigma = mid, xi_hat = xi_mid, method = "bisection", iters = i))
      if (xi_mid > xi_target) lo <- mid else hi <- mid # xi decreasing in sigma
    }
    return(list(sigma = mid, xi_hat = xi_mid, method = "bisection_maxiter", iters = max_iter))
  }

  # non-monotone: global grid search over a denser mesh
  grid <- exp(seq(log(sigma_range[1]), log(sigma_range[2]), length.out = 60))
  xi_grid <- vapply(grid, xi_of_sigma, numeric(1))
  best <- which.min(abs(xi_grid - xi_target))
  list(sigma = grid[best], xi_hat = xi_grid[best], method = "grid_search", iters = length(grid))
}
