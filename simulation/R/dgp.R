# Data-generating mechanisms (manuscript Section 3.2)
#
# Five marginal families for the latent variable Z, and five
# dependence structures f(). For "circular" and "heterocedastica",
# the observed X is Z bounded into (-1,1) via its own probability
# integral transform, X = 2*F_Z(Z) - 1 (always exactly Uniform(-1,1)
# regardless of Z's family; see paper Section 3.2). For "lineal",
# "cuadratica", "senoidal", X = Z directly (no bounding needed).

FAMILIES <- c("normal", "uniform", "exponential", "t", "mixture")
STRUCTURES <- c("lineal", "cuadratica", "senoidal", "circular", "heterocedastica")

#' Draw n i.i.d. values from one of the five latent marginal families.
rlatent <- function(n, family, df = 3) {
  family <- match.arg(family, FAMILIES)
  switch(family,
    normal      = rnorm(n),
    uniform     = runif(n, -1, 1),
    exponential = rexp(n, rate = 1),
    t           = rt(n, df = df),
    mixture     = {
      comp <- rbinom(n, 1, 0.5)
      ifelse(comp == 1, rnorm(n, -2, 0.5), rnorm(n, 2, 0.5))
    }
  )
}

#' Theoretical CDF of the latent family at points z (vectorized).
plandau_family <- function(z, family, df = 3) {
  family <- match.arg(family, FAMILIES)
  switch(family,
    normal      = pnorm(z),
    uniform     = punif(z, -1, 1),
    exponential = pexp(z, rate = 1),
    t           = pt(z, df = df),
    mixture     = 0.5 * pnorm(z, -2, 0.5) + 0.5 * pnorm(z, 2, 0.5)
  )
}

#' Bound a latent variable into (-1,1) via its own probability
#' integral transform. Always exactly Uniform(-1,1), by construction.
bound_pm1 <- function(z, family, df = 3) {
  2 * plandau_family(z, family, df = df) - 1
}

#' Generate one (X,Y) sample under a given structure, latent family,
#' and noise level sigma. lambda is only used by "heterocedastica".
#'
#' Returns a list(X=, Y=, Z=) -- Z is the pre-bounding latent variable
#' for circular/heterocedastica (NULL otherwise), kept because H3's
#' "forma marginal" for those two structures refers to Z, not to the
#' already-bounded X (paper Section 3.2).
generate_xy <- function(n, structure, family, sigma, lambda = 0.5, df = 3) {
  structure <- match.arg(structure, STRUCTURES)
  E <- rlatent(n, family, df = df)

  if (structure %in% c("lineal", "cuadratica", "senoidal")) {
    X <- rlatent(n, family, df = df)
    f <- switch(structure,
      lineal     = X,
      cuadratica = X^2,
      senoidal   = sin(2 * pi * X)
    )
    Y <- f + sigma * E
    return(list(X = X, Y = Y, Z = NULL))
  }

  # circular / heterocedastica: need bounded support
  Z <- rlatent(n, family, df = df)
  X <- bound_pm1(Z, family, df = df)

  if (structure == "circular") {
    Zsign <- sample(c(-1, 1), n, replace = TRUE)
    Y <- Zsign * sqrt(pmax(1 - X^2, 0)) + sigma * E
  } else { # heterocedastica
    s <- as.numeric(abs(X) <= 0.5)
    Y <- 3 * (s * (1 - lambda) + lambda) * E
  }
  list(X = X, Y = Y, Z = Z)
}
