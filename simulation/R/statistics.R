# Statistic wrappers (manuscript Section 2, Table "tab:comparison")
#
# xicor(x, y) computes xi_n in the sense of "how well y is determined
# by x" -- i.e. xicor(X, Y) == xi_n(X,Y) in the manuscript's notation
# (verified empirically: xicor(X, X^2) ~ 1, xicor(X^2, X) ~ 0.26 for
# X ~ Unif(-1,1) noiseless, matching the directionality expected from
# Chatterjee 2021).

suppressPackageStartupMessages({
  library(XICOR)
  library(energy)
  library(TauStar)
})
