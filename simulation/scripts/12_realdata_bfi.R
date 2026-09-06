# Real-data illustration: psych::bfi (manuscript Section 5, Table
# "tab:realdata")
#
# Three illustrative pairs, one per variable-type comparison in the
# factorial design:
#   1. Continua-continua: age vs. neuroticism composite (mean of N1-N5,
#      no reverse-scoring needed per the standard bfi scoring key).
#   2. Ordinal-ordinal: two raw Likert items of the same facet (N1, N2)
#      -- reports xi_n under B independent tie-break replications to
#      illustrate H2b's variance-from-tie-breaking finding directly on
#      real data.
#   3. Binaria-continua: gender vs. extraversion composite (mean of
#      E1,E2 reverse-scored + E3,E4,E5 as-is, per the standard bfi key).
#
# NOTE: run from the simulation/ directory, e.g.
#   cd simulation && Rscript scripts/12_realdata_bfi.R
source("R/statistics.R")

suppressPackageStartupMessages(library(psych))
data(bfi)

reverse6 <- function(x) 7 - x # 1-6 Likert reversal

## Pair 1: age (continuous) vs. neuroticism composite (continuous-ish)
d1 <- na.omit(bfi[, c("age", paste0("N", 1:5))])
neuro <- rowMeans(d1[, paste0("N", 1:5)])
age <- d1$age

## Pair 2: two raw Likert items, same facet (N1, N2)
d2 <- na.omit(bfi[, c("N1", "N2")])
B_ties <- 20
xi_N1N2 <- replicate(B_ties, xicor(d2$N1, d2$N2))
xi_N2N1 <- replicate(B_ties, xicor(d2$N2, d2$N1))

## Pair 3: gender (binary) vs. extraversion composite
d3 <- na.omit(bfi[, c("gender", paste0("E", 1:5))])
extra <- rowMeans(cbind(reverse6(d3$E1), reverse6(d3$E2), d3$E3, d3$E4, d3$E5))
gender <- d3$gender

compute_row <- function(X, Y) {
  c(xi_XY = xicor(X, Y), xi_YX = xicor(Y, X), dcor = dcor(X, Y),
    taustar = tStar(X, Y), pearson = cor(X, Y), spearman = cor(X, Y, method = "spearman"))
}

tab <- data.frame(
  `Edad-Neuroticismo` = compute_row(age, neuro),
  `Item1-Item2 (Likert)` = c(xi_XY = mean(xi_N1N2), xi_YX = mean(xi_N2N1),
                              dcor = dcor(d2$N1, d2$N2), taustar = tStar(d2$N1, d2$N2),
                              pearson = cor(d2$N1, d2$N2), spearman = cor(d2$N1, d2$N2, method = "spearman")),
  `Genero-Extraversion` = compute_row(gender, extra),
  check.names = FALSE
)

cat("=== Real-data illustration: bfi ===\n")
cat("n (pair 1, age-neuroticism):", nrow(d1), "\n")
cat("n (pair 2, N1-N2 Likert):", nrow(d2), "\n")
cat("n (pair 3, gender-extraversion):", nrow(d3), "\n")
cat("SD of xi_n(N1,N2) across", B_ties, "tie-breaks:", round(sd(xi_N1N2), 4), "\n\n")
print(round(tab, 4))

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(tab, file.path(out_dir, "12_realdata_bfi.csv"))
