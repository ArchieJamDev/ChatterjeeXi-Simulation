# H3 follow-up: Type I error and power per family (manuscript Table
# "tab:h3", added per external review -- H3 promises to evaluate
# "sesgo, varianza, y potencia").
#
# Runs ONE family per invocation (pass as a command-line argument) to
# keep memory use low on this machine (previous all-in-one run was
# OOM-killed by the kernel). R_power and nperm kept deliberately small.
# Appends its one-row result to results/07b_h3_typeI_power.csv.
#
# NOTE: run from the simulation/ directory, once per family, e.g.
#   cd simulation
#   for f in normal uniform exponential t mixture; do
#     Rscript scripts/07b_h3_typeI_power.R $f
#   done
source("R/dgp.R")
source("R/statistics.R")
source("R/calibration.R")

args <- commandArgs(trailingOnly = TRUE)
fam <- args[1]
stopifnot(fam %in% FAMILIES)

set.seed(20260906)

n_ref <- 100
R_power <- 50
nperm <- 49
alpha <- 0.05
xi_target <- 0.5

mono <- check_monotonicity("lineal", fam, seed = 1)
cal <- calibrate_sigma("lineal", fam, xi_target = xi_target, mono_check = mono, seed = 1)
sigma <- cal$sigma
rm(mono); gc(verbose = FALSE)

rej_xi_h0 <- rej_dcor_h0 <- rej_xi_h1 <- rej_dcor_h1 <- logical(R_power)
for (r in seq_len(R_power)) {
  X0 <- rlatent(n_ref, fam); Y0 <- rlatent(n_ref, fam)
  rej_xi_h0[r] <- xicor(X0, Y0, pvalue = TRUE)$pval < alpha
  rej_dcor_h0[r] <- dcor.test(X0, Y0, R = nperm)$p.value < alpha

  xy1 <- generate_xy(n_ref, "lineal", fam, sigma)
  rej_xi_h1[r] <- xicor(xy1$X, xy1$Y, pvalue = TRUE)$pval < alpha
  rej_dcor_h1[r] <- dcor.test(xy1$X, xy1$Y, R = nperm)$p.value < alpha

  if (r %% 10 == 0) gc(verbose = FALSE)
}

row <- data.frame(
  familia = fam, sigma = sigma,
  xi_n_tipo1 = mean(rej_xi_h0), xi_n_potencia = mean(rej_xi_h1),
  dcor_tipo1 = mean(rej_dcor_h0), dcor_potencia = mean(rej_dcor_h1)
)
cat(fam, ": xi_tipo1=", round(row$xi_n_tipo1,3), " xi_pot=", round(row$xi_n_potencia,3),
    " dcor_tipo1=", round(row$dcor_tipo1,3), " dcor_pot=", round(row$dcor_potencia,3), "\n")

out_file <- "results/07b_h3_typeI_power.csv"
dir.create("results", showWarnings = FALSE, recursive = TRUE)
write.table(row, out_file, sep = ",", row.names = FALSE,
            col.names = !file.exists(out_file), append = file.exists(out_file))
