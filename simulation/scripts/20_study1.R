# Study 1 (BRM revision): type I error and power of xi_n's tests with discrete
# margins (binary, ordinal, composite scores) and realistic category counts.
#
# One invocation = one block + one shard, so the GitHub Actions workflow
# (.github/workflows/study1.yml) can run the shards as parallel jobs:
#
#   block        shard        what it computes
#   typeI_asym   n            tie-aware and naive asymptotic type I error under exact
#                             independence, 3 X types x 8 Y types, R = 5000
#   typeI_perm   n (100/250/500)  same datasets + permutation test (B = 499),
#                             X in {cont, bin50} x Y in {cont, bin10, bin05, ord7f}, R = 1000
#   power_asym   structure    asymptotic power under dependence calibrated on the LATENT
#                             variable (xi = 0.15), 3 X types x 8 Y types x 5 n, R = 2000
#   power_perm   n (100/250)  asymptotic + permutation power (B = 199), 3 structures x
#                             X in {cont, bin50} x Y in {cont, bin10, bin05, ord7f}, R = 500
#
# NOTE: run from the simulation/ directory, e.g.
#   cd simulation && Rscript scripts/20_study1.R typeI_asym 100
# Set STUDY1_SMOKE=1 for a tiny run (R = 20) that writes to a temp directory.
suppressPackageStartupMessages(library(XICOR))
source("R/dgp.R")
source("R/calibration.R")
source("R/study1_lib.R")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
block <- args[1]; shard <- args[2]
smoke <- identical(Sys.getenv("STUDY1_SMOKE"), "1")
out_dir <- Sys.getenv("STUDY1_OUT", unset = if (smoke) tempdir() else "results/study1")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

N_ALL <- c(50, 100, 250, 500, 1000)
STRUCTS <- c("lineal", "cuadratica", "senoidal")
XI_TARGET <- 0.15
PERM_X <- c("cont", "bin50")
PERM_Y <- c("cont", "bin10", "bin05", "ord7f")

cfg <- list(typeI_asym = list(R = 5000, B = 0), typeI_perm = list(R = 1000, B = 499),
            power_asym = list(R = 2000, B = 0), power_perm = list(R = 500, B = 199))
stopifnot(block %in% names(cfg))
if (smoke) cfg <- lapply(cfg, function(z) list(R = 20, B = if (z$B > 0) 29 else 0))
R <- cfg[[block]]$R; B <- cfg[[block]]$B
if (nzchar(Sys.getenv("STUDY1_R"))) R <- as.integer(Sys.getenv("STUDY1_R"))  # validation runs

shard_idx <- if (shard %in% STRUCTS) match(shard, STRUCTS) else as.integer(shard)
seed <- 20261100L + 10000L * match(block, names(cfg)) + shard_idx
log_cell <- function(...) { cat(sprintf(...), "\n"); flush.console() }

build_margins <- function(types, pop_latent, factory = make_margin) {
  set.seed(1)  # fixed stream for the cut points, independent of the replication stream
  setNames(lapply(types, function(t) factory(t, pop_latent)), types)
}

rows <- list()
t0 <- Sys.time()
add_row <- function(extra, cell) rows[[length(rows) + 1]] <<- cbind(extra, cell)
elapsed <- function() as.numeric(difftime(Sys.time(), t0, units = "mins"))

if (block %in% c("typeI_asym", "typeI_perm")) {
  n <- as.integer(shard)
  stopifnot(n %in% (if (block == "typeI_asym") N_ALL else c(100, 250, 500)))
  xt_set <- if (block == "typeI_asym") X_TYPES else PERM_X
  yt_set <- if (block == "typeI_asym") Y_TYPES else PERM_Y
  pop <- { set.seed(1); rnorm(2e5) }
  mx <- build_margins(xt_set, pop); my <- build_margins(yt_set, pop)
  set.seed(seed)
  for (xt in xt_set) for (yt in yt_set) {
    gen <- function() list(x = mx[[xt]](rnorm(n)), y = my[[yt]](rnorm(n)))
    add_row(data.frame(block = block, shard = shard, structure = "indep", sigma = NA_real_,
                       xi_latent = 0, xtype = xt, ytype = yt, n = n), run_cell(gen, R, B))
    log_cell("[%s %s] X=%s Y=%s done (%.1f min)", block, shard, xt, yt, elapsed())
  }
} else {
  structs <- if (block == "power_asym") shard else STRUCTS
  ns <- if (block == "power_asym") N_ALL else as.integer(shard)
  if (block == "power_perm") stopifnot(ns %in% c(100, 250))
  xt_set <- if (block == "power_asym") X_TYPES else PERM_X
  yt_set <- if (block == "power_asym") Y_TYPES else PERM_Y
  for (st in structs) {
    cal <- calibrate_latent(st, XI_TARGET)
    log_cell("[%s %s] %s calibrated: sigma=%.4f latent xi=%.4f (%s)", block, shard, st,
             cal$sigma, cal$xi_hat, cal$method)
    set.seed(1); big <- generate_xy(2e5, st, "normal", cal$sigma)
    mx <- build_margins(xt_set, big$X); my <- build_margins(yt_set, big$Y)
    set.seed(seed)
    for (n in ns) for (xt in xt_set) for (yt in yt_set) {
      gen <- function() {
        xy <- generate_xy(n, st, "normal", cal$sigma)
        list(x = mx[[xt]](xy$X), y = my[[yt]](xy$Y))
      }
      add_row(data.frame(block = block, shard = shard, structure = st, sigma = cal$sigma,
                         xi_latent = cal$xi_hat, xtype = xt, ytype = yt, n = n),
              run_cell(gen, R, B))
      log_cell("[%s %s] %s n=%d X=%s Y=%s done (%.1f min)", block, shard, st, n, xt, yt, elapsed())
    }
  }
}

res <- do.call(rbind, rows)
res$seed <- seed
write.csv(res, file.path(out_dir, sprintf("%s_%s.csv", block, shard)), row.names = FALSE)
writeLines(c(R.version.string, paste("XICOR", as.character(packageVersion("XICOR"))),
             paste("seed", seed)), file.path(out_dir, sprintf("%s_%s_env.txt", block, shard)))
cat("DONE", block, shard, "->", out_dir, "\n")
