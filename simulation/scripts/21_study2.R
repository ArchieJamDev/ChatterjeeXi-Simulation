# Study 2 (BRM revision): xi_n against Spearman, Pearson, distance correlation and
# tau* -- type I error and power with discrete responses, several sample sizes and
# a NEUTRAL calibration of the dependence (signal share s, see R/study2_lib.R).
#
# One invocation = one block + one shard (GitHub Actions: .github/workflows/study2.yml):
#
#   block        shard         what it computes
#   typeI        n (50..400)   type I error of all five tests under independence;
#                              X continuous, 5 Y types; R = 1000, dcor B = 199
#   power_cheap  structure     power of xi, Spearman, Pearson; s in {.05,.13,.26,.50},
#                              n in {50,100,200,400,800}, 5 Y types; R = 2000
#   power_full   structure     power of all five tests (adds dCor and tau*);
#                              s in {.13,.26}, n in {50,100,200,400}, 4 Y types; R = 500
#   power_circ_hi circular     same as power_full for the circular structure only, at the
#                              much higher signal shares {.70,.85,.95} (at s <= .50 no test has
#                              power for a circular relationship, so that range is uninformative)
#
# X is continuous throughout (discretizing X is studied in Study 1).
#
# NOTE: run from the simulation/ directory, e.g.
#   cd simulation && Rscript scripts/21_study2.R typeI 100
# STUDY2_SMOKE=1 -> tiny run to a temp dir; STUDY2_R overrides R; STUDY2_OUT sets the output dir.
suppressPackageStartupMessages({library(XICOR); library(energy); library(TauStar)})
source("R/dgp.R")
source("R/study1_lib.R")
source("R/study2_lib.R")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
block <- args[1]; shard <- args[2]
smoke <- identical(Sys.getenv("STUDY2_SMOKE"), "1")
out_dir <- Sys.getenv("STUDY2_OUT", unset = if (smoke) tempdir() else "results/study2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cfg <- list(typeI = list(R = 1000, B = 199), power_cheap = list(R = 2000, B = 0),
            power_full = list(R = 500, B = 199), power_circ_hi = list(R = 500, B = 199))
stopifnot(block %in% names(cfg))
if (smoke) cfg <- lapply(cfg, function(z) list(R = 10, B = if (z$B > 0) 19 else 0))
R <- cfg[[block]]$R; B <- cfg[[block]]$B
if (nzchar(Sys.getenv("STUDY2_R"))) R <- as.integer(Sys.getenv("STUDY2_R"))

shard_idx <- if (shard %in% STRUCTS2) match(shard, STRUCTS2) else as.integer(shard)
seed <- 20262100L + 10000L * match(block, names(cfg)) + shard_idx
log_cell <- function(...) { cat(sprintf(...), "\n"); flush.console() }
t0 <- Sys.time(); elapsed <- function() as.numeric(difftime(Sys.time(), t0, units = "mins"))

build_margins <- function(types, pop_latent) {
  set.seed(1)
  setNames(lapply(types, function(t) make_margin(t, pop_latent)), types)
}

rows <- list()
add_row <- function(extra, cell) rows[[length(rows) + 1]] <<- cbind(extra, cell)

if (block == "typeI") {
  n <- as.integer(shard); stopifnot(n %in% c(50, 100, 200, 400))
  pop <- { set.seed(1); rnorm(2e5) }
  my <- build_margins(Y_TYPES2, pop)
  set.seed(seed)
  for (yt in Y_TYPES2) {
    gen <- function() list(x = rnorm(n), y = my[[yt]](rnorm(n)))
    add_row(data.frame(block = block, shard = shard, structure = "indep", share = 0, sigma = NA_real_,
                       ytype = yt, n = n), run_cell2(gen, R, full = TRUE, B = B))
    log_cell("[typeI %s] Y=%s done (%.1f min)", shard, yt, elapsed())
  }
} else {
  stopifnot(shard %in% STRUCTS2)
  if (block == "power_circ_hi") stopifnot(shard == "circular")
  st <- shard
  shares <- switch(block, power_cheap = S_ALL, power_full = c(0.13, 0.26),
                   power_circ_hi = c(0.70, 0.85, 0.95))
  ns <- if (block == "power_cheap") c(50, 100, 200, 400, 800) else c(50, 100, 200, 400)
  yts <- if (block == "power_cheap") Y_TYPES2 else Y_FULL
  full <- block != "power_cheap"
  for (s in shares) {
    sig <- sigma_for_share(st, s)
    set.seed(1); big <- generate_xy(2e5, st, "normal", sig)
    my <- build_margins(yts, big$Y)
    set.seed(seed + as.integer(round(1000 * s)))
    for (n in ns) for (yt in yts) {
      gen <- function() { xy <- generate_xy(n, st, "normal", sig); list(x = xy$X, y = my[[yt]](xy$Y)) }
      add_row(data.frame(block = block, shard = shard, structure = st, share = s, sigma = sig,
                         ytype = yt, n = n), run_cell2(gen, R, full = full, B = B))
      log_cell("[%s %s] s=%.2f n=%d Y=%s done (%.1f min)", block, shard, s, n, yt, elapsed())
    }
  }
}

res <- do.call(rbind, rows)
res$seed <- seed
write.csv(res, file.path(out_dir, sprintf("%s_%s.csv", block, shard)), row.names = FALSE)
writeLines(c(R.version.string, paste("XICOR", as.character(packageVersion("XICOR"))),
             paste("energy", as.character(packageVersion("energy"))),
             paste("TauStar", as.character(packageVersion("TauStar"))), paste("seed", seed)),
           file.path(out_dir, sprintf("%s_%s_env.txt", block, shard)))
cat("DONE", block, shard, "->", out_dir, "\n")
