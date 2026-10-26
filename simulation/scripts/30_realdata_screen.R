# Real-data screening (BRM revision): every pair of variables of two public psych
# datasets, analysed with xi_n, Spearman, Pearson, distance correlation and tau*.
#
# PRE-SPECIFIED RULE (no pair is chosen by its result): all pairs among
#   * psych::bfi: 25 items (A1-O5, as stored, 1-6), 5 Big Five composite scores
#     (mean of the 5 items after reverse-keying, 7 - x for A1, C4, C5, E1, E2, O2, O5;
#     NA if any item is missing), gender, education, age = 33 variables, 528 pairs, minus
#     the 25 pairs (item, own composite) -> 503 pairs;
#   * psych::sat.act: gender, education, age, ACT, SATV, SATQ -> 15 pairs.
# Each pair uses its complete cases.
#
# PRE-SPECIFIED OUTCOMES (computed here, summarized later):
#   * estimates: xi(X,Y) and xi(Y,X) averaged over B_TB random tie-breaks with their SD,
#     Spearman rho with Fisher-type CI (variance 1.06/(n-3)), Pearson r with CI, dCor,
#     tau* (TauStar), n, number of distinct values;
#   * p-values: xi asymptotic tie-aware test per direction (median over tie-breaks and share
#     of tie-breaks with p < .05), Spearman and Pearson t tests, dCor permutation test
#     (B = 499), tau* test (mode chosen from the data). No CI is given for xi: the bootstrap
#     fails for xi_n (Dalitz et al., 2024) and no meaningful CI exists in general.
#   * decision rule used in the summary: association "detected" by a measure if its
#     Benjamini-Hochberg-adjusted p < .05 across all pairs (also reported unadjusted);
#     xi's p-value = min(1, 2 * min over the two directions of the median p).
#
# Each pair is seeded from its own index (20263000 + index), so the results do not depend
# on how pairs are split into shards.
#
# NOTE: run from the simulation/ directory, e.g.
#   cd simulation && Rscript scripts/30_realdata_screen.R 1 12   # shard 1 of 12
# REALDATA_SMOKE=1 -> tiny run (few pairs, B = 19, 3 tie-breaks) to a temp dir;
# REALDATA_OUT sets the output directory; REALDATA_MAXPAIRS limits pairs per shard.
suppressPackageStartupMessages({library(XICOR); library(energy); library(TauStar); library(psych)})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
shard <- as.integer(args[1]); nsh <- as.integer(args[2])
stopifnot(shard >= 1, shard <= nsh)
smoke <- identical(Sys.getenv("REALDATA_SMOKE"), "1")
B_PERM <- if (smoke) 19L else 499L
B_TB <- if (smoke) 3L else 20L
out_dir <- Sys.getenv("REALDATA_OUT", unset = if (smoke) tempdir() else "results/realdata")
max_pairs <- if (nzchar(Sys.getenv("REALDATA_MAXPAIRS"))) as.integer(Sys.getenv("REALDATA_MAXPAIRS")) else
  if (smoke) 4L else Inf
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- data --------------------------------------------------------------------
data(bfi, package = "psych"); data(sat.act, package = "psych")
KEYS <- list(A = c("-A1", "A2", "A3", "A4", "A5"), C = c("C1", "C2", "C3", "-C4", "-C5"),
             E = c("-E1", "-E2", "E3", "E4", "E5"), N = c("N1", "N2", "N3", "N4", "N5"),
             O = c("O1", "-O2", "O3", "O4", "-O5"))
d_bfi <- bfi
for (sc in names(KEYS)) {
  items <- sub("^-", "", KEYS[[sc]]); neg <- grepl("^-", KEYS[[sc]])
  m <- sapply(seq_along(items), function(i) if (neg[i]) 7 - bfi[[items[i]]] else bfi[[items[i]]])
  d_bfi[[paste0(sc, "_comp")]] <- rowMeans(m)       # NA if any item is missing
}
bfi_items <- paste0(rep(c("A", "C", "E", "N", "O"), each = 5), 1:5)
bfi_vars <- c(bfi_items, paste0(c("A", "C", "E", "N", "O"), "_comp"), "gender", "education", "age")
is_own <- function(a, b) {
  f <- function(comp, item) grepl("_comp$", comp) && item %in% bfi_items && substr(item, 1, 1) == substr(comp, 1, 1)
  f(a, b) || f(b, a)
}
pb <- t(combn(bfi_vars, 2)); pb <- pb[!apply(pb, 1, function(p) is_own(p[1], p[2])), , drop = FALSE]
ps <- t(combn(c("gender", "education", "age", "ACT", "SATV", "SATQ"), 2))
pairs <- rbind(data.frame(dataset = "bfi", a = pb[, 1], b = pb[, 2]),
               data.frame(dataset = "sat.act", a = ps[, 1], b = ps[, 2]))
pairs$index <- seq_len(nrow(pairs))
stopifnot(nrow(pairs) == 503 + 15)
data_of <- list(bfi = d_bfi, sat.act = sat.act)

## ---- one pair -----------------------------------------------------------------
analyze_pair <- function(x, y) {
  n <- length(x)
  sp <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  rho <- unname(sp$estimate); sp_ci <- tanh(atanh(rho) + c(-1, 1) * qnorm(0.975) * sqrt(1.06 / (n - 3)))
  pe <- cor.test(x, y)
  xi_dir <- function(a, b) {
    r <- replicate(B_TB, { o <- XICOR::xicor(a, b, pvalue = TRUE); c(o$xi, o$pval) })
    c(mean = mean(r[1, ]), sd = sd(r[1, ]), p_med = median(r[2, ]), flip = mean(r[2, ] < 0.05))
  }
  xy <- xi_dir(x, y); yx <- xi_dir(y, x)
  dc <- energy::dcor.test(x, y, R = B_PERM)
  nw <- 0
  tau <- tryCatch(withCallingHandlers(TauStar::tauStarTest(x, y),
                  warning = function(w) { nw <<- nw + 1; invokeRestart("muffleWarning") }),
                  error = function(e) NULL)
  data.frame(
    n = n, kx = length(unique(x)), ky = length(unique(y)),
    spearman = rho, spearman_lo = sp_ci[1], spearman_hi = sp_ci[2], spearman_p = sp$p.value,
    pearson = unname(pe$estimate), pearson_lo = pe$conf.int[1], pearson_hi = pe$conf.int[2],
    pearson_p = pe$p.value,
    xi_xy = xy["mean"], xi_xy_sd = xy["sd"], xi_xy_p = xy["p_med"], xi_xy_flip = xy["flip"],
    xi_yx = yx["mean"], xi_yx_sd = yx["sd"], xi_yx_p = yx["p_med"], xi_yx_flip = yx["flip"],
    dcor = unname(dc$estimates["dCor"]), dcor_p = dc$p.value,
    tau = if (is.null(tau)) NA_real_ else tau$tStar,
    tau_p = if (is.null(tau)) NA_real_ else as.numeric(tau$pVal),
    tau_mode = if (is.null(tau)) NA_character_ else tau$mode, tau_warn = nw > 0,
    row.names = NULL)
}

mine <- pairs[pairs$index %% nsh == shard %% nsh, ]
mine <- head(mine, max_pairs)
rows <- list(); t0 <- Sys.time()
for (k in seq_len(nrow(mine))) {
  p <- mine[k, ]
  d <- data_of[[p$dataset]][, c(p$a, p$b)]; d <- d[complete.cases(d), ]
  set.seed(20263000L + p$index)
  rows[[k]] <- cbind(p, analyze_pair(d[[1]], d[[2]]))
  cat(sprintf("[shard %d/%d] pair %d/%d (%s: %s ~ %s) done (%.1f min)\n", shard, nsh, k, nrow(mine),
              p$dataset, p$a, p$b, as.numeric(difftime(Sys.time(), t0, units = "mins")))); flush.console()
}
res <- do.call(rbind, rows)
write.csv(res, file.path(out_dir, sprintf("screen_%02d_of_%02d.csv", shard, nsh)), row.names = FALSE)
writeLines(c(R.version.string, paste("XICOR", as.character(packageVersion("XICOR"))),
             paste("energy", as.character(packageVersion("energy"))),
             paste("TauStar", as.character(packageVersion("TauStar"))),
             paste("psych", as.character(packageVersion("psych"))),
             paste("B_PERM", B_PERM, "B_TB", B_TB, "seed_base 20263000")),
           file.path(out_dir, sprintf("screen_%02d_of_%02d_env.txt", shard, nsh)))
cat("DONE shard", shard, "of", nsh, "->", out_dir, "\n")
