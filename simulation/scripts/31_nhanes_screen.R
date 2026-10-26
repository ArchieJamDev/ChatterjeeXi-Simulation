# Real-data screening, second dataset: NHANES 2017-2018 (US CDC / NCHS), public-use files,
# downloaded anonymously (no registration) from
#   https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/{DEMO_J,DPQ_J,SLQ_J,BMX_J}.xpt
# Same per-pair analysis as 30_realdata_screen.R (xi_n averaged over random tie-breaks, Spearman,
# Pearson, dCor permutation test with B = 499, tau*), applied to EVERY pair of a variable set
# fixed in advance (adults aged >= 18, complete cases per pair):
#   PHQ-9 items 1-9 (DPQ010-DPQ090, 0-3; codes 7 and 9 -> NA), PHQ-9 total (0-27; NA if any item
#   is missing), weekday and weekend sleep hours (SLD012, SLD013), age (RIDAGEYR), gender
#   (RIAGENDR), body-mass index (BMXBMI), family income-to-poverty ratio (INDFMPIR), education
#   (DMDEDUC2, codes 7 and 9 -> NA) = 17 variables, 136 pairs, minus the 9 pairs (item, PHQ-9
#   total) -> 127 pairs.
# The sleep-depression pair was expected a priori to be non-monotone (a U shape: both short and
# long sleep go with more depressive symptoms); it is one of the 127 pairs, not selected by result.
# CAVEAT: the analysis is UNWEIGHTED (NHANES survey weights and design are ignored), so estimates
# describe the sample, not the US population.
#
# Decision rules, outputs and the treatment of xi's confidence interval are the same as in
# 30_realdata_screen.R. Each pair is seeded from its own index (20264000 + index).
#
# NOTE: run from the simulation/ directory, e.g.
#   cd simulation && Rscript scripts/31_nhanes_screen.R 1 12     # shard 1 of 12
# NHANES_SMOKE=1 -> tiny run (4 pairs, B = 19, 3 tie-breaks); NHANES_OUT sets the output
# directory; NHANES_DATA sets the download cache (default data/nhanes, git-ignored).
suppressPackageStartupMessages({library(XICOR); library(energy); library(TauStar); library(foreign)})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
shard <- as.integer(args[1]); nsh <- as.integer(args[2])
stopifnot(shard >= 1, shard <= nsh)
smoke <- identical(Sys.getenv("NHANES_SMOKE"), "1")
B_PERM <- if (smoke) 19L else 499L
B_TB <- if (smoke) 3L else 20L
out_dir <- Sys.getenv("NHANES_OUT", unset = if (smoke) tempdir() else "results/nhanes")
data_dir <- Sys.getenv("NHANES_DATA", unset = "data/nhanes")
max_pairs <- if (smoke) 4L else Inf
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

## ---- data (anonymous download from the CDC) ---------------------------------------
BASE <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/"
get_xpt <- function(name) {
  f <- file.path(data_dir, paste0(name, ".xpt"))
  if (!file.exists(f) || file.size(f) < 1e4) {
    for (try in 1:4) {
      ok <- tryCatch({ download.file(paste0(BASE, name, ".xpt"), f, mode = "wb", quiet = TRUE); TRUE },
                     error = function(e) FALSE)
      if (ok && file.size(f) > 1e4) break
      Sys.sleep(5 * try)
    }
  }
  stopifnot(file.exists(f), file.size(f) > 1e4)
  foreign::read.xport(f)
}
demo <- get_xpt("DEMO_J"); dpq <- get_xpt("DPQ_J"); slq <- get_xpt("SLQ_J"); bmx <- get_xpt("BMX_J")
d <- Reduce(function(a, b) merge(a, b, by = "SEQN"),
            list(demo[, c("SEQN", "RIDAGEYR", "RIAGENDR", "INDFMPIR", "DMDEDUC2")], dpq,
                 slq[, c("SEQN", "SLD012", "SLD013")], bmx[, c("SEQN", "BMXBMI")]))
d <- d[d$RIDAGEYR >= 18, ]
items <- paste0("DPQ0", c(10, 20, 30, 40, 50, 60, 70, 80, 90))
for (v in items) d[[v]][d[[v]] %in% c(7, 9)] <- NA
d$DMDEDUC2[d$DMDEDUC2 %in% c(7, 9)] <- NA
vv <- data.frame(age = d$RIDAGEYR, gender = d$RIAGENDR, bmi = d$BMXBMI, income_ratio = d$INDFMPIR,
                 education = d$DMDEDUC2, sleep_weekday = d$SLD012, sleep_weekend = d$SLD013)
for (i in 1:9) vv[[paste0("phq", i)]] <- d[[items[i]]]
vv$phq_total <- ifelse(rowSums(is.na(d[, items])) == 0, rowSums(d[, items]), NA)
vars <- c(paste0("phq", 1:9), "phq_total", "sleep_weekday", "sleep_weekend", "age", "gender", "bmi",
          "income_ratio", "education")
pr <- t(combn(vars, 2))
pr <- pr[!(grepl("^phq[1-9]$", pr[, 1]) & pr[, 2] == "phq_total") &
           !(grepl("^phq[1-9]$", pr[, 2]) & pr[, 1] == "phq_total"), , drop = FALSE]
pairs <- data.frame(dataset = "nhanes_2017_2018", a = pr[, 1], b = pr[, 2], index = seq_len(nrow(pr)))
stopifnot(nrow(pairs) == 127)

## ---- one pair (identical to 30_realdata_screen.R) ---------------------------------
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
  dd <- vv[, c(p$a, p$b)]; dd <- dd[complete.cases(dd), ]
  set.seed(20264000L + p$index)
  rows[[k]] <- cbind(p, analyze_pair(dd[[1]], dd[[2]]))
  cat(sprintf("[shard %d/%d] pair %d/%d (%s ~ %s) done (%.1f min)\n", shard, nsh, k, nrow(mine),
              p$a, p$b, as.numeric(difftime(Sys.time(), t0, units = "mins")))); flush.console()
}
res <- do.call(rbind, rows)
write.csv(res, file.path(out_dir, sprintf("screen_%02d_of_%02d.csv", shard, nsh)), row.names = FALSE)
fs <- file.path(data_dir, paste0(c("DEMO_J", "DPQ_J", "SLQ_J", "BMX_J"), ".xpt"))
writeLines(c(R.version.string, paste("XICOR", as.character(packageVersion("XICOR"))),
             paste("energy", as.character(packageVersion("energy"))),
             paste("TauStar", as.character(packageVersion("TauStar"))),
             paste("B_PERM", B_PERM, "B_TB", B_TB, "seed_base 20264000"),
             paste(basename(fs), file.size(fs), "bytes md5", tools::md5sum(fs))),
           file.path(out_dir, sprintf("screen_%02d_of_%02d_env.txt", shard, nsh)))
cat("DONE shard", shard, "of", nsh, "->", out_dir, "\n")
