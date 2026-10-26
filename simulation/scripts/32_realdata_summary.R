# Summary of the real-data screening (scripts 30 and 31) under the PRE-SPECIFIED decision rule
# (see the header of 30_realdata_screen.R):
#   * a measure "detects" an association in a pair if its Benjamini-Hochberg-adjusted p-value,
#     computed within each dataset and measure across all its pairs, is < .05 (unadjusted .05 also
#     reported);
#   * xi's p-value = min(1, 2 * min(median p of xi(X,Y), median p of xi(Y,X)));
#   * "non-monotone-only" pair: detected by at least one of xi, dCor, tau* and by neither
#     Spearman nor Pearson; "monotone-only": detected by Spearman or Pearson and by none of
#     xi, dCor, tau*;
#   * tie-break instability: share of pairs in which xi's decision (p < .05) differs across the 20
#     random tie-breaks (strictly between 0 and 1 of the tie-breaks, or between .25 and .75).
# Writes CSV tables to results/realdata_summary/.
#
# NOTE: run from the simulation/ directory:  cd simulation && Rscript scripts/32_realdata_summary.R
rd <- do.call(rbind, lapply(list.files("results/realdata", "^screen_.*[0-9]\\.csv$", full.names = TRUE), read.csv))
nh <- do.call(rbind, lapply(list.files("results/nhanes", "^screen_.*[0-9]\\.csv$", full.names = TRUE), read.csv))
d <- rbind(rd, nh)
d <- d[order(d$dataset, d$index), ]
stopifnot(nrow(d) == 518 + 127, !anyDuplicated(d[, c("dataset", "index")]))

M <- c("spearman", "pearson", "xi", "dcor", "tau")
d$p_spearman <- d$spearman_p; d$p_pearson <- d$pearson_p
d$p_xi <- pmin(1, 2 * pmin(d$xi_xy_p, d$xi_yx_p)); d$p_dcor <- d$dcor_p; d$p_tau <- d$tau_p
for (m in M) {
  d[[paste0("padj_", m)]] <- ave(d[[paste0("p_", m)]], d$dataset, FUN = function(p) p.adjust(p, "BH"))
  d[[paste0("det_", m)]] <- d[[paste0("padj_", m)]] < 0.05
  d[[paste0("raw_", m)]] <- d[[paste0("p_", m)]] < 0.05
}
mono <- d$det_spearman | d$det_pearson; gen <- d$det_xi | d$det_dcor | d$det_tau
d$cat <- ifelse(gen & !mono, "non-monotone-only", ifelse(mono & !gen, "monotone-only",
        ifelse(mono & gen, "both", "none")))
dir.create("results/realdata_summary", showWarnings = FALSE, recursive = TRUE)

## 1. detection counts
cnt <- do.call(rbind, lapply(split(d, d$dataset), function(x)
  data.frame(dataset = x$dataset[1], pairs = nrow(x), measure = M,
             detected_bh = sapply(M, function(m) sum(x[[paste0("det_", m)]])),
             detected_unadjusted = sapply(M, function(m) sum(x[[paste0("raw_", m)]])))))
write.csv(cnt, "results/realdata_summary/detection_counts.csv", row.names = FALSE)

## 2. concordance between monotone and general measures
conc <- do.call(rbind, lapply(split(d, d$dataset), function(x) data.frame(
  dataset = x$dataset[1], pairs = nrow(x), table(factor(x$cat, c("both", "non-monotone-only", "monotone-only", "none"))) |> as.list(),
  xi_detects_but_spearman_pearson_do_not = sum(x$det_xi & !mono[d$dataset == x$dataset[1]]),
  spearman_detects_but_xi_does_not = sum(x$det_spearman & !x$det_xi), check.names = FALSE)))
write.csv(conc, "results/realdata_summary/concordance.csv", row.names = FALSE)

## 3. pairs detected only by the general measures
nm <- d[d$cat == "non-monotone-only",
        c("dataset", "a", "b", "n", "spearman", "pearson", "xi_xy", "xi_yx", "dcor", "tau",
          "padj_xi", "padj_dcor", "padj_tau", "det_xi", "det_dcor", "det_tau")]
nm <- nm[order(nm$dataset, -nm$dcor), ]
write.csv(nm, "results/realdata_summary/non_monotone_only_pairs.csv", row.names = FALSE)

## 4. tie-break instability of xi's decision
ti <- do.call(rbind, lapply(split(d, d$dataset), function(x) {
  fl <- pmax(pmin(x$xi_xy_flip, 1 - x$xi_xy_flip), pmin(x$xi_yx_flip, 1 - x$xi_yx_flip))  # distance from a stable 0/1
  sd_max <- pmax(x$xi_xy_sd, x$xi_yx_sd)
  data.frame(dataset = x$dataset[1], pairs = nrow(x),
             share_decision_not_constant = mean(fl > 0), share_decision_unstable_25_75 = mean(fl >= 0.25),
             median_sd_xi = median(sd_max), p90_sd_xi = quantile(sd_max, .9), max_sd_xi = max(sd_max),
             median_abs_xi = median(pmax(abs(x$xi_xy), abs(x$xi_yx))))
}))
write.csv(ti, "results/realdata_summary/tiebreak_instability.csv", row.names = FALSE)

## 5. key NHANES pairs (sleep, BMI, age, income, gender vs PHQ-9 total)
key <- d[d$dataset == "nhanes_2017_2018" &
           ((d$a %in% c("sleep_weekday", "sleep_weekend", "bmi", "age", "income_ratio", "gender") & d$b == "phq_total") |
              (d$b %in% c("sleep_weekday", "sleep_weekend", "bmi", "age", "income_ratio", "gender") & d$a == "phq_total")),
         c("a", "b", "n", "spearman", "spearman_p", "pearson", "pearson_p", "xi_xy", "xi_yx", "p_xi",
           "dcor", "dcor_p", "tau", "tau_p", "det_spearman", "det_pearson", "det_xi", "det_dcor", "det_tau")]
write.csv(key, "results/realdata_summary/nhanes_key_pairs.csv", row.names = FALSE)

options(width = 140)
cat("\n== detection counts (BH FDR 5% / unadjusted .05) ==\n"); print(cnt, row.names = FALSE)
cat("\n== concordance ==\n"); print(conc, row.names = FALSE)
cat("\n== tie-break instability of xi ==\n"); print(ti, row.names = FALSE, digits = 3)
cat("\n== NHANES key pairs ==\n"); print(key, row.names = FALSE, digits = 3)
