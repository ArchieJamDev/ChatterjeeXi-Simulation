# Tables and figures of the BRM manuscript (paper_brm/), generated from the result CSV files only,
# so that no number is typed by hand. Tables are written as bare `tabular` environments to
# ../paper_brm/tables/ (the captions and notes live in the manuscript); figures as PDF to
# ../paper_brm/figures/.
#
# TABLE_LANG=es (default) or en selects the language of labels.
# NOTE: run from the simulation/ directory:  cd simulation && Rscript scripts/40_make_tables.R
# Figure 3 needs the NHANES files; they are downloaded (anonymously) into data/nhanes/ if missing.
LANG <- Sys.getenv("TABLE_LANG", "es")
L <- function(es, en) if (LANG == "es") es else en
out_tab <- "../paper_brm/tables"; out_fig <- "../paper_brm/figures"
dir.create(out_tab, FALSE, TRUE); dir.create(out_fig, FALSE, TRUE)

## ---- helpers -----------------------------------------------------------------
fmt <- function(x, d = 3) { s <- sprintf(paste0("%.", d, "f"), x); s <- sub("^(-?)0\\.", "\\1.", s); sub("^-", "$-$", s) }
pfmt <- function(p) ifelse(is.na(p), "--", ifelse(p < 0.001, "$<$ .001", fmt(p)))
read_dir <- function(dir, pattern) do.call(rbind, lapply(list.files(dir, pattern, full.names = TRUE), read.csv))
write_tab <- function(name, lines) writeLines(lines, file.path(out_tab, paste0(name, ".tex")))
row_tex <- function(cells) paste0(paste(cells, collapse = " & "), " \\\\")
n80 <- function(ns, ps, top) {
  if (ps[1] >= 0.8) return(paste0("$<$", ns[1]))
  for (i in 2:length(ns)) if (ps[i] >= 0.8) {
    x0 <- log(ns[i - 1]); x1 <- log(ns[i])
    return(as.character(round(exp(x0 + (0.8 - ps[i - 1]) * (x1 - x0) / (ps[i] - ps[i - 1])))))
  }
  paste0("$>$", top)
}
YLAB <- c(cont = L("Continua", "Continuous"), bin50 = L("Binaria 50/50", "Binary 50/50"),
          bin20 = L("Binaria 80/20", "Binary 80/20"), bin10 = L("Binaria 90/10", "Binary 90/10"),
          bin05 = L("Binaria 95/5", "Binary 95/5"),
          ord5s = L("Ordinal, 5 cat. (asim\\'etrica)", "Ordinal, 5 categories (skewed)"),
          ord7f = L("Ordinal, 7 cat. (efecto suelo)", "Ordinal, 7 categories (floor effect)"),
          comp5 = L("Compuesto de 5 \\'items", "5-item composite"))
XLAB <- c(cont = L("Continua", "Continuous"), bin50 = L("Binaria 50/50", "Binary 50/50"),
          ord5s = L("Ordinal, 5 cat.", "Ordinal, 5 cat."))
SLAB <- c(lineal = L("Lineal", "Linear"), cuadratica = L("Cuadr\\'atica", "Quadratic"),
          senoidal = L("Sinusoidal", "Sinusoidal"), circular = L("Circular", "Circular"))
YPLAIN <- c(cont = L("Continua", "Continuous"), bin50 = L("Binaria 50/50", "Binary 50/50"), bin10 = L("Binaria 90/10", "Binary 90/10"),
            ord5s = L("Ordinal 5 cat.", "Ordinal 5 cat."), comp5 = L("Compuesto 5 \u00edtems", "5-item composite"))
SPLAIN <- c(lineal = L("Lineal", "Linear"), cuadratica = L("Cuadr\u00e1tica", "Quadratic"), senoidal = "Sinusoidal", circular = "Circular")
MEAS <- c("xi", "spearman", "pearson", "dcor", "tau")
MLAB <- c(xi = "$\\xi_n$", spearman = "Spearman", pearson = "Pearson", dcor = "dCor", tau = "$\\tau^*$")
MCOL <- c(xi = "#D55E00", spearman = "#0072B2", pearson = "#56B4E9", dcor = "#009E73", tau = "#CC79A7")
MPCH <- c(xi = 16, spearman = 15, pearson = 0, dcor = 17, tau = 2)
MPLAIN <- c(xi = "xi_n", spearman = "Spearman", pearson = "Pearson", dcor = "dCor", tau = "tau*")

s1_I <- read_dir("results/study1", "^typeI_asym_.*\\.csv$")
s1_P <- read_dir("results/study1", "^power_asym_.*\\.csv$")
s2_I <- read_dir("results/study2", "^typeI_[0-9]+\\.csv$")
s2_c <- read_dir("results/study2", "^power_cheap_.*\\.csv$")
s2_f <- read_dir("results/study2", "^power_full_.*\\.csv$")
s2_h <- read_dir("results/study2", "^power_circ_hi_.*\\.csv$")
stopifnot(nrow(s1_I) == 120, nrow(s1_P) == 360, nrow(s2_I) == 20, nrow(s2_c) == 400,
          nrow(s2_f) == 128, nrow(s2_h) == 48)

## ---- Table: Study 1, type I error of the xi test ---------------------------------
NS1 <- c(50, 100, 250, 500, 1000); YORD <- c("cont", "bin50", "bin20", "bin10", "bin05", "ord5s", "ord7f", "comp5")
d <- s1_I[s1_I$xtype == "cont", ]
lines <- c("\\begin{tabular}{lccccccc}", "\\toprule",
  paste0(" & \\multicolumn{5}{c}{", L("Prueba con ajuste por empates, $\\alpha = .05$", "Tie-aware test, $\\alpha = .05$"),
         "} & ", L("Ingenua", "Naive"), " & ", L("Con ajuste", "Tie-aware"), " \\\\"),
  "\\cmidrule(lr){2-6}",
  row_tex(c(L("Respuesta $Y$", "Response $Y$"), paste0("$n=", NS1, "$"), L("$\\alpha = .05$", "$\\alpha = .05$"), "$\\alpha = .01$")),
  "\\midrule")
for (y in YORD) {
  r <- d[d$ytype == y, ]; r <- r[order(r$n), ]
  lines <- c(lines, row_tex(c(YLAB[y], fmt(r$rej05), fmt(mean(r$rej05_naive)), fmt(mean(r$rej01)))))
}
write_tab("s1_typeI", c(lines, "\\bottomrule", "\\end{tabular}"))
# check used in the text: is the pattern the same when X is discretized?
chk <- aggregate(rej05 ~ ytype + xtype, s1_I, mean)
cat("\n[check] Study 1 type I, pooled over n, by Y and X type:\n"); print(reshape(chk, idvar = "ytype", timevar = "xtype", direction = "wide"), digits = 3)

## ---- Table: Study 1, power when X is discretized ---------------------------------
d <- s1_P[s1_P$ytype == "cont" & s1_P$n %in% c(100, 250), ]
lines <- c("\\begin{tabular}{lcccccc}", "\\toprule",
  paste0(" & \\multicolumn{3}{c}{$n = 100$} & \\multicolumn{3}{c}{$n = 250$} \\\\"),
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
  row_tex(c(L("Estructura", "Structure"), rep(XLAB[c("cont", "bin50", "ord5s")], 2))), "\\midrule")
for (st in c("lineal", "cuadratica", "senoidal")) {
  cells <- c()
  for (n in c(100, 250)) for (x in c("cont", "bin50", "ord5s"))
    cells <- c(cells, fmt(d$rej05[d$structure == st & d$n == n & d$xtype == x]))
  lines <- c(lines, row_tex(c(SLAB[st], cells)))
}
write_tab("s1_Xdisc", c(lines, "\\bottomrule", "\\end{tabular}"))

## ---- Table: tie-breaking variance (earlier modular simulations, script 03) --------
tb <- read.csv("results/03_h2b.csv")
lab <- c(Continua = L("Continua", "Continuous"), Binaria = L("Binaria", "Binary"), `Likert-3` = "Likert-3",
         `Likert-5` = "Likert-5", `Likert-7` = "Likert-7")
lines <- c("\\begin{tabular}{lccc}", "\\toprule",
  row_tex(c(L("$X$ discretizada en", "$X$ discretized into"), L("Entre conjuntos de datos", "Between datasets"),
            L("Por desempate", "Tie-breaking"), L("\\% por desempate", "\\% tie-breaking"))), "\\midrule")
for (i in seq_len(nrow(tb)))
  lines <- c(lines, row_tex(c(lab[tb$condicion[i]], fmt(tb$var_entre[i], 4), fmt(tb$var_dentro[i], 4),
                              sprintf("%.0f", 100 * tb$var_dentro[i] / (tb$var_entre[i] + tb$var_dentro[i])))))
write_tab("s1_tiebreak", c(lines, "\\bottomrule", "\\end{tabular}"))
kv <- read.csv("results/06_h2d_kvarying.csv"); h2d <- read.csv("results/05_h2d.csv")
cat("\n[check] composite scores: variance by k:\n"); print(kv[, c("k", "varianza_empirica")], digits = 3)
cat("[check] single item vs composite k=10:\n"); print(h2d, digits = 3)

## ---- Table: Study 2, type I error of the five tests -----------------------------
Y2 <- c("cont", "bin50", "bin10", "ord5s", "comp5")
lines <- c("\\begin{tabular}{lcccccc}", "\\toprule",
  row_tex(c(L("Respuesta $Y$", "Response $Y$"), MLAB[MEAS], L("Avisos de $\\tau^*$", "$\\tau^*$ warnings"))), "\\midrule",
  paste0("\\multicolumn{7}{l}{\\emph{", L("Nivel de significaci\\'on .05", "Significance level .05"), "}} \\\\"))
for (y in Y2) { r <- s2_I[s2_I$ytype == y, ]
  lines <- c(lines, row_tex(c(YLAB[y], fmt(sapply(MEAS, function(m) mean(r[[paste0("rej05_", m)]]))),
                              sprintf("%.0f\\%%", 100 * mean(r$warn_rate_tau)))))}
lines <- c(lines, "\\midrule", paste0("\\multicolumn{7}{l}{\\emph{", L("Nivel de significaci\\'on .01", "Significance level .01"), "}} \\\\"))
for (y in Y2) { r <- s2_I[s2_I$ytype == y, ]
  lines <- c(lines, row_tex(c(YLAB[y], fmt(sapply(MEAS, function(m) mean(r[[paste0("rej01_", m)]]))), "")))}
write_tab("s2_typeI", c(lines, "\\bottomrule", "\\end{tabular}"))

## ---- Table: Study 2, n for 80% power ---------------------------------------------
lines <- c("\\begin{tabular}{lllccccc}", "\\toprule",
  row_tex(c(L("Estructura", "Structure"), "$R^2$", L("Respuesta $Y$", "Response $Y$"), MLAB[MEAS])), "\\midrule")
first <- TRUE
for (st in c("lineal", "cuadratica", "senoidal")) {
  if (!first) lines <- c(lines, "\\addlinespace"); first <- FALSE
  for (s in c(0.13, 0.26)) for (y in c("cont", "ord5s")) {
    cells <- sapply(MEAS, function(m) {
      if (m %in% c("xi", "spearman", "pearson")) {
        r <- s2_c[s2_c$structure == st & s2_c$ytype == y & abs(s2_c$share - s) < 1e-9, ]; r <- r[order(r$n), ]
        n80(r$n, r[[paste0("rej05_", m)]], 800)
      } else {
        r <- s2_f[s2_f$structure == st & s2_f$ytype == y & abs(s2_f$share - s) < 1e-9, ]; r <- r[order(r$n), ]
        n80(r$n, r[[paste0("rej05_", m)]], 400)
      }
    })
    lines <- c(lines, row_tex(c(if (s == 0.13 && y == "cont") SLAB[st] else "", fmt(s, 2), YLAB[y], cells)))
  }
}
write_tab("s2_n80", c(lines, "\\bottomrule", "\\end{tabular}"))

## ---- Table: Study 2, circular structure at high signal ----------------------------
lines <- c("\\begin{tabular}{llccccc}", "\\toprule", row_tex(c("$R^2$", "$n$", MLAB[MEAS])), "\\midrule")
for (s in c(0.70, 0.85, 0.95)) for (n in c(100, 200, 400)) {
  r <- s2_h[s2_h$ytype == "cont" & abs(s2_h$share - s) < 1e-9 & s2_h$n == n, ]
  lines <- c(lines, row_tex(c(if (n == 100) fmt(s, 2) else "", n, fmt(sapply(MEAS, function(m) r[[paste0("rej05_", m)]])))))
}
write_tab("s2_circular", c(lines, "\\bottomrule", "\\end{tabular}"))

## ---- Study 3 tables ------------------------------------------------------------------
dc <- read.csv("results/realdata_summary/detection_counts.csv"); cn <- read.csv("results/realdata_summary/concordance.csv", check.names = FALSE)
DS <- c("bfi", "nhanes_2017_2018", "sat.act"); DSL <- c(bfi = "\\texttt{bfi}", nhanes_2017_2018 = "NHANES", sat.act = "\\texttt{sat.act}")
cellv <- function(ds, m) dc$detected_bh[dc$dataset == ds & dc$measure == m]
lines <- c("\\begin{tabular}{lccc}", "\\toprule", row_tex(c("", DSL[DS])), "\\midrule",
  row_tex(c(L("Pares analizados", "Pairs analyzed"), sapply(DS, function(ds) cn$pairs[cn$dataset == ds]))), "\\midrule",
  paste0("\\multicolumn{4}{l}{\\emph{", L("Pares detectados por cada medida", "Pairs detected by each measure"), "}} \\\\"))
for (m in c("spearman", "pearson", "xi", "dcor", "tau")) lines <- c(lines, row_tex(c(MLAB[m], sapply(DS, function(ds) cellv(ds, m)))))
lines <- c(lines, "\\midrule",
  row_tex(c(L("Solo por medidas generales ($\\xi_n$, dCor o $\\tau^*$)", "Only by general measures ($\\xi_n$, dCor or $\\tau^*$)"),
            sapply(DS, function(ds) cn[["non-monotone-only"]][cn$dataset == ds]))),
  row_tex(c(L("Solo por Spearman o Pearson", "Only by Spearman or Pearson"), sapply(DS, function(ds) cn[["monotone-only"]][cn$dataset == ds]))),
  row_tex(c(L("$\\xi_n$ detecta y Spearman y Pearson no", "$\\xi_n$ detects and Spearman and Pearson do not"),
            sapply(DS, function(ds) cn$xi_detects_but_spearman_pearson_do_not[cn$dataset == ds]))))
write_tab("s3_detection", c(lines, "\\bottomrule", "\\end{tabular}"))

kp <- read.csv("results/realdata_summary/nhanes_key_pairs.csv")
other <- ifelse(kp$a == "phq_total", kp$b, kp$a)
VLAB <- c(sleep_weekday = L("Sue\\~no entre semana (h)", "Weekday sleep (h)"), sleep_weekend = L("Sue\\~no en fin de semana (h)", "Weekend sleep (h)"),
          bmi = L("\\'Indice de masa corporal", "Body-mass index"), age = L("Edad", "Age"),
          income_ratio = L("Raz\\'on ingreso/pobreza", "Income-to-poverty ratio"), gender = L("G\\'enero", "Gender"))
ord <- c("sleep_weekday", "sleep_weekend", "bmi", "age", "income_ratio", "gender")
dag <- function(det) ifelse(det, "$^{\\dagger}$", "")
lines <- c("\\begin{tabular}{lrcccc}", "\\toprule",
  row_tex(c(paste0(L("Variable (frente a PHQ-9 total)", "Variable (vs. total PHQ-9)")), "$n$", "Spearman $\\rho$ ($p$)", "$\\xi_n$ ($p$)", "dCor ($p$)", "$\\tau^*$ ($p$)")), "\\midrule")
for (v in ord) { r <- kp[other == v, ]
  xi_val <- max(r$xi_xy, r$xi_yx)
  lines <- c(lines, row_tex(c(VLAB[v], r$n,
    paste0(fmt(r$spearman), " (", pfmt(r$spearman_p), ")", dag(r$det_spearman)),
    paste0(fmt(xi_val), " (", pfmt(r$p_xi), ")", dag(r$det_xi)),
    paste0(fmt(r$dcor), " (", pfmt(r$dcor_p), ")", dag(r$det_dcor)),
    paste0(fmt(r$tau, 4), " (", pfmt(r$tau_p), ")", dag(r$det_tau)))))
}
write_tab("s3_nhanes_key", c(lines, "\\bottomrule", "\\end{tabular}"))

ti <- read.csv("results/realdata_summary/tiebreak_instability.csv")
lines <- c("\\begin{tabular}{lcccc}", "\\toprule",
  row_tex(c("", L("Pares", "Pairs"), L("DE de $\\xi_n$ entre desempates (mediana)", "SD of $\\xi_n$ across tie-breaks (median)"),
             L("$|\\xi_n|$ (mediana)", "$|\\xi_n|$ (median)"), L("Decisi\\'on inestable (\\%)", "Unstable decision (\\%)"))), "\\midrule")
for (ds in DS) { r <- ti[ti$dataset == ds, ]
  lines <- c(lines, row_tex(c(DSL[ds], r$pairs, fmt(r$median_sd_xi, 4), fmt(r$median_abs_xi, 4), sprintf("%.0f", 100 * r$share_decision_unstable_25_75))))}
write_tab("s3_tiebreak", c(lines, "\\bottomrule", "\\end{tabular}"))

## ---- Figures ----------------------------------------------------------------------------
pdf_fig <- function(name, w, h, expr) { pdf(file.path(out_fig, paste0(name, ".pdf")), width = w, height = h, family = "Helvetica"); on.exit(dev.off()); expr() }

## Figure 1: Study 1, power by structure and type of response (X continuous)
pdf_fig(paste0("fig1_study1_", LANG), 6.6, 3.0, function() {
  par(mfrow = c(1, 3), mar = c(3.4, 3.4, 2, 0.6), oma = c(2.4, 0, 0, 0), mgp = c(2.1, 0.6, 0), cex.axis = 0.8, cex.lab = 0.9)
  yt <- c("cont", "bin50", "bin10", "ord5s", "comp5")
  col <- c(cont = "black", bin50 = "#D55E00", bin10 = "#CC79A7", ord5s = "#0072B2", comp5 = "#009E73")
  pch <- c(cont = 16, bin50 = 15, bin10 = 17, ord5s = 1, comp5 = 2)
  for (st in c("lineal", "cuadratica", "senoidal")) {
    plot(NA, xlim = range(NS1), ylim = c(0, 1), log = "x", xlab = "n", ylab = L("Potencia", "Power"), main = SPLAIN[st], cex.main = 0.95, xaxt = "n")
    axis(1, NS1, NS1); abline(h = 0.05, lty = 3, col = "grey60")
    for (y in yt) { r <- s1_P[s1_P$structure == st & s1_P$xtype == "cont" & s1_P$ytype == y, ]; r <- r[order(r$n), ]
      lines(r$n, r$rej05, col = col[y], type = "b", pch = pch[y], lty = 1) }
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend("bottom", legend = YPLAIN[yt], col = col, pch = pch, lty = 1, horiz = TRUE, bty = "n", cex = 0.7, xpd = NA)
})

## Figure 2: Study 2, power of the five measures (R^2 = .26; circular: R^2 = .95)
pdf_fig(paste0("fig2_study2_", LANG), 6.8, 4.6, function() {
  par(mfrow = c(2, 4), mar = c(3.2, 4.0, 2.0, 0.5), oma = c(2.6, 0, 0, 0), mgp = c(2.0, 0.6, 0), cex.axis = 0.75, cex.lab = 0.85)
  for (y in c("cont", "ord5s")) for (st in c("lineal", "cuadratica", "senoidal", "circular")) {
    src <- if (st == "circular") s2_h[abs(s2_h$share - 0.95) < 1e-9, ] else s2_f[abs(s2_f$share - 0.26) < 1e-9, ]
    src <- src[src$structure == st & src$ytype == y, ]
    ttl <- paste0(SPLAIN[st], if (st == "circular") " (R² = .95)" else " (R² = .26)")
    rowlab <- if (y == "cont") L("Y continua", "Y continuous") else L("Y ordinal, 5 cat.", "Y ordinal, 5 cat.")
    plot(NA, xlim = c(50, 400), ylim = c(0, 1), log = "x", xlab = "n",
         ylab = if (st == "lineal") paste0(rowlab, "\n", L("Potencia", "Power")) else "",
         main = if (y == "cont") ttl else "", cex.main = 0.85, xaxt = "n")
    axis(1, c(50, 100, 200, 400), c(50, 100, 200, 400)); abline(h = 0.05, lty = 3, col = "grey60")
    for (m in MEAS) { r <- src[order(src$n), ]; lines(r$n, r[[paste0("rej05_", m)]], col = MCOL[m], pch = MPCH[m], type = "b", lty = 1) }
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend("bottom", legend = c(expression(xi[n]), "Spearman", "Pearson", "dCor", expression(tau^"*")), col = MCOL[MEAS], pch = MPCH[MEAS], lty = 1, horiz = TRUE, bty = "n", cex = 0.8, xpd = NA)
})

## Figure 3: NHANES, mean PHQ-9 total by weekday sleep hours and by BMI
suppressPackageStartupMessages(library(foreign))
data_dir <- "data/nhanes"; dir.create(data_dir, FALSE, TRUE)
get_xpt <- function(nm) { f <- file.path(data_dir, paste0(nm, ".xpt"))
  if (!file.exists(f) || file.size(f) < 1e4) download.file(paste0("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/", nm, ".xpt"), f, mode = "wb", quiet = TRUE)
  foreign::read.xport(f) }
dm <- get_xpt("DEMO_J"); dq <- get_xpt("DPQ_J"); sl <- get_xpt("SLQ_J"); bm <- get_xpt("BMX_J")
nh <- Reduce(function(a, b) merge(a, b, by = "SEQN"), list(dm[, c("SEQN", "RIDAGEYR")], dq, sl[, c("SEQN", "SLD012")], bm[, c("SEQN", "BMXBMI")]))
nh <- nh[nh$RIDAGEYR >= 18, ]; items <- paste0("DPQ0", c(10, 20, 30, 40, 50, 60, 70, 80, 90))
for (v in items) nh[[v]][nh[[v]] %in% c(7, 9)] <- NA
nh$phq <- ifelse(rowSums(is.na(nh[, items])) == 0, rowSums(nh[, items]), NA)
band_stats <- function(x, y, breaks, labels) {
  g <- cut(x, breaks, right = FALSE, labels = labels); ok <- !is.na(g) & !is.na(y)
  m <- tapply(y[ok], g[ok], mean); s <- tapply(y[ok], g[ok], sd); k <- tapply(y[ok], g[ok], length)
  data.frame(band = names(m), mean = as.numeric(m), lo = as.numeric(m - 1.96 * s / sqrt(k)), hi = as.numeric(m + 1.96 * s / sqrt(k)), n = as.numeric(k)) }
sl_b <- band_stats(nh$SLD012, nh$phq, c(0, 5, 6, 7, 8, 9, 10, 20), c("<5", "5-6", "6-7", "7-8", "8-9", "9-10", ">=10"))
bm_b <- band_stats(nh$BMXBMI, nh$phq, c(0, 20, 25, 30, 35, 40, 200), c("<20", "20-25", "25-30", "30-35", "35-40", ">=40"))
write.csv(sl_b, file.path(out_tab, "fig3_sleep_bands.csv"), row.names = FALSE); write.csv(bm_b, file.path(out_tab, "fig3_bmi_bands.csv"), row.names = FALSE)
pdf_fig(paste0("fig3_nhanes_", LANG), 6.4, 2.9, function() {
  par(mfrow = c(1, 2), mar = c(3.6, 3.6, 2, 0.8), mgp = c(2.2, 0.6, 0), cex.axis = 0.8, cex.lab = 0.9)
  for (k in 1:2) { b <- if (k == 1) sl_b else bm_b
    plot(seq_len(nrow(b)), b$mean, ylim = c(2, 6), xaxt = "n", pch = 16, type = "b", xlab = if (k == 1) L("Horas de sue\u00f1o entre semana", "Weekday sleep (hours)") else L("\u00cdndice de masa corporal", "Body-mass index"),
         ylab = L("PHQ-9 total (media)", "Total PHQ-9 (mean)"), main = if (k == 1) L("a) Sue\u00f1o", "a) Sleep") else "b) IMC / BMI", cex.main = 0.95)
    axis(1, seq_len(nrow(b)), b$band, cex.axis = 0.75, gap.axis = -1); arrows(seq_len(nrow(b)), b$lo, seq_len(nrow(b)), b$hi, angle = 90, code = 3, length = 0.03) }
})
cat("\nDone. Tables:", paste(list.files(out_tab, "\\.tex$"), collapse = ", "), "\nFigures:", paste(list.files(out_fig, "\\.pdf$"), collapse = ", "), "\n")
