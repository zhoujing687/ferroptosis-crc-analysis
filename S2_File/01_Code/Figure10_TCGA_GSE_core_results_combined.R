###############################################################################
# Figure 10: combined TCGA-COAD and GSE42387 core results
#
# A-C: TCGA-COAD survival and stage-related analyses
# D:   GSE42387 11-gene average resistant-versus-parental effects
# E:   GSE42387 P/B/P-B fixed-effects estimates
#
# Usage:
# Rscript Figure10_TCGA_GSE_core_results_combined.R <analysis_root> <outdir>
###############################################################################

set.seed(42)
options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: Rscript Figure10_TCGA_GSE_core_results_combined.R <analysis_root> <outdir>",
    call. = FALSE
  )
}

analysis_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
out_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("The 'survival' package is required.", call. = FALSE)
}
suppressPackageStartupMessages(library(survival))

tcga_dir <- file.path(analysis_root, "S1_File", "03_Source_Data")
gse_dir <- file.path(analysis_root, "S2_File", "03_Source_Data")
if (!dir.exists(tcga_dir)) tcga_dir <- file.path(analysis_root, "02_TCGA_COAD")
if (!dir.exists(gse_dir)) gse_dir <- file.path(analysis_root, "03_GSE42387")

read_source <- function(directory, filename) {
  path <- file.path(directory, filename)
  if (!file.exists(path)) stop("Missing source-data file: ", path, call. = FALSE)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

tcga_km <- read_source(tcga_dir, "TCGA_Figure10_panelA_KM_source_data.csv")
tcga_cox <- read_source(tcga_dir, "TCGA_Figure10_panelB_primary_Cox_source_data.csv")
tcga_stage <- read_source(tcga_dir, "TCGA_Figure10_panelC_P_B_PminusB_source_data.csv")
tcga_stage_tests <- read_source(tcga_dir, "TCGA_Figure10_panelC_stage_comparisons.csv")
gse_genes <- read_source(gse_dir, "GSE42387_11_genes_limma_main_contrasts.csv")
gse_scores <- read_source(gse_dir, "GSE42387_PBI_fixed_effects_contrasts.csv")

tcga_km$OS_time <- as.numeric(tcga_km$OS_time)
tcga_km$OS <- as.integer(tcga_km$OS)
tcga_km$index_group <- factor(tcga_km$index_group, levels = c("High", "Low"))
tcga_stage$stage <- factor(tcga_stage$stage, levels = c("II", "III"))

orange <- "#D4743C"
blue <- "#1874A6"
red <- "#F04B3A"
cyan <- "#43B8D3"
grid_col <- "#E8E8E8"
line_col <- "#555555"
pressure_genes <- c("ACSL4", "LPCAT3", "TFRC", "NCOA4")
target_genes <- c(
  pressure_genes, "SLC7A11", "GPX4", "AIFM2", "GCH1", "DHODH", "NQO1", "NFE2L2"
)

fmt_p <- function(p) {
  if (!length(p) || is.na(p[[1L]])) return("NA")
  if (p[[1L]] < 0.001) return("<0.001")
  sprintf("%.3f", p[[1L]])
}

panel_tag <- function(label, line = 0.4, adj = -0.02) {
  mtext(label, side = 3, line = line, adj = adj, font = 2, cex = 1.32)
}

step_coordinates <- function(time, value, start = 1) {
  if (!length(time)) return(list(x = 0, y = start))
  list(
    x = c(0, as.vector(rbind(time, time))),
    y = c(start, as.vector(rbind(c(start, value[-length(value)]), value)))
  )
}

group_fit <- function(group) {
  survival::survfit(
    survival::Surv(OS_time, OS) ~ 1,
    data = tcga_km[tcga_km$index_group == group, , drop = FALSE],
    conf.type = "log-log"
  )
}

km_fits <- setNames(lapply(c("High", "Low"), group_fit), c("High", "Low"))
logrank <- survival::survdiff(
  survival::Surv(OS_time, OS) ~ index_group,
  data = tcga_km
)
logrank_p <- stats::pchisq(logrank$chisq, df = length(logrank$n) - 1L, lower.tail = FALSE)

original_x_max <- 4500
original_x_ticks <- seq(0, original_x_max, by = 900)
km_palette <- c(High = red, Low = cyan)

draw_km_curve <- function() {
  par(mar = c(2.7, 4.8, 2.8, 0.5), family = "sans", lend = "round")
  plot(
    NA, xlim = c(0, original_x_max), ylim = c(0, 105), xaxs = "i",
    axes = FALSE, xlab = "", ylab = ""
  )
  abline(v = original_x_ticks, h = seq(0, 100, 25), col = grid_col, lwd = 0.7)
  axis(1, at = original_x_ticks, labels = original_x_ticks, cex.axis = 0.71)
  axis(2, at = seq(0, 100, 25), las = 1, cex.axis = 0.71)
  box(col = line_col, lwd = 0.7)
  mtext("Overall survival probability (%)", side = 2, line = 3.3, cex = 0.78)

  for (group in c("High", "Low")) {
    fit <- km_fits[[group]]
    lower <- step_coordinates(fit$time, 100 * fit$lower, 100)
    upper <- step_coordinates(fit$time, 100 * fit$upper, 100)
    polygon(
      c(lower$x, rev(upper$x)), c(lower$y, rev(upper$y)),
      col = adjustcolor(km_palette[[group]], alpha.f = 0.24), border = NA
    )
  }
  for (group in c("High", "Low")) {
    fit <- km_fits[[group]]
    curve <- step_coordinates(fit$time, 100 * fit$surv, 100)
    lines(curve$x, curve$y, col = km_palette[[group]], lwd = 1.55)
    censor_times <- sort(unique(tcga_km$OS_time[
      tcga_km$index_group == group & tcga_km$OS == 0 & tcga_km$OS_time <= original_x_max
    ]))
    if (length(censor_times)) {
      censor_surv <- summary(fit, times = censor_times, extend = TRUE)$surv * 100
      points(censor_times, censor_surv, pch = 3, cex = 0.55, lwd = 0.8,
             col = km_palette[[group]])
    }
  }
  text(250, 9, paste0("log-rank p = ", fmt_p(logrank_p)), adj = 0, cex = 0.75)
  legend(
    "topright", inset = c(0.01, 0.01), bty = "n", horiz = TRUE,
    legend = c("High Pᵀ−Bᵀ", "Low Pᵀ−Bᵀ"), col = km_palette[c("High", "Low")],
    lwd = 1.55, pch = 3, pt.cex = 0.65, cex = 0.68
  )
  panel_tag("A", line = 1.1)
}

draw_km_risk <- function() {
  par(mar = c(2.7, 4.8, 0.3, 0.5), family = "sans")
  plot(
    NA, xlim = c(0, original_x_max), ylim = c(0.45, 2.55), xaxs = "i",
    axes = FALSE, xlab = "", ylab = ""
  )
  rect(0, 0.55, original_x_max, 2.45, border = line_col, lwd = 0.65)
  abline(v = original_x_ticks, col = grid_col, lwd = 0.65)
  segments(0, 1.5, original_x_max, 1.5, col = grid_col, lwd = 0.65)
  axis(1, at = original_x_ticks, labels = original_x_ticks, cex.axis = 0.66)
  mtext("Time (days)", side = 1, line = 2.0, cex = 0.70)
  text(0, 2.67, "Number at risk", adj = 0, cex = 0.74, xpd = NA, font = 2)
  y_positions <- c(High = 1.98, Low = 1.02)
  for (group in names(y_positions)) {
    n_risk <- vapply(original_x_ticks, function(t) {
      sum(tcga_km$OS_time[tcga_km$index_group == group] >= t)
    }, integer(1L))
    text(-0.02 * original_x_max, y_positions[[group]], group, adj = 1,
         col = km_palette[[group]], cex = 0.67, xpd = NA)
    text(original_x_ticks, rep(y_positions[[group]], length(original_x_ticks)),
         n_risk, cex = 0.67, xpd = NA)
  }
}

draw_km_censor <- function() {
  par(mar = c(3.0, 4.8, 0.3, 0.5), family = "sans")
  censor_tables <- lapply(c("High", "Low"), function(group) {
    times <- tcga_km$OS_time[
      tcga_km$index_group == group & tcga_km$OS == 0 & tcga_km$OS_time <= original_x_max
    ]
    if (!length(times)) return(data.frame(time = numeric(), count = integer()))
    tab <- as.data.frame(table(times), stringsAsFactors = FALSE)
    data.frame(time = as.numeric(as.character(tab$times)), count = tab$Freq)
  })
  names(censor_tables) <- c("High", "Low")
  max_count <- max(1L, vapply(censor_tables, function(x) {
    if (nrow(x)) max(x$count) else 0L
  }, integer(1L)))
  plot(
    NA, xlim = c(0, original_x_max), ylim = c(0, max_count + 0.25), xaxs = "i",
    axes = FALSE, xlab = "", ylab = ""
  )
  abline(v = original_x_ticks, h = seq(0, max_count, 1), col = grid_col, lwd = 0.65)
  axis(1, at = original_x_ticks, labels = original_x_ticks, cex.axis = 0.64)
  axis(2, at = seq(0, max_count, 1), las = 1, cex.axis = 0.62)
  box(col = line_col, lwd = 0.65)
  mtext("Time (days)", side = 1, line = 2.0, cex = 0.69)
  mtext("Number censored", side = 2, line = 3.2, cex = 0.71)
  for (group in names(censor_tables)) {
    z <- censor_tables[[group]]
    if (nrow(z)) {
      offset <- if (group == "High") -5 else 5
      segments(z$time + offset, 0, z$time + offset, z$count,
               col = km_palette[[group]], lwd = 0.8)
    }
  }
}

draw_cox <- function() {
  par(mar = c(2.2, 0.4, 2.0, 0.3), family = "sans", xpd = NA)
  n_rows <- nrow(tcga_cox)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0.15, n_rows + 2.35))
  row_y <- rev(seq_len(n_rows))
  for (i in seq_len(n_rows)) {
    if (i %% 2L == 1L) {
      rect(0.01, row_y[[i]] - 0.42, 0.99, row_y[[i]] + 0.42,
           col = "#F1F1F1", border = NA)
    }
  }
  header_low <- n_rows + 1.48
  header_high <- n_rows + 2.22
  title_low <- n_rows + 0.68
  rect(0.01, header_low, 0.99, header_high, border = line_col, lwd = 0.65)
  rect(0.01, title_low, 0.99, header_low, border = line_col, lwd = 0.65)
  rect(0.01, 0.20, 0.99, title_low, border = line_col, lwd = 0.65)
  segments(c(0.24, 0.31, 0.67), header_low, c(0.24, 0.31, 0.67), header_high,
           col = line_col, lwd = 0.65)
  segments(c(0.24, 0.31, 0.67), 0.20, c(0.24, 0.31, 0.67), title_low,
           col = "#BBBBBB", lwd = 0.55)
  text(0.035, mean(c(header_low, header_high)), "Variable", adj = 0, cex = 0.66, font = 2)
  text(0.275, mean(c(header_low, header_high)), "N", cex = 0.66, font = 2)
  text(0.50, mean(c(header_low, header_high)), "Hazard ratio", cex = 0.66, font = 2)
  text(0.965, mean(c(header_low, header_high)), "p", adj = 1, cex = 0.66, font = 2)
  text(0.50, mean(c(title_low, header_low)), "Primary multivariable Cox model",
       cex = 0.66, font = 2)

  map_hr <- function(x) {
    0.325 + (log(x) - log(0.5)) / (log(5) - log(0.5)) * (0.635 - 0.325)
  }
  segments(map_hr(1), 0.34, map_hr(1), title_low - 0.10,
           lty = 2, col = line_col, lwd = 0.65)
  labels <- c(
    "Pressure score\n(per 1 SD)", "Buffering score\n(per 1 SD)",
    "Age\n(per 10 years)", "Stage\nIII vs II"
  )
  for (i in seq_len(n_rows)) {
    text(0.035, row_y[[i]], labels[[i]], adj = 0, cex = 0.59, font = 2)
    text(0.275, row_y[[i]], tcga_cox$n[[i]], cex = 0.61)
    segments(map_hr(tcga_cox$CI95_low[[i]]), row_y[[i]],
             map_hr(tcga_cox$CI95_high[[i]]), row_y[[i]], lwd = 0.95)
    points(map_hr(tcga_cox$HR[[i]]), row_y[[i]], pch = 15, cex = 0.62)
    text(0.685, row_y[[i]], sprintf(
      "%.2f (%.2f, %.2f)", tcga_cox$HR[[i]], tcga_cox$CI95_low[[i]], tcga_cox$CI95_high[[i]]
    ), adj = 0, cex = 0.55)
    text(0.965, row_y[[i]], fmt_p(tcga_cox$p_value[[i]]), adj = 1, cex = 0.57)
  }
  axis_y <- 0.20
  segments(0.325, axis_y, 0.635, axis_y, col = line_col, lwd = 0.55)
  for (tick in c(0.5, 1, 2, 5)) {
    x <- map_hr(tick)
    segments(x, axis_y, x, axis_y - 0.06, col = line_col, lwd = 0.55)
    text(x, axis_y - 0.12, tick, cex = 0.50)
  }
  panel_tag("B", line = 0.6)
}

draw_stage <- function() {
  par(mar = c(3.1, 3.6, 2.4, 0.3), family = "sans", xpd = NA)
  score_order <- c("P", "B_NQO1", "I_NQO1")
  score_labels <- c(P = "Pᵀ", B_NQO1 = "Bᵀ", I_NQO1 = "Pᵀ−Bᵀ")
  positions <- matrix(c(1, 2, 4, 5, 7, 8), nrow = 3L, byrow = TRUE,
                      dimnames = list(score_order, c("II", "III")))
  value_min <- min(tcga_stage$value, na.rm = TRUE)
  value_max <- max(tcga_stage$value, na.rm = TRUE)
  value_range <- value_max - value_min
  y_min <- floor(value_min * 2) / 2
  bracket_y <- value_max + max(0.32, value_range * 0.065)
  y_top <- bracket_y + max(0.52, value_range * 0.11)
  plot(NA, xlim = c(0.25, 8.75), ylim = c(y_min, y_top), axes = FALSE, xlab = "", ylab = "")
  y_ticks <- pretty(c(y_min, value_max), n = 5)
  y_ticks <- y_ticks[y_ticks >= y_min & y_ticks <= value_max]
  abline(h = y_ticks, col = "#E7EEF1", lwd = 0.65)
  abline(h = 0, col = "#B8C2C8", lwd = 0.65)
  axis(1, at = rowMeans(positions), labels = score_labels, cex.axis = 0.58,
       lwd = 0, lwd.ticks = 0)
  axis(2, at = y_ticks, las = 1, cex.axis = 0.61)
  mtext("Composite score (gene-wise Z units)", side = 2, line = 2.6, cex = 0.65)

  stage_palette <- c(II = red, III = cyan)
  for (i in seq_along(score_order)) {
    score <- score_order[[i]]
    for (stage_name in c("II", "III")) {
      values <- tcga_stage$value[tcga_stage$score == score & tcga_stage$stage == stage_name]
      position <- positions[score, stage_name]
      boxplot(values, at = position, add = TRUE, axes = FALSE, outline = FALSE,
              width = 0.66, border = stage_palette[[stage_name]],
              col = adjustcolor(stage_palette[[stage_name]], alpha.f = 0.04),
              lwd = 0.85, pars = list(medlwd = 1.25, whisklty = 1, staplewex = 0.45))
      set.seed(420 + i * 10 + match(stage_name, c("II", "III")))
      points(position + runif(length(values), -0.22, 0.22), values,
             pch = 16, cex = 0.21, col = adjustcolor(stage_palette[[stage_name]], alpha.f = 0.70))
    }
    test <- tcga_stage_tests[tcga_stage_tests$score == score, , drop = FALSE]
    x_left <- positions[score, "II"]
    x_right <- positions[score, "III"]
    segments(x_left, bracket_y, x_right, bracket_y, col = line_col, lwd = 0.55)
    segments(c(x_left, x_right), bracket_y, c(x_left, x_right), bracket_y - 0.12,
             col = line_col, lwd = 0.55)
    text(mean(c(x_left, x_right)), bracket_y + 0.13,
         paste0("p = ", fmt_p(test$p_value)), cex = 0.50)
  }
  legend(
    "topright", inset = c(0, -0.03), xpd = NA, horiz = TRUE, bty = "n",
    legend = c(
      sprintf("Stage II (n=%d)", length(unique(tcga_stage$patient_id[tcga_stage$stage == "II"]))),
      sprintf("Stage III (n=%d)", length(unique(tcga_stage$patient_id[tcga_stage$stage == "III"])))
    ), col = stage_palette, pch = 16, pt.cex = 0.60, cex = 0.51
  )
  panel_tag("C", line = 0.8)
}

gene_average_contrasts <- c("OxR_average_3_cell_lines", "SN38R_average_3_cell_lines")
gene_average <- gse_genes[gse_genes$contrast %in% gene_average_contrasts, , drop = FALSE]
gene_limit <- max(abs(c(gene_average$ci_lower_log2, gene_average$ci_upper_log2)), na.rm = TRUE)
gene_limit <- ceiling(gene_limit * 5) / 5

score_order <- c("P", "B_NQO1", "I_NQO1")
forest_contrasts <- c(
  "HCT116_OxR_vs_Parental", "HT29_OxR_vs_Parental", "LoVo_OxR_vs_Parental",
  "OxR_average_3_cell_lines", "HCT116_SN38R_vs_Parental",
  "HT29_SN38R_vs_Parental", "LoVo_SN38R_vs_Parental", "SN38R_average_3_cell_lines"
)
forest_labels <- c(
  "HCT116 OxR", "HT29 OxR", "LoVo OxR", "OxR average",
  "HCT116 SN38R", "HT29 SN38R", "LoVo SN38R", "SN38R average"
)
score_primary <- gse_scores[
  gse_scores$analysis_role == "Primary framework" &
    gse_scores$score %in% score_order & gse_scores$contrast %in% forest_contrasts,
  , drop = FALSE
]
score_limit <- max(abs(c(score_primary$ci_lower, score_primary$ci_upper)), na.rm = TRUE)
score_limit <- ceiling(score_limit * 2) / 2

oxr_col <- "#C86B3C"
sn38r_col <- "#3C91B7"
neutral_band_1 <- "#F7F7F7"
neutral_band_2 <- "#F1F3F4"
sensitivity_band <- "#F2F3F3"

draw_gene_forest_final <- function() {
  par(mar = c(3.25, 4.95, 1.95, 0.35), family = "sans", xpd = NA)
  gene_y <- rev(seq_along(target_genes))
  gene_axis_labels <- c(target_genes[-length(target_genes)], "NFE2L2†")
  plot(
    NA, xlim = c(-gene_limit, gene_limit), ylim = c(0.40, 12.65), axes = FALSE,
    xlab = "", ylab = ""
  )
  rect(-gene_limit, 7.50, gene_limit, 11.55, col = neutral_band_1, border = NA)
  rect(-gene_limit, 1.50, gene_limit, 7.50, col = neutral_band_2, border = NA)
  rect(-gene_limit, 0.45, gene_limit, 1.50, col = sensitivity_band, border = NA)
  segments(0, 0.45, 0, 11.55, lty = 2, lwd = 0.75, col = "#777777")
  text(-gene_limit + 0.04, 11.72, "Pressure-domain genes", adj = 0,
       col = line_col, font = 2, cex = 0.61)
  text(-gene_limit + 0.04, 7.22, "Buffering-domain genes", adj = 0,
       col = line_col, font = 2, cex = 0.61)
  text(-gene_limit + 0.04, 1.39, "Sensitivity definition only", adj = 0,
       col = "#666666", font = 2, cex = 0.55)

  for (i in seq_along(target_genes)) {
    d <- gene_average[gene_average$target_gene == target_genes[[i]], , drop = FALSE]
    d <- d[match(gene_average_contrasts, d$contrast), , drop = FALSE]
    segments(d$ci_lower_log2[[1L]], gene_y[[i]] + 0.15,
             d$ci_upper_log2[[1L]], gene_y[[i]] + 0.15, col = oxr_col, lwd = 1.20)
    points(d$logFC[[1L]], gene_y[[i]] + 0.15, pch = 15, cex = 0.70, col = oxr_col)
    segments(d$ci_lower_log2[[2L]], gene_y[[i]] - 0.15,
             d$ci_upper_log2[[2L]], gene_y[[i]] - 0.15, col = sn38r_col, lwd = 0.95)
    points(d$logFC[[2L]], gene_y[[i]] - 0.15, pch = 16, cex = 0.53, col = sn38r_col)
  }
  axis(1, at = seq(-1, 1, by = 0.5), cex.axis = 0.58)
  axis(2, at = gene_y, labels = gene_axis_labels, las = 1, tick = FALSE, cex.axis = 0.60)
  mtext("Difference in log2 expression (resistant - parental)", side = 1, line = 2.10, cex = 0.61)
  box(bty = "l", lwd = 0.55)
  legend("topright", inset = c(0, 0.005), bty = "n", horiz = TRUE,
         legend = c("OxR average", "SN38R average"), pch = c(15, 16),
         pt.cex = c(0.84, 0.64), col = c(oxr_col, sn38r_col), cex = 0.58, x.intersp = 0.65)
  panel_tag("D", line = 0.55, adj = -0.12)
}

draw_score_forest_final <- function() {
  par(mar = c(3.25, 0.65, 1.95, 0.35), family = "sans", xpd = NA)
  score_y <- c(8.05, 6.92, 5.79, 4.66, 3.50, 2.37, 1.24, 0.36)
  score_offsets <- c(-0.19, 0, 0.19)
  score_pch <- c(P = 16, B_NQO1 = 17, I_NQO1 = 18)
  plot(NA, xlim = c(-score_limit, score_limit), ylim = c(-0.25, 9.75), axes = FALSE,
       xlab = "", ylab = "")
  rect(-score_limit, 4.12, score_limit, 8.62, col = neutral_band_1, border = NA)
  rect(-score_limit, 0.02, score_limit, 3.90, col = neutral_band_2, border = NA)
  segments(0, 4.12, 0, 8.62, lty = 2, lwd = 0.75, col = "#777777")
  segments(0, 0.02, 0, 3.90, lty = 2, lwd = 0.75, col = "#777777")
  text(-score_limit + 0.05, 8.80, "OxR", adj = 0, col = oxr_col, font = 2, cex = 0.66)
  text(-score_limit + 0.05, 4.00, "SN38R", adj = 0, col = sn38r_col, font = 2, cex = 0.66)

  for (row_i in seq_along(forest_contrasts)) {
    treatment_col <- if (row_i <= 4L) oxr_col else sn38r_col
    is_average <- row_i %in% c(4L, 8L)
    text(-score_limit + 0.06, score_y[[row_i]], forest_labels[[row_i]], adj = 0,
         cex = 0.57, font = if (is_average) 2 else 1, col = line_col)
    for (score_i in seq_along(score_order)) {
      d <- score_primary[
        score_primary$score == score_order[[score_i]] &
          score_primary$contrast == forest_contrasts[[row_i]], , drop = FALSE
      ]
      segments(d$ci_lower, score_y[[row_i]] + score_offsets[[score_i]],
               d$ci_upper, score_y[[row_i]] + score_offsets[[score_i]],
               col = treatment_col, lwd = if (is_average) 1.18 else 0.84)
      points(d$estimate, score_y[[row_i]] + score_offsets[[score_i]],
             pch = score_pch[[score_order[[score_i]]]],
             cex = if (is_average) 0.76 else 0.50, col = treatment_col)
    }
  }
  axis(1, at = seq(-2, 2, by = 1), cex.axis = 0.58)
  mtext("Difference in mean gene-wise Z score (resistant - parental)", side = 1, line = 2.10, cex = 0.61)
  box(bty = "l", lwd = 0.55)
  legend("top", inset = c(0, 0.005), bty = "n", horiz = TRUE,
         legend = c(expression(P^T), expression(B^T), expression(P^T - B^T)),
         pch = c(16, 17, 18), pt.cex = c(0.90, 0.98, 0.98),
         col = line_col, cex = 0.68, x.intersp = 0.75)
  panel_tag("E", line = 0.55, adj = -0.08)
}

draw_combined_figure <- function() {
  layout_matrix <- rbind(
    c(rep(1, 9), rep(4, 6)),
    c(rep(1, 9), rep(4, 6)),
    c(rep(1, 9), rep(4, 6)),
    c(rep(1, 9), rep(4, 6)),
    c(rep(2, 9), rep(5, 6)),
    c(rep(2, 9), rep(5, 6)),
    c(rep(3, 9), rep(5, 6)),
    c(rep(3, 9), rep(5, 6)),
    c(rep(3, 9), rep(5, 6)),
    c(rep(6, 7), rep(7, 8)),
    c(rep(6, 7), rep(7, 8)),
    c(rep(6, 7), rep(7, 8)),
    c(rep(6, 7), rep(7, 8)),
    c(rep(6, 7), rep(7, 8)),
    c(rep(6, 7), rep(7, 8)),
    c(rep(6, 7), rep(7, 8))
  )
  layout(
    layout_matrix,
    widths = rep(1, 15),
    heights = c(rep(0.85, 4), rep(0.55, 2), rep(0.62, 3), rep(0.70, 7))
  )
  par(oma = c(0.15, 0.15, 0.15, 0.15))
  draw_km_curve()
  draw_km_risk()
  draw_km_censor()
  draw_cox()
  draw_stage()
  draw_gene_forest_final()
  draw_score_forest_final()
}

save_figure_formats <- function(stem, width_mm = 183, height_mm = 220, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  grDevices::cairo_pdf(
    file.path(out_dir, paste0(stem, ".pdf")),
    width = width_in, height = height_in, family = "Arial", pointsize = 8
  )
  draw_combined_figure()
  dev.off()

  png(
    file.path(out_dir, paste0(stem, ".png")),
    width = width_in, height = height_in, units = "in", res = 300,
    pointsize = 8, type = if (capabilities("cairo")) "cairo" else getOption("bitmapType")
  )
  draw_combined_figure()
  dev.off()

  tiff(
    file.path(out_dir, paste0(stem, ".tiff")),
    width = width_in, height = height_in, units = "in", res = dpi,
    pointsize = 8, compression = "lzw"
  )
  draw_combined_figure()
  dev.off()
}

stem <- "Figure10_TCGA_COAD_GSE42387_core_results_combined"
save_figure_formats(stem)

caption <- paste(
  "Figure 10. TCGA-COAD clinical associations and GSE42387 resistance-related transcriptional analyses.",
  "(A) Descriptive overall-survival curves for 289 stage II-III TCGA-COAD tumors stratified at the median exploratory pressure-to-buffering transcriptional index (Pᵀ−Bᵀ); the median split was used for visualization only.",
  "(B) Primary multivariable Cox model including continuous ferroptotic pressure-related transcriptional (Pᵀ) and antioxidant-buffering-related transcriptional (Bᵀ) scores, age per 10 years, and pathological stage. Hazard ratios and 95% confidence intervals are reported for 55 deaths.",
  "(C) Distributions of Pᵀ, Bᵀ, and Pᵀ−Bᵀ by pathological stage; p values are from two-sided Wilcoxon rank-sum tests.",
  "(D) Limma empirical-Bayes estimates of equally weighted average resistant-minus-parental effects across HCT116, HT29, and LoVo cells for 11 prespecified genes. Squares and circles denote OxR and SN38R averages, respectively; orange and blue indicate the pressure and buffering domains, respectively; horizontal lines indicate 95% confidence intervals. Pressure- and buffering-domain genes are separated by labeled bands. NFE2L2† was used only in the prespecified sensitivity definition and was not a constituent of the primary Bᵀ score.",
  "(E) Ordinary fixed-effects estimates for Pᵀ, Bᵀ, and Pᵀ−Bᵀ in the same nine-group GSE42387 design. Orange and blue denote oxaliplatin-resistant (OxR) and SN-38-resistant (SN38R) contrasts, respectively; circle, triangle, and diamond symbols denote Pᵀ, Bᵀ, and Pᵀ−Bᵀ. Composite scores are expressed in mean gene-wise Z-score units.",
  "GSE42387 included 27 arrays, comprising three independent passage replicates for each cell-line-by-status group.",
  sep = " "
)
writeLines(caption, file.path(out_dir, paste0(stem, "_caption.txt")), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(out_dir, paste0(stem, "_R_sessionInfo.txt")))
