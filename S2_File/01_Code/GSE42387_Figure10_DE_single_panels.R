###############################################################################
# Figure 10D and Figure 10E: standalone forest plots
#
# This script reads the locked GSE42387 result tables and only redraws the
# figures. It does not rerun limma, refit the models, or change any estimates.
#
# Usage:
#   Rscript GSE42387_Figure10_DE_single_panels.R <analysis_dir> <output_dir>
#
# Outputs:
#   Figure10D_GSE42387_11_gene_effects.png/.pdf/.tiff
#   Figure10E_GSE42387_PBI_fixed_effects.png/.pdf/.tiff
###############################################################################

options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "Usage: Rscript GSE42387_Figure10_DE_single_panels.R ",
    "<analysis_dir> <output_dir>"
  )
}

analysis_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

genes_file <- file.path(analysis_dir, "GSE42387_11_genes_limma_main_contrasts.csv")
scores_file <- file.path(analysis_dir, "GSE42387_PBI_fixed_effects_contrasts.csv")
if (!file.exists(genes_file) || !file.exists(scores_file)) {
  stop("Required GSE42387 result CSV files were not found in analysis_dir.")
}

genes <- read.csv(genes_file, check.names = FALSE)
scores <- read.csv(scores_file, check.names = FALSE)

required_gene_cols <- c(
  "target_gene", "contrast", "logFC", "ci_lower_log2", "ci_upper_log2"
)
required_score_cols <- c(
  "score", "contrast", "estimate", "ci_lower", "ci_upper", "analysis_role"
)
if (!all(required_gene_cols %in% names(genes))) {
  stop("The gene result table is missing one or more required columns.")
}
if (!all(required_score_cols %in% names(scores))) {
  stop("The score result table is missing one or more required columns.")
}

# Figure 10D uses domain-specific colors and lighter secondary treatment marks;
# Figure 10E uses the darker treatment colors for all score estimates.
pressure_gene_col <- "#D4743C"
pressure_gene_light <- "#E2A47F"
buffering_gene_col <- "#1874A6"
buffering_gene_light <- "#69A6C6"
oxr_col <- "#C86B3C"
sn38r_col <- "#3C91B7"
ink <- "#252525"
reference_col <- "#737373"
pressure_band <- "#FBF1EC"
buffering_band <- "#EEF6FA"

pressure_genes <- c("ACSL4", "LPCAT3", "TFRC", "NCOA4")
target_genes <- c(
  pressure_genes, "SLC7A11", "GPX4", "AIFM2", "GCH1", "DHODH", "NQO1", "NFE2L2"
)
gene_contrasts <- c("OxR_average_3_cell_lines", "SN38R_average_3_cell_lines")

score_order <- c("P", "B_NQO1", "I_NQO1")
score_contrasts <- c(
  "HCT116_OxR_vs_Parental", "HT29_OxR_vs_Parental", "LoVo_OxR_vs_Parental",
  "OxR_average_3_cell_lines", "HCT116_SN38R_vs_Parental",
  "HT29_SN38R_vs_Parental", "LoVo_SN38R_vs_Parental", "SN38R_average_3_cell_lines"
)
score_rows <- c(
  "HCT116 OxR", "HT29 OxR", "LoVo OxR", "OxR average",
  "HCT116 SN38R", "HT29 SN38R", "LoVo SN38R", "SN38R average"
)

get_one <- function(data, condition, label) {
  out <- data[condition, , drop = FALSE]
  if (nrow(out) != 1L) stop("Expected one row for ", label, "; found ", nrow(out), ".")
  out
}

gene_data <- do.call(rbind, lapply(target_genes, function(gene) {
  d <- genes[genes$target_gene == gene & genes$contrast %in% gene_contrasts, , drop = FALSE]
  d <- d[match(gene_contrasts, d$contrast), , drop = FALSE]
  if (nrow(d) != 2L || anyNA(d$contrast)) stop("Missing gene contrast for ", gene, ".")
  d
}))
gene_data$domain <- ifelse(gene_data$target_gene %in% pressure_genes, "Pressure", "Buffering")

score_data <- scores[
  scores$analysis_role == "Primary framework" &
    scores$score %in% score_order & scores$contrast %in% score_contrasts,
  , drop = FALSE
]
score_data <- do.call(rbind, lapply(score_order, function(score) {
  d <- score_data[score_data$score == score, , drop = FALSE]
  d <- d[match(score_contrasts, d$contrast), , drop = FALSE]
  if (nrow(d) != length(score_contrasts) || anyNA(d$contrast)) {
    stop("Missing primary score contrast for ", score, ".")
  }
  d
}))

save_figure_formats <- function(stem, width, height, draw) {
  pdf(file.path(output_dir, paste0(stem, ".pdf")), width = width, height = height,
      family = "sans", useDingbats = FALSE)
  draw()
  dev.off()

  png(file.path(output_dir, paste0(stem, ".png")), width = width, height = height,
      units = "in", res = 300,
      type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
  draw()
  dev.off()

  tiff(file.path(output_dir, paste0(stem, ".tiff")), width = width, height = height,
       units = "in", res = 600, compression = "lzw")
  draw()
  dev.off()
}

draw_figure10d <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old))
  par(mar = c(4.75, 5.25, 1.40, 0.35), family = "sans", xpd = NA)

  gene_y <- rev(seq_along(target_genes))
  gene_axis_labels <- c(target_genes[-length(target_genes)], "NFE2L2\u2020")
  gene_limit <- 1.25
  plot(NA, xlim = c(-gene_limit, gene_limit), ylim = c(0.40, 12.70),
       axes = FALSE, xlab = "", ylab = "")
  rect(-gene_limit, 7.50, gene_limit, 11.55, col = pressure_band, border = NA)
  rect(-gene_limit, 1.50, gene_limit, 7.50, col = buffering_band, border = NA)
  rect(-gene_limit, 0.45, gene_limit, 1.50, col = buffering_band, border = NA)
  segments(0, 0.45, 0, 11.55, lty = 2, lwd = 0.85, col = reference_col)
  text(-gene_limit + 0.04, 11.72, "Pressure domain", adj = 0,
       col = pressure_gene_col, font = 2, cex = 0.64)
  text(-gene_limit + 0.04, 7.22, "Buffering domain", adj = 0,
       col = buffering_gene_col, font = 2, cex = 0.64)
  text(-gene_limit + 0.04, 1.39, "Sensitivity definition only", adj = 0,
       col = "#666666", font = 2, cex = 0.58)

  for (i in seq_along(target_genes)) {
    d <- gene_data[gene_data$target_gene == target_genes[[i]], , drop = FALSE]
    d <- d[match(gene_contrasts, d$contrast), , drop = FALSE]
    in_pressure_domain <- target_genes[[i]] %in% pressure_genes
    oxr_gene_col <- if (in_pressure_domain) pressure_gene_col else buffering_gene_col
    sn38r_gene_col <- if (in_pressure_domain) pressure_gene_light else buffering_gene_light
    segments(d$ci_lower_log2[[1]], gene_y[[i]] + 0.15,
             d$ci_upper_log2[[1]], gene_y[[i]] + 0.15,
             col = oxr_gene_col, lwd = 1.35)
    points(d$logFC[[1]], gene_y[[i]] + 0.15,
           pch = 15, cex = 0.78, col = oxr_gene_col)
    segments(d$ci_lower_log2[[2]], gene_y[[i]] - 0.15,
             d$ci_upper_log2[[2]], gene_y[[i]] - 0.15,
             col = sn38r_gene_col, lwd = 1.05)
    points(d$logFC[[2]], gene_y[[i]] - 0.15,
           pch = 16, cex = 0.58, col = sn38r_gene_col)
  }
  axis(1, at = seq(-1, 1, by = 0.5), cex.axis = 0.65)
  axis(2, at = gene_y, labels = gene_axis_labels, las = 1,
       tick = FALSE, cex.axis = 0.68)
  mtext("Difference in log2 expression (resistant - parental)",
        side = 1, line = 2.65, cex = 0.70)
  box(bty = "l", lwd = 0.60)
  legend("topright", inset = c(0, 0.005), bty = "n", horiz = TRUE,
         legend = c("OxR average", "SN38R average"), pch = c(15, 16),
         pt.cex = c(1.02, 0.78), col = ink,
         cex = 0.70, x.intersp = 0.75)
  mtext("D", side = 3, line = 0.30, adj = -0.08, font = 2, cex = 1.10)
}

draw_figure10e <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old))
  par(mar = c(4.75, 4.80, 1.40, 0.35), family = "sans", xpd = NA)

  score_y <- c(8.05, 6.92, 5.79, 4.66, 3.50, 2.37, 1.24, 0.36)
  score_limit <- 2.15
  score_offsets <- c(-0.19, 0, 0.19)
  score_pch <- c(P = 16, B_NQO1 = 17, I_NQO1 = 18)
  plot(NA, xlim = c(-score_limit, score_limit), ylim = c(-0.25, 9.75),
       axes = FALSE, xlab = "", ylab = "")
  rect(-score_limit, 4.12, score_limit, 8.62,
       col = pressure_band, border = NA)
  rect(-score_limit, 0.02, score_limit, 3.90,
       col = buffering_band, border = NA)
  segments(0, 4.12, 0, 8.62, lty = 2, lwd = 0.85, col = reference_col)
  segments(0, 0.02, 0, 3.90, lty = 2, lwd = 0.85, col = reference_col)
  text(-score_limit + 0.05, 8.80, "OxR", adj = 0,
       col = oxr_col, font = 2, cex = 0.72)
  text(-score_limit + 0.05, 4.00, "SN38R", adj = 0,
       col = sn38r_col, font = 2, cex = 0.72)

  for (row_i in seq_along(score_contrasts)) {
    treatment_col <- if (row_i <= 4) oxr_col else sn38r_col
    is_average <- row_i %in% c(4, 8)
    text(-score_limit + 0.06, score_y[[row_i]], score_rows[[row_i]], adj = 0,
         cex = 0.61, font = if (is_average) 2 else 1, col = ink)
    for (score_i in seq_along(score_order)) {
      d <- get_one(
        score_data,
        score_data$score == score_order[[score_i]] &
          score_data$contrast == score_contrasts[[row_i]],
        paste(score_order[[score_i]], score_contrasts[[row_i]])
      )
      line_width <- if (is_average) 1.30 else 0.92
      point_size <- if (is_average) 0.85 else 0.54
      segments(d$ci_lower, score_y[[row_i]] + score_offsets[[score_i]],
               d$ci_upper, score_y[[row_i]] + score_offsets[[score_i]],
               col = treatment_col, lwd = line_width)
      points(d$estimate, score_y[[row_i]] + score_offsets[[score_i]],
             pch = score_pch[[score_order[[score_i]]]], cex = point_size,
             col = treatment_col)
    }
  }
  axis(1, at = seq(-2, 2, by = 1), cex.axis = 0.65)
  axis(2, at = numeric(), labels = FALSE, tick = FALSE)
  mtext("Difference in mean gene-wise Z score (resistant - parental)",
        side = 1, line = 2.65, cex = 0.70)
  box(bty = "l", lwd = 0.60)
  legend("top", inset = c(0, 0.005), bty = "n", horiz = TRUE,
         legend = c(expression(P^T), expression(B^T), expression(P^T - B^T)),
         pch = c(16, 17, 18), pt.cex = c(1.08, 1.16, 1.16),
         col = ink, cex = 0.82, x.intersp = 0.90)
  mtext("E", side = 3, line = 0.30, adj = -0.08, font = 2, cex = 1.10)
}

save_figure_formats("Figure10D_GSE42387_11_gene_effects", 4.80, 4.45, draw_figure10d)
save_figure_formats("Figure10E_GSE42387_PBI_fixed_effects", 4.80, 4.45, draw_figure10e)

writeLines(
  c(
    "Figure 10D and Figure 10E were redrawn from the existing locked GSE42387 result tables.",
    "No statistical model was refit and no estimates were changed.",
    "Figure 10D: squares and circles denote OxR and SN38R averages, respectively; orange and blue indicate the pressure and buffering domains, respectively.",
    "Figure 10E: circle = P^T, triangle = B^T, diamond = P^T-B^T; color identifies OxR or SN38R.",
    "The dagger after NFE2L2 identifies the sensitivity-definition-only gene."
  ),
    file.path(output_dir, "Figure10_DE_single_panels_README.txt"),
  useBytes = TRUE
)

message("Finished. Standalone Figure 10D and Figure 10E files are in: ", output_dir)
