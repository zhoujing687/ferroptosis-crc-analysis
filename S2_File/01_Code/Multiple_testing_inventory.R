###############################################################################
# Inventory of multiplicity families as implemented in archived analysis code.
# Added 2026-11-21. This audit does not establish prospective preregistration.
# Usage: Rscript Multiple_testing_inventory.R <TCGA_03_Source_Data> <GSE_03_Source_Data> <outdir>
###############################################################################
set.seed(42)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Supply TCGA/GSE source-data directories and output directory.")
tcga <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
gse <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
out <- args[[3L]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
read <- function(dir, name) read.csv(file.path(dir, name), check.names = FALSE)
rows <- list()
add <- function(analysis, family, hypotheses, families, method, p_column, scope, source) {
  rows[[length(rows) + 1L]] <<- data.frame(Analysis = analysis, Family = family,
    Hypotheses_per_family = hypotheses, Number_of_families = families, Adjustment = method,
    Adjusted_p_column = p_column, Scope = scope, Evidence = source)
}
de <- read(gse, "GSE42387_full_matrix_limma_summary.csv")
stopifnot(nrow(de) == 8L, length(unique(de$features_tested)) == 1L)
add("GSE42387 full-feature limma", "Each main contrast", de$features_tested[[1L]], 8L, "BH",
    "adj.P.Val", "Eight separate genome-wide contrast families", "GSE analysis: topTable(adjust.method = BH)")
add("GSE42387 full-feature limma", "Each omnibus test", de$features_tested[[1L]], 4L, "BH",
    "adj.P.Val", "Four separate genome-wide omnibus families", "GSE analysis: topTable(adjust.method = BH)")
target <- read(gse, "GSE42387_11_genes_limma_main_contrasts.csv")
target_omni <- read(gse, "GSE42387_11_genes_limma_omnibus_tests.csv")
stopifnot(all(table(target$contrast) == 11L), all(table(target_omni$test) == 11L))
add("GSE42387 target-gene limma", "Each main contrast", 11L, length(unique(target$contrast)), "BH",
    "p_adjusted_BH_11_genes_within_contrast", "Eleven targets, including the NFE2L2 sensitivity gene", "GSE analysis: ave(P.Value, contrast, p.adjust)")
add("GSE42387 target-gene limma", "Each omnibus test", 11L, length(unique(target_omni$test)), "BH",
    "p_adjusted_BH_11_genes_within_test", "Eleven targets within each omnibus test", "GSE analysis: ave(P.Value, test, p.adjust)")
for (hc3 in c(FALSE, TRUE)) {
  name <- if (hc3) "GSE42387_PBI_HC3_robust_sensitivity_contrasts.csv" else "GSE42387_PBI_fixed_effects_contrasts.csv"
  d <- read(gse, name)
  for (role in c("Primary framework", "Sensitivity")) {
    for (average in c(TRUE, FALSE)) {
      idx <- d$analysis_role == role & grepl("average", d$contrast) == average
      p_col <- if (hc3) "p_value_HC3" else "p_value"
      q_col <- if (hc3) "p_adjusted_BH_HC3_prespecified_family" else "p_adjusted_BH_prespecified_family"
      stopifnot(isTRUE(all.equal(p.adjust(d[[p_col]][idx], "BH"), d[[q_col]][idx], tolerance = 1e-12)))
      add(if (hc3) "GSE42387 composite HC3 sensitivity" else "GSE42387 ordinary composite contrasts",
          paste(role, if (average) "average effects" else "within-cell-line effects"), sum(idx), 1L,
          "BH", q_col, if (average) "Two drugs jointly across included scores" else "Six drug-by-cell-line contrasts jointly across included scores",
          "GSE analysis: four explicit index masks and p.adjust")
    }
  }
}
d <- read(gse, "GSE42387_PBI_fixed_effects_omnibus_tests.csv")
for (role in c("Primary framework", "Sensitivity")) {
  idx <- d$analysis_role == role
  stopifnot(isTRUE(all.equal(p.adjust(d$p_value[idx], "BH"), d$p_adjusted_BH_prespecified_family[idx], tolerance = 1e-12)))
  add("GSE42387 ordinary composite omnibus", role, sum(idx), 1L, "BH",
      "p_adjusted_BH_prespecified_family", "Four omnibus tests jointly across included scores", "GSE analysis: p.adjust by analysis_role")
}
stage <- read(tcga, "TCGA_Figure10_panelC_stage_comparisons.csv")
stopifnot(nrow(stage) == 3L, isTRUE(all.equal(p.adjust(stage$p_value, "BH"), stage$p_BH, tolerance = 1e-12)))
add("TCGA-COAD stage associations", "Primary P, B, and P-B stage comparisons", 3L, 1L, "BH", "p_BH",
    "Three two-sided Wilcoxon tests jointly", "TCGA analysis: stage_score_tests$p_BH <- p.adjust")
cox <- read(tcga, "TCGA_PBI_Cox_models.csv")
for (model in unique(cox$model)) {
  count <- sum(cox$model == model)
  add("TCGA-COAD Cox model", model, count, 1L, "None", "Not applicable; raw p_value",
      "Model coefficients include exposure score(s), age, and stage; exploratory, not multiplicity-adjusted", "Archived TCGA_PBI_Cox_models.csv and coxph model definitions")
}
add("TCGA-COAD descriptive survival curves", "Median-split P-B log-rank", 1L, 1L, "None",
    "Not applicable", "Descriptive visualization; not a validated cutpoint or a confirmatory biomarker test", "TCGA analysis: survdiff and unadjusted logrank_p")
inventory <- do.call(rbind, rows)
write.csv(inventory, file.path(out, "S7_Table.csv"), row.names = FALSE, na = "")
writeLines(c("BH means Benjamini-Hochberg. Families are documented as implemented in the archived code.",
  "This inventory does not establish that analyses were preregistered before outcomes were inspected.",
  "Diagnostic Shapiro-Wilk, Brown-Forsythe, and proportional-hazards tests are diagnostic screens; they were not included in efficacy families.",
  "No new p-values or correction scheme have been substituted in the archived analyses."),
  file.path(out, "S7_Table_notes.txt"))
capture.output(sessionInfo(), file = file.path(out, "S7_Table_R_sessionInfo.txt"))
print(inventory[, c("Analysis", "Family", "Hypotheses_per_family", "Number_of_families", "Adjustment")], row.names = FALSE)
