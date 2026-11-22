###############################################################################
# Analysis: GSE42387 pressure-buffering analysis in 27 arrays
# Date: 2026-08-12
# Random seed: 42
# R: >= 4.5.0
# Packages: limma, statmod
#
# Single-gene/full-feature inference: nine-group limma empirical Bayes.
# Composite P/B/I inference: ordinary fixed-effects linear models using the
# same nine-group design and prespecified contrasts (no empirical Bayes).
###############################################################################

set.seed(42)
options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: Rscript GSE42387_PBI_fixed_effects_analysis.R <series_matrix.gz> <GPL16297.txt> <outdir>")
}
series_file <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
platform_file <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
out_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)

if (!requireNamespace("limma", quietly = TRUE)) stop("Bioconductor package limma is required.")
if (!requireNamespace("statmod", quietly = TRUE)) stop("CRAN package statmod is required.")

pressure_genes <- c("ACSL4", "LPCAT3", "TFRC", "NCOA4")
buffer_genes_nqo1 <- c("SLC7A11", "GPX4", "AIFM2", "GCH1", "DHODH", "NQO1")
buffer_genes_nfe2l2 <- c("SLC7A11", "GPX4", "AIFM2", "GCH1", "DHODH", "NFE2L2")
target_genes <- unique(c(pressure_genes, buffer_genes_nqo1, buffer_genes_nfe2l2))

write_csv <- function(x, name, row.names = FALSE) {
  write.csv(x, file.path(out_dir, name), row.names = row.names, na = "", fileEncoding = "UTF-8")
}

write_csv_gz <- function(x, name) {
  con <- gzfile(file.path(out_dir, name), open = "wt")
  on.exit(close(con), add = TRUE)
  write.csv(x, con, row.names = FALSE, na = "")
}

z_score <- function(x) {
  s <- sd(x)
  if (!is.finite(s) || s <= 0) stop("Zero variance prevents gene-wise Z standardization.")
  as.numeric((x - mean(x)) / s)
}

save_figure_formats <- function(stem, width, height, draw) {
  pdf(file.path(out_dir, paste0(stem, ".pdf")), width = width, height = height,
      family = "sans", useDingbats = FALSE)
  draw(); dev.off()
  png(file.path(out_dir, paste0(stem, ".png")), width = width, height = height,
      units = "in", res = 300, type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
  draw(); dev.off()
  tiff(file.path(out_dir, paste0(stem, ".tiff")), width = width, height = height,
       units = "in", res = 600, compression = "lzw")
  draw(); dev.off()
}

extract_metadata_line <- function(prefix, lines) {
  hit <- grep(paste0("^", prefix, "\t"), lines, value = TRUE)
  if (length(hit) != 1L) stop("Expected one metadata line for ", prefix)
  gsub(intToUtf8(34), "", strsplit(hit, "\t")[[1L]][-1L], fixed = TRUE)
}

# Read the processed series matrix and all 27 arrays.
geo_lines <- readLines(gzfile(series_file), warn = FALSE)
sample_titles <- extract_metadata_line("!Sample_title", geo_lines)
sample_ids <- extract_metadata_line("!Sample_geo_accession", geo_lines)
table_start <- grep("^!series_matrix_table_begin", geo_lines) + 1L
table_end <- grep("^!series_matrix_table_end", geo_lines) - 1L
if (length(table_start) != 1L || length(table_end) != 1L) stop("Series table boundaries not found.")
matrix_table <- read.delim(
  text = paste(geo_lines[table_start:table_end], collapse = "\n"), header = TRUE,
  sep = "\t", quote = "\"", check.names = FALSE, comment.char = ""
)
rm(geo_lines); invisible(gc())
if (!"ID_REF" %in% names(matrix_table) || anyDuplicated(matrix_table$ID_REF)) {
  stop("Invalid or duplicated ID_REF values.")
}
row.names(matrix_table) <- matrix_table$ID_REF
expression_all <- as.matrix(matrix_table[, setdiff(names(matrix_table), "ID_REF"), drop = FALSE])
storage.mode(expression_all) <- "double"
rm(matrix_table); invisible(gc())

metadata <- data.frame(
  sample = sample_ids,
  title = sample_titles,
  cell_line = sub("[-_].*", "", sample_titles),
  status = ifelse(grepl("-OxPt", sample_titles, fixed = TRUE), "OxR",
                  ifelse(grepl("-SN38", sample_titles, fixed = TRUE), "SN38R", "Parental")),
  stringsAsFactors = FALSE
)
metadata$cell_line <- factor(metadata$cell_line, levels = c("HCT116", "HT29", "LoVo"))
metadata$status <- factor(metadata$status, levels = c("Parental", "OxR", "SN38R"))
metadata$group <- factor(paste(metadata$cell_line, metadata$status, sep = "_"), levels = c(
  "HCT116_Parental", "HCT116_OxR", "HCT116_SN38R",
  "HT29_Parental", "HT29_OxR", "HT29_SN38R",
  "LoVo_Parental", "LoVo_OxR", "LoVo_SN38R"
))
cell_counts <- table(metadata$cell_line, metadata$status)
if (nrow(metadata) != 27L || anyNA(metadata$group) || any(cell_counts != 3L)) {
  stop("Expected 3 cell lines x 3 statuses x 3 independent passage replicates.")
}
if (!setequal(metadata$sample, colnames(expression_all))) stop("Sample identifiers do not match.")
expression_selected <- expression_all[, metadata$sample, drop = FALSE]
if (anyNA(expression_selected) || any(!is.finite(expression_selected))) stop("Expression matrix has invalid values.")

# Read and audit GPL16297 before any score calculation.
platform_lines <- readLines(platform_file, warn = FALSE)
platform_start <- grep("^!platform_table_begin", platform_lines) + 1L
platform_end <- grep("^!platform_table_end", platform_lines) - 1L
if (length(platform_start) != 1L || length(platform_end) != 1L) stop("Platform table boundaries not found.")
platform <- read.delim(
  text = paste(platform_lines[platform_start:platform_end], collapse = "\n"),
  header = TRUE, sep = "\t", quote = "", fill = TRUE, check.names = FALSE, comment.char = ""
)
rm(platform_lines); invisible(gc())
required_annotation <- c("ID", "GB_ACC", "Entrez Gene ID", "Gene Symbol", "Gene Name")
if (!all(required_annotation %in% names(platform))) stop("GPL16297 annotation columns are missing.")
if (anyDuplicated(platform$ID)) stop("GPL16297 contains duplicated feature IDs.")

annotation_index <- match(row.names(expression_selected), platform$ID)
feature_annotation <- data.frame(
  feature_id = row.names(expression_selected),
  platform_id_matched = !is.na(annotation_index),
  accession = platform$GB_ACC[annotation_index],
  entrez_gene_id = platform[["Entrez Gene ID"]][annotation_index],
  gene_symbol_raw = platform[["Gene Symbol"]][annotation_index],
  gene_name = platform[["Gene Name"]][annotation_index],
  stringsAsFactors = FALSE
)
feature_annotation$gene_symbol_clean <- trimws(feature_annotation$gene_symbol_raw)
feature_annotation$ambiguous_multi_gene <- grepl("///|//|;|,|\\|", feature_annotation$gene_symbol_clean)

target_probe_rows <- do.call(rbind, lapply(target_genes, function(gene) {
  rows <- feature_annotation[
    !is.na(feature_annotation$gene_symbol_clean) &
      feature_annotation$gene_symbol_clean == gene &
      !feature_annotation$ambiguous_multi_gene &
      feature_annotation$platform_id_matched,
    , drop = FALSE
  ]
  if (nrow(rows) == 0L) {
    return(data.frame(feature_id = NA_character_, platform_id_matched = FALSE,
                      accession = NA_character_, entrez_gene_id = NA_character_,
                      gene_symbol_raw = NA_character_, gene_name = NA_character_,
                      gene_symbol_clean = gene, ambiguous_multi_gene = NA,
                      target_gene = gene, stringsAsFactors = FALSE))
  }
  rows$target_gene <- gene
  rows
}))
row.names(target_probe_rows) <- NULL
probe_coverage <- do.call(rbind, lapply(target_genes, function(gene) {
  rows <- target_probe_rows[target_probe_rows$target_gene == gene & !is.na(target_probe_rows$feature_id), ]
  data.frame(
    gene = gene, n_valid_unique_probes = nrow(rows),
    feature_ids = paste(rows$feature_id, collapse = ";"),
    aggregation_rule = if (nrow(rows) > 1L) "arithmetic mean before gene-wise Z standardization" else "single valid probe",
    status = if (nrow(rows) >= 1L) "PASS" else "FAIL",
    stringsAsFactors = FALSE
  )
}))
write_csv(target_probe_rows, "GSE42387_target_probe_mapping_audit.csv")
write_csv(probe_coverage, "GSE42387_target_probe_coverage_summary.csv")
if (any(probe_coverage$status != "PASS")) stop("At least one target gene lacks a valid unique probe.")

# Prelocked collapse rule: average all valid uniquely mapped probes per gene/sample.
gene_expression <- sapply(target_genes, function(gene) {
  ids <- target_probe_rows$feature_id[target_probe_rows$target_gene == gene & !is.na(target_probe_rows$feature_id)]
  colMeans(expression_selected[ids, , drop = FALSE])
})
gene_expression <- as.data.frame(gene_expression, check.names = FALSE)
gene_expression$sample <- row.names(gene_expression)
gene_expression <- gene_expression[match(metadata$sample, gene_expression$sample), , drop = FALSE]
if (anyNA(gene_expression[, target_genes, drop = FALSE])) stop("Gene-level collapse produced missing values.")

z_parameters <- do.call(rbind, lapply(target_genes, function(gene) {
  data.frame(gene = gene, mean_log2_expression = mean(gene_expression[[gene]]),
             sd_log2_expression = sd(gene_expression[[gene]]), n = nrow(gene_expression))
}))
write_csv(z_parameters, "GSE42387_gene_Z_standardization_parameters.csv")
score_data <- cbind(metadata, gene_expression[, target_genes, drop = FALSE])
for (gene in target_genes) score_data[[paste0("Z_", gene)]] <- z_score(score_data[[gene]])
score_data$P <- rowMeans(score_data[, paste0("Z_", pressure_genes), drop = FALSE])
score_data$B_NQO1 <- rowMeans(score_data[, paste0("Z_", buffer_genes_nqo1), drop = FALSE])
score_data$B_NFE2L2 <- rowMeans(score_data[, paste0("Z_", buffer_genes_nfe2l2), drop = FALSE])
score_data$I_NQO1 <- score_data$P - score_data$B_NQO1
score_data$I_NFE2L2 <- score_data$P - score_data$B_NFE2L2
write_csv(score_data, "GSE42387_PBI_analysis_dataset.csv")

# Nine-group design shared by limma and ordinary composite-score models.
design <- model.matrix(~ 0 + group, data = metadata)
colnames(design) <- levels(metadata$group)
row.names(design) <- metadata$sample
if (qr(design)$rank != ncol(design)) stop("Nine-group design is not full rank.")

contrast_vector <- function(...) {
  x <- setNames(rep(0, ncol(design)), colnames(design))
  values <- list(...)
  for (item in values) x[names(item)] <- x[names(item)] + as.numeric(item)
  x
}
v <- function(group, weight) setNames(weight, group)

main_contrasts <- cbind(
  HCT116_OxR_vs_Parental = contrast_vector(v("HCT116_OxR", 1), v("HCT116_Parental", -1)),
  HT29_OxR_vs_Parental = contrast_vector(v("HT29_OxR", 1), v("HT29_Parental", -1)),
  LoVo_OxR_vs_Parental = contrast_vector(v("LoVo_OxR", 1), v("LoVo_Parental", -1)),
  OxR_average_3_cell_lines = contrast_vector(
    v("HCT116_OxR", 1/3), v("HCT116_Parental", -1/3),
    v("HT29_OxR", 1/3), v("HT29_Parental", -1/3),
    v("LoVo_OxR", 1/3), v("LoVo_Parental", -1/3)),
  HCT116_SN38R_vs_Parental = contrast_vector(v("HCT116_SN38R", 1), v("HCT116_Parental", -1)),
  HT29_SN38R_vs_Parental = contrast_vector(v("HT29_SN38R", 1), v("HT29_Parental", -1)),
  LoVo_SN38R_vs_Parental = contrast_vector(v("LoVo_SN38R", 1), v("LoVo_Parental", -1)),
  SN38R_average_3_cell_lines = contrast_vector(
    v("HCT116_SN38R", 1/3), v("HCT116_Parental", -1/3),
    v("HT29_SN38R", 1/3), v("HT29_Parental", -1/3),
    v("LoVo_SN38R", 1/3), v("LoVo_Parental", -1/3))
)
row.names(main_contrasts) <- colnames(design)

interaction_contrasts <- cbind(
  OxR_HCT116_vs_HT29 = contrast_vector(v("HCT116_OxR", 1), v("HCT116_Parental", -1), v("HT29_OxR", -1), v("HT29_Parental", 1)),
  OxR_HCT116_vs_LoVo = contrast_vector(v("HCT116_OxR", 1), v("HCT116_Parental", -1), v("LoVo_OxR", -1), v("LoVo_Parental", 1)),
  SN38R_HCT116_vs_HT29 = contrast_vector(v("HCT116_SN38R", 1), v("HCT116_Parental", -1), v("HT29_SN38R", -1), v("HT29_Parental", 1)),
  SN38R_HCT116_vs_LoVo = contrast_vector(v("HCT116_SN38R", 1), v("HCT116_Parental", -1), v("LoVo_SN38R", -1), v("LoVo_Parental", 1)),
  DrugSpecificity_HCT116_vs_HT29 = contrast_vector(v("HCT116_OxR", 1), v("HCT116_SN38R", -1), v("HT29_OxR", -1), v("HT29_SN38R", 1)),
  DrugSpecificity_HCT116_vs_LoVo = contrast_vector(v("HCT116_OxR", 1), v("HCT116_SN38R", -1), v("LoVo_OxR", -1), v("LoVo_SN38R", 1))
)
row.names(interaction_contrasts) <- colnames(design)
omnibus_sets <- list(
  OxR_cell_line_heterogeneity = c("OxR_HCT116_vs_HT29", "OxR_HCT116_vs_LoVo"),
  SN38R_cell_line_heterogeneity = c("SN38R_HCT116_vs_HT29", "SN38R_HCT116_vs_LoVo"),
  Drug_specificity_cell_line_heterogeneity = c("DrugSpecificity_HCT116_vs_HT29", "DrugSpecificity_HCT116_vs_LoVo"),
  Global_cell_line_by_treatment_interaction = c("OxR_HCT116_vs_HT29", "OxR_HCT116_vs_LoVo", "SN38R_HCT116_vs_HT29", "SN38R_HCT116_vs_LoVo")
)

write_csv(metadata, "GSE42387_sample_metadata.csv")
write.csv(design, file.path(out_dir, "GSE42387_nine_group_design_matrix.csv"), row.names = TRUE)
write.csv(main_contrasts, file.path(out_dir, "GSE42387_main_contrast_matrix.csv"), row.names = TRUE)
write.csv(interaction_contrasts, file.path(out_dir, "GSE42387_interaction_contrast_matrix.csv"), row.names = TRUE)

# Full-feature and target-level limma empirical-Bayes analysis.
fit <- limma::lmFit(expression_selected, design)
fit_main <- limma::eBayes(limma::contrasts.fit(fit, main_contrasts), trend = TRUE, robust = TRUE)
fit_interaction <- limma::eBayes(limma::contrasts.fit(fit, interaction_contrasts), trend = TRUE, robust = TRUE)

add_annotation <- function(tab) {
  tab$feature_id <- row.names(tab)
  ann <- feature_annotation[match(tab$feature_id, feature_annotation$feature_id), , drop = FALSE]
  cbind(ann, tab[, setdiff(names(tab), "feature_id"), drop = FALSE])
}

main_tables <- lapply(colnames(main_contrasts), function(name) {
  tab <- limma::topTable(fit_main, coef = name, number = Inf, sort.by = "none", adjust.method = "BH", confint = TRUE)
  tab <- add_annotation(tab)
  tab$contrast <- name
  tab$ci_lower_log2 <- tab$CI.L
  tab$ci_upper_log2 <- tab$CI.R
  write_csv_gz(tab, paste0("GSE42387_limma_", name, "_all_features.csv.gz"))
  tab
})
names(main_tables) <- colnames(main_contrasts)

target_main <- do.call(rbind, lapply(names(main_tables), function(contrast_name) {
  tab <- main_tables[[contrast_name]]
  do.call(rbind, lapply(target_genes, function(gene) {
    ids <- target_probe_rows$feature_id[target_probe_rows$target_gene == gene & !is.na(target_probe_rows$feature_id)]
    selected <- tab[tab$feature_id %in% ids, , drop = FALSE]
    selected$target_gene <- gene
    selected$contrast <- contrast_name
    selected
  }))
}))
row.names(target_main) <- NULL
target_main$p_adjusted_BH_11_genes_within_contrast <- ave(
  target_main$P.Value, target_main$contrast, FUN = function(x) p.adjust(x, method = "BH")
)
write_csv(target_main, "GSE42387_11_genes_limma_main_contrasts.csv")

target_omnibus <- do.call(rbind, lapply(names(omnibus_sets), function(test_name) {
  coef_names <- omnibus_sets[[test_name]]
  tab <- limma::topTable(fit_interaction, coef = match(coef_names, colnames(interaction_contrasts)),
                         number = Inf, sort.by = "none", adjust.method = "BH")
  tab <- add_annotation(tab)
  write_csv_gz(tab, paste0("GSE42387_limma_", test_name, "_all_features.csv.gz"))
  selected <- do.call(rbind, lapply(target_genes, function(gene) {
    ids <- target_probe_rows$feature_id[target_probe_rows$target_gene == gene & !is.na(target_probe_rows$feature_id)]
    out <- tab[tab$feature_id %in% ids, , drop = FALSE]
    out$target_gene <- gene
    common_columns <- c(
      "feature_id", "platform_id_matched", "accession", "entrez_gene_id",
      "gene_symbol_raw", "gene_name", "gene_symbol_clean", "ambiguous_multi_gene",
      "AveExpr", "F", "P.Value", "adj.P.Val", "target_gene"
    )
    out[, intersect(common_columns, names(out)), drop = FALSE]
  }))
  selected$test <- test_name
  selected
}))
row.names(target_omnibus) <- NULL
target_omnibus$p_adjusted_BH_11_genes_within_test <- ave(
  target_omnibus$P.Value, target_omnibus$test, FUN = function(x) p.adjust(x, method = "BH")
)
write_csv(target_omnibus, "GSE42387_11_genes_limma_omnibus_tests.csv")

de_summary <- do.call(rbind, lapply(names(main_tables), function(name) {
  tab <- main_tables[[name]]
  data.frame(contrast = name, features_tested = nrow(tab),
             FDR_below_0_05 = sum(tab$adj.P.Val < 0.05, na.rm = TRUE),
             positive_FDR_below_0_05 = sum(tab$adj.P.Val < 0.05 & tab$logFC > 0, na.rm = TRUE),
             negative_FDR_below_0_05 = sum(tab$adj.P.Val < 0.05 & tab$logFC < 0, na.rm = TRUE))
}))
write_csv(de_summary, "GSE42387_full_matrix_limma_summary.csv")

# Ordinary fixed-effects models for composite scores.
estimate_lm_contrast <- function(fit, contrast, name) {
  beta <- coef(fit)
  contrast <- contrast[names(beta)]
  estimate <- as.numeric(sum(contrast * beta))
  se <- as.numeric(sqrt(t(contrast) %*% vcov(fit) %*% contrast))
  df <- df.residual(fit)
  critical <- qt(0.975, df)
  data.frame(contrast = name, estimate = estimate, standard_error = se,
             ci_lower = estimate - critical * se, ci_upper = estimate + critical * se,
             t_statistic = estimate / se, degrees_freedom = df,
             p_value = 2 * pt(abs(estimate / se), df, lower.tail = FALSE))
}

estimate_lm_contrast_hc3 <- function(fit, contrast, name) {
  x <- model.matrix(fit)
  residual <- residuals(fit)
  leverage <- hatvalues(fit)
  bread <- solve(crossprod(x))
  adjusted_squared_residual <- residual^2 / (1 - leverage)^2
  meat <- crossprod(x, x * adjusted_squared_residual)
  covariance_hc3 <- bread %*% meat %*% bread
  beta <- coef(fit)
  contrast <- contrast[names(beta)]
  estimate <- as.numeric(sum(contrast * beta))
  se <- as.numeric(sqrt(t(contrast) %*% covariance_hc3 %*% contrast))
  critical <- qnorm(0.975)
  data.frame(
    contrast = name, estimate = estimate, robust_standard_error_HC3 = se,
    ci_lower_HC3 = estimate - critical * se, ci_upper_HC3 = estimate + critical * se,
    z_statistic_HC3 = estimate / se,
    p_value_HC3 = 2 * pnorm(abs(estimate / se), lower.tail = FALSE)
  )
}

brown_forsythe <- function(values, groups) {
  med <- tapply(values, groups, median)
  dev <- abs(values - med[as.character(groups)])
  tab <- anova(lm(dev ~ groups))
  c(F = tab$`F value`[[1L]], p = tab$`Pr(>F)`[[1L]])
}

omnibus_test <- function(values, name) {
  keep <- switch(name,
    OxR_cell_line_heterogeneity = metadata$status %in% c("Parental", "OxR"),
    SN38R_cell_line_heterogeneity = metadata$status %in% c("Parental", "SN38R"),
    Drug_specificity_cell_line_heterogeneity = metadata$status %in% c("OxR", "SN38R"),
    Global_cell_line_by_treatment_interaction = rep(TRUE, length(values)))
  d <- data.frame(value = values[keep], cell_line = droplevels(metadata$cell_line[keep]),
                  status = droplevels(metadata$status[keep]))
  additive <- lm(value ~ cell_line + status, data = d)
  interaction <- lm(value ~ cell_line * status, data = d)
  tab <- anova(additive, interaction)
  data.frame(test = name, numerator_df = tab$Df[[2L]], denominator_df = df.residual(interaction),
             F_statistic = tab$F[[2L]], p_value = tab$`Pr(>F)`[[2L]], n = nrow(d))
}

analyze_composite <- function(score) {
  d <- data.frame(value = score_data[[score]], design, check.names = FALSE)
  fit_score <- lm(reformulate(colnames(design), response = "value", intercept = FALSE), data = d)
  contrasts <- do.call(rbind, lapply(colnames(main_contrasts), function(name) {
    estimate_lm_contrast(fit_score, main_contrasts[, name], name)
  }))
  contrasts$score <- score
  robust_contrasts <- do.call(rbind, lapply(colnames(main_contrasts), function(name) {
    estimate_lm_contrast_hc3(fit_score, main_contrasts[, name], name)
  }))
  robust_contrasts$score <- score
  omnibus <- do.call(rbind, lapply(names(omnibus_sets), function(name) omnibus_test(score_data[[score]], name)))
  omnibus$score <- score
  residual <- residuals(fit_score)
  bf <- brown_forsythe(score_data[[score]], metadata$group)
  diagnostics <- data.frame(
    score = score, residual_shapiro_W = unname(shapiro.test(residual)$statistic),
    residual_shapiro_p = shapiro.test(residual)$p.value,
    brown_forsythe_F = unname(bf[["F"]]), brown_forsythe_p = unname(bf[["p"]]),
    maximum_cooks_distance = max(cooks.distance(fit_score)), cooks_threshold_4_over_n = 4 / nrow(d),
    R_squared = summary(fit_score)$r.squared, adjusted_R_squared = summary(fit_score)$adj.r.squared
  )
  list(
    model = fit_score, contrasts = contrasts, robust_contrasts = robust_contrasts,
    omnibus = omnibus, diagnostics = diagnostics
  )
}

score_names <- c("P", "B_NQO1", "I_NQO1", "B_NFE2L2", "I_NFE2L2")
score_models <- setNames(lapply(score_names, analyze_composite), score_names)
score_contrasts <- do.call(rbind, lapply(score_models, `[[`, "contrasts"))
score_robust_contrasts <- do.call(rbind, lapply(score_models, `[[`, "robust_contrasts"))
score_omnibus <- do.call(rbind, lapply(score_models, `[[`, "omnibus"))
score_diagnostics <- do.call(rbind, lapply(score_models, `[[`, "diagnostics"))
row.names(score_contrasts) <- row.names(score_robust_contrasts) <- NULL
row.names(score_omnibus) <- row.names(score_diagnostics) <- NULL

score_contrasts$analysis_role <- ifelse(score_contrasts$score %in% c("P", "B_NQO1", "I_NQO1"),
                                        "Primary framework", "Sensitivity")
score_contrasts$p_adjusted_BH_prespecified_family <- NA_real_
primary_average <- score_contrasts$analysis_role == "Primary framework" & score_contrasts$contrast %in% c("OxR_average_3_cell_lines", "SN38R_average_3_cell_lines")
primary_cell <- score_contrasts$analysis_role == "Primary framework" & !grepl("average", score_contrasts$contrast)
sensitivity_average <- score_contrasts$analysis_role == "Sensitivity" & score_contrasts$contrast %in% c("OxR_average_3_cell_lines", "SN38R_average_3_cell_lines")
sensitivity_cell <- score_contrasts$analysis_role == "Sensitivity" & !grepl("average", score_contrasts$contrast)
for (idx in list(primary_average, primary_cell, sensitivity_average, sensitivity_cell)) {
  score_contrasts$p_adjusted_BH_prespecified_family[idx] <- p.adjust(score_contrasts$p_value[idx], method = "BH")
}

score_omnibus$analysis_role <- ifelse(score_omnibus$score %in% c("P", "B_NQO1", "I_NQO1"),
                                      "Primary framework", "Sensitivity")
score_omnibus$p_adjusted_BH_prespecified_family <- NA_real_
for (role in unique(score_omnibus$analysis_role)) {
  idx <- score_omnibus$analysis_role == role
  score_omnibus$p_adjusted_BH_prespecified_family[idx] <- p.adjust(score_omnibus$p_value[idx], method = "BH")
}

write_csv(score_contrasts, "GSE42387_PBI_fixed_effects_contrasts.csv")
score_robust_contrasts$analysis_role <- ifelse(
  score_robust_contrasts$score %in% c("P", "B_NQO1", "I_NQO1"),
  "Primary framework", "Sensitivity"
)
score_robust_contrasts$p_adjusted_BH_HC3_prespecified_family <- NA_real_
primary_average_hc3 <- score_robust_contrasts$analysis_role == "Primary framework" &
  score_robust_contrasts$contrast %in% c("OxR_average_3_cell_lines", "SN38R_average_3_cell_lines")
primary_cell_hc3 <- score_robust_contrasts$analysis_role == "Primary framework" &
  !grepl("average", score_robust_contrasts$contrast)
sensitivity_average_hc3 <- score_robust_contrasts$analysis_role == "Sensitivity" &
  score_robust_contrasts$contrast %in% c("OxR_average_3_cell_lines", "SN38R_average_3_cell_lines")
sensitivity_cell_hc3 <- score_robust_contrasts$analysis_role == "Sensitivity" &
  !grepl("average", score_robust_contrasts$contrast)
for (idx in list(primary_average_hc3, primary_cell_hc3, sensitivity_average_hc3, sensitivity_cell_hc3)) {
  score_robust_contrasts$p_adjusted_BH_HC3_prespecified_family[idx] <- p.adjust(
    score_robust_contrasts$p_value_HC3[idx], method = "BH"
  )
}
write_csv(score_robust_contrasts, "GSE42387_PBI_HC3_robust_sensitivity_contrasts.csv")
write_csv(score_omnibus, "GSE42387_PBI_fixed_effects_omnibus_tests.csv")
write_csv(score_diagnostics, "GSE42387_PBI_model_diagnostics.csv")

# Descriptive group summaries.
score_summary <- do.call(rbind, lapply(score_names, function(score) {
  do.call(rbind, lapply(levels(metadata$group), function(group_name) {
    x <- score_data[[score]][metadata$group == group_name]
    data.frame(score = score, group = group_name, n = length(x), mean = mean(x), sd = sd(x),
               median = median(x), min = min(x), max = max(x))
  }))
}))
write_csv(score_summary, "GSE42387_PBI_group_descriptive_summary.csv")

# Publication figures for primary P/B/I scores and target-gene average effects.
forest_contrasts <- c(
  "HCT116_OxR_vs_Parental", "HT29_OxR_vs_Parental", "LoVo_OxR_vs_Parental", "OxR_average_3_cell_lines",
  "HCT116_SN38R_vs_Parental", "HT29_SN38R_vs_Parental", "LoVo_SN38R_vs_Parental", "SN38R_average_3_cell_lines"
)
forest_labels <- c("HCT116 OxR", "HT29 OxR", "LoVo OxR", "OxR average",
                   "HCT116 SN38R", "HT29 SN38R", "LoVo SN38R", "SN38R average")

draw_score_forest <- function() {
  old <- par(no.readonly = TRUE); on.exit(par(old))
  par(mfrow = c(1, 3), mar = c(4.3, 6.7, 2.2, 0.7), oma = c(0, 0, 1.0, 0), family = "sans")
  palette <- c(rep("#D4743C", 4), rep("#1874A6", 4))
  score_labels <- c(P = "Pressure (P)", B_NQO1 = "Buffering (B)", I_NQO1 = "P-B index")
  for (score in c("P", "B_NQO1", "I_NQO1")) {
    d <- score_contrasts[score_contrasts$score == score, , drop = FALSE]
    d <- d[match(forest_contrasts, d$contrast), , drop = FALSE]
    y <- rev(seq_len(nrow(d)))
    xlim <- range(c(d$ci_lower, d$ci_upper, 0), finite = TRUE)
    plot(d$estimate, y, xlim = xlim, ylim = c(0.5, 8.5), xlab = "Resistant - parental\n(mean-Z-score units)",
         ylab = "", yaxt = "n", pch = 15, col = palette, bty = "n", main = score_labels[[score]])
    abline(v = 0, lty = 2, col = "#777777")
    segments(d$ci_lower, y, d$ci_upper, y, lwd = 1.4)
    axis(2, at = y, labels = forest_labels, las = 1, tick = FALSE, cex.axis = 0.75)
  }
}

draw_gene_average_forest <- function() {
  old <- par(no.readonly = TRUE); on.exit(par(old))
  par(mfrow = c(1, 2), mar = c(4.3, 5.8, 2.2, 0.8), family = "sans")
  for (contrast_name in c("OxR_average_3_cell_lines", "SN38R_average_3_cell_lines")) {
    d <- target_main[target_main$contrast == contrast_name, , drop = FALSE]
    d <- d[match(target_genes, d$target_gene), , drop = FALSE]
    y <- rev(seq_len(nrow(d)))
    xlim <- range(c(d$ci_lower_log2, d$ci_upper_log2, 0), finite = TRUE)
    plot(d$logFC, y, xlim = xlim, ylim = c(0.5, length(y) + 0.5),
         xlab = "Average resistant - parental (log2)", ylab = "", yaxt = "n", pch = 15,
         col = ifelse(d$target_gene %in% pressure_genes, "#D4743C", "#1874A6"),
         bty = "n", main = sub("_average_3_cell_lines", " average", contrast_name))
    abline(v = 0, lty = 2, col = "#777777")
    segments(d$ci_lower_log2, y, d$ci_upper_log2, y, lwd = 1.4)
    axis(2, at = y, labels = d$target_gene, las = 1, tick = FALSE, cex.axis = 0.8)
  }
}

# Per-gene distribution figures retain the balanced nine-group structure.
# The 10 primary genes are shown separately from the NFE2L2 sensitivity gene.
primary_target_genes <- c(pressure_genes, buffer_genes_nqo1)
gene_panel_long <- do.call(rbind, lapply(target_genes, function(gene) {
  data.frame(
    sample = score_data$sample,
    sample_title = score_data$title,
    cell_line = as.character(score_data$cell_line),
    status = as.character(score_data$status),
    gene = gene,
    analysis_role = if (gene == "NFE2L2") "Sensitivity" else "Primary 10-gene framework",
    expression_log2 = score_data[[gene]],
    stringsAsFactors = FALSE
  )
}))
row.names(gene_panel_long) <- NULL
write_csv(gene_panel_long, "GSE42387_11_gene_distribution_source_data.csv")

cell_line_contrasts <- c(
  HCT116_OxR_vs_Parental = "HCT116",
  HT29_OxR_vs_Parental = "HT29",
  LoVo_OxR_vs_Parental = "LoVo",
  HCT116_SN38R_vs_Parental = "HCT116",
  HT29_SN38R_vs_Parental = "HT29",
  LoVo_SN38R_vs_Parental = "LoVo"
)
gene_panel_stats <- target_main[
  target_main$contrast %in% names(cell_line_contrasts),
  c("target_gene", "contrast", "logFC", "ci_lower_log2", "ci_upper_log2",
    "P.Value", "p_adjusted_BH_11_genes_within_contrast"),
  drop = FALSE
]
gene_panel_stats$cell_line <- unname(cell_line_contrasts[gene_panel_stats$contrast])
gene_panel_stats$resistant_status <- ifelse(grepl("_OxR_", gene_panel_stats$contrast),
                                             "OxR", "SN38R")
gene_panel_stats$contrast_definition <- paste0(
  gene_panel_stats$cell_line, " ", gene_panel_stats$resistant_status, " minus parental"
)
gene_panel_stats$analysis_role <- ifelse(
  gene_panel_stats$target_gene == "NFE2L2", "Sensitivity", "Primary 10-gene framework"
)
gene_panel_stats <- gene_panel_stats[, c(
  "target_gene", "analysis_role", "cell_line", "resistant_status", "contrast",
  "contrast_definition", "logFC", "ci_lower_log2", "ci_upper_log2", "P.Value",
  "p_adjusted_BH_11_genes_within_contrast"
)]
write_csv(gene_panel_stats, "GSE42387_11_gene_cell_line_limma_statistics.csv")
write_csv(gene_panel_stats[gene_panel_stats$analysis_role == "Primary 10-gene framework", ],
          "GSE42387_10_primary_genes_cell_line_limma_statistics.csv")
write_csv(gene_panel_stats[gene_panel_stats$analysis_role == "Sensitivity", ],
          "GSE42387_NFE2L2_sensitivity_cell_line_limma_statistics.csv")

gene_panel_palette <- c(Parental = "#F04B3A", OxR = "#43B8D3", SN38R = "#43B8D3")

fmt_q <- function(q) {
  if (!is.finite(q)) return("q=NA")
  if (q < 0.001) return("q<0.001")
  paste0("q=", sprintf("%.3f", q))
}

draw_gene_distribution_panel <- function(gene, resistant_status, show_y_label = TRUE) {
  cell_lines <- c("HCT116", "HT29", "LoVo")
  d <- gene_panel_long[
    gene_panel_long$gene == gene &
      gene_panel_long$status %in% c("Parental", resistant_status),
    , drop = FALSE
  ]
  stats <- gene_panel_stats[
    gene_panel_stats$target_gene == gene &
      gene_panel_stats$resistant_status == resistant_status,
    , drop = FALSE
  ]
  stats <- stats[match(cell_lines, stats$cell_line), , drop = FALSE]
  if (nrow(d) != 18L || nrow(stats) != 3L || anyNA(stats$cell_line)) {
    stop("Per-gene distribution panel is missing balanced samples or limma statistics.")
  }

  values <- d$expression_log2
  y_range <- diff(range(values))
  if (!is.finite(y_range) || y_range <= 0) y_range <- 1
  y_min <- min(values) - 0.08 * y_range
  bracket_y <- max(values) + 0.09 * y_range
  y_max <- max(values) + 0.22 * y_range
  plot(NA, xlim = c(0.50, 3.50), ylim = c(y_min, y_max), axes = FALSE,
       xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  y_ticks <- pretty(c(y_min, max(values)), n = 4)
  y_ticks <- y_ticks[y_ticks >= y_min & y_ticks <= max(values)]
  abline(h = y_ticks, col = "#E7EEF1", lwd = 0.7)
  abline(v = c(1.5, 2.5), col = "#D0D0D0", lwd = 0.65)
  axis(1, at = 1:3, labels = cell_lines, line = -0.65, padj = 0.15,
       cex.axis = 0.70, col.axis = "#111111", lwd = 0, lwd.ticks = 0)
  axis(2, at = y_ticks, las = 1, cex.axis = 0.70)
  if (show_y_label) mtext("Normalized log2 expression", side = 2, line = 2.7, cex = 0.74)
  box(bty = "l", col = "#555555", lwd = 0.7)

  for (i in seq_along(cell_lines)) {
    cell_line <- cell_lines[[i]]
    positions <- c(Parental = i - 0.16, resistant = i + 0.16)
    parental <- d$expression_log2[d$cell_line == cell_line & d$status == "Parental"]
    resistant <- d$expression_log2[d$cell_line == cell_line & d$status == resistant_status]
    for (status_name in c("Parental", resistant_status)) {
      vals <- if (status_name == "Parental") parental else resistant
      position <- if (status_name == "Parental") positions[["Parental"]] else positions[["resistant"]]
      boxplot(vals, at = position, add = TRUE, axes = FALSE, outline = FALSE,
              boxwex = 0.30, border = gene_panel_palette[[status_name]],
              col = adjustcolor(gene_panel_palette[[status_name]], alpha.f = 0.04),
              lwd = 0.9, pars = list(medlwd = 1.3, whisklty = 1, staplewex = 0.45))
      set.seed(4200 + match(gene, target_genes) * 100 + i * 10 +
                 match(status_name, c("Parental", resistant_status)))
      points(position + runif(length(vals), -0.042, 0.042), vals, pch = 16, cex = 0.48,
             col = adjustcolor(gene_panel_palette[[status_name]], alpha.f = 0.88))
    }
    bracket_drop <- 0.035 * y_range
    segments(positions[["Parental"]], bracket_y, positions[["resistant"]], bracket_y,
             col = "#555555", lwd = 0.65)
    segments(c(positions[["Parental"]], positions[["resistant"]]), bracket_y,
             c(positions[["Parental"]], positions[["resistant"]]) - bracket_drop,
             col = "#555555", lwd = 0.65)
    text(i, bracket_y + 0.025 * y_range,
         fmt_q(stats$p_adjusted_BH_11_genes_within_contrast[[i]]), cex = 0.59)
  }

  usr <- par("usr")
  strip_height <- 0.12 * diff(usr[3:4])
  strip_fill <- if (gene %in% pressure_genes) "#F7E4D8" else "#DCEBF1"
  rect(usr[[1L]], usr[[4L]], usr[[2L]], usr[[4L]] + strip_height,
       col = strip_fill, border = NA, xpd = NA)
  text(mean(usr[1:2]), usr[[4L]] + strip_height / 2, gene,
       font = 2, cex = 0.76, xpd = NA)
}

draw_10_gene_distribution_grid <- function(resistant_status) {
  old <- par(no.readonly = TRUE)
  on.exit({layout(1); par(old)})
  layout(matrix(seq_len(10L), nrow = 5L, ncol = 2L, byrow = TRUE))
  par(mar = c(1.8, 3.7, 2.2, 0.7), oma = c(0.3, 0.3, 3.2, 0.3), family = "sans")
  for (i in seq_along(primary_target_genes)) {
    draw_gene_distribution_panel(primary_target_genes[[i]], resistant_status,
                                 show_y_label = i %% 2L == 1L)
  }
  mtext(paste0("GSE42387: parental versus ", resistant_status,
               " across three colorectal cancer cell lines"),
        outer = TRUE, side = 3, line = 1.5, cex = 1.02, font = 2)
  mtext("Parental", outer = TRUE, side = 3, line = 0.20, adj = 0.83,
        col = gene_panel_palette[["Parental"]], cex = 0.72, font = 2)
  mtext(resistant_status, outer = TRUE, side = 3, line = 0.20, adj = 0.94,
        col = gene_panel_palette[[resistant_status]], cex = 0.72, font = 2)
}

draw_nfe2l2_sensitivity <- function() {
  old <- par(no.readonly = TRUE)
  on.exit({layout(1); par(old)})
  layout(matrix(1:2, nrow = 1L))
  par(mar = c(3.1, 3.8, 2.2, 0.7), oma = c(2.0, 0.3, 3.0, 0.3), family = "sans")
  draw_gene_distribution_panel("NFE2L2", "OxR", show_y_label = TRUE)
  draw_gene_distribution_panel("NFE2L2", "SN38R", show_y_label = TRUE)
  mtext("NFE2L2 sensitivity analysis", outer = TRUE, side = 3, line = 1.3,
        cex = 1.02, font = 2)
  mtext("Left: OxR versus parental; right: SN38R versus parental. q = limma BH-adjusted p value within each 11-gene contrast.",
        outer = TRUE, side = 1, line = 0.5, cex = 0.68)
}

draw_score_distributions <- function() {
  old <- par(no.readonly = TRUE); on.exit(par(old))
  par(mfrow = c(1, 3), mar = c(7.5, 4.4, 2.0, 0.7), oma = c(0, 0, 1.0, 0), family = "sans")
  colors <- rep(c("#B9B9B9", "#D4743C", "#1874A6"), 3)
  score_labels <- c(P = "Pressure (P)", B_NQO1 = "Buffering (B)", I_NQO1 = "P-B index")
  for (score in c("P", "B_NQO1", "I_NQO1")) {
    boxplot(score_data[[score]] ~ metadata$group, names = levels(metadata$group), las = 2,
            ylab = "Composite score (mean gene Z)", xlab = "", main = score_labels[[score]],
            col = colors, border = "#555555", outline = FALSE, cex.axis = 0.72)
    stripchart(score_data[[score]] ~ metadata$group, vertical = TRUE, method = "jitter",
               pch = 16, cex = 0.65, col = "#222222", add = TRUE)
    abline(h = 0, lty = 3, col = "#777777")
  }
}

diagnostic_score_labels <- c(
  P = "Pressure (P)",
  B_NQO1 = "Buffering (B)",
  I_NQO1 = "P-B index"
)

diagnostic_flagged_samples <- do.call(rbind, lapply(names(diagnostic_score_labels), function(score) {
  fitted_model <- score_models[[score]]$model
  sample_names <- row.names(model.frame(fitted_model))
  cook <- cooks.distance(fitted_model)
  selected <- head(order(cook, decreasing = TRUE), 3L)
  data.frame(
    score = score,
    score_label = diagnostic_score_labels[[score]],
    label_number = seq_along(selected),
    sample = sample_names[selected],
    sample_title = metadata$title[match(sample_names[selected], metadata$sample)],
    cooks_distance = unname(cook[selected]),
    stringsAsFactors = FALSE
  )
}))
row.names(diagnostic_flagged_samples) <- NULL
write_csv(diagnostic_flagged_samples, "GSE42387_PBI_flagged_diagnostic_samples.csv")

diagnostic_data <- function(score) {
  fitted_model <- score_models[[score]]$model
  sample_names <- row.names(model.frame(fitted_model))
  selected_table <- diagnostic_flagged_samples[
    diagnostic_flagged_samples$score == score,
    ,
    drop = FALSE
  ]
  list(
    model = fitted_model,
    sample_names = sample_names,
    fitted = fitted(fitted_model),
    residual = residuals(fitted_model),
    standardized_residual = rstandard(fitted_model),
    cooks_distance = cooks.distance(fitted_model),
    selected = match(selected_table$sample, sample_names),
    selected_table = selected_table
  )
}

label_flagged_points <- function(x, y, selected) {
  positions <- rep(c(2, 3, 4), length.out = length(selected))
  text(
    x[selected], y[selected], labels = seq_along(selected),
    pos = positions, offset = 0.55, cex = 0.80, font = 2,
    col = "#0072B2", xpd = NA
  )
}

draw_diagnostic_panel <- function(diagnostic, panel) {
  if (panel == "residuals") {
    x <- diagnostic$fitted
    y <- diagnostic$residual
    plot(
      x, y, pch = 1, cex = 0.82, col = "#333333",
      xlab = "Fitted values", ylab = "Residuals",
      main = "Residuals vs fitted", bty = "l"
    )
    abline(h = 0, lty = 3, col = "#777777")
    lines(lowess(x, y, f = 2 / 3), col = "#D55E00", lwd = 1.4)
    label_flagged_points(x, y, diagnostic$selected)
  } else if (panel == "qq") {
    qq <- qqnorm(diagnostic$standardized_residual, plot.it = FALSE)
    plot(
      qq$x, qq$y, pch = 1, cex = 0.82, col = "#333333",
      xlab = "Theoretical quantiles", ylab = "Standardized residuals",
      main = "Normal Q-Q", bty = "l"
    )
    qqline(diagnostic$standardized_residual, lty = 3, col = "#777777", lwd = 1.2)
    # qqnorm() keeps y in observation order and ranks x to match each observation.
    label_flagged_points(qq$x, qq$y, diagnostic$selected)
  } else if (panel == "scale") {
    x <- diagnostic$fitted
    y <- sqrt(abs(diagnostic$standardized_residual))
    plot(
      x, y, pch = 1, cex = 0.82, col = "#333333",
      xlab = "Fitted values", ylab = expression(sqrt("|Standardized residuals|")),
      main = "Scale-location", bty = "l"
    )
    lines(lowess(x, y, f = 2 / 3), col = "#D55E00", lwd = 1.4)
    label_flagged_points(x, y, diagnostic$selected)
  } else if (panel == "cook") {
    x <- seq_along(diagnostic$cooks_distance)
    y <- diagnostic$cooks_distance
    y_limit <- max(c(y * 1.17, 4 / length(y)), na.rm = TRUE)
    plot(
      x, y, type = "h", lwd = 1.2, col = "#333333",
      xlab = "Observation index", ylab = "Cook distance",
      main = "Influence", ylim = c(0, y_limit), bty = "l"
    )
    abline(h = 4 / length(y), lty = 3, col = "#D55E00")
    text(
      x[diagnostic$selected], y[diagnostic$selected],
      labels = seq_along(diagnostic$selected), pos = 3, offset = 0.35,
      cex = 0.80, font = 2, col = "#0072B2", xpd = NA
    )
  }
}

draw_diagnostic_key <- function(diagnostic, include_score_title = FALSE) {
  plot.new()
  if (include_score_title) {
    text(
      0, 0.92, diagnostic_score_labels[[diagnostic$selected_table$score[[1L]]]],
      adj = c(0, 1), font = 2, cex = 0.92
    )
  }
  text(0, 0.76, "Flagged by Cook distance", adj = c(0, 1), cex = 0.73, col = "#555555")
  y_positions <- c(0.59, 0.40, 0.21)
  for (i in seq_len(nrow(diagnostic$selected_table))) {
    row <- diagnostic$selected_table[i, , drop = FALSE]
    text(
      0, y_positions[[i]],
      sprintf("%d = %s", row$label_number, row$sample),
      adj = c(0, 0.5), cex = 0.74, font = 2, col = "#0072B2"
    )
    text(
      0, y_positions[[i]] - 0.09, row$sample_title,
      adj = c(0, 0.5), cex = 0.64, col = "#555555"
    )
  }
}

draw_diagnostic_key_wide <- function(diagnostic) {
  par(mar = c(0.2, 0.8, 0.2, 0.8), family = "sans")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  text(
    0.02, 0.83, "Flagged by Cook distance",
    adj = c(0, 0.5), cex = 0.72, col = "#555555"
  )
  x_positions <- c(0.02, 0.35, 0.68)
  for (i in seq_len(nrow(diagnostic$selected_table))) {
    row <- diagnostic$selected_table[i, , drop = FALSE]
    text(
      x_positions[[i]], 0.54,
      sprintf("%d = %s", row$label_number, row$sample),
      adj = c(0, 0.5), cex = 0.77, font = 2, col = "#0072B2"
    )
    text(
      x_positions[[i]], 0.25, row$sample_title,
      adj = c(0, 0.5), cex = 0.67, col = "#555555"
    )
  }
}

draw_diagnostics_overview <- function() {
  old <- par(no.readonly = TRUE)
  on.exit({layout(1); par(old)})
  layout(matrix(seq_len(15L), nrow = 3L, byrow = TRUE), widths = c(1, 1, 1, 1, 0.82))
  par(mar = c(3.8, 4.0, 2.0, 0.8), oma = c(0, 0, 1.8, 0), family = "sans")
  for (score in names(diagnostic_score_labels)) {
    diagnostic <- diagnostic_data(score)
    for (panel in c("residuals", "qq", "scale", "cook")) {
      draw_diagnostic_panel(diagnostic, panel)
    }
    draw_diagnostic_key(diagnostic, include_score_title = TRUE)
  }
  mtext("Ordinary fixed-effects model diagnostics", outer = TRUE, line = 0.3, cex = 1.05, font = 2)
}

draw_single_score_diagnostics <- function(score) {
  old <- par(no.readonly = TRUE)
  on.exit({layout(1); par(old)})
  layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3L, byrow = TRUE), heights = c(1, 1, 0.38))
  par(mar = c(4.2, 4.5, 2.2, 1.0), oma = c(0, 0, 2.0, 0), family = "sans")
  diagnostic <- diagnostic_data(score)
  for (panel in c("residuals", "qq", "scale", "cook")) {
    draw_diagnostic_panel(diagnostic, panel)
  }
  draw_diagnostic_key_wide(diagnostic)
  mtext(
    paste0(diagnostic_score_labels[[score]], ": ordinary fixed-effects diagnostics"),
    outer = TRUE, line = 0.3, cex = 1.05, font = 2
  )
}

save_figure_formats("GSE42387_PBI_fixed_effects_forest", 7.0, 4.8, draw_score_forest)
save_figure_formats("GSE42387_11_genes_average_effect_forest", 7.0, 5.8, draw_gene_average_forest)
save_figure_formats("GSE42387_10_primary_genes_OxR_distributions", 8.0, 13.0,
                    function() draw_10_gene_distribution_grid("OxR"))
save_figure_formats("GSE42387_10_primary_genes_SN38R_distributions", 8.0, 13.0,
                    function() draw_10_gene_distribution_grid("SN38R"))
save_figure_formats("GSE42387_NFE2L2_sensitivity_distributions", 8.0, 3.8,
                    draw_nfe2l2_sensitivity)
save_figure_formats("GSE42387_PBI_group_distributions", 7.0, 5.6, draw_score_distributions)
save_figure_formats("GSE42387_PBI_model_diagnostics", 12.0, 8.5, draw_diagnostics_overview)
for (score in names(diagnostic_score_labels)) {
  output_score_name <- c(P = "P", B_NQO1 = "B", I_NQO1 = "P_minus_B")[[score]]
  drawing_function <- local({
    score_name <- score
    function() draw_single_score_diagnostics(score_name)
  })
  save_figure_formats(
    paste0("GSE42387_", output_score_name, "_model_diagnostics"),
    7.0, 7.4, drawing_function
  )
}

validation <- data.frame(
  metric = c("features", "arrays", "cell_lines", "statuses", "minimum_replicates_per_cell",
             "missing_expression_values", "target_genes", "target_genes_with_valid_probe",
             "target_genes_with_multiple_valid_probes", "design_rank", "design_columns"),
  value = c(nrow(expression_selected), ncol(expression_selected), nlevels(metadata$cell_line),
            nlevels(metadata$status), min(cell_counts), sum(is.na(expression_selected)),
            length(target_genes), sum(probe_coverage$status == "PASS"),
            sum(probe_coverage$n_valid_unique_probes > 1L), qr(design)$rank, ncol(design))
)
write_csv(validation, "GSE42387_structural_validation.csv")

input_files <- c(series_file, platform_file)
write_csv(data.frame(file = basename(input_files), md5 = unname(tools::md5sum(input_files))),
          "GSE42387_input_file_md5.csv")
write_csv(data.frame(
  package = c("R", "limma", "statmod"),
  version = c(paste(R.version$major, R.version$minor, sep = "."),
              as.character(packageVersion("limma")), as.character(packageVersion("statmod")))
), "GSE42387_package_versions.csv")
capture.output(sessionInfo(), file = file.path(out_dir, "GSE42387_R_sessionInfo.txt"))

methods_text <- c(
  "GSE42387 pressure-buffering analysis (locked 2026-08-12)",
  "All 27 processed arrays were retained in a balanced nine-group design comprising parental, oxaliplatin-resistant, and SN-38-resistant states in HCT116, HT29, and LoVo cells.",
  "Probes without gene annotation or with ambiguous multi-gene annotation were excluded. If multiple valid probes uniquely mapped to one target gene, probe-level log2 expression values were averaged by sample before gene-wise standardization.",
  "Target genes were Z-standardized across all 27 arrays. P, B, and P-B were calculated as unweighted mean-Z-score composites.",
  "Full-feature and single-gene analyses used limma empirical-Bayes moderation. Composite scores used ordinary fixed-effects linear models with the same nine-group design and prespecified contrasts.",
  "Raw P values and BH-adjusted P values were reported. GSE42387 was interpreted as exploratory cell-line support, not clinical or functional validation."
)
writeLines(methods_text, file.path(out_dir, "GSE42387_manuscript_methods.txt"), useBytes = TRUE)

manifest_path <- file.path(out_dir, "_analysis_outputs.md")
manifest_connection <- file(manifest_path, open = "wt", encoding = "UTF-8")
writeLines(c(
  "# GSE42387 analysis outputs", "", "Generated: 2026-08-12", "",
  "## Core composite-score tables",
  "- `GSE42387_PBI_fixed_effects_contrasts.csv`: P/B/I estimates and 95% CIs.",
  "- `GSE42387_PBI_HC3_robust_sensitivity_contrasts.csv`: HC3 robust-SE sensitivity estimates.",
  "- `GSE42387_PBI_fixed_effects_omnibus_tests.csv`: interaction and heterogeneity tests.",
  "- `GSE42387_PBI_model_diagnostics.csv`: residual, variance, influence, and fit checks.", ""
), manifest_connection, useBytes = TRUE)
writeLines(c(
  "## Mapping and expression",
  "- `GSE42387_target_probe_coverage_summary.csv`: mandatory pre-score coverage audit.",
  "- `GSE42387_target_probe_mapping_audit.csv`: target probe details.",
  "- `GSE42387_PBI_analysis_dataset.csv`: gene expression, gene Z scores, and composites.", "",
  "## Limma outputs",
  "- Eleven-gene main-contrast and omnibus tables are provided as CSV.",
  "- Full-feature contrast and omnibus results are provided as compressed CSV files.", ""
), manifest_connection, useBytes = TRUE)
writeLines(c(
  "## Figures",
  "- Main P/B/I forest, 11-gene average-effect forest, group distributions, and diagnostics are provided in PDF, PNG, and TIFF.",
  "- `GSE42387_10_primary_genes_OxR_distributions.*`: 10-gene parental-versus-OxR distributions stratified by cell line.",
  "- `GSE42387_10_primary_genes_SN38R_distributions.*`: 10-gene parental-versus-SN38R distributions stratified by cell line.",
  "- `GSE42387_NFE2L2_sensitivity_distributions.*`: separate NFE2L2 sensitivity distributions.",
  "- `GSE42387_11_gene_distribution_source_data.csv` and cell-line limma statistics tables provide figure traceability.", ""
), manifest_connection, useBytes = TRUE)
writeLines(c(
  "- Model diagnostics are provided as one overview and as separate 2 x 2 figures for P, B, and P-B. Full GSM identifiers are listed outside the plotting regions and linked to numbered points selected by Cook distance.", "",
  "## Interpretation",
  "- Composite scores use ordinary fixed-effects models; limma moderation is reserved for feature-level inference."
), manifest_connection, useBytes = TRUE)
close(manifest_connection)

cat("GSE42387 P/B/I analysis complete.\n")
cat("Arrays:", ncol(expression_selected), " Features:", nrow(expression_selected), "\n")
cat("Probe coverage:\n"); print(probe_coverage, row.names = FALSE)
cat("Primary composite average effects:\n")
print(score_contrasts[score_contrasts$analysis_role == "Primary framework" &
                        score_contrasts$contrast %in% c("OxR_average_3_cell_lines", "SN38R_average_3_cell_lines"),
                      c("score", "contrast", "estimate", "ci_lower", "ci_upper", "p_value", "p_adjusted_BH_prespecified_family")],
      row.names = FALSE)
cat("Output directory:", out_dir, "\n")
