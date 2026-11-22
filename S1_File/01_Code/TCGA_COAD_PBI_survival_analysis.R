###############################################################################
# Analysis: TCGA-COAD pressure-buffering composite survival analysis
# Date: 2026-08-12
# Random seed: 42
# R: >= 4.5.0
# Required package: survival
#
# Inputs are read-only. All derived files are written to --outdir.
###############################################################################

set.seed(42)
options(stringsAsFactors = FALSE, scipen = 999)

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("The 'survival' package is required.", call. = FALSE)
}
suppressPackageStartupMessages(library(survival))

parse_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[[i]], "--") || i == length(x)) {
      stop("Arguments must be supplied as --name value pairs.", call. = FALSE)
    }
    out[[substring(x[[i]], 3L)]] <- x[[i + 1L]]
    i <- i + 2L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("expression", "probemap", "clinical-dir", "outdir")
missing_args <- required[vapply(required, function(x) is.null(args[[x]]), logical(1L))]
if (length(missing_args) > 0L) {
  stop("Missing argument(s): --", paste(missing_args, collapse = ", --"), call. = FALSE)
}

expression_path <- normalizePath(args$expression, winslash = "/", mustWork = TRUE)
probemap_path <- normalizePath(args$probemap, winslash = "/", mustWork = TRUE)
clinical_dir <- normalizePath(args$`clinical-dir`, winslash = "/", mustWork = TRUE)
out_dir <- normalizePath(args$outdir, winslash = "/", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)

access_date <- if (is.null(args$`data-access-date`)) "2026-02-22" else args$`data-access-date`
expression_version <- if (is.null(args$`expression-version`)) "05-09-2024" else args$`expression-version`
min_followup <- if (is.null(args$`min-followup-days`)) 1 else as.numeric(args$`min-followup-days`)
if (!is.finite(min_followup) || min_followup < 0) stop("Invalid minimum follow-up.")

pressure_genes <- c("ACSL4", "LPCAT3", "TFRC", "NCOA4")
buffer_genes_nqo1 <- c("SLC7A11", "GPX4", "AIFM2", "GCH1", "DHODH", "NQO1")
buffer_genes_nfe2l2 <- c("SLC7A11", "GPX4", "AIFM2", "GCH1", "DHODH", "NFE2L2")
target_genes <- unique(c(pressure_genes, buffer_genes_nqo1, buffer_genes_nfe2l2))
na_strings <- c("", "--", "'--", "not reported", "Not Reported", "unknown", "Unknown")

write_csv <- function(x, name) {
  write.csv(x, file.path(out_dir, name), row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

require_columns <- function(x, columns, label) {
  absent <- setdiff(columns, names(x))
  if (length(absent) > 0L) stop(label, " missing: ", paste(absent, collapse = ", "))
}

numeric_or_na <- function(x) suppressWarnings(as.numeric(x))

max_or_na <- function(x) {
  x <- numeric_or_na(x)
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

row_max_or_na <- function(a, b) {
  out <- pmax(a, b, na.rm = TRUE)
  out[!is.finite(out)] <- NA_real_
  out
}

single_consistent <- function(x, variable, patient) {
  values <- unique(trimws(as.character(x[!is.na(x)])))
  values <- values[values != ""]
  if (length(values) > 1L) stop("Inconsistent ", variable, " for ", patient)
  if (length(values) == 0L) return(NA_character_)
  values[[1L]]
}

extract_patient_id <- function(x) toupper(substr(trimws(as.character(x)), 1L, 12L))

extract_sample_type <- function(x) {
  vapply(strsplit(toupper(as.character(x)), "-", fixed = TRUE), function(parts) {
    if (length(parts) < 4L) return(NA_character_)
    substr(parts[[4L]], 1L, 2L)
  }, character(1L))
}

collapse_stage <- function(x) {
  x <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  out[grepl("^STAGE II([^I]|$)", x)] <- "II"
  out[grepl("^STAGE III", x)] <- "III"
  out
}

z_score <- function(x) {
  s <- sd(x)
  if (!is.finite(s) || s <= 0) stop("A target gene or composite has zero variance.")
  as.numeric((x - mean(x)) / s)
}

save_figure_formats <- function(stem, width, height, draw) {
  pdf(file.path(out_dir, paste0(stem, ".pdf")), width = width, height = height,
      family = "sans", useDingbats = FALSE)
  draw()
  dev.off()
  png(file.path(out_dir, paste0(stem, ".png")), width = width, height = height,
      units = "in", res = 300, type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
  draw()
  dev.off()
  tiff(file.path(out_dir, paste0(stem, ".tiff")), width = width, height = height,
       units = "in", res = 600, compression = "lzw")
  draw()
  dev.off()
}

read_target_expression <- function(path, id_to_gene) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  header <- strsplit(readLines(con, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]
  if (length(header) < 2L || header[[1L]] != "Ensembl_ID") {
    stop("Unexpected expression-matrix header.")
  }
  samples <- header[-1L]
  found <- list()
  repeat {
    lines <- readLines(con, n = 5000L, warn = FALSE)
    if (length(lines) == 0L) break
    ids <- sub("\t.*$", "", lines)
    for (j in which(ids %in% names(id_to_gene))) {
      parts <- strsplit(lines[[j]], "\t", fixed = TRUE)[[1L]]
      values <- numeric_or_na(parts[-1L])
      if (length(values) != length(samples) || any(!is.finite(values))) {
        stop("Invalid expression row: ", parts[[1L]])
      }
      found[[parts[[1L]]]] <- values
    }
  }
  absent <- setdiff(names(id_to_gene), names(found))
  if (length(absent) > 0L) stop("Target Ensembl IDs absent: ", paste(absent, collapse = ", "))
  out <- data.frame(sample = samples, check.names = FALSE)
  for (id in names(id_to_gene)) out[[id_to_gene[[id]]]] <- found[[id]]
  out
}

# Target mapping is a hard pre-analysis gate.
probemap <- read.delim(probemap_path, stringsAsFactors = FALSE, check.names = FALSE,
                       quote = "", comment.char = "")
require_columns(probemap, c("id", "gene"), "GENCODE v36 probemap")
mapping_all <- probemap[probemap$gene %in% target_genes, c("id", "gene"), drop = FALSE]
coverage <- data.frame(
  gene = target_genes,
  n_probemap_rows = vapply(target_genes, function(g) sum(mapping_all$gene == g), integer(1L)),
  stringsAsFactors = FALSE
)
coverage$status <- ifelse(coverage$n_probemap_rows == 1L, "PASS: unique", "FAIL")
write_csv(coverage, "TCGA_target_gene_coverage_audit.csv")
if (any(coverage$n_probemap_rows != 1L)) stop("TCGA target-gene mapping audit failed.")
mapping <- mapping_all[match(target_genes, mapping_all$gene), , drop = FALSE]
mapping$role <- ifelse(mapping$gene %in% pressure_genes, "Pressure",
                       ifelse(mapping$gene == "NFE2L2", "Buffering sensitivity", "Buffering primary"))
write_csv(mapping, "TCGA_target_gene_mapping.csv")

expression <- read_target_expression(expression_path, setNames(mapping$gene, mapping$id))
expression$sample <- toupper(expression$sample)
expression$patient_id <- extract_patient_id(expression$sample)
expression$sample_type_code <- extract_sample_type(expression$sample)
if (any(expression[, target_genes, drop = FALSE] < 0)) {
  stop("Negative target-gene values conflict with expected log2(TPM + 1) input.")
}

sample_type_counts <- as.data.frame(table(expression$sample_type_code, useNA = "ifany"))
names(sample_type_counts) <- c("sample_type_code", "n_samples")
write_csv(sample_type_counts, "TCGA_sample_type_counts.csv")

primary_expression <- expression[expression$sample_type_code == "01", , drop = FALSE]
patient_expression <- do.call(rbind, lapply(split(primary_expression, primary_expression$patient_id), function(d) {
  out <- data.frame(
    patient_id = d$patient_id[[1L]], representative_sample = sort(d$sample)[[1L]],
    primary_sample_count = nrow(d), stringsAsFactors = FALSE
  )
  for (gene in target_genes) out[[gene]] <- mean(d[[gene]])
  out
}))
row.names(patient_expression) <- NULL
write_csv(patient_expression, "TCGA_patient_level_target_expression_log2TPM.csv")

# Collapse GDC clinical and follow-up records to one row per patient.
clinical_path <- file.path(clinical_dir, "clinical.tsv")
follow_path <- file.path(clinical_dir, "follow_up.tsv")
clinical <- read.delim(clinical_path, stringsAsFactors = FALSE, check.names = FALSE,
                       quote = "", comment.char = "", na.strings = na_strings)
follow_up <- read.delim(follow_path, stringsAsFactors = FALSE, check.names = FALSE,
                        quote = "", comment.char = "", na.strings = na_strings)
required_clinical <- c(
  "cases.submitter_id", "demographic.vital_status", "demographic.days_to_death",
  "demographic.age_at_index", "demographic.gender", "diagnoses.age_at_diagnosis",
  "diagnoses.ajcc_pathologic_stage", "diagnoses.diagnosis_is_primary_disease",
  "diagnoses.days_to_last_follow_up"
)
require_columns(clinical, required_clinical, "clinical.tsv")
require_columns(follow_up, c("cases.submitter_id", "follow_ups.days_to_follow_up"), "follow_up.tsv")
clinical$patient_id <- extract_patient_id(clinical$cases.submitter_id)
follow_up$patient_id <- extract_patient_id(follow_up$cases.submitter_id)
clinical <- clinical[nchar(clinical$patient_id) == 12L, , drop = FALSE]
follow_up <- follow_up[nchar(follow_up$patient_id) == 12L, , drop = FALSE]

clinical_case <- do.call(rbind, lapply(split(clinical, clinical$patient_id), function(d) {
  patient <- d$patient_id[[1L]]
  is_primary <- tolower(trimws(as.character(d$diagnoses.diagnosis_is_primary_disease))) %in% c("true", "yes", "1")
  stage_rows <- d[is_primary & !is.na(d$diagnoses.ajcc_pathologic_stage), , drop = FALSE]
  if (nrow(stage_rows) == 0L) stage_rows <- d[!is.na(d$diagnoses.ajcc_pathologic_stage), , drop = FALSE]
  raw_stage <- single_consistent(stage_rows$diagnoses.ajcc_pathologic_stage, "pathological stage", patient)
  primary_age_days <- max_or_na(d$diagnoses.age_at_diagnosis[is_primary])
  age_index <- max_or_na(d$demographic.age_at_index)
  data.frame(
    patient_id = patient,
    vital_status = single_consistent(d$demographic.vital_status, "vital status", patient),
    death_days = max_or_na(d$demographic.days_to_death),
    age = if (is.finite(age_index)) age_index else primary_age_days / 365.25,
    gender = single_consistent(d$demographic.gender, "gender", patient),
    raw_pathological_stage = raw_stage,
    stage = collapse_stage(raw_stage),
    diagnosis_last_followup_days = max_or_na(d$diagnoses.days_to_last_follow_up),
    stringsAsFactors = FALSE
  )
}))
row.names(clinical_case) <- NULL

follow_case <- do.call(rbind, lapply(split(follow_up, follow_up$patient_id), function(d) {
  data.frame(patient_id = d$patient_id[[1L]],
             followup_file_days = max_or_na(d$follow_ups.days_to_follow_up))
}))
row.names(follow_case) <- NULL
clinical_case <- merge(clinical_case, follow_case, by = "patient_id", all.x = TRUE)
clinical_case$last_followup_days <- row_max_or_na(
  clinical_case$diagnosis_last_followup_days, clinical_case$followup_file_days
)
clinical_case$OS <- ifelse(tolower(clinical_case$vital_status) == "dead", 1,
                           ifelse(tolower(clinical_case$vital_status) == "alive", 0, NA_real_))
clinical_case$OS_time <- ifelse(clinical_case$OS == 1,
                                clinical_case$death_days, clinical_case$last_followup_days)
write_csv(clinical_case, "TCGA_GDC_clinical_by_patient.csv")

linked <- merge(patient_expression, clinical_case, by = "patient_id", all.x = TRUE)
stage_cohort <- linked[linked$stage %in% c("II", "III"), , drop = FALSE]
complete_variables <- c("OS", "OS_time", "age", "stage", target_genes)
analysis_data <- stage_cohort[complete.cases(stage_cohort[, complete_variables, drop = FALSE]), , drop = FALSE]
analysis_data <- analysis_data[analysis_data$OS_time >= min_followup, , drop = FALSE]
analysis_data$stage <- factor(analysis_data$stage, levels = c("II", "III"))
analysis_data$age10 <- analysis_data$age / 10
if (nrow(analysis_data) < 20L || sum(analysis_data$OS == 1) < 5L) {
  stop("Insufficient complete cases or events.")
}

# Gene-wise standardization is calculated only after the final cohort is fixed.
z_parameters <- do.call(rbind, lapply(target_genes, function(gene) {
  data.frame(gene = gene, mean_log2TPM = mean(analysis_data[[gene]]),
             sd_log2TPM = sd(analysis_data[[gene]]), n = nrow(analysis_data))
}))
write_csv(z_parameters, "TCGA_gene_Z_standardization_parameters.csv")
for (gene in target_genes) analysis_data[[paste0("Z_", gene)]] <- z_score(analysis_data[[gene]])

analysis_data$P <- rowMeans(analysis_data[, paste0("Z_", pressure_genes), drop = FALSE])
analysis_data$B_NQO1 <- rowMeans(analysis_data[, paste0("Z_", buffer_genes_nqo1), drop = FALSE])
analysis_data$B_NFE2L2 <- rowMeans(analysis_data[, paste0("Z_", buffer_genes_nfe2l2), drop = FALSE])
analysis_data$I_NQO1 <- analysis_data$P - analysis_data$B_NQO1
analysis_data$I_NFE2L2 <- analysis_data$P - analysis_data$B_NFE2L2
for (score in c("P", "B_NQO1", "B_NFE2L2", "I_NQO1", "I_NFE2L2")) {
  analysis_data[[paste0(score, "_perSD")]] <- z_score(analysis_data[[score]])
}
write_csv(analysis_data, "TCGA_PBI_survival_analysis_dataset.csv")

flow <- data.frame(
  step = c("STAR-TPM sample columns", "Primary tumor samples", "Unique primary-tumor patients",
           "Stage II-III patients", "Complete cases with follow-up >= minimum"),
  unit = c("samples", "samples", "patients", "patients", "patients"),
  n = c(nrow(expression), nrow(primary_expression), nrow(patient_expression),
        nrow(stage_cohort), nrow(analysis_data)),
  stringsAsFactors = FALSE
)
flow$excluded_from_previous <- c(NA_integer_, pmax(0, head(flow$n, -1L) - tail(flow$n, -1L)))
write_csv(flow, "TCGA_patient_flow.csv")

score_summary <- do.call(rbind, lapply(c("P", "B_NQO1", "B_NFE2L2", "I_NQO1", "I_NFE2L2"), function(score) {
  x <- analysis_data[[score]]
  data.frame(score = score, n = length(x), mean = mean(x), sd = sd(x), median = median(x),
             q1 = unname(quantile(x, 0.25)), q3 = unname(quantile(x, 0.75)), min = min(x), max = max(x))
}))
write_csv(score_summary, "TCGA_PBI_descriptive_summary.csv")

score_correlations <- do.call(rbind, lapply(c("B_NQO1", "B_NFE2L2"), function(buffer) {
  pearson <- cor.test(analysis_data$P, analysis_data[[buffer]], method = "pearson")
  spearman <- cor.test(analysis_data$P, analysis_data[[buffer]], method = "spearman", exact = FALSE)
  data.frame(buffer_definition = buffer, pearson_r = unname(pearson$estimate),
             pearson_ci_low = pearson$conf.int[[1L]], pearson_ci_high = pearson$conf.int[[2L]],
             pearson_p = pearson$p.value, spearman_rho = unname(spearman$estimate),
             spearman_p = spearman$p.value)
}))
write_csv(score_correlations, "TCGA_P_B_correlations.csv")

tidy_cox <- function(fit, model, analysis_role) {
  s <- summary(fit)
  data.frame(
    model = model, analysis_role = analysis_role, variable = row.names(s$coefficients),
    logHR = s$coefficients[, "coef"], SE = s$coefficients[, "se(coef)"],
    HR = s$conf.int[, "exp(coef)"], CI95_low = s$conf.int[, "lower .95"],
    CI95_high = s$conf.int[, "upper .95"], z = s$coefficients[, "z"],
    p_value = s$coefficients[, "Pr(>|z|)"], n = fit$n, events = fit$nevent,
    concordance = unname(s$concordance[[1L]]), stringsAsFactors = FALSE, row.names = NULL
  )
}

models <- list(
  primary_joint_NQO1 = coxph(Surv(OS_time, OS) ~ P_perSD + B_NQO1_perSD + age10 + stage,
                              data = analysis_data, ties = "efron", x = TRUE, y = TRUE),
  secondary_index_NQO1 = coxph(Surv(OS_time, OS) ~ I_NQO1_perSD + age10 + stage,
                                data = analysis_data, ties = "efron", x = TRUE, y = TRUE),
  sensitivity_joint_NFE2L2 = coxph(Surv(OS_time, OS) ~ P_perSD + B_NFE2L2_perSD + age10 + stage,
                                    data = analysis_data, ties = "efron", x = TRUE, y = TRUE),
  sensitivity_index_NFE2L2 = coxph(Surv(OS_time, OS) ~ I_NFE2L2_perSD + age10 + stage,
                                    data = analysis_data, ties = "efron", x = TRUE, y = TRUE)
)
model_roles <- c("Primary", "Secondary exploratory", "Sensitivity", "Sensitivity")
cox_results <- do.call(rbind, Map(function(fit, name, role) tidy_cox(fit, name, role),
                                  models, names(models), model_roles))
write_csv(cox_results, "TCGA_PBI_Cox_models.csv")

ph_results <- do.call(rbind, lapply(names(models), function(name) {
  test <- cox.zph(models[[name]], transform = "km")
  data.frame(model = name, variable = row.names(test$table), chisq = test$table[, "chisq"],
             df = test$table[, "df"], p_value = test$table[, "p"], row.names = NULL)
}))
write_csv(ph_results, "TCGA_PBI_proportional_hazards_tests.csv")

vif_for_model <- function(fit, model_name) {
  x <- as.data.frame(fit$x)
  do.call(rbind, lapply(names(x), function(variable) {
    others <- setdiff(names(x), variable)
    r2 <- if (length(others) == 0L) 0 else summary(lm(reformulate(others, variable), data = x))$r.squared
    data.frame(model = model_name, variable = variable, VIF = 1 / (1 - r2))
  }))
}
vif_results <- do.call(rbind, Map(vif_for_model, models, names(models)))
write_csv(vif_results, "TCGA_PBI_collinearity_VIF.csv")

reverse_km <- survfit(Surv(OS_time, 1 - OS) ~ 1, data = analysis_data, conf.type = "log-log")
reverse_table <- summary(reverse_km)$table
followup <- data.frame(
  n = nrow(analysis_data), deaths = sum(analysis_data$OS == 1), censored = sum(analysis_data$OS == 0),
  median_followup_days = unname(reverse_table["median"]),
  ci95_low_days = unname(reverse_table["0.95LCL"]),
  ci95_high_days = unname(reverse_table["0.95UCL"]), method = "Reverse Kaplan-Meier"
)
write_csv(followup, "TCGA_followup_summary.csv")

# Descriptive median-split KM analysis for the primary index.
index_cutoff <- median(analysis_data$I_NQO1)
analysis_data$index_group <- factor(ifelse(analysis_data$I_NQO1 > index_cutoff, "High", "Low"),
                                    levels = c("Low", "High"))
km_fit <- survfit(Surv(OS_time, OS) ~ index_group, data = analysis_data, conf.type = "log-log")
logrank <- survdiff(Surv(OS_time, OS) ~ index_group, data = analysis_data)
logrank_p <- pchisq(logrank$chisq, df = 1L, lower.tail = FALSE)
km_summary <- do.call(rbind, lapply(levels(analysis_data$index_group), function(g) {
  d <- analysis_data[analysis_data$index_group == g, , drop = FALSE]
  data.frame(group = g, n = nrow(d), deaths = sum(d$OS == 1), censored = sum(d$OS == 0),
             median_cutoff = index_cutoff, logrank_chisq = unname(logrank$chisq), logrank_p = logrank_p)
}))
write_csv(km_summary, "TCGA_I_NQO1_descriptive_KM_summary.csv")
write_csv(analysis_data, "TCGA_PBI_survival_analysis_dataset.csv")

draw_km <- function() {
  old <- par(no.readonly = TRUE)
  on.exit({layout(1); par(old)})
  layout(matrix(c(1, 2), ncol = 1), heights = c(4.2, 1.15))
  par(mar = c(1.2, 4.8, 1.0, 1.0), family = "sans")
  cols <- c("#1874A6", "#C44E52")
  line_types <- c(1, 2)
  max_time <- max(analysis_data$OS_time)
  tick_years <- seq(0, floor(max_time / 365.25), by = 2)
  tick_days <- tick_years * 365.25
  plot(km_fit, col = cols, lty = line_types, lwd = 2.2, mark.time = TRUE, xscale = 365.25,
       xlab = "", ylab = "Overall survival probability", bty = "l", xaxt = "n",
       conf.int = FALSE, cex.axis = 0.9, cex.lab = 1.0)
  axis(1, at = tick_days, labels = FALSE)
  legend("bottomleft", legend = c(
    sprintf("Low index (n = %d)", sum(analysis_data$index_group == "Low")),
    sprintf("High index (n = %d)", sum(analysis_data$index_group == "High"))
  ), col = cols, lty = line_types, lwd = 2.2, bty = "n", cex = 0.85)
  legend("topright", legend = sprintf("Log-rank p = %.3f", logrank_p), bty = "n", cex = 0.9)

  risk_low <- summary(survfit(Surv(OS_time, OS) ~ 1,
                              data = analysis_data[analysis_data$index_group == "Low", ]),
                      times = tick_days, extend = TRUE)$n.risk
  risk_high <- summary(survfit(Surv(OS_time, OS) ~ 1,
                               data = analysis_data[analysis_data$index_group == "High", ]),
                       times = tick_days, extend = TRUE)$n.risk
  par(mar = c(3.7, 4.8, 0.2, 1.0), family = "sans")
  plot.new()
  plot.window(xlim = c(0, max_time), ylim = c(0.4, 2.8), xaxs = "i")
  axis(1, at = tick_days, labels = tick_years, cex.axis = 0.85)
  mtext("Time (years)", side = 1, line = 2.4, cex = 0.95)
  mtext("Number at risk", side = 2, line = 3.2, cex = 0.82)
  text(tick_days, rep(2, length(tick_days)), labels = risk_low, cex = 0.78, xpd = NA)
  text(tick_days, rep(1, length(tick_days)), labels = risk_high, cex = 0.78, xpd = NA)
  label_x <- -0.35 * 365.25
  text(label_x, 2, labels = "Low", pos = 2, xpd = NA, col = cols[[1L]], cex = 0.78)
  text(label_x, 1, labels = "High", pos = 2, xpd = NA, col = cols[[2L]], cex = 0.78)
}

exposure_rows <- c("P_perSD", "B_NQO1_perSD", "I_NQO1_perSD", "B_NFE2L2_perSD", "I_NFE2L2_perSD")
forest <- cox_results[cox_results$variable %in% exposure_rows, , drop = FALSE]
forest$label <- ifelse(
  forest$model == "primary_joint_NQO1" & forest$variable == "P_perSD", "P: primary joint model",
  ifelse(forest$model == "primary_joint_NQO1" & forest$variable == "B_NQO1_perSD", "B (NQO1): primary joint model",
  ifelse(forest$model == "secondary_index_NQO1", "P-B (NQO1): secondary model",
  ifelse(forest$model == "sensitivity_joint_NFE2L2" & forest$variable == "P_perSD", "P: NFE2L2 sensitivity model",
  ifelse(forest$model == "sensitivity_joint_NFE2L2", "B (NFE2L2): sensitivity model",
         "P-B (NFE2L2): sensitivity model")))))
forest$point_color <- c("#D4743C", "#1874A6", "#333333", "#D4743C", "#5B8C5A", "#6A4C93")
forest$point_shape <- c(15, 15, 18, 0, 0, 1)

draw_forest <- function() {
  old <- par(no.readonly = TRUE); on.exit(par(old))
  par(mar = c(4.5, 11.5, 2.0, 1.0), family = "sans")
  y <- rev(seq_len(nrow(forest)))
  xlim <- range(c(forest$CI95_low, forest$CI95_high, 1), finite = TRUE)
  plot(forest$HR, y, log = "x", xlim = xlim, ylim = c(0.5, nrow(forest) + 0.5),
       xlab = "Hazard ratio per 1-SD increase (95% CI)", ylab = "", yaxt = "n",
       pch = forest$point_shape, col = forest$point_color,
       bty = "n", cex = 1.0)
  abline(v = 1, lty = 2, col = "#777777")
  segments(forest$CI95_low, y, forest$CI95_high, y, lwd = 1.5)
  axis(2, at = y, labels = forest$label, las = 1, tick = FALSE, cex.axis = 0.82)
}

draw_scatter <- function() {
  old <- par(no.readonly = TRUE); on.exit(par(old))
  par(mar = c(4.6, 4.8, 2.0, 1.0), family = "sans")
  cols <- c(II = "#1874A6", III = "#C44E52")
  plot(analysis_data$P, analysis_data$B_NQO1,
       col = adjustcolor(cols[as.character(analysis_data$stage)], alpha.f = 0.65),
       pch = 16, xlab = "Pressure score (mean gene Z)",
       ylab = "Buffering score (mean gene Z)", bty = "l")
  abline(lm(B_NQO1 ~ P, data = analysis_data), col = "#333333", lwd = 1.5)
  legend("topleft", legend = c("Stage II", "Stage III"), col = cols, pch = 16, bty = "n")
  legend("bottomright", legend = sprintf("Pearson r = %.3f", score_correlations$pearson_r[1]), bty = "n")
}

# Original-manuscript visual style, updated for the locked P/B/I analysis.
original_palette <- c(High = "#F04B3A", Low = "#43B8D3")
stage_palette <- c(II = "#F04B3A", III = "#43B8D3")
original_x_max <- 4500
original_x_ticks <- seq(0, original_x_max, by = 900)
original_left_margin <- 5.4
original_right_margin <- 0.9

fmt_p <- function(p) {
  if (length(p) == 0L || is.na(p[[1L]])) return("NA")
  if (p[[1L]] < 0.001) return("<0.001")
  sprintf("%.3f", p[[1L]])
}

step_coordinates <- function(time, value, start = 1) {
  if (length(time) == 0L) return(list(x = 0, y = start))
  x <- c(0, as.vector(rbind(time, time)))
  y_before <- c(start, value[-length(value)])
  y <- c(start, as.vector(rbind(y_before, value)))
  list(x = x, y = y)
}

index_group_fit <- function(group) {
  survfit(
    Surv(OS_time, OS) ~ 1,
    data = analysis_data[analysis_data$index_group == group, , drop = FALSE],
    conf.type = "log-log"
  )
}

draw_original_km_curve <- function(panel_label = NULL) {
  par(mar = c(3.1, original_left_margin, 3.8, original_right_margin),
      family = "sans", lend = "round")
  plot(NA, xlim = c(0, original_x_max), ylim = c(0, 105), xaxs = "i",
       axes = FALSE, xlab = "", ylab = "")
  abline(v = original_x_ticks, h = seq(0, 100, 25), col = "#E8E8E8", lwd = 0.8)
  axis(1, at = original_x_ticks, labels = original_x_ticks, cex.axis = 0.82)
  axis(2, at = seq(0, 100, 25), las = 1, cex.axis = 0.82)
  box(col = "#555555", lwd = 0.8)
  mtext("Overall survival probability (%)", side = 2, line = 3.8, cex = 1.0)

  fits <- setNames(lapply(c("High", "Low"), index_group_fit), c("High", "Low"))
  for (group in c("High", "Low")) {
    fit <- fits[[group]]
    lower <- step_coordinates(fit$time, 100 * fit$lower, 100)
    upper <- step_coordinates(fit$time, 100 * fit$upper, 100)
    polygon(c(lower$x, rev(upper$x)), c(lower$y, rev(upper$y)),
            col = adjustcolor(original_palette[[group]], alpha.f = 0.25), border = NA)
  }
  for (group in c("High", "Low")) {
    fit <- fits[[group]]
    curve <- step_coordinates(fit$time, 100 * fit$surv, 100)
    lines(curve$x, curve$y, col = original_palette[[group]], lwd = 1.8)
    censor_times <- sort(unique(analysis_data$OS_time[
      analysis_data$index_group == group & analysis_data$OS == 0 &
        analysis_data$OS_time <= original_x_max
    ]))
    if (length(censor_times) > 0L) {
      censor_surv <- summary(fit, times = censor_times, extend = TRUE)$surv * 100
      points(censor_times, censor_surv, pch = 3, cex = 0.70, lwd = 1.0,
             col = original_palette[[group]])
    }
  }

  text(320, 9, paste0("log-rank p = ", fmt_p(logrank_p)), adj = 0, cex = 0.95)
  text(1650, 113.5, "P-B index group", font = 2, cex = 0.82, xpd = NA)
  legend(x = 2550, y = 118, xjust = 0, yjust = 1, xpd = NA, horiz = TRUE, bty = "n",
         legend = c("High", "Low"), col = original_palette[c("High", "Low")],
         lwd = 1.8, pch = 3, pt.cex = 0.8, cex = 0.82)
  if (!is.null(panel_label)) {
    mtext(panel_label, side = 3, line = 1.6, adj = -0.03, font = 2, cex = 1.7)
  }
}

draw_original_risk <- function() {
  par(mar = c(3.6, original_left_margin, 2.0, original_right_margin), family = "sans")
  plot(NA, xlim = c(0, original_x_max), ylim = c(0.45, 2.55), xaxs = "i",
       axes = FALSE, xlab = "", ylab = "")
  rect(0, 0.55, original_x_max, 2.45, border = "#666666", lwd = 0.8)
  abline(v = original_x_ticks, col = "#E8E8E8", lwd = 0.8)
  segments(0, 1.5, original_x_max, 1.5, col = "#E8E8E8", lwd = 0.8)
  axis(1, at = original_x_ticks, labels = original_x_ticks, cex.axis = 0.80)
  mtext("Time (days)", side = 1, line = 2.5, cex = 0.90)
  text(0, 2.66, "Number at risk", adj = 0, cex = 1.0, xpd = NA)
  y_positions <- c(High = 1.98, Low = 1.02)
  for (group in names(y_positions)) {
    n_risk <- vapply(original_x_ticks, function(t) {
      sum(analysis_data$OS_time[analysis_data$index_group == group] >= t)
    }, integer(1L))
    text(-0.025 * original_x_max, y_positions[[group]], group, adj = 1,
         col = original_palette[[group]], cex = 0.84, xpd = NA)
    text(original_x_ticks, rep(y_positions[[group]], length(original_x_ticks)),
         n_risk, cex = 0.90, xpd = NA)
  }
}

draw_original_censor <- function() {
  par(mar = c(4.0, original_left_margin, 2.0, original_right_margin), family = "sans")
  censor_tables <- lapply(c("High", "Low"), function(group) {
    times <- analysis_data$OS_time[
      analysis_data$index_group == group & analysis_data$OS == 0 &
        analysis_data$OS_time <= original_x_max
    ]
    if (length(times) == 0L) return(data.frame(time = numeric(), count = integer()))
    tab <- as.data.frame(table(times), stringsAsFactors = FALSE)
    data.frame(time = as.numeric(as.character(tab$times)), count = tab$Freq)
  })
  names(censor_tables) <- c("High", "Low")
  max_count <- max(1L, vapply(censor_tables, function(x) {
    if (nrow(x)) max(x$count) else 0L
  }, integer(1L)))
  plot(NA, xlim = c(0, original_x_max), ylim = c(0, max_count + 0.25), xaxs = "i",
       axes = FALSE, xlab = "", ylab = "")
  abline(v = original_x_ticks, h = seq(0, max_count, 1), col = "#E8E8E8", lwd = 0.8)
  axis(1, at = original_x_ticks, labels = original_x_ticks, cex.axis = 0.78)
  axis(2, at = seq(0, max_count, 1), las = 1, cex.axis = 0.75)
  box(col = "#666666", lwd = 0.8)
  mtext("Time (days)", side = 1, line = 2.5, cex = 0.88)
  mtext("Number censored", side = 2, line = 3.7, cex = 0.92)
  text(0, max_count + 0.42, "Number of censoring", adj = 0, cex = 1.0, xpd = NA)
  for (group in names(censor_tables)) {
    z <- censor_tables[[group]]
    if (nrow(z) > 0L) {
      offset <- if (group == "High") -5 else 5
      segments(z$time + offset, 0, z$time + offset, z$count,
               col = original_palette[[group]], lwd = 1.0)
    }
  }
}

prepare_cox_panel <- function(model_name) {
  z <- cox_results[cox_results$model == model_name, , drop = FALSE]
  expected <- if (model_name == "primary_joint_NQO1") {
    c("P_perSD", "B_NQO1_perSD", "age10", "stageIII")
  } else {
    c("I_NQO1_perSD", "age10", "stageIII")
  }
  z <- z[match(expected, z$variable), , drop = FALSE]
  if (any(is.na(z$variable))) stop("Cox figure source is missing an expected model term.")
  z$display_label <- if (model_name == "primary_joint_NQO1") {
    c("Pressure score\n(per 1 SD)", "Buffering score\n(per 1 SD)",
      "Age\n(per 10 years)", "Stage\nIII vs II")
  } else {
    c("P-B index\n(per 1 SD)", "Age\n(per 10 years)", "Stage\nIII vs II")
  }
  z
}

primary_forest <- prepare_cox_panel("primary_joint_NQO1")
secondary_forest <- prepare_cox_panel("secondary_index_NQO1")

draw_original_cox <- function(forest_data = primary_forest,
                              title = "Primary multivariable Cox regression",
                              panel_label = NULL) {
  par(mar = c(2.7, 0.8, 2.2, 0.7), family = "sans", xpd = NA)
  n_rows <- nrow(forest_data)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0.2, n_rows + 2.45))
  row_y <- rev(seq_len(n_rows))
  row_low <- row_y - 0.43
  row_high <- row_y + 0.43

  for (i in seq_len(n_rows)) {
    if (i %% 2L == 1L) rect(0.01, row_low[[i]], 0.99, row_high[[i]],
                            col = "#F0F0F0", border = NA)
  }
  header_low <- n_rows + 1.55
  header_high <- n_rows + 2.30
  title_low <- n_rows + 0.72
  rect(0.01, header_low, 0.99, header_high, border = "#555555", lwd = 0.8)
  rect(0.01, title_low, 0.99, header_low, border = "#555555", lwd = 0.8)
  rect(0.01, 0.22, 0.99, title_low, border = "#555555", lwd = 0.8)
  segments(c(0.22, 0.29, 0.66), header_low, c(0.22, 0.29, 0.66), header_high,
           col = "#555555", lwd = 0.8)
  segments(c(0.22, 0.29, 0.66), 0.22, c(0.22, 0.29, 0.66), title_low,
           col = "#BBBBBB", lwd = 0.7)

  text(0.04, mean(c(header_low, header_high)), "Variable", adj = 0, cex = 0.88, font = 2)
  text(0.255, mean(c(header_low, header_high)), "N", cex = 0.88, font = 2)
  text(0.50, mean(c(header_low, header_high)), "Hazard ratio", cex = 0.88, font = 2)
  text(0.965, mean(c(header_low, header_high)), "p", adj = 1, cex = 0.88, font = 2)
  text(0.50, mean(c(title_low, header_low)), title, cex = 0.86, font = 2)

  hr_min <- 0.5
  hr_max <- 5.0
  forest_left <- 0.305
  forest_right <- 0.625
  map_hr <- function(x) {
    forest_left + (log(x) - log(hr_min)) / (log(hr_max) - log(hr_min)) *
      (forest_right - forest_left)
  }
  segments(map_hr(1), 0.36, map_hr(1), title_low - 0.12,
           lty = 2, col = "#555555", lwd = 0.8)
  for (i in seq_len(n_rows)) {
    text(0.04, row_y[[i]], forest_data$display_label[[i]], adj = 0, cex = 0.79, font = 2)
    text(0.255, row_y[[i]], as.character(forest_data$n[[i]]), cex = 0.80)
    segments(map_hr(forest_data$CI95_low[[i]]), row_y[[i]],
             map_hr(forest_data$CI95_high[[i]]), row_y[[i]], lwd = 1.15)
    segments(map_hr(forest_data$CI95_low[[i]]), row_y[[i]] - 0.07,
             map_hr(forest_data$CI95_low[[i]]), row_y[[i]] + 0.07, lwd = 1.0)
    segments(map_hr(forest_data$CI95_high[[i]]), row_y[[i]] - 0.07,
             map_hr(forest_data$CI95_high[[i]]), row_y[[i]] + 0.07, lwd = 1.0)
    points(map_hr(forest_data$HR[[i]]), row_y[[i]], pch = 15, cex = 0.78)
    text(0.675, row_y[[i]],
         sprintf("%.2f (%.2f, %.2f)", forest_data$HR[[i]],
                 forest_data$CI95_low[[i]], forest_data$CI95_high[[i]]),
         adj = 0, cex = 0.76)
    text(0.965, row_y[[i]], fmt_p(forest_data$p_value[[i]]), adj = 1, cex = 0.76)
  }
  tick_values <- c(0.5, 1, 2, 5)
  axis_y <- 0.22
  segments(forest_left, axis_y, forest_right, axis_y, col = "#555555", lwd = 0.7)
  for (tick in tick_values) {
    x <- map_hr(tick)
    segments(x, axis_y, x, axis_y - 0.08, col = "#555555", lwd = 0.7)
    text(x, axis_y - 0.14, tick, cex = 0.67)
  }
  if (!is.null(panel_label)) {
    mtext(panel_label, side = 3, line = 0.2, adj = -0.04, font = 2, cex = 1.7)
  }
}

stage_test_B <- wilcox.test(B_NQO1 ~ stage, data = analysis_data, exact = FALSE)

stage_score_specs <- data.frame(
  score = c("P", "B_NQO1", "I_NQO1"),
  label = c("Pressure (P)", "Buffering (B)", "P-B index"),
  stringsAsFactors = FALSE
)

stage_score_long <- do.call(rbind, lapply(seq_len(nrow(stage_score_specs)), function(i) {
  score <- stage_score_specs$score[[i]]
  data.frame(
    patient_id = analysis_data$patient_id,
    stage = as.character(analysis_data$stage),
    score = score,
    score_label = stage_score_specs$label[[i]],
    value = analysis_data[[score]],
    stringsAsFactors = FALSE
  )
}))
row.names(stage_score_long) <- NULL

stage_score_tests <- do.call(rbind, lapply(seq_len(nrow(stage_score_specs)), function(i) {
  score <- stage_score_specs$score[[i]]
  stage_ii <- analysis_data[[score]][analysis_data$stage == "II"]
  stage_iii <- analysis_data[[score]][analysis_data$stage == "III"]
  test <- suppressWarnings(wilcox.test(
    x = stage_iii, y = stage_ii, alternative = "two.sided",
    exact = FALSE, conf.int = TRUE, conf.level = 0.95
  ))
  data.frame(
    score = score,
    score_label = stage_score_specs$label[[i]],
    contrast = "Stage III minus Stage II",
    n_stage_ii = length(stage_ii),
    n_stage_iii = length(stage_iii),
    median_stage_ii = median(stage_ii),
    q1_stage_ii = unname(quantile(stage_ii, 0.25)),
    q3_stage_ii = unname(quantile(stage_ii, 0.75)),
    median_stage_iii = median(stage_iii),
    q1_stage_iii = unname(quantile(stage_iii, 0.25)),
    q3_stage_iii = unname(quantile(stage_iii, 0.75)),
    hodges_lehmann_difference = unname(test$estimate),
    CI95_low = unname(test$conf.int[[1L]]),
    CI95_high = unname(test$conf.int[[2L]]),
    wilcoxon_W = unname(test$statistic),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}))
stage_score_tests$p_BH <- p.adjust(stage_score_tests$p_value, method = "BH")
write_csv(stage_score_long, "TCGA_Figure10_panelC_P_B_PminusB_source_data.csv")
write_csv(stage_score_tests, "TCGA_Figure10_panelC_stage_comparisons.csv")

draw_original_stage <- function(panel_label = NULL) {
  par(mar = c(4.3, 5.3, 3.4, 0.8), family = "sans", xpd = NA)
  score <- analysis_data$B_NQO1
  score_min <- floor(min(score) * 2) / 2
  score_max <- ceiling(max(score) * 2) / 2
  score_range <- score_max - score_min
  y_top <- score_max + max(0.55, score_range * 0.14)
  plot(NA, xlim = c(0.4, 2.6), ylim = c(score_min, y_top),
       axes = FALSE, xlab = "", ylab = "")
  y_ticks <- pretty(c(score_min, score_max), n = 5)
  y_ticks <- y_ticks[y_ticks >= score_min & y_ticks <= score_max]
  abline(h = y_ticks, col = "#E7EEF1", lwd = 0.8)
  axis(1, at = c(1, 2), labels = c("Stage II", "Stage III"),
       cex.axis = 0.85, lwd = 0, lwd.ticks = 0)
  axis(2, at = y_ticks, las = 1, cex.axis = 0.88)
  mtext("Buffering score (mean gene Z)", side = 2, line = 3.8, cex = 1.0)

  boxplot(B_NQO1 ~ stage, data = analysis_data, at = c(1, 2), add = TRUE,
          axes = FALSE, outline = FALSE, width = c(0.52, 0.52),
          border = unname(stage_palette[c("II", "III")]),
          col = adjustcolor(unname(stage_palette[c("II", "III")]), alpha.f = 0.04),
          lwd = 1.25, pars = list(medlwd = 1.8, whisklty = 1, staplewex = 0.45))
  for (i in seq_along(levels(analysis_data$stage))) {
    group <- levels(analysis_data$stage)[[i]]
    values <- analysis_data$B_NQO1[analysis_data$stage == group]
    set.seed(42 + i)
    points(i + runif(length(values), -0.16, 0.16), values,
           pch = 16, cex = 0.55, col = stage_palette[[group]])
    segments(i - 0.26, median(values), i + 0.26, median(values),
             col = stage_palette[[group]], lwd = 2.0)
  }

  bracket_y <- score_max + max(0.23, score_range * 0.055)
  bracket_drop <- max(0.08, score_range * 0.018)
  p_label <- paste0("p = ", fmt_p(stage_test_B$p.value))
  p_gap_half <- max(0.14, strwidth(p_label, units = "user", cex = 0.88) / 2 + 0.04)
  segments(1, bracket_y, 1.5 - p_gap_half, bracket_y, col = "#555555", lwd = 0.9)
  segments(1.5 + p_gap_half, bracket_y, 2, bracket_y, col = "#555555", lwd = 0.9)
  segments(c(1, 2), bracket_y, c(1, 2), bracket_y - bracket_drop,
           col = "#555555", lwd = 0.8)
  text(1.5, bracket_y, p_label, cex = 0.88)

  usr <- par("usr")
  strip_bottom <- usr[[4L]] + score_range * 0.035
  strip_top <- usr[[4L]] + score_range * 0.13
  rect(usr[[1L]], strip_bottom, usr[[2L]], strip_top, col = "#DCEBF1", border = NA)
  text(1.5, mean(c(strip_bottom, strip_top)), "Stage II-III TCGA-COAD",
       col = "#244A68", cex = 0.78)
  if (!is.null(panel_label)) {
    mtext(panel_label, side = 3, line = 0.2, adj = -0.04, font = 2, cex = 1.7)
  }
}

draw_stage_all_scores <- function(panel_label = NULL) {
  par(mar = c(5.0, 5.3, 4.3, 0.8), family = "sans", xpd = NA)
  group_positions <- matrix(c(1, 2, 4, 5, 7, 8), nrow = 3L, byrow = TRUE,
                            dimnames = list(stage_score_specs$score, c("II", "III")))
  value_min <- min(stage_score_long$value)
  value_max <- max(stage_score_long$value)
  value_range <- value_max - value_min
  y_min <- floor(value_min * 2) / 2
  bracket_y <- value_max + max(0.38, value_range * 0.075)
  y_top <- bracket_y + max(0.72, value_range * 0.14)
  plot(NA, xlim = c(0.25, 8.75), ylim = c(y_min, y_top),
       axes = FALSE, xlab = "", ylab = "")
  y_ticks <- pretty(c(y_min, value_max), n = 6)
  y_ticks <- y_ticks[y_ticks >= y_min & y_ticks <= value_max]
  abline(h = y_ticks, col = "#E7EEF1", lwd = 0.8)
  abline(h = 0, col = "#B8C2C8", lwd = 0.8)
  axis(1, at = rowMeans(group_positions), labels = stage_score_specs$label,
       cex.axis = 0.78, lwd = 0, lwd.ticks = 0)
  axis(2, at = y_ticks, las = 1, cex.axis = 0.82)
  mtext("Composite score (gene-wise Z units)", side = 2, line = 3.8, cex = 0.94)

  for (i in seq_len(nrow(stage_score_specs))) {
    score <- stage_score_specs$score[[i]]
    for (stage_name in c("II", "III")) {
      position <- group_positions[score, stage_name]
      values <- analysis_data[[score]][analysis_data$stage == stage_name]
      boxplot(values, at = position, add = TRUE, axes = FALSE, outline = FALSE,
              width = 0.72, border = stage_palette[[stage_name]],
              col = adjustcolor(stage_palette[[stage_name]], alpha.f = 0.04),
              lwd = 1.1, pars = list(medlwd = 1.7, whisklty = 1, staplewex = 0.45))
      set.seed(420 + i * 10 + match(stage_name, c("II", "III")))
      points(position + runif(length(values), -0.24, 0.24), values,
             pch = 16, cex = 0.36, col = adjustcolor(stage_palette[[stage_name]], alpha.f = 0.82))
      segments(position - 0.36, median(values), position + 0.36, median(values),
               col = stage_palette[[stage_name]], lwd = 1.7)
    }

    x_left <- group_positions[i, "II"]
    x_right <- group_positions[i, "III"]
    label <- paste0("p = ", fmt_p(stage_score_tests$p_value[[i]]),
                    "\nBH p = ", fmt_p(stage_score_tests$p_BH[[i]]))
    gap_half <- max(0.24, strwidth(label, units = "user", cex = 0.58) / 2 + 0.06)
    x_mid <- mean(c(x_left, x_right))
    segments(x_left, bracket_y, x_mid - gap_half, bracket_y,
             col = "#555555", lwd = 0.8)
    segments(x_mid + gap_half, bracket_y, x_right, bracket_y,
             col = "#555555", lwd = 0.8)
    segments(c(x_left, x_right), bracket_y,
             c(x_left, x_right), bracket_y - max(0.10, value_range * 0.018),
             col = "#555555", lwd = 0.8)
    text(x_mid, bracket_y, label, cex = 0.58)
  }

  usr <- par("usr")
  strip_bottom <- usr[[4L]] + value_range * 0.025
  strip_top <- usr[[4L]] + value_range * 0.10
  rect(usr[[1L]], strip_bottom, usr[[2L]], strip_top, col = "#DCEBF1", border = NA)
  text(mean(usr[1:2]), mean(c(strip_bottom, strip_top)), "Stage II-III TCGA-COAD",
       col = "#244A68", cex = 0.72)
  legend("topright", inset = c(0.01, -0.105), xpd = NA, horiz = TRUE, bty = "n",
         legend = c(
           sprintf("Stage II (n=%d)", sum(analysis_data$stage == "II")),
           sprintf("Stage III (n=%d)", sum(analysis_data$stage == "III"))
         ), col = stage_palette[c("II", "III")],
         pch = 16, pt.cex = 0.85, cex = 0.67)
  if (!is.null(panel_label)) {
    mtext(panel_label, side = 3, line = 0.2, adj = -0.04, font = 2, cex = 1.7)
  }
}

draw_original_km <- function() {
  layout(matrix(c(1, 2, 3), ncol = 1L), heights = c(3.4, 1.35, 1.35))
  draw_original_km_curve()
  draw_original_risk()
  draw_original_censor()
}

draw_original_combined <- function() {
  layout_matrix <- rbind(
    c(rep(1, 7), rep(4, 5)),
    c(rep(1, 7), rep(4, 5)),
    c(rep(1, 7), rep(4, 5)),
    c(rep(1, 7), rep(4, 5)),
    c(rep(2, 7), rep(5, 5)),
    c(rep(2, 7), rep(5, 5)),
    c(rep(3, 7), rep(5, 5)),
    c(rep(3, 7), rep(5, 5)),
    c(rep(3, 7), rep(5, 5))
  )
  layout(layout_matrix, heights = c(rep(0.86, 4), rep(0.65, 2), rep(0.75, 3)))
  draw_original_km_curve("A")
  draw_original_risk()
  draw_original_censor()
  draw_original_cox(primary_forest, "Primary multivariable Cox regression", "B")
  draw_stage_all_scores("C")
}

panel_a_source <- analysis_data[, c("patient_id", "OS_time", "OS", "index_group", "I_NQO1")]
names(panel_a_source)[names(panel_a_source) == "I_NQO1"] <- "P_minus_B_index"
panel_b_source <- primary_forest[, c("variable", "display_label", "HR", "CI95_low",
                                     "CI95_high", "p_value", "n", "events")]
write_csv(panel_a_source, "TCGA_Figure10_panelA_KM_source_data.csv")
write_csv(panel_b_source, "TCGA_Figure10_panelB_primary_Cox_source_data.csv")
write_csv(secondary_forest, "TCGA_secondary_index_Cox_forest_source_data.csv")

save_figure_formats("TCGA_I_NQO1_descriptive_KM", 7.0, 5.8, draw_km)
save_figure_formats("TCGA_PBI_Cox_forest", 7.0, 4.8, draw_forest)
save_figure_formats("TCGA_P_B_correlation", 5.8, 5.2, draw_scatter)
save_figure_formats("Figure10A_P_minus_B_KM_original_style", 8.8, 7.5, draw_original_km)
save_figure_formats("Figure10B_primary_Cox_original_style", 7.5, 5.3,
                    function() draw_original_cox(primary_forest))
save_figure_formats("Figure10C_buffering_by_stage_original_style", 7.2, 6.2,
                    function() draw_original_stage())
save_figure_formats("Figure10C_P_B_PminusB_by_stage_original_style", 9.0, 6.2,
                    function() draw_stage_all_scores())
save_figure_formats("Figure10_ABC_original_style", 15.0, 9.0, draw_original_combined)
save_figure_formats("TCGA_secondary_P_minus_B_Cox_complete", 7.5, 4.8,
                    function() draw_original_cox(
                      secondary_forest, "Secondary multivariable Cox regression"
                    ))

input_files <- c(expression_path, probemap_path, clinical_path, follow_path)
write_csv(data.frame(file = basename(input_files), md5 = unname(tools::md5sum(input_files))),
          "TCGA_input_file_md5.csv")
capture.output(sessionInfo(), file = file.path(out_dir, "TCGA_R_sessionInfo.txt"))

methods_text <- c(
  "TCGA-COAD pressure-buffering analysis (locked 2026-08-12)",
  "Patient-level primary-tumor expression was derived from the UCSC Xena GDC STAR-TPM matrix (already log2[TPM+1]).",
  "Within the final stage II-III complete-case cohort, each target gene was Z-standardized before unweighted domain means were calculated.",
  "P comprised ACSL4, LPCAT3, TFRC, and NCOA4. Primary B comprised SLC7A11, GPX4, AIFM2, GCH1, DHODH, and NQO1.",
  "The primary Cox model included P, B, age per 10 years, and pathological stage; the secondary model used P-B with the same covariates.",
  "Sensitivity analyses replaced NQO1 with NFE2L2. Composite HRs are reported per 1-SD increase. Median-split KM analysis is descriptive."
)
writeLines(methods_text, file.path(out_dir, "TCGA_manuscript_methods.txt"), useBytes = TRUE)

manifest <- c(
  "# TCGA-COAD analysis outputs", "", "Generated: 2026-08-12", "",
  "## Core tables",
  "- `TCGA_PBI_Cox_models.csv`: primary, secondary, and sensitivity Cox estimates with 95% CIs.",
  "- `TCGA_PBI_proportional_hazards_tests.csv`: Schoenfeld tests.",
  "- `TCGA_PBI_collinearity_VIF.csv`: predictor collinearity audit.",
  "- `TCGA_target_gene_coverage_audit.csv`: prespecified mapping gate.",
  "- `TCGA_gene_Z_standardization_parameters.csv`: cohort-specific means and SDs.",
  "- `TCGA_PBI_survival_analysis_dataset.csv`: patient-level derived analysis data.", "",
  "## Manuscript-facing figures",
  "- `Figure10_ABC_original_style.*`: original-style A/B/C composite in PDF, PNG, and 600-dpi TIFF.",
  "- `Figure10A_P_minus_B_KM_original_style.*`: descriptive median-split P-B Kaplan-Meier panel.",
  "- `Figure10B_primary_Cox_original_style.*`: primary Cox table/forest including P, B, age, and stage.",
  "- `Figure10C_P_B_PminusB_by_stage_original_style.*`: P, B, and P-B distributions by pathological stage.",
  "- `TCGA_Figure10_panelC_stage_comparisons.csv`: Stage III-minus-II Hodges-Lehmann differences, 95% CIs, raw p values, and BH-adjusted p values.",
  "- `TCGA_secondary_P_minus_B_Cox_complete.*`: secondary P-B Cox model including age and stage.",
  "- Panel-specific source-data CSV files accompany all manuscript-facing figures.", "",
  "## Earlier analytical figures",
  "- The compact Kaplan-Meier, exposure-only Cox overview, and P-B correlation figures are retained for audit history.", "",
  "## Interpretation", "- The P-B measure is an exploratory transcriptional index, not a direct ferroptosis assay or validated clinical score."
)
writeLines(manifest, file.path(out_dir, "_analysis_outputs.md"), useBytes = TRUE)

cat("TCGA-COAD P/B/I analysis complete.\n")
cat("Patients:", nrow(analysis_data), " Events:", sum(analysis_data$OS == 1), "\n")
print(cox_results[cox_results$variable %in% exposure_rows,
                  c("model", "variable", "HR", "CI95_low", "CI95_high", "p_value")], row.names = FALSE)
cat("\nStage III versus Stage II score comparisons:\n")
print(stage_score_tests[, c("score_label", "hodges_lehmann_difference", "CI95_low",
                            "CI95_high", "p_value", "p_BH")], row.names = FALSE)
cat("Output directory:", out_dir, "\n")
