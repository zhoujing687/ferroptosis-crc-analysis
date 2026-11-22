###############################################################################
# Descriptive baseline table from the archived final 289-patient dataset.
# Added 2026-11-21. No cohort selection, score, or association model is changed.
# Usage: Rscript TCGA_baseline_characteristics.R <03_Source_Data> <output_dir>
###############################################################################
set.seed(42)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply the archived source-data directory and output directory.")
source_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
out <- args[[2L]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
d <- read.csv(file.path(source_dir, "TCGA_PBI_survival_analysis_dataset.csv"), check.names = FALSE)
followup <- read.csv(file.path(source_dir, "TCGA_followup_summary.csv"), check.names = FALSE)
if (nrow(d) != 289L || anyDuplicated(d$patient_id) || sum(d$OS == 1) != 55L) {
  stop("Unexpected archived analysis cohort.")
}
summary_rows <- list()
add <- function(characteristic, summary, n_available, n_missing, definition) {
  summary_rows[[length(summary_rows) + 1L]] <<- data.frame(
    Characteristic = characteristic, Summary = summary,
    N_available = n_available, N_missing = n_missing, Definition = definition)
}
n <- nrow(d)
add("Patients", as.character(n), n, 0L, "Final archived stage II-III analysis cohort; one row per patient")
age <- d$age[!is.na(d$age)]
q <- quantile(age, c(0.25, 0.75), names = FALSE, type = 7)
add("Age, years", sprintf("%.1f (%.1f-%.1f)", median(age), q[[1L]], q[[2L]]),
    length(age), sum(is.na(d$age)), "Median (interquartile range)")
for (variable in c("gender", "stage", "OS")) {
  available <- !is.na(d[[variable]]) & nzchar(as.character(d[[variable]]))
  levels <- switch(variable, gender = c("female", "male"), stage = c("II", "III"), OS = c("1", "0"))
  for (level in levels) {
    count <- sum(as.character(d[[variable]]) == level, na.rm = TRUE)
    label <- switch(variable, gender = paste("Recorded gender:", level),
                    stage = paste("Pathological stage", level),
                    OS = if (level == "1") "Deaths" else "Censored")
    add(label, sprintf("%d (%.1f%%)", count, 100 * count / n),
        sum(available), sum(!available), "Count (percentage of all 289 patients)")
  }
}
add("Potential follow-up, days", sprintf("%.0f (%.0f-%.0f)",
    followup$median_followup_days, followup$ci95_low_days, followup$ci95_high_days),
    n, sum(is.na(d$OS_time)), "Reverse Kaplan-Meier median (95% confidence interval), copied from archived summary")
add("OS event and time complete", sprintf("%d (%.1f%%)", sum(complete.cases(d[, c("OS", "OS_time")])),
    100 * mean(complete.cases(d[, c("OS", "OS_time")]))),
    sum(complete.cases(d[, c("OS", "OS_time")])), sum(!complete.cases(d[, c("OS", "OS_time")])),
    "Complete in the final selected cohort; not the original linked TCGA population")
table <- do.call(rbind, summary_rows)
write.csv(table, file.path(out, "S6_Table.csv"), row.names = FALSE, na = "")
capture.output(sessionInfo(), file = file.path(out, "S6_Table_R_sessionInfo.txt"))
print(table[, c("Characteristic", "Summary", "N_missing")], row.names = FALSE)
