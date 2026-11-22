# Reproduce the aggregate bibliometric figure panels without database access.
target_r_version <- "4.5.3"
if (as.character(getRversion()) != target_r_version) {
  warning("This package was validated with R ", target_r_version,
          "; current runtime is R ", getRversion(), ".", call. = FALSE)
}
if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) {
  invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
}
script_file <- if (sys.nframe() > 0) sys.frame(1)$ofile else NULL
if (is.null(script_file)) {
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(script_arg)) script_file <- sub("^--file=", "", script_arg[1])
}
if (is.null(script_file)) stop("Source this saved run_all.R file or run Rscript run_all.R.")
package_root <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))

required_packages <- c("ggplot2", "patchwork", "readxl", "ragg")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "),
       ". See README.md for installation instructions.")
}
for (folder in c("Figures_PNG", "Figures_TIFF", "Source_Data", "Session_Info")) {
  dir.create(file.path(package_root, folder), recursive = TRUE, showWarnings = FALSE)
}

scripts <- c(
  "01_Fig3B_Annual_Publication_Trend.R",
  "02_Fig3CD_Publication_Lifecycle.R",
  "03_Fig5D_Journal_Publication_Trends.R",
  "04_Fig6C_Author_Production_Over_Time.R"
)
for (script in scripts) {
  cat("\nRunning ", script, "\n", sep = "")
  source(file.path(package_root, "Code", script),
         local = new.env(parent = environment()), encoding = "UTF-8")
}

capture.output(sessionInfo(),
               file = file.path(package_root, "Session_Info", "R_sessionInfo.txt"))
write.csv(
  data.frame(
    software = c("R", required_packages),
    version = c(as.character(getRversion()),
                vapply(required_packages, function(x) as.character(packageVersion(x)), character(1)))
  ),
  file.path(package_root, "Session_Info", "Software_Versions.csv"), row.names = FALSE
)
cat("\nCompleted: five figure panels in PNG and TIFF, with source-data tables.\n")
