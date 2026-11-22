###############################################################################
# Descriptive PCA reconstructed from the archived GSE42387 processed matrix.
# Added 2026-11-21; the historical PCA implementation was not archived.
# This equivalent implementation reproduces the reported PC1/PC2 percentages.
# No differential-expression model, composite score, or sample is changed.
###############################################################################
set.seed(42)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript GSE42387_PCA_reconstructed.R <series_matrix.gz> <outdir>")
}
input <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
out <- args[[2L]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
lines <- readLines(gzfile(input), warn = FALSE)
start <- grep("^!series_matrix_table_begin", lines) + 1L
end <- grep("^!series_matrix_table_end", lines) - 1L
if (length(start) != 1L || length(end) != 1L) stop("Invalid series table boundaries.")
tab <- read.delim(text = paste(lines[start:end], collapse = "\n"),
                  header = TRUE, quote = "\"", check.names = FALSE)
expr <- as.matrix(tab[, -1L, drop = FALSE])
storage.mode(expr) <- "double"
row.names(expr) <- tab[[1L]]
if (ncol(expr) != 27L || any(!is.finite(expr))) stop("Expected 27 finite arrays.")
variance <- apply(expr, 1L, var)
selected <- order(variance, decreasing = TRUE)[seq_len(500L)]
fit <- prcomp(t(expr[selected, , drop = FALSE]), center = TRUE, scale. = FALSE)
pct <- 100 * fit$sdev^2 / sum(fit$sdev^2)
scores <- data.frame(sample = row.names(fit$x), fit$x[, 1:2, drop = FALSE])
write.csv(scores, file.path(out, "GSE42387_PCA_scores_reconstructed.csv"), row.names = FALSE)
write.csv(data.frame(component = paste0("PC", seq_along(pct)), variance_percent = pct),
          file.path(out, "GSE42387_PCA_variance_reconstructed.csv"), row.names = FALSE)
write.csv(data.frame(rank = seq_along(selected), feature_id = row.names(expr)[selected],
                     sample_variance = variance[selected]),
          file.path(out, "GSE42387_PCA_top500_features_reconstructed.csv"), row.names = FALSE)
writeLines(c(
  "Reconstructed equivalent PCA; the historical PCA code was not archived.",
  "Input: archived processed log2-expression matrix, all 27 arrays and 32,705 features.",
  "Select 500 features with greatest sample variance across arrays (R var, denominator n-1).",
  "Transpose to samples by features, center each feature, do not variance-scale.",
  "PCA uses R prcomp singular value decomposition; signs of components are arbitrary.",
  sprintf("PC1 = %.6f%%; PC2 = %.6f%%.", pct[[1L]], pct[[2L]]),
  "This describes global sample structure; it does not remove cell-line or batch effects.",
  "No sample was excluded and no inferential result was changed."
), file.path(out, "GSE42387_PCA_reconstruction_methods.txt"))
capture.output(sessionInfo(), file = file.path(out, "GSE42387_PCA_sessionInfo.txt"))
cat(sprintf("PC1 %.6f%%; PC2 %.6f%%\n", pct[[1L]], pct[[2L]]))
