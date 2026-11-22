# GSE42387 analysis outputs

Generated: 2026-08-12

## Input data
- `../03_Input_Data/GSE42387_series_matrix.txt.gz`: processed GEO series matrix for GSE42387, containing all 27 arrays.
- `../03_Input_Data/GPL16297.txt`: GPL16297 platform annotation used to map probe IDs to gene symbols.
- Source, size, MD5, and SHA-256 details are recorded in `../03_Input_Data/README.txt`; the legacy MD5 table is `GSE42387_input_file_md5.csv`.
- Complete rerun from the package root: `Rscript 01_Code/GSE42387_PBI_fixed_effects_analysis.R 03_Input_Data/GSE42387_series_matrix.txt.gz 03_Input_Data/GPL16297.txt 04_Recomputed_Outputs`.

## Core composite-score tables
- `GSE42387_PBI_fixed_effects_contrasts.csv`: P/B/I estimates and 95% CIs.
- `GSE42387_PBI_HC3_robust_sensitivity_contrasts.csv`: HC3 robust-SE sensitivity estimates.
- `GSE42387_PBI_fixed_effects_omnibus_tests.csv`: interaction and heterogeneity tests.
- `GSE42387_PBI_model_diagnostics.csv`: residual, variance, influence, and fit checks.

## Mapping and expression
- `GSE42387_target_probe_coverage_summary.csv`: mandatory pre-score coverage audit.
- `GSE42387_target_probe_mapping_audit.csv`: target probe details.
- `GSE42387_PBI_analysis_dataset.csv`: gene expression, gene Z scores, and composites.

## Limma outputs
- Eleven-gene main-contrast and omnibus tables are provided as CSV.
- Full-feature contrast and omnibus results are provided as compressed CSV files.

## Figures
- Main P/B/I forest, 11-gene average-effect forest, group distributions, and diagnostics are provided in PDF, PNG, and TIFF.
- Standalone manuscript-style Figure 10D/E redraw is provided by `01_Code/GSE42387_Figure10_DE_single_panels.R`; it reads the locked result tables above and does not refit models.
- `GSE42387_10_primary_genes_OxR_distributions.*`: 10-gene parental-versus-OxR distributions stratified by cell line.
- `GSE42387_10_primary_genes_SN38R_distributions.*`: 10-gene parental-versus-SN38R distributions stratified by cell line.
- `GSE42387_NFE2L2_sensitivity_distributions.*`: separate NFE2L2 sensitivity distributions.
- `GSE42387_11_gene_distribution_source_data.csv` and cell-line limma statistics tables provide figure traceability.

- Model diagnostics are provided as one overview and as separate 2 x 2 figures for P, B, and P-B. Full GSM identifiers are listed outside the plotting regions and linked to numbered points selected by Cook distance.

## Interpretation
- Composite scores use ordinary fixed-effects models; limma moderation is reserved for feature-level inference.
