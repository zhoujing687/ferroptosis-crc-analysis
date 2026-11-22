# TCGA-COAD analysis outputs

Generated: 2026-08-12

## Core tables
- `TCGA_PBI_Cox_models.csv`: primary, secondary, and sensitivity Cox estimates with 95% CIs.
- `TCGA_PBI_proportional_hazards_tests.csv`: Schoenfeld tests.
- `TCGA_PBI_collinearity_VIF.csv`: predictor collinearity audit.
- `TCGA_target_gene_coverage_audit.csv`: prespecified mapping gate.
- `TCGA_gene_Z_standardization_parameters.csv`: cohort-specific means and SDs.
- `TCGA_PBI_survival_analysis_dataset.csv`: patient-level derived analysis data.

## Manuscript-facing figures
- `Figure10_ABC_original_style.*`: original-style A/B/C composite in PDF, PNG, and 600-dpi TIFF.
- `Figure10A_P_minus_B_KM_original_style.*`: descriptive median-split P-B Kaplan-Meier panel.
- `Figure10B_primary_Cox_original_style.*`: primary Cox table/forest including P, B, age, and stage.
- `Figure10C_P_B_PminusB_by_stage_original_style.*`: P, B, and P-B distributions by pathological stage.
- `TCGA_Figure10_panelC_stage_comparisons.csv`: Stage III-minus-II Hodges-Lehmann differences, 95% CIs, raw p values, and BH-adjusted p values.
- `TCGA_secondary_P_minus_B_Cox_complete.*`: secondary P-B Cox model including age and stage.
- Panel-specific source-data CSV files accompany all manuscript-facing figures.

## Earlier analytical figures
- The compact Kaplan-Meier, exposure-only Cox overview, and P-B correlation figures are retained for audit history.

## Interpretation
- The P-B measure is an exploratory transcriptional index, not a direct ferroptosis assay or validated clinical score.
