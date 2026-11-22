# README
## TCGA-COAD Reproducibility Package

### 1. Overview

This package provides the R script and supporting source-data and audit materials for reproducing the exploratory TCGA-COAD survival analysis of the Pᵀ, Bᵀ, and Pᵀ−Bᵀ transcriptional scores.

The current manuscript labels the combined TCGA-COAD/GSE42387 core results as Figure 10, with panels A–C showing the TCGA-COAD analyses. This archive standardizes historical Figure 11 script filenames and internal references to Figure 10 to match the current manuscript. This is a file-naming and cross-reference change only; the analytical methods, data, and numerical results are unchanged.

### 2. Primary Score Definitions

- **Pᵀ (ferroptotic pressure-related transcriptional score):** ACSL4, LPCAT3, TFRC, and NCOA4.
- **Bᵀ (antioxidant-buffering-related transcriptional score):** SLC7A11, GPX4, AIFM2, GCH1, DHODH, and NQO1.
- **Sensitivity analysis:** NFE2L2 replaces NQO1; NFE2L2 is not added to the primary Bᵀ score.

### 3. Software Requirements

The analysis was performed using:

- R version 4.5.2
- survival version 3.8-3

The recorded environment is documented in `03_Source_Data/TCGA_R_sessionInfo.txt`.

### 4. Input Data

#### Included Clinical Inputs

The original clinical files from the 22 February 2026 [GDC TCGA-COAD export](https://portal.gdc.cancer.gov/projects/TCGA-COAD) are included in the package:

- `02_Input_Data/clinical.tsv`
- `02_Input_Data/follow_up.tsv`

Both files match the MD5 checksums recorded for the original analysis in `03_Source_Data/TCGA_input_file_md5.csv`. Use the included files for reproduction, as a current GDC export may differ from this historical snapshot.

#### Required External Inputs

Only the expression matrix and gene annotation file need to be obtained separately:

1. [TCGA-COAD.star_tpm.tsv.gz](https://gdc.xenahubs.net/download/TCGA-COAD.star_tpm.tsv.gz): UCSC Xena GDC hub dataset version 05-09-2024, GENCODE v36 annotation, supplied as log2(TPM + 1) values.
2. [gencode.v36.annotation.gtf.gene.probemap](https://gdc.xenahubs.net/download/gencode.v36.annotation.gtf.gene.probemap): GENCODE v36 gene mapping file.

### 5. Preparing the Input Directory

After extracting the archive, set the working directory to the package root (the directory containing `01_Code`). Key files already included in the package are:

```text
01_Code/
|-- TCGA_COAD_PBI_survival_analysis.R
`-- Figure10_TCGA_GSE_core_results_combined.R
02_Input_Data/
|-- clinical.tsv
`-- follow_up.tsv
03_Source_Data/
`-- TCGA_input_file_md5.csv
```

Place the downloaded expression matrix and gene mapping file directly in `02_Input_Data`. Verify their MD5 checksums against `03_Source_Data/TCGA_input_file_md5.csv` before running the analysis. Keep the included clinical files unchanged.

### 6. Reproducing the TCGA-COAD Analysis

Run the following command from the package root to generate fresh analysis outputs in `04_Recomputed_Outputs`:

```bash
Rscript "01_Code/TCGA_COAD_PBI_survival_analysis.R" --expression "02_Input_Data/TCGA-COAD.star_tpm.tsv.gz" --probemap "02_Input_Data/gencode.v36.annotation.gtf.gene.probemap" --clinical-dir "02_Input_Data" --outdir "04_Recomputed_Outputs" --data-access-date "2026-02-22" --expression-version "05-09-2024"
```

### 7. Combined Figure Script

The combined Figure 10 script redraws from both archives. Extract S1 File and
S2 File into sibling `S1_File` and `S2_File` directories, and run from their parent:

```bash
Rscript "S1_File/01_Code/Figure10_TCGA_GSE_core_results_combined.R" "." "Figure10_Redrawn"
```

It reads `S1_File/03_Source_Data` and `S2_File/03_Source_Data`. The historical
`02_TCGA_COAD`/`03_GSE42387` source-table layout is also accepted. The generated
panel-D legend now correctly separates the drug symbols from domain colors;
no estimate or plotted data was changed.

### 8. Clinical baseline table

The added descriptive script uses the archived final analysis dataset and the
archived reverse-Kaplan-Meier follow-up summary, without changing patient
selection or fitting new association models:

```bash
Rscript "01_Code/TCGA_baseline_characteristics.R" "03_Source_Data" "04_Descriptive_Tables"
```

`S6_Table.csv` records the denominator, missingness and summary definition.

### 9. Reproducibility Note

The package includes the historical clinical inputs and records checksums for all four required input files. Exact rerunning of the analysis also requires obtaining the STAR-TPM expression matrix and GENCODE v36 gene mapping file described above. Recomputed outputs are written to a separate directory and do not replace the archived source data.
