# README
## GSE42387 Reproducibility Package

### 1. Overview

This package provides the materials required to reproduce the exploratory transcriptomic analyses of the GSE42387 colorectal cancer cell-line dataset.

All 27 arrays are retained. The dataset comprises HCT116, HT29, and LoVo cells, each represented under parental, oxaliplatin-resistant (OxR), and SN-38-resistant (SN38R) conditions, with three independent passage replicates per cell-line-by-condition group.

The current manuscript labels the combined TCGA-COAD/GSE42387 core results as Figure 10, with panels A–C showing the TCGA-COAD analyses and panels D–E showing the GSE42387 analyses. This archive standardizes historical Figure 11 script filenames and internal references to Figure 10 to match the current manuscript. This is a file-naming and cross-reference change only; the analytical methods, data, and numerical results are unchanged.

### 2. Primary Score Definitions

- **Pᵀ (ferroptotic pressure-related transcriptional score):** ACSL4, LPCAT3, TFRC, and NCOA4.
- **Bᵀ (antioxidant-buffering-related transcriptional score):** SLC7A11, GPX4, AIFM2, GCH1, DHODH, and NQO1.
- **Sensitivity analysis:** NFE2L2 replaces NQO1; NFE2L2 is not added to the primary Bᵀ score.

### 3. Package Contents and Input Data

The archive contains the fixed-effects/limma analysis script, the standalone script for reconstructing Figure 10D and Figure 10E, source-data and audit outputs, and the public GEO input files used in the analysis. Its top-level directories are `01_Code`, `02_Input_Data`, and `03_Source_Data`. The analysis script filename now exactly matches the command below. The ten primary targets and the NFE2L2 replacement sensitivity definition are unchanged.

The twelve archived full-feature result tables are supplied separately in
S4 File (six within-cell-line resistant-versus-parental contrast tables) and
S5 File (two equally weighted average contrasts and four omnibus tables).
Splitting the original 41.6-MB archive preserves every original compressed
table and keeps each submitted ZIP below 20 MB. These results can also be
regenerated from the included inputs using the complete analysis command.

The required input files are:

- `GSE42387_series_matrix.txt.gz` — processed GEO series matrix containing all 27 arrays.
- `GPL16297.txt` — corresponding platform annotation.

Both input files are included in `02_Input_Data`. Their provenance and MD5/SHA-256 values are documented in `02_Input_Data/README.txt`; the MD5 values are also listed in `03_Source_Data/GSE42387_input_file_md5.csv`.

### 4. Software Requirements

The analysis was performed using:

- R version 4.5.2
- limma version 3.66.0
- statmod version 1.5.2

The recorded environment is documented in `03_Source_Data/GSE42387_R_sessionInfo.txt`.

### 5. Reproducing the Complete GSE42387 Analysis

The required GEO input files are included in the archive, so no additional GEO download is required.

After extracting the archive, set the working directory to the package root (the directory containing `01_Code`). Run:

```bash
Rscript "01_Code/GSE42387_PBI_fixed_effects_analysis.R" "02_Input_Data/GSE42387_series_matrix.txt.gz" "02_Input_Data/GPL16297.txt" "04_Recomputed_Outputs"
```

The command reruns the fixed-effects and limma analyses and writes the recomputed analysis outputs to `04_Recomputed_Outputs`.

### 6. Reconstructing Figure 10D and Figure 10E

The package also includes a standalone script for independently reconstructing the manuscript-style Figure 10D and Figure 10E panels from the provided statistical result tables.

From the package root, run:

```bash
Rscript "01_Code/GSE42387_Figure10_DE_single_panels.R" "03_Source_Data" "04_Figure10_DE_single_panels"
```

This script reads the archived result tables `GSE42387_11_genes_limma_main_contrasts.csv` and `GSE42387_PBI_fixed_effects_contrasts.csv` from `03_Source_Data`. It writes the redrawn panels to the separate `04_Figure10_DE_single_panels` directory without refitting the statistical models.

For the combined Figure 10, extract S1 File and S2 File into sibling `S1_File`
and `S2_File` directories and run from their parent:

```bash
Rscript "S2_File/01_Code/Figure10_TCGA_GSE_core_results_combined.R" "." "Figure10_Redrawn"
```

The script reads both `03_Source_Data` subdirectories, or the historical
`02_TCGA_COAD`/`03_GSE42387` source-table layout.

### 7. Descriptive PCA reconstruction and testing inventory

The original PCA script was not archived. The added equivalent implementation
uses the 500 features with greatest sample variance across all 27 arrays,
transposes to samples by features, centers each feature, and does not scale
features. It reproduces PC1 = 44.562486% and PC2 = 40.931249% from the archived
processed matrix. Component signs are arbitrary. This retrospective
reconstruction is descriptive and does not change any sample or inference:

```bash
Rscript "01_Code/GSE42387_PCA_reconstructed.R" "02_Input_Data/GSE42387_series_matrix.txt.gz" "04_Reconstructed_PCA"
```

The distributed `04_Reconstructed_PCA` directory contains the reconstructed
scores, selected feature IDs and variances, component variance percentages,
method note and session information.

After extracting both packages into sibling directories, S7 Table can be
rebuilt from the archived results without refitting models:

```bash
Rscript "S2_File/01_Code/Multiple_testing_inventory.R" "S1_File/03_Source_Data" "S2_File/03_Source_Data" "Testing_Inventory"
```

The inventory documents the correction families as implemented; it does not
establish that they were prospectively registered.

### 8. Platform-Specific Note

On Windows, if R reports a path-encoding error, run the commands from a package path that contains only ASCII characters.
