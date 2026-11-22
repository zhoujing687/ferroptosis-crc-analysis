# Figure manifest

Generated: 2026-08-12
Target: Frontiers in Pharmacology

## Input data

- `../03_Input_Data/GSE42387_series_matrix.txt.gz`: processed GEO series matrix used to read the 27 expression arrays.
- `../03_Input_Data/GPL16297.txt`: GPL16297 platform annotation used for probe-to-gene mapping.
- Input provenance and SHA-256/MD5 checksums are documented in `../03_Input_Data/README.txt`.

## Gene-level expression distributions

- Primary panel: 10 genes (ACSL4, LPCAT3, TFRC, NCOA4, SLC7A11, GPX4, AIFM2, GCH1, DHODH, and NQO1).
- Oxaliplatin resistance: one 5 x 2 grid comparing parental and OxR samples within HCT116, HT29, and LoVo cells.
- SN-38 resistance: one 5 x 2 grid comparing parental and SN38R samples within HCT116, HT29, and LoVo cells.
- Sensitivity panel: NFE2L2 shown separately for both drug-resistance states.
- Replicates: all three independent passages per cell line and condition are displayed.
- Statistics: cell-line-specific limma contrasts with BH adjustment across the 11 tested genes within each contrast.
- Source data: expression values and matching limma results are supplied as CSV files.
- Formats: vector PDF, 300-dpi PNG, and 600-dpi LZW-compressed TIFF.

## Figure 10D/E standalone redraw

- `01_Code/GSE42387_Figure10_DE_single_panels.R` redraws the manuscript-style standalone D/E panels from the locked gene-contrast and P/B/I contrast tables.
- This is a plotting-only script: no limma model or composite-score model is refit.
- Outputs: PDF, 300-dpi PNG, and 600-dpi LZW-compressed TIFF.

## Model diagnostics

- Primary readable versions: separate P, B, and P-B 2 x 2 diagnostic figures.
- Overview version: three score rows, with a dedicated key column for each row.
- Data-region labels: numeric labels 1-3 only.
- Identifier key: complete GSM identifiers and sample titles outside the data regions.
- Selection rule: three observations with the largest Cook's distances in each model.
- Reference line in the influence panel: Cook's distance = 4/n.
- Formats: vector PDF, 300-dpi PNG, and 600-dpi LZW-compressed TIFF.
- Statistical outputs: unchanged; figure revision affected presentation only.

critic_pass: yes
refine_rounds: 3
