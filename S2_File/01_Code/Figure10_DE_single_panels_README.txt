Figure 10D/E standalone redraw

Script: GSE42387_Figure10_DE_single_panels.R

Inputs (from 03_Source_Data):
- GSE42387_11_genes_limma_main_contrasts.csv
- GSE42387_PBI_fixed_effects_contrasts.csv

Run from the package root with base R:

Rscript 01_Code/GSE42387_Figure10_DE_single_panels.R 03_Source_Data 04_Figure10_DE_single_panels

The script redraws the manuscript-style standalone Figure 10D and Figure 10E panels in PDF, 300-dpi PNG, and 600-dpi LZW TIFF formats. It does not rerun limma, refit the composite-score models, or change any estimates.

Windows note: if R reports a file-name conversion error under a non-ASCII extraction path, copy the package to an ASCII path such as C:\RQA_GSE42387 and rerun the command.
