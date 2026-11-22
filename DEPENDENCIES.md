# Recorded software environments

These versions are transcribed from the existing packages. They describe
archived runs, not newly installed dependencies or a tested lockfile.

| Component | Recorded R | Additional packages |
| --- | --- | --- |
| TCGA-COAD analysis and combined Figure 10 | 4.5.2 | survival 3.8-3 |
| GSE42387 fixed-effects/limma analysis | 4.5.2 | limma 3.66.0; statmod 1.5.2 |
| S3 aggregate bibliometric reconstruction | 4.5.3 | ggplot2 4.0.3; patchwork 1.3.2; readxl 1.5.0; ragg 1.5.2 |

The baseline-table, multiplicity-inventory, retrospective PCA, and standalone
Figure 10D/E scripts use packages supplied with R. See their archived session
records for the runtime used for each derived artifact. Base/recommended R
packages and graphics devices also form part of the environment.

Authoritative environment records:

- `S1_File/03_Source_Data/TCGA_R_sessionInfo.txt`
- `S1_File/04_Descriptive_Tables/S6_Table_R_sessionInfo.txt`
- `S2_File/03_Source_Data/GSE42387_R_sessionInfo.txt`
- `S2_File/03_Source_Data/GSE42387_package_versions.csv`
- `S2_File/04_Reconstructed_PCA/GSE42387_PCA_sessionInfo.txt`
- `S2_File/05_Multiple_Testing_Inventory/S7_Table_R_sessionInfo.txt`
- `S3_File/Session_Info/R_sessionInfo.txt`
- `S3_File/Session_Info/Software_Versions.csv`

If packages are missing, install CRAN packages using `install.packages()` and
Bioconductor packages such as limma using `BiocManager::install()` in an
appropriate R/Bioconductor environment. A current default installation may
select versions different from the archived versions. Exact-version restoration
has not been prepared or tested as part of this folder assembly.

VOSviewer, CiteSpace, bibliometrix/Biblioshiny, Pajek, and other historical
desktop analysis steps are outside these standalone R scripts. Their full
historical project configuration is not recreated by installing the packages
listed above. Do not interpret this dependency table as a complete environment
for every bibliometric analysis in the manuscript.
