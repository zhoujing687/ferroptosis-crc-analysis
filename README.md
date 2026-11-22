# Ferroptosis and colorectal cancer: analysis code

This repository brings together the existing S1-S3 supporting packages.
It contains the archived scripts, their supporting
inputs, and reference outputs. It is not a reconstruction of the complete
historical bibliometric workflow. No statistical model was rerun during this
folder assembly, and no manuscript figure was replaced.

## Layout

| Directory | Contents |
| --- | --- |
| `S1_File/` | TCGA-COAD survival analysis, baseline table, Figure 10 code, clinical inputs, and archived results |
| `S2_File/` | GSE42387 analysis, public GEO inputs, Figure 10 code, retrospective PCA reconstruction, and multiplicity inventory |
| `S3_File/` | Code, aggregate inputs, and reference outputs for five bibliometric figure panels |
| `Documentation/` | Current submission versions of S6 and the expanded S7 CSV tables |
| `Provenance/` | Source archive hashes, file origins, and assembly verification report |

There are 13 R script files. The combined Figure 10 script occurs in both S1
and S2; these two copies are identical. Original package paths and bytes are
preserved so existing relative references continue to work. The root
`MANIFEST_SHA256.csv` lists the distributed files except itself; the original manifests
inside S1-S3 retain their original scope.

Read [DEPENDENCIES.md](DEPENDENCIES.md) before running code.

## Run commands

Commands below use paths relative to this folder. Run each block in a fresh
terminal starting at this folder. Quote paths containing spaces. On Windows,
if an R build cannot handle the local non-ASCII path, copy the entire folder
to an ASCII-only path before running it.

### TCGA-COAD

The clinical and follow-up TSV files are included. Obtain these two external
inputs separately and put them in `S1_File/02_Input_Data/`:

- [TCGA-COAD.star_tpm.tsv.gz](https://gdc.xenahubs.net/download/TCGA-COAD.star_tpm.tsv.gz)
- [gencode.v36.annotation.gtf.gene.probemap](https://gdc.xenahubs.net/download/gencode.v36.annotation.gtf.gene.probemap)

The archived analysis used Xena dataset version `05-09-2024`, GENCODE v36,
with supplied expression values already on the log2(TPM + 1) scale. Compare
downloads with `S1_File/03_Source_Data/TCGA_input_file_md5.csv` before use;
the live files may differ from the historical snapshot. Keep the included
clinical snapshot unchanged.

```sh
cd S1_File
Rscript "01_Code/TCGA_COAD_PBI_survival_analysis.R" --expression "02_Input_Data/TCGA-COAD.star_tpm.tsv.gz" --probemap "02_Input_Data/gencode.v36.annotation.gtf.gene.probemap" --clinical-dir "02_Input_Data" --outdir "04_Recomputed_Outputs" --data-access-date "2026-02-22" --expression-version "05-09-2024"
```

The metadata arguments above describe the archived inputs. Update them if
you intentionally analyze another snapshot. The script writes to the supplied
output directory. To rebuild the baseline table from the archived final cohort:

```sh
cd S1_File
Rscript "01_Code/TCGA_baseline_characteristics.R" "03_Source_Data" "04_Recomputed_Baseline"
```

### GSE42387

The processed series matrix and platform annotation are included, with
provenance in `S2_File/02_Input_Data/README.txt`. All 27 arrays are retained.
They represent three cell lines and three treatment states, with three passage
replicates per group. The historical description of independent triplicates
does not establish three independently generated resistant models.

```sh
cd S2_File
Rscript "01_Code/GSE42387_PBI_fixed_effects_analysis.R" "02_Input_Data/GSE42387_series_matrix.txt.gz" "02_Input_Data/GPL16297.txt" "04_Recomputed_Outputs"
```

The separate S4/S5 submission archives hold twelve full-feature result tables.
They are not duplicated here; the main GSE42387 script can regenerate them.

Reconstruct descriptive PCA from the included processed matrix:

```sh
cd S2_File
Rscript "01_Code/GSE42387_PCA_reconstructed.R" "02_Input_Data/GSE42387_series_matrix.txt.gz" "04_Recomputed_PCA"
```

This is an explicitly retrospective implementation; the historical PCA script
was not archived. It uses the 500 highest-variance features, centers them,
and does not scale them.

### Figure 10 and multiplicity inventory

The following commands read the archived result tables. They do not require
the external TCGA expression matrix. Run from this folder:

```sh
Rscript "S1_File/01_Code/Figure10_TCGA_GSE_core_results_combined.R" "." "runs/Figure10_Redrawn"
Rscript "S2_File/01_Code/GSE42387_Figure10_DE_single_panels.R" "S2_File/03_Source_Data" "runs/Figure10_DE"
Rscript "S2_File/01_Code/Multiple_testing_inventory.R" "S1_File/03_Source_Data" "S2_File/03_Source_Data" "runs/Testing_Inventory"
```

`Multiple_testing_inventory.R` produces the original eight-column inventory,
as archived in `S2_File/05_Multiple_Testing_Inventory/`. The current submission
table `Documentation/S7_Table.csv` has six additional explanatory columns on
analysis units, sample size, degrees of freedom, reference distributions,
confidence intervals, and family boundaries. The script does not generate
these added columns. Both versions are retained with their distinct origins.

### Five bibliometric panels

S3 covers Fig 3B, Fig 3C, Fig 3D, Fig 5D, and Fig 6C. Its entry point writes
within its own folder and replaces some reference outputs and session files.
Run a working copy to preserve the archive. For example, in PowerShell from
this folder, choose a new destination for each run:

```powershell
New-Item -ItemType Directory -Path "runs" -Force | Out-Null
Copy-Item -LiteralPath "S3_File" -Destination "runs/S3_work" -Recurse -ErrorAction Stop
Rscript --vanilla "runs/S3_work/run_all.R"
```

Use `runs/S3_work` only when it does not already exist. Consult
`S3_File/README.md` for the aggregate data provenance, lifecycle equations,
fixed parameters, and the limitations of the extrapolations.

## Scope and known limitations

- The historical Scopus conversion/Python deduplication code was not archived.
  Reason-specific exclusion logs, identity crosswalks, and some software
  project settings also remain unavailable. This folder does not recreate them.
- The available bibliometric R code operates on aggregate inputs. It does not
  reproduce every search, screening, collaboration, or network-analysis step.
  The current screening description is one initial reviewer followed by a
  second reviewer checking all records, not independent dual screening.
- The original main Figure 5 remains the author's retained figure. Its PLOS ONE
  cumulative count for 2017 is 5, whereas the archived aggregate input and S3
  redraw use 10. Copying S3 does not resolve that discrepancy or replace the
  main figure. The original S5 Table local-citation ranking is also retained.
- Lifecycle projections are conditional descriptive extrapolations, not
  validated forecasts. Transcriptomic analyses are exploratory; their
  inclusion does not establish clinical utility or ferroptosis sensitivity.
- Original package READMEs and revision notes are preserved as historical
  records. This root README clarifies the assembly scope and the current S7
  table; it does not silently edit those records.
- Web of Science and Scopus record-level exports are not included. Data-source
  terms and third-party rights must be considered separately from licensing
  author-written code.

## Verification

`Provenance/ASSEMBLY_REPORT.json` records the checks performed when assembling
this folder. Every extracted file was compared with its ZIP member, the two
current CSV tables were copied byte for byte, and source archive hashes were
checked before and after assembly. These are packaging checks, not a fresh
statistical validation.
