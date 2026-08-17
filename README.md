# Repository-name-GWAS-MetaAnalysis-QC
Quality control, harmonization and METAL meta-analysis of GWAS summary statistics.
# GWAS Meta-Analysis QC and Workflow

This repository contains the quality-control, allele harmonization, and meta-analysis workflow used for a three-study GWAS meta-analysis involving the DGI, MAGIC, and Sardinia datasets.

The workflow includes:

- Initial assessment of GWAS summary statistics
- Identification and removal of ambiguous palindromic SNPs
- Allele harmonization across datasets
- Conversion of DGI allele coding
- Identification of SNPs shared across studies
- Cross-study allele consistency checks
- Selection of the final common SNP set
- METAL-based GWAS meta-analysis

A total of 2,124 SNPs were retained for the final three-study meta-analysis.

The repository contains the analysis scripts, input datasets, quality-control files, workflow documentation, and meta-analysis results to facilitate reproducibility of the analysis.
## Analysis Workflow

DGI GWAS
   ↓
MAGIC GWAS
   ↓
Sardinia GWAS
   ↓
SNP matching
   ↓
Allele harmonization
   ↓
Common SNP filtering
   ↓
2,124 common SNPs
   ↓
METAL meta-analysis
   ↓
METAANALYSIS_1.tbl
