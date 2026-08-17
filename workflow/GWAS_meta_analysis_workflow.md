# GWAS Meta-analysis Workflow

## Overview

This workflow describes the quality control, allele harmonization,
SNP matching, and meta-analysis of GWAS summary statistics from
three studies:

- DGI
- MAGIC
- Sardinia

## 1. Input GWAS datasets

The original GWAS summary statistics were obtained for the three
study datasets.

Input datasets:

- DGI
- MAGIC
- Sardinia

## 2. Initial quality control

The datasets were inspected for:

- SNP identifiers
- Chromosomal position
- Effect allele
- Non-effect allele
- Effect estimates
- Standard errors
- P-values
- Effect allele frequencies
- Sample size

## 3. Allele coding assessment

The DGI dataset used numerical allele coding:

1 = A
2 = C
3 = G
4 = T

Allele coding was converted to nucleotide representation to allow
comparison with MAGIC and Sardinia.

## 4. SNP matching

SNPs were matched using the rsID.

The datasets were compared to identify SNPs shared across studies.

## 5. Allele harmonization

Effect and non-effect alleles were compared between datasets.

SNPs were classified according to allele orientation, including:

- Same orientation
- Reversed orientation
- Complementary orientation
- Complementary + reversed orientation
- Incompatible alleles

Where necessary, effect estimates were reversed so that all studies
used the same effect allele.

## 6. Common SNP selection

SNPs present and compatible across all three datasets were retained.

A final set of 2,124 SNPs was used for the three-study meta-analysis.

## 7. Final QC datasets

The harmonized datasets were:

- DGI_FINAL.txt
- MAGIC_FINAL.txt
- SARDINIA_FINAL3.txt

Each dataset contained:

SNP, CHR, BP, EA, NEA, BETA, SE, P, EAF, N

## 8. Meta-analysis

Meta-analysis was performed using METAL.

The STDERR scheme was used with:

- Effect size: BETA
- Standard error: SE
- Sample size: N
- Effect allele: EA
- Non-effect allele: NEA

## 9. Final results

The final meta-analysis included 2,124 SNPs.

The results were written to:

METAANALYSIS_1.tbl

The METAL output contains:

- MarkerName
- Allele1
- Allele2
- Effect
- StdErr
- P-value
- Direction
