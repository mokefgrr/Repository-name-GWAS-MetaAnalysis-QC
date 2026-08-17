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

## Analysis Workflow

The GWAS meta-analysis was performed using the following workflow:

1. **Input GWAS datasets**
   - DGI
   - MAGIC
   - Sardinia

2. **Initial dataset assessment**
   - Examined SNP identifiers, chromosome, position, effect allele, non-effect allele, beta, standard error, P-value, effect allele frequency, and sample size.
   - Assessed allele coding and dataset structure.

3. **Allele coding conversion**
   - The DGI dataset contained numerical allele coding (1–4).
   - The numerical alleles were converted to nucleotide alleles (A/C/G/T) to allow harmonization with the MAGIC and Sardinia datasets.

4. **SNP matching**
   - SNPs were matched across datasets using rsID.
   - 4,592 SNPs were shared between the relevant datasets.

5. **Allele compatibility assessment**
   - Effect and non-effect alleles were compared across datasets.
   - SNPs were assessed for:
     - Same orientation
     - Reversed orientation
     - Complementary orientation
     - Complementary + reversed orientation
     - Incompatible alleles

6. **Allele harmonization**
   - Effect alleles were aligned across studies.
   - Where the effect allele orientation was reversed, the beta coefficient was reversed so that all studies represented effects for the same allele.

7. **Common SNP selection**
   - SNPs compatible across all three studies were retained.
   - 2,124 SNPs were retained for the final three-study analysis.

8. **Final QC datasets**
   - DGI_FINAL.txt
   - MAGIC_FINAL.txt
   - SARDINIA_FINAL3.txt

   All three final datasets contained 2,124 SNPs.

9. **GWAS meta-analysis**
   - Meta-analysis was performed using METAL.
   - The `SCHEME STDERR` approach was used.
   - Effect estimates were combined using beta coefficients and standard errors.
   - Sample size was supplied as the analysis weight.

10. **Final meta-analysis output**
    - The final analysis included 2,124 SNPs.
    - Results were generated in:
      `METAANALYSIS_1.tbl`
    - The METAL output contains:
      - MarkerName
      - Allele1
      - Allele2
      - Effect
      - StdErr
      - P-value
      - Direction
