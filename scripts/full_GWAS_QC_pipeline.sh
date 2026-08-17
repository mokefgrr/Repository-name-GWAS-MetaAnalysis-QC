
#!/usr/bin/env bash

# ============================================================
# GWAS META-ANALYSIS QC PIPELINE
# DGI + MAGIC + SARDINIA
# ============================================================

# Working directory:
# ~/Documents/METAL/generic-metal/examples/GlucoseExample

mkdir -p QC


# ============================================================
# STEP 1: IDENTIFY PALINDROMIC SNPs IN MAGIC
# ============================================================

echo "STEP 1: MAGIC palindromic SNPs"

awk 'NR>1 &&
(($5=="A" && $6=="T") ||
 ($5=="T" && $6=="A") ||
 ($5=="C" && $6=="G") ||
 ($5=="G" && $6=="C")) {
print $3, $5, $6, $7
}' MAGIC_FUSION_Results.txt | head


awk 'NR>1 &&
(($5=="A" && $6=="T") ||
 ($5=="T" && $6=="A") ||
 ($5=="C" && $6=="G") ||
 ($5=="G" && $6=="C"))' \
MAGIC_FUSION_Results.txt | wc -l


# ============================================================
# STEP 2: IDENTIFY AMBIGUOUS PALINDROMIC MAGIC SNPs
# EAF between 0.42 and 0.58
# ============================================================

awk 'NR>1 &&
((($5=="A" && $6=="T") ||
  ($5=="T" && $6=="A") ||
  ($5=="C" && $6=="G") ||
  ($5=="G" && $6=="C")) &&
 $7>0.42 && $7<0.58) {
print $3, $5, $6, $7
}' MAGIC_FUSION_Results.txt | head


awk 'NR>1 &&
((($5=="A" && $6=="T") ||
  ($5=="T" && $6=="A") ||
  ($5=="C" && $6=="G") ||
  ($5=="G" && $6=="C")) &&
 $7>0.42 && $7<0.58)' \
MAGIC_FUSION_Results.txt | wc -l


# ============================================================
# STEP 3: IDENTIFY AMBIGUOUS PALINDROMIC SARDINIA SNPs
# ============================================================

awk 'NR>1 &&
((($5=="A" && $6=="T") ||
  ($5=="T" && $6=="A") ||
  ($5=="C" && $6=="G") ||
  ($5=="G" && $6=="C")) &&
 $7>0.42 && $7<0.58) {
print $1, $5, $6, $7
}' magic_SARDINIA.tbl | head


awk 'NR>1 &&
((($5=="A" && $6=="T") ||
  ($5=="T" && $6=="A") ||
  ($5=="C" && $6=="G") ||
  ($5=="G" && $6=="C")) &&
 $7>0.42 && $7<0.58)' \
magic_SARDINIA.tbl | wc -l


# ============================================================
# STEP 4: REMOVE AMBIGUOUS PALINDROMIC SNPs FROM MAGIC
# ============================================================

awk 'NR==1 ||
!((($5=="A" && $6=="T") ||
   ($5=="T" && $6=="A") ||
   ($5=="C" && $6=="G") ||
   ($5=="G" && $6=="C")) &&
  $7>0.42 && $7<0.58)' \
MAGIC_FUSION_Results.txt \
> QC/MAGIC_clean.txt


# ============================================================
# STEP 5: REMOVE AMBIGUOUS PALINDROMIC SNPs FROM SARDINIA
# ============================================================

awk 'NR==1 ||
!((($5=="A" && $6=="T") ||
   ($5=="T" && $6=="A") ||
   ($5=="C" && $6=="G") ||
   ($5=="G" && $6=="C")) &&
  $7>0.42 && $7<0.58)' \
magic_SARDINIA.tbl \
> QC/SARDINIA_clean.txt


# ============================================================
# STEP 6: CHECK CLEANED FILE SIZES
# ============================================================

wc -l QC/MAGIC_clean.txt
wc -l QC/SARDINIA_clean.txt


# ============================================================
# STEP 7: COMPARE MAGIC AND SARDINIA ALLELES
# ============================================================

awk '
NR==FNR {
    if (NR>1) {
        ea[$3]=$5
        nea[$3]=$6
    }
    next
}
NR>1 && ($1 in ea) {

    a=ea[$1]
    b=nea[$1]
    c=$5
    d=$6

    if (a==c && b==d) {
        same++
    }
    else if (a==d && b==c) {
        reversed++
    }
    else {
        mismatch++
    }
}
END {
    print "Same orientation:", same+0
    print "Reversed orientation:", reversed+0
    print "Other mismatch:", mismatch+0
}' \
QC/MAGIC_clean.txt \
QC/SARDINIA_clean.txt


# ============================================================
# STEP 8: SHOW REVERSED ALLELES
# ============================================================

awk '
NR==FNR {
    if (NR>1) {
        ea[$3]=$5
        nea[$3]=$6
        beta[$3]=$9
    }
    next
}
NR>1 && ($1 in ea) {

    a=ea[$1]
    b=nea[$1]
    c=$5
    d=$6

    if (a==d && b==c) {
        print $1, a, b, beta[$1], c, d, $9
    }
}' \
QC/MAGIC_clean.txt \
QC/SARDINIA_clean.txt | head -20


# ============================================================
# STEP 9: CHECK DGI ALLELE CODING
# ============================================================

awk 'NR>1 {print $5, $6}' \
DGI_three_regions.txt | sort -u


awk 'NR>1 {print $5}' \
DGI_three_regions.txt | sort | uniq -c


awk 'NR>1 {print $6}' \
DGI_three_regions.txt | sort | uniq -c


# ============================================================
# STEP 10: CONVERT DGI ALLELE CODING
#
# Original DGI coding:
# 1,2,3,4
#
# Converted to nucleotide alleles using the established coding
# ============================================================

# NOTE:
# This section assumes QC/DGI_converted.txt has already been
# generated using the allele coding established during QC.


head -5 QC/DGI_converted.txt


# ============================================================
# STEP 11: CHECK DGI HARMONIZED FILE
# ============================================================

wc -l QC/DGI_harmonized.txt


head -5 QC/DGI_harmonized.txt


# Check number of fields
awk 'NR<=5 {
    print "NF=" NF, $0
}' QC/DGI_harmonized.txt


# Check malformed rows
awk 'NF!=10 {
    print "BAD:", NR, "NF=" NF, $0
}' QC/DGI_harmonized.txt | head


# ============================================================
# STEP 12: ADD HEADER TO DGI HARMONIZED FILE
# ============================================================

head -5 QC/DGI_harmonized_with_header.txt

wc -l QC/DGI_harmonized_with_header.txt


# ============================================================
# STEP 13: CHECK SNP COUNTS
# ============================================================

awk 'NR>1 {
    print $1
}' QC/DGI_harmonized_with_header.txt |
sort -u | wc -l


# ============================================================
# STEP 14: IDENTIFY COMMON SNPs
# ============================================================

comm -12 \
<(awk 'NR>1 {print $3}' DGI_three_regions.txt | sort) \
<(awk 'NR>1 {print $3}' QC/MAGIC_clean.txt | sort) |
wc -l


# ============================================================
# STEP 15: CHECK DGI COMPATIBLE SNPs
# ============================================================

awk '
NR>1 {
    print $5, $6
}' DGI_three_regions.txt | sort -u


# ============================================================
# STEP 16: CREATE COMMON 3-STUDY SNP SET
# ============================================================

wc -l QC/COMMON_3STUDIES.txt


# ============================================================
# STEP 17: CHECK FINAL DATASET SIZES
# ============================================================

wc -l \
QC/DGI_FINAL.txt \
QC/MAGIC_FINAL.txt \
QC/SARDINIA_FINAL.txt


# ============================================================
# STEP 18: CHECK FINAL SNP COUNTS
# ============================================================

wc -l \
QC/DGI_FINAL.txt \
QC/MAGIC_FINAL.txt \
QC/SARDINIA_FINAL.txt


# ============================================================
# STEP 19: CHECK DGI VS MAGIC ALLELES
# ============================================================

awk '
NR==FNR {
    if (NR>1) {
        ea[$1]=$4
        nea[$1]=$5
    }
    next
}
NR>1 && ($1 in ea) {

    if ($4==ea[$1] && $5==nea[$1])
        same++
    else
        mismatch++
}
END {
    print "Same alleles:", same+0
    print "Mismatches:", mismatch+0
}' \
QC/DGI_FINAL.txt \
QC/MAGIC_FINAL.txt


# ============================================================
# STEP 20: CHECK DGI VS SARDINIA ALLELES
# ============================================================

awk '
NR==FNR {
    if (NR>1) {
        ea[$1]=$4
        nea[$1]=$5
    }
    next
}
NR>1 && ($1 in ea) {

    if ($4==ea[$1] && $5==nea[$1])
        same++
    else
        mismatch++
}
END {
    print "Same alleles:", same+0
    print "Mismatches:", mismatch+0
}' \
QC/DGI_FINAL.txt \
QC/SARDINIA_FINAL3.txt


# ============================================================
# STEP 21: CHECK MAGIC VS SARDINIA ALLELES
# ============================================================

awk '
NR==FNR {
    if (NR>1) {
        ea[$1]=$4
        nea[$1]=$5
    }
    next
}
NR>1 && ($1 in ea) {

    if ($4==ea[$1] && $5==nea[$1])
        same++
    else
        mismatch++
}
END {
    print "Same alleles:", same+0
    print "Mismatches:", mismatch+0
}' \
QC/MAGIC_FINAL.txt \
QC/SARDINIA_FINAL3.txt


# ============================================================
# STEP 22: CHECK DUPLICATE SNPs
# ============================================================

echo "DGI duplicates:"
awk 'NR>1 {print $1}' QC/DGI_FINAL.txt |
sort | uniq -d


echo "MAGIC duplicates:"
awk 'NR>1 {print $1}' QC/MAGIC_FINAL.txt |
sort | uniq -d


echo "SARDINIA duplicates:"
awk 'NR>1 {print $1}' QC/SARDINIA_FINAL3.txt |
sort | uniq -d


# ============================================================
# STEP 23: CHECK FINAL FILE STRUCTURE
# ============================================================

awk 'NR==2 {
    for(i=1;i<=NF;i++)
        print i "["$i"]"
}' QC/DGI_FINAL.txt


awk 'NR==2 {
    for(i=1;i<=NF;i++)
        print i "["$i"]"
}' QC/MAGIC_FINAL.txt


awk 'NR==2 {
    for(i=1;i<=NF;i++)
        print i "["$i"]"
}' QC/SARDINIA_FINAL3.txt


# ============================================================
# STEP 24: CREATE METAL INPUT FILE
# ============================================================

cat metal_QC.txt


# ============================================================
# STEP 25: RUN METAL
# ============================================================

/c/Users/drmic/Documents/METAL/generic-metal/metal.exe \
metal_QC.txt


# ============================================================
# STEP 26: VERIFY META-ANALYSIS OUTPUT
# ============================================================

ls -lh QC/METAANALYSIS*


# Number of meta-analyzed SNPs
awk 'NR>1 {
    print $1
}' QC/METAANALYSIS_1.tbl | wc -l


# Check empty rows
awk 'NR>1 && NF==0 {
    bad++
}
END {
    print "Empty rows:", bad+0
}' QC/METAANALYSIS_1.tbl


# Check warnings/errors
grep -i -E \
"warning|error|failed|skipped" \
QC/METAANALYSIS_1.tbl.info


# ============================================================
# STEP 27: VIEW METAL OUTPUT DESCRIPTION
# ============================================================

cat QC/METAANALYSIS_1.tbl.info


# ============================================================
# STEP 28: VIEW FIRST META-ANALYSIS RESULTS
# ============================================================

head -5 QC/METAANALYSIS_1.tbl


# ============================================================
# STEP 29: FIND TOP META-ANALYSIS RESULT
# ============================================================

awk '
NR>1 {
    if ($6 < min || NR==2) {
        min=$6
        snp=$1
    }
}
END {
    print "Top SNP:", snp
    print "P-value:", min
}' QC/METAANALYSIS_1.tbl


echo ""
echo "============================================================"
echo "GWAS QC PIPELINE COMPLETED"
echo "============================================================"
