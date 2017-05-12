# Colon and Rectal Cancer Data

## Outline

From: https://tcga-data.nci.nih.gov/docs/publications/coadread_2012/

Comprehensive Molecular Characterization of Human Colon and Rectal Cancer

Nature, Volume 487 Number 7407, July 19, 2012 [doi:10.1038/nature11252]
Abstract

To characterize somatic alterations in colorectal carcinoma (CRC), we conducted
genome-scale analysis of 276 samples, analyzing exome sequence, DNA copy
number, promoter methylation, mRNA and microRNA expression. A subset (97)
underwent low-depth-of-coverage whole-genome sequencing. 16% of CRC have
hypermutation, three quarters of which have the expected high microsatellite
instability (MSI), usually with hypermethylation and MLH1 silencing, but one
quarter has somatic mismatch repair gene mutations. Excluding hypermutated
cancers, colon and rectum cancers have remarkably similar patterns of genomic
alteration. Twenty-four genes are significantly mutated. In addition to the
expected APC, TP53, SMAD4, PIK3CA and KRAS mutations, we found frequent
mutations in ARID1A, SOX9, and FAM123B/WTX. Recurrent copy number alterations
include potentially drug-targetable amplifications of ERBB2 and newly
discovered amplification of IGF2. Recurrent chromosomal translocations include
fusion of NAV2 and WNT pathway member TCF7L1. Integrative analyses suggest new
markers for aggressive CRC and important role for MYC-directed transcriptional
activation and repression. 

## Data pre-processing

(from the Supplementary Material PDF S1, "Microarray Expression Profiling")

Briefly, 2 ug of total RNA of sample (n=220) and Stratagene Universal Human
Reference were amplified and labeled using Agilent’s Low RNA Input Linear
Amplification Kit. Sample and reference were co-hybridized on a Custom Agilent
244K Gene Expression Microarray (AMDID019760). The expression data was Lowess
normalized and the ratio of the Cy5 channel (sample) and Cy3 channel
(reference) were log2 transformed to create gene expression values for 23,199
probesets. Probesets without gene annotations and genes with missing data in ≥
20% of the samples were removed, resulting in 13,994 genes available for
further analysis. Missing values in the remaining genes were imputed with the
mean value across all samples. PCA (JMP Genomics, v.4.0) analysis indicated
that the source of the RNA (BCR) was responsible for 23% of the variance of the
microarray data, which was normalized out using JMP Genomics’ Batch Correction
procedure.
