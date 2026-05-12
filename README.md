# Integrative single-cell profiling of melanoma reveals a tumor microenvironment signature predictive of immunotherapy response

This repository contains the R code, analytical pipelines, and model development scripts for the study.

## Project Overview
Immunotherapy has transformed cancer treatment, yet many patients fail to respond. This study presents a stepwise integrative analysis of pre-treatment single-cell and bulk tissue datasets to identify predictive molecular features across diverse patient populations. 

---

## Data Resources
The analysis integrates multiple publicly available datasets as described in the manuscript:
Specifically, the datasets analyzed during the current study are available in the Gene Expression Omnibus (GEO), ArrayExpress, and the KU Leuven Research Data Repository. Specifically, single-cell data were retrieved under accessions E-MTAB-13770, GSE120575 and DOI: 10.48804/GSAXBN . The Bulk RNA-seq data were retrieved from GSE78220, GSE213145, E-MTAB-11729 and https://github.com/riazn/bms038_analysis. 
Competing interests


---

## Analysis Scripts & Reproducibility

### `1.Discovery_ScRNAseq_Cohort.R`
This script handles the processing of the primary melanoma scRNA-seq discovery cohort.
**Workflow:** Quality control (gene counts, mitochondrial/hemoglobin content), log-normalization, and Harmony-based integration to mitigate study-specific batch effects.

**Annotation:** Utilizes *Tirosh et al.* markers for non-immune cells and *SingleR* with the *Blueprint Encode* database for immune populations. Functional T-cell states are annotated via *ProjectTILs*.

**Statistics:** Implements a pseudobulk approach by aggregating UMI counts per patient and cell type using `AggregateExpression`. This treats the patient as the unit of observation to prevent false-positive inflation. Differential expression is performed via `DESeq2`, including 'Study' as a covariate to account for technical variability.

 **Key Results Generated:**
**Figure 2A:** UMAP of cellular composition in the discovery cohort.
     
**Figure 3A:** UMAP of T-cell functional states.
     
 **Figure 3B:** Violin plots showing DESeq2-normalized expression of key regulatory genes in T-cell functional states (TNFRSFs, BATF, VCAM1 etc.).
     
**Supplementary Figure 1A,B:** Study distribution and PTPRC (CD45+) expression UMAPs.
     
  **Supplementary Figure 2A:** Heatmap showcasing the top differentially expressed marker genes used to validate the identity of each assigned cell type in the discovery cohort.
  
**Supplementary Table 1:** Pseudobulk DE results (Rs vs NRs) per cell type.
    
**Supplementary Table 4:** Markers of T-cell functional states, focusing on the 13BM and Table 2 (original publication) genes.


### `2.Validation_ScRNAseq_Cohort.R`
This script validates findings in the independent Sade-Feldman et al. scRNA-seq Transcripts Per Million dataset.

**Workflow:** Processing of TPM values and patient-level pseudobulk aggregation.

**Statistics:** Employs `limma` with sample-specific array weights and empirical Bayes moderation to account for continuous TPM data noise for the Rs vs NRs comparison.

**Key Visuals Generated:**

 **Figure 2B:** UMAP of the single-cell validation cohort.
    
**Supplementary Figure 2B:** Heatmap showcasing the top differentially expressed marker genes used to validate the identity of each assigned cell type in the validation cohort.
    
 **Supplementary Table 2:** `limma` pseudobulk results for Rs vs NRs.

### `3. Bulk_RNAseq_Model_Development.R`

This script integrates the bulk transcriptomics data and builds the 13BM predictive model and does gene set overrepresentation analysis based on its correlating genes.

**Workflow:** Batch effect correction using `ComBat_seq` while preserving biological variation and then DESeq2 analysis to finid the DE genes (Rs vs NRs) in bulk setting. 

**Model:** Logistic regression trained on 13 prioritized genes (13BM) with 10-fold cross-validation.

**Key Visuals Generated:**

**Figure 4A-C:** ROC curves for the Melanoma training cohort, Bladder cancer validation, and low-TMB subgroup.
    
 **Figure 5:** Dot plots of Gene Set Overrepresentation Analysis (ORA) based on 13BM score correlations.
    
 **Supplementary Table 5:** DESeq2 results for the integrated bulk melanoma dataset.
    
 **Supplementary Table 6:** GO: Biological Processes gene sets overrepresentation results of the 13BM signature score correlator genes (positive and negative).

---

## Citation
Please cite the following paper if you use this code or the 13BM signature in your research:


