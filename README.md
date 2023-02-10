# Single-Cell-RNA-seq Analysis

I have replicated two data analyses based on the single cell RNA-seq data from the study "Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq" (https://www.science.org/doi/10.1126/science.aad0501). The study, published on April 8, 2016, aimed to explore the diverse genetic and physical properties of melanoma tumors using single-cell RNA sequencing (RNA-seq) on 4645 individual cells obtained from 17 patients. The sequencing analyzed malignant, immune, stromal, and endothelial cells and revealed transcriptional heterogeneity among the malignant cells within the same tumor, which was related to the cell cycle, local environment, and drug resistance.

I began by downloading the single cell RNA-seq data "GSE72056" (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056) with 4645 cells and 23689 genes, and preprocessing it. I then divided the data into two sets: one set of 1061 malignant cells and another set of 2716 non-malignant cells. The code for the malignant cells can be found at "https://github.com/balqees-mansour/Single-Cell-RNA-Analysis/blob/main/malignant_melanoma_tumors/Malignant.code.md" and the code for the non-malignant cells at "https://github.com/balqees-mansour/Single-Cell-RNA-Analysis/blob/main/non.malignant/Non.malignant.celltypes.md". My next step was to analyze the data using the Seurat package in R.


 


