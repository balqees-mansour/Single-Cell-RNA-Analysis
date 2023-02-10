# Single-Cell-RNA-seq Analysis
I Reproduce two data analyses using the single cell RNA-seq data of the following study: https://www.science.org/doi/10.1126/science.aad0501

In order to reproduce the experiment published on 8 April 2016 (Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq, "https://www.science.org/doi/10.1126/science.aad0501") I started with preparing the data using R and bash, my first step,is to download the single cell RNA seq data " GSE72056" from here  ( https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056 ) with 4645 cells and 23689 genes. then preprocess the data and subset our data into two distinct sets, malignant data with 1061 cells, which you can find here"https://github.com/balqees-mansour/Single-Cell-RNA-Analysis/blob/main/malignant_melanoma_tumors/Malignant.code.md", and Non-malignant data with 2716 cells, which you can find it here"https://github.com/balqees-mansour/Single-Cell-RNA-Analysis/blob/main/non.malignant/Non.malignant.celltypes.md
My second step is to conduct analysis using Suerat Package in R then I.

The authors investigated the diverse genetic and physical characteristics of melanoma tumors, they utilized single-cell RNA sequencing (RNA-seq) on 4645 individual cells taken from 17 patients. This sequencing analyzed malignant cells, immune cells, stromal cells, and endothelial cells and uncovered transcriptional diversity among the malignant cells within the same tumor, which was linked to factors such as the cell cycle, local environment, and drug resistance.
 


