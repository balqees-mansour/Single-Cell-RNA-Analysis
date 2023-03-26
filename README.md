# Single-Cell-RNA-seq Analysis

I have replicated two data analyses based on the single cell RNA-seq data from the study "Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq" (https://www.science.org/doi/10.1126/science.aad0501). The study, published on April 8, 2016, aimed to explore the diverse genetic and physical properties of melanoma tumors using single-cell RNA sequencing (RNA-seq) on 4645 individual cells obtained from 17 patients. The sequencing analyzed malignant, immune, stromal, and endothelial cells and revealed transcriptional heterogeneity among the malignant cells within the same tumor, which was related to the cell cycle, local environment, and drug resistance.

I began by downloading the single cell RNA-seq data "GSE72056" (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056) with 4645 cells and 23689 genes, and preprocessing it, I then divided the data into two sets: one set of 1061 malignant cells and another set of 2716 non-malignant cells, The code for the malignant cells can be found at:
#### https://github.com/balqees-mansour/Single-Cell-RNA-Analysis/tree/main/Malignant
and the code for the non-malignant cells at :
#### https://github.com/balqees-mansour/Single-Cell-RNA-Analysis/blob/main/non.malignant/Non.malignant.celltypes.md

My next step was to analyze the data using the Seurat package in R.

I selected malignant data from 6 tumors (78,79,88,80,81,89) with more than 50 cells, as specified in the paper. However, I did not include mel84 because it only had 14 cells, so I substituted it with mel89.
Next, I selected non-malignant data from 13 tumors with more than 100 cells, according to the paper. This included Mel (67, 75,79 ,88, 80, 89, 84,94 ,53 ,58, 60, 72, 74). I replaced mel81 and mel78 with others as they had a low cell count less than 100.

![image](https://user-images.githubusercontent.com/87857777/227800001-ac30904e-47bf-485a-bd53-ccea37e7b6d0.png)




 
## My analysis was performed using Seurat:

### Step 0: Import the Seurat package.

### Step 1: Create a Seurat object.

### Step 2: Quality Control.
Following the creation of the Seurat object, I conducted quality control on the data. This typically involved filtering:

Cells with too few genes detected, as they often represent cells that were not sequenced in enough depth to be reliably characterized.
Cells with too many genes detected, as they may represent doublets or multiplets.
Cells with a high mitochondrial transcript percentage.
### Step 3: Normalization and Log Transformation.
To make gene expression levels between different cells comparable, a normalization step was necessary. The most commonly used normalization in scRNA-seq data analysis is similar to TPM (Transcripts Per Million reads). This involves normalizing the feature expression measurements for each cell to the total expression, then multiplying by a scale factor (default of 10000). The expression values were then log-transformed to better fit a normal distribution. It's worth noting that a pseudocount was added to each value before the log-transformation to ensure genes with zero transcripts detected in a cell still had zero values after the log-transform. I used the normalization method "LogNormalize" which combines these two steps into one parameter.
### Step 4: Feature Selection for Heterogeneity Analysis
I identified highly variable features/genes, which are genes with the highest expression variation across cells, for further heterogeneity analysis.
 
### Step 5. Data Scaling
Since different genes have varying base expression levels and distributions, their contribution to the analysis may be unequal without any data transformation. To prevent our analysis from relying solely on highly expressed genes, data scaling is applied using the selected features, as is standard in data science.

### Step 6. Linear Dimensionality Reduction using Principal Component Analysis (PCA)
The benefits of dimension reduction include:

Data becomes more compact, resulting in faster computation.
In scRNA-seq data, which is inherently sparse, summarizing measurements of related features significantly improves signal robustness.
I used 15 PCs, as specified in the paper.
Step 7. Non-linear Dimension Reduction for Visualization
I ran t-SNE on the two sets of data with resolutions of 0.1 for the malignant subset and 0.03 for the non-malignant subset.

### Step 8. Identification of Markers and Cell Clustering
I used the Human Cell Atlas (https://www.humancellatlas.org/) to determine cell types based on their markers.

### Non-malignant cell types clusters 
Cluster 0: 👍 T cells CD2,  CD3D,  CD3E,  IL32,  NKG7

Cluster1 : 👍 B cells CD79A BANK1 CD19 IGLL5.

Cluster 2: 👍 Macrophages FCER1G CD14 TYROBO  CST3  C1QB. 

Cluster 3: 👍 CAF DCN LUM COL3A1 COL1A2 COL1A1 

Cluster 4: 👍 Endothelial cells  CLDN5 CCL21 EFEMP1 IGFBP7 TFPI.

Cluster 5: 👍 NK  AQP1 PLVAP RAMP3 SPARCL1 IGFBP7.

## Conclusions : I will continue with the remaining analysis until I confirm these results.

### 1.importance of in-tercellular communication for tumor phenotype.

The interactions between different cell types within a tumor play a crucial role in determining the overall behavior and characteristics of the tumor. This can involve the transfer of signals, growth factors, and other molecules between cells, which can affect the gene expression profiles and other properties of different cell types and contribute to the development and progression of the tumor. The importance of intercellular communication in tumor phenotype highlights the complexity of the tumor microenvironment and the need for a systems-level understanding of the underlying processes that shape tumor biology.

### 2.interactions between CAFs and immune cell profile

The abundance of immune cells in melanoma core biopsies is influenced by factors derived from the tumor stroma. This implies that, instead of considering the overall expression level of genes in the bulk tissue sample, the composition of different cell types within the tumor should be taken into account when developing diagnostic and therapeutic strategies for melanoma. The authors are suggesting that a more nuanced understanding of the cellular interactions within the tumor will lead to more effective treatments.

### 3.a subset of genes expressed by one cell type (CAFs) may influence the proportion of other cell types (T cells).

Single-cell profiles from a few tumors to understand the composition of a large collection of bulk profiles from the Cancer Genome Atlas (TCGA). By analyzing the single-cell profiles, the researchers discovered various microenvironments that were associated with different malignant cell types. They also identified a subset of genes that were expressed by one specific cell type, such as cancer-associated fibroblasts (CAFs), that may have an impact on the proportion of other cell types in the tumors, such as T cells.

### 4.leveraging single-cell profiles for a few tumors to deconvolve a large collection of bulk

The use of single-cell profiles to analyze a large set of bulk tissue samples from The Cancer Genome Atlas (TCGA). By breaking down bulk tissue into its individual cell components and analyzing the unique profiles of these individual cells, researchers can gain more detailed insights into the composition of the tissue and the types of cells that are present. This information can be used to better understand the underlying biology of cancer and to inform the development of new treatments.

### 5.potential bio-markers for distinguishing exhausted and cyto-toxic T cells that may aid in selecting patients for immune checkpoint blockade

The presence of certain markers that could differentiate between exhausted and cytotoxic T cells. These markers could then be used to determine which patients would be most suitable for treatment with immune checkpoint blockade therapy. In other words, by identifying these biomarkers, researchers and clinicians can potentially identify which patients are more likely to respond well to this type of treatment.

## (Cell-cell communications and Bulk RNA Analysis).
******To be continued ****** 
