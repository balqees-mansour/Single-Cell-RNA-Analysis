# Single-Cell-RNA-seq Analysis

I have replicated two data analyses based on the single cell RNA-seq data from the study "Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq" (https://www.science.org/doi/10.1126/science.aad0501). The study, published on April 8, 2016, aimed to explore the diverse genetic and physical properties of melanoma tumors using single-cell RNA sequencing (RNA-seq) on 4645 individual cells obtained from 17 patients. The sequencing analyzed malignant, immune, stromal, and endothelial cells and revealed transcriptional heterogeneity among the malignant cells within the same tumor, which was related to the cell cycle, local environment, and drug resistance.

I began by downloading the single cell RNA-seq data "GSE72056" (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056) with 4645 cells and 23689 genes, and preprocessing it. I then divided the data into two sets: one set of 1061 malignant cells and another set of 2716 non-malignant cells. The code for the malignant cells can be found at "https://github.com/balqees-mansour/Single-Cell-RNA-Analysis/blob/main/malignant_melanoma_tumors/Malignant.code.md" and the code for the non-malignant cells at "https://github.com/balqees-mansour/Single-Cell-RNA-Analysis/blob/main/non.malignant/Non.malignant.celltypes.md". My next step was to analyze the data using the Seurat package in R.

I selected malignant data from 6 tumors (78,79,88,80,81,89) with more than 50 cells, as specified in the paper. However, I did not include mel84 because it only had 14 cells, so I substituted it with mel78.
Next, I selected non-malignant data from 13 tumors with more than 100 cells, according to the paper. This included Mel (53, 58, 80, 60, 72, 74, 75) instead of (78,79), 67 instead of 81, 84, 89, 94,88. I replaced mel81 and mel78 with others as they had a low cell count less than 100.

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
