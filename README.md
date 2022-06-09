# HSC_TCL1A_overexpression_scRNAseq
TCL1A overexpression scRNAseq analysis in R for Weinstock &amp; Gopakumar paper

The BCL files were demultiplexed using 8 base pair 10X sample indexes and cellranger mkfastq to generate paired-end FASTQ. We ran cellranger count to align the reads to the hg38 reference genome from GenBank using STAR aligner as well as perform filtering, barcode counting, and UMI counting. The alignment results were used to quantify the expression level of human genes and generation of gene-barcode matrix.

Each sample's cellranger matrix was then loaded in a SeuratObject_4.1.0 using Seurat_4.1.1 (https://github.com/satijalab/seurat). Low quality cells, doublets and potential dead cells were removed according to the percentage of mitochondrial genes and number of genes and UMIs expressed in each cell (nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA > 2500 & percent.mt < 10). Clean count matrixes from each sample were then combined using Seurat’s merge function. The merged gene expression data was normalized using sctransform based normalization while removing confounding variables, percentage of mitochondrial genes and sample origin. Then, cell cycle scores were assigned using Seurat's CellCycleScoring function. The difference between the G2M and S phase scores was then calculated and regressed out using sctransform based normalization to minimize differences due to differences in cell cycle phase among proliferating cells. The cell surface feature output was normalized using centered log-ratio (CLR) normalization, computed independently for each feature.

The 4 datasets were integrated using Harmony (https://github.com/immunogenomics/harmony) on SCTranform normalized gene counts to group cells by cell type while correcting for sample origin. Dimensionality reduction via PCA and UMAP embedding was performed on the integrated dataset. Identities of the cell clusters were determined using canonical RNA cell type markers and cell surface feature expression patterns. HSC/MPP clusters were identified by staining positively for CD34 and CD49f, and negatively for CD38, CD45RA, and CD11a. The common myeloid progenitor cluster was identified by staining positively for CD34 and CD38, and negatively for CD45RA and CD49f. The granulocyte macrophage progenitor cluster was identified by staining positively for CD34, CD38, and CD45RA, and negatively for CD49f. The difference between the proportion of cells in HSC/MPP 1-4 clusters between control-eGFP and TCL1A-eGFP transduced cells was calculated by a proportion test using the Single Cell Proportion Test R package (https://github.com/rpolicastro/scProportionTest). To reconstruct the pseudotime trajectory of the HSC/MPP and CMP clusters, Monocle 3 pseudotime analysis was performed using the central node of the HSC/MPP1 Cluster as the root node.
  
Comments describing steps of analysis and generated figures are included in .R file.
