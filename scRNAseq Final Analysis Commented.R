library(scCATCH)
library(Seurat)
library(CytoTRACE)
library(ggplot2)
library(sctransform)
library(tidyverse)
library(RColorBrewer)
library(gprofiler2)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(scProportionTest)
library(harmony)
library(cowplot)
library(data.table)
library(RColorBrewer)
library(colorspace)
library(scater)
library(TSCAN)
library(SingleCellExperiment)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
####################################################################################################################
#create merged object
wt_control_data <- Read10X(data.dir = "JG186_control/outs/filtered_feature_bc_matrix/")
wt_control <- CreateSeuratObject(counts = wt_control_data$`Gene Expression`, project = "wt_wt-control")
wt_control_adt <- CreateAssayObject(counts = wt_control_data$`Antibody Capture`)
wt_control[["ADT"]] <- wt_control_adt
wt_control@meta.data$snp <- "wt_wt"
wt_control@meta.data$edit <- "control"

wt_tcl1a_data <- Read10X(data.dir = "JG186_tcl1a/outs/filtered_feature_bc_matrix/")
wt_tcl1a <- CreateSeuratObject(counts = wt_tcl1a_data$`Gene Expression`, project = "wt_wt-tcl1a_overexpression")
wt_tcl1a_adt <- CreateAssayObject(counts = wt_tcl1a_data$`Antibody Capture`)
wt_tcl1a[["ADT"]] <- wt_tcl1a_adt
wt_tcl1a@meta.data$snp <- "wt_wt"
wt_tcl1a@meta.data$edit <- "tcl1a_overexpression"

alt_control_data <- Read10X(data.dir = "JG196_control/outs/filtered_feature_bc_matrix/")
alt_control <- CreateSeuratObject(counts = alt_control_data$`Gene Expression`, project = "alt_alt-control")
alt_control_adt <- CreateAssayObject(counts = alt_control_data$`Antibody Capture`)
alt_control[["ADT"]] <- alt_control_adt
alt_control@meta.data$snp <- "alt_alt"
alt_control@meta.data$edit <- "control"

alt_tcl1a_data <- Read10X(data.dir = "JG196_tcl1a/outs/filtered_feature_bc_matrix/")
alt_tcl1a <- CreateSeuratObject(counts = alt_tcl1a_data$`Gene Expression`, project = "alt_alt-tcl1a_overexpression")
alt_tcl1a_adt <- CreateAssayObject(counts = alt_tcl1a_data$`Antibody Capture`)
alt_tcl1a[["ADT"]] <- alt_tcl1a_adt
alt_tcl1a@meta.data$snp <- "alt_alt"
alt_tcl1a@meta.data$edit <- "tcl1a_overexpression"

merged_object <- merge(wt_control, y = c(wt_tcl1a, alt_control, alt_tcl1a), add.cell.ids = c("WC", "WT", "AC", "AT"), project = "merged")

#sctransform based normalization, regressing out cell cycle, percent mito, and sample origin variability
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
DefaultAssay(merged_object) <- "RNA"
merged_object <- PercentageFeatureSet(merged_object, pattern = "^MT-", col.name = "percent.mt")
merged_object_subset <- subset(merged_object, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA > 2500 & percent.mt < 20)
DefaultAssay(merged_object_subset) <- "RNA"
merged_object_subset <- SCTransform(merged_object_subset, assay = 'RNA',
                                    new.assay.name = 'SCT',
                                    vars.to.regress = c('percent.mt', 'orig.ident'))
merged_object_subset <- CellCycleScoring(
  merged_object_subset,
  s.features = s.genes,
  g2m.features = g2m.genes,
  assay = 'SCT',
  set.ident = TRUE
)
merged_object_subset$CC.Difference <- merged_object_subset$S.Score - merged_object_subset$G2M.Score
merged_object_subset <- SCTransform(
  merged_object_subset,
  assay = 'RNA',
  new.assay.name = 'SCT',
  vars.to.regress = c('percent.mt', 'orig.ident', 'CC.Difference'))

#centered log-ratio (CLR) normalization of ADT
DefaultAssay(merged_object_subset) <- "ADT"
merged_object_subset <- NormalizeData(merged_object_subset, normalization.method = "CLR", margin = 2, assay = "ADT")
saveRDS(merged_object_subset, file = "20mt_merge_then_sct_pre_harmony.rds")

#harmony integration and dimensionality reduction
seurat_obj <- readRDS("20mt_merge_then_sct_pre_harmony.rds")
set.seed(2)
seurat_obj <- seurat_obj %>%
  RunHarmony("orig.ident", assay.use="SCT")
DefaultAssay(seurat_obj) <- "SCT"
seurat_obj <- RunUMAP(seurat_obj, reduction='harmony', dims = 1:30)
seurat_obj <- FindNeighbors(seurat_obj, reduction='harmony')
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
saveRDS(seurat_obj, file = "20mt_merge_then_sct_then_harmony_cell_cycle_regressed_seed2.rds")

#formatting and labeling of seurat object
seurat_obj <- readRDS("20mt_merge_then_sct_then_harmony_cell_cycle_regressed_seed2.rds")
Idents(seurat_obj) <- seurat_obj@meta.data$seurat_clusters
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
new.cluster.ids <- c("Eosinophil", "Neutrophil", "Mast Cell", "GMP",
                     "HSC/MPP 2", "HSC/MPP 1", "CMP", "HSC/MPP 4",
                     "HSC/MPP 3", "Macrophage", "undefined", "undefined", "undefined")
Idents(seurat_obj) <- plyr::mapvalues(x = Idents(seurat_obj), from = current.cluster.ids, to = new.cluster.ids)
seurat_obj@meta.data$new.id <- Idents(seurat_obj)
table(seurat_obj@meta.data$new.id, seurat_obj@meta.data$SCT_snn_res.0.5)
desired_levels <- c('HSC/MPP 1', 'HSC/MPP 2', 'HSC/MPP 3', 'HSC/MPP 4', 
                    'CMP', 'GMP', 'Macrophage','Neutrophil', 'Eosinophil',
                    'Mast Cell','undefined')
Idents(seurat_obj) <- factor(x = Idents(seurat_obj), levels = desired_levels)

#Figures
#Supplemental Table  14###############################################################################################
seurat_obj <- PrepSCTFindMarkers(seurat_obj)
cell_markers <- FindAllMarkers(seurat_obj, assay = "SCT", min.pct =  0.25, only.pos = TRUE)
cell_markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20_cell_markers
write.csv(top20_cell_markers,"top20_genes.csv")

#Figure 5D###########################################################################################################
Idents(seurat_obj) <- seurat_obj$orig.ident
seurat_obj$orig.ident <- factor(x = seurat_obj$orig.ident, levels = c('wt_wt-control', 'wt_wt-tcl1a_overexpression',
                                                                      'alt_alt-control', 'alt_alt-tcl1a_overexpression'))
my_cols <- c("#12477C","#67001E","#A1CDE1","#EC9274")
SplitPlot <- DimPlot(seurat_obj, split.by = 'orig.ident', 
                     ncol = 2, cols = my_cols) & NoLegend() & NoAxes()
desired_levels <- c('HSC/MPP 1', 'HSC/MPP 2', 'HSC/MPP 3', 'HSC/MPP 4', 
                    'CMP', 'GMP', 'Macrophage','Neutrophil', 'Eosinophil',
                    'Mast Cell','undefined')
Idents(seurat_obj) <- factor(x = Idents(seurat_obj), levels = desired_levels)
ClusterPlot <- DimPlot(seurat_obj, label.size = 4) + NoLegend() + NoAxes()
ClusterPlotBold <- LabelClusters(ClusterPlot, id = "ident",  fontface = "bold")
SplitPlot + ClusterPlotBold

#Extended Figure 15C###########################################################################################################
seurat_obj_subset <- subset(x = seurat_obj, idents = c('HSC/MPP 1', 'HSC/MPP 2', 'HSC/MPP 3', 'HSC/MPP 4', 
                                                       'CMP', 'GMP', 'Macrophage','Neutrophil', 'Eosinophil',
                                                       'Mast Cell'))
DimPlot(seurat_obj_subset, label = TRUE)
phase_count_table <- data.frame(unclass(table(Idents(seurat_obj_subset), seurat_obj_subset@meta.data$Phase)))
phase_count_table_percentage <- apply(phase_count_table, 1, function(x){x*100/sum(x,na.rm=T)})
phase_count_tibble <- phase_count_table_percentage %>% as_tibble()
phase_count_tibble$Phase <- c("G1", "G2M", "S")
data_long <- gather(phase_count_tibble, CellType, Percent, "HSC/MPP 1":"Mast Cell", factor_key=TRUE)
data_long$CellType <- factor(data_long$CellType,levels = c('HSC/MPP 1', 'HSC/MPP 2', 'HSC/MPP 3', 'HSC/MPP 4', 
                                                           'CMP', 'GMP', 'Macrophage','Neutrophil', 'Eosinophil',
                                                           'Mast Cell'))
CellCycleByCluster <- ggplot(data_long, aes(fill=Phase, y=Percent, x=CellType)) + 
  geom_bar(position="stack", stat="identity") + 
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Cell Cluster")
CellCycleByCluster

#Extended Figure 16A###########################################################################################################
celltype_table <- data.frame(unclass(table(Idents(seurat_obj_subset), seurat_obj_subset@meta.data$orig.ident)))
celltype_table_percentage <- apply(celltype_table, 2, function(x){x*100/sum(x,na.rm=T)})
celltype_tibble <- celltype_table_percentage %>% as_tibble()
celltype_tibble$CellType <- c('HSC/MPP 1', 'HSC/MPP 2', 'HSC/MPP 3', 'HSC/MPP 4', 
                              'CMP', 'GMP', 'Macrophage','Neutrophil', 'Eosinophil',
                              'Mast Cell')
celltype_data_long <- gather(celltype_tibble, Condition, Percent, "alt_alt.control":"wt_wt.tcl1a_overexpression", factor_key=TRUE)
celltype_data_long$Condition <- factor(celltype_data_long$Condition,levels = c('wt_wt.control', 'wt_wt.tcl1a_overexpression',
                                                                               'alt_alt.control', 'alt_alt.tcl1a_overexpression'))
celltype_data_long$CellType <- factor(celltype_data_long$CellType,levels = c('HSC/MPP 1', 'HSC/MPP 2', 'HSC/MPP 3', 'HSC/MPP 4', 
                                                                             'CMP', 'GMP', 'Macrophage','Neutrophil', 'Eosinophil',
                                                                             'Mast Cell'))
CellTypesByCondition <- ggplot(celltype_data_long, aes(fill=CellType, y=Percent, x=Condition)) + 
  geom_bar(position="stack", stat="identity") + 
  theme_cowplot() +
  background_grid() +
  theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust=1)) +
  scale_x_discrete(labels=c("wt_wt.control" = "Donor 1 (G/G) control eGFP",
                            "wt_wt.tcl1a_overexpression" = "Donor 1 (G/G) TCL1A eGFP",
                            "alt_alt.control" = "Donor 2 (T/T) control eGFP",
                            "alt_alt.tcl1a_overexpression" = "Donor 2 (T/T) TCL1A eGFP")) +
  xlab("Genotype and Lentivirus") + scale_fill_discrete(name = "Cell Clusters")
CellTypesByCondition

#Figure 5F###########################################################################################################
HSC_MPP <- subset(x = seurat_obj, idents = c("HSC/MPP 1", "HSC/MPP 2", "HSC/MPP 3", "HSC/MPP 4"))
table(HSC_MPP@meta.data$orig.ident, HSC_MPP@meta.data$new.id)
RR_prop_test <- sc_utils(HSC_MPP)
RR_prop_test <- permutation_test(
  RR_prop_test, cluster_identity = "new.id",
  sample_1 = "wt_wt-control", sample_2 = "wt_wt-tcl1a_overexpression",
  sample_identity = "orig.ident"
)
RR_plot_data <- copy(RR_prop_test@results$permutation)
RR_plot_data[, clusters := fct_reorder(factor(clusters), obs_log2FD)]
AA_prop_test <- sc_utils(HSC_MPP)
AA_prop_test <- permutation_test(
  RR_prop_test, cluster_identity = "new.id",
  sample_1 = "alt_alt-control", sample_2 = "alt_alt-tcl1a_overexpression",
  sample_identity = "orig.ident"
)
AA_plot_data <- copy(AA_prop_test@results$permutation)
AA_plot_data[, clusters := fct_reorder(factor(clusters), obs_log2FD)]
RR_plot_data$Genotype <- "rs2887399 G/G"
AA_plot_data$Genotype <- "rs2887399 T/T"
combo_prop_plot <- AA_plot_data %>% 
  bind_rows(RR_plot_data)
combo_permutation_plot <- ggplot(combo_prop_plot, aes(x = clusters, y = obs_log2FD)) +
  geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = Genotype),
                  position = position_dodge(width = -0.2)) +
  theme_bw() +
  ylab("Observed Log2FD in TCL1A Overexpression") +
  xlab("Clusters") +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c("salmon", "grey")) +
  theme(axis.title.x = element_text(face="bold", size=10)) +
  theme(axis.title.y = element_text(face="bold", size=10)) +
  coord_flip()
combo_permutation_plot

#Extended Figure 15A###########################################################################################################
DefaultAssay(seurat_obj_subset) <- "ADT"
ADT_plots <- FeaturePlot(seurat_obj_subset, features = c("CD34-TotalSeqB","CD38-TotalSeqB","CD45RA-TotalSeqB",
                                                         "CD49f-TotalSeqB","CD11a-TotalSeqB","CD117-TotalSeqB"),
                         ncol=2,cols = c("yellow", "blue"),order=T, min.cutoff = 'q5', max.cutoff = 'q95')
ADT_plots

#Extended Figure 15D###########################################################################################################
DefaultAssay(seurat_obj_subset) <- "SCT"
stress_UMAP_plots <- FeaturePlot(seurat_obj_subset, 
                                 features = c("GADD45B","GDF15","DDIT3","ATF4","TXNIP","BBC3","CDKN1A",
                                              "PPP1R15A"),
                                 ncol=2,cols = c("yellow", "blue"),order=T)
stress_UMAP_plots

#Figure 5E###########################################################################################################
DefaultAssay(seurat_obj_subset) <- "SCT"
Feature_Dot_Plot <- DotPlot(seurat_obj_subset, features = c("BAD","TXNIP","BBC3","GADD45B","CCNG2","CDKN1B","CDKN1A","BCL2L11","SOD2",
                                                            "NFE2L2","GDF15","ATF4","DDIT3","PPP1R15A",
                                                            "CENPA","PCNA","TOP2A","MKI67",
                                                            "TPSAB1","CLC","CSF2RB","MPO","ELANE","SPI1","MRC1","HBB",
                                                            "ITGA2B","MEF2C","IGSF10","HEMGN","MECOM","FAM30A"),
                            cols = "RdYlBu") +
  coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  cowplot::theme_cowplot() + RotatedAxis()
Feature_Dot_Plot

#Extended Figure 16B###########################################################################################################
complete_celltype <- data.frame(unclass(table(Idents(seurat_obj), seurat_obj@meta.data$orig.ident)))
complete_celltype_table_percentage <- as.data.frame(apply(complete_celltype, 2, function(x){x/sum(x,na.rm=T)}))
complete_celltype_table_percentage$WWC_counts <- complete_celltype_table_percentage$wt_wt.control*5443
complete_celltype_table_percentage$WWT_counts <- complete_celltype_table_percentage$wt_wt.tcl1a_overexpression*11882
complete_celltype_table_percentage$AAC_counts <- complete_celltype_table_percentage$alt_alt.control*15403
complete_celltype_table_percentage$AAT_counts <- complete_celltype_table_percentage$alt_alt.tcl1a_overexpression*43705
cell_counts <- complete_celltype_table_percentage %>% select(contains('count'))
hsc_cell_counts <- filter(cell_counts, grepl("HSC", row.names(cell_counts))) %>% rownames_to_column()
hsc_cell_counts_long <- gather(hsc_cell_counts, celltype, measurement, WWC_counts:AAT_counts, factor_key=TRUE)
cfu_y_label <- "Absolute Counts Per 1000 HSCs at Day 0"
hsc_cell_counts_plot <- hsc_cell_counts_long %>%
  ggplot(aes(x=celltype, y=measurement, fill=rowname)) +
  geom_bar(stat="identity")+
  theme_cowplot() +
  background_grid() + 
  labs(x = NULL, y = cfu_y_label) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_x_discrete(labels=c("WWC_counts" = "Donor 1 (G/G) control eGFP",
                            "WWT_counts" = "Donor 1 (G/G) TCL1A eGFP",
                            "AAC_counts" = "Donor 2 (T/T) control eGFP",
                            "AAT_counts" = "Donor 2 (T/T) TCL1A eGFP")) +
  xlab("Genotype and Lentivirus") + scale_fill_discrete(name = "Cell Clusters")
hsc_cell_counts_plot

#Extended Figure 15B###########################################################################################################
set.seed(1234)
HSC_MPP_CMP <- subset(x = seurat_obj, idents = c("HSC/MPP 1", "HSC/MPP 2", "HSC/MPP 3", "HSC/MPP 4", "CMP"))
hsc_mpp.cds <- as.cell_data_set(HSC_MPP_CMP)
hsc_mpp.cds <- cluster_cells(cds = hsc_mpp.cds, reduction_method = "UMAP")
hsc_mpp.cds <- learn_graph(hsc_mpp.cds, use_partition = TRUE)
hsc_mpp.cds <- order_cells(hsc_mpp.cds, reduction_method = "UMAP") #choose root node in middle of HSC/MPP 1 cluster
HSC_MPP_CMP <- AddMetaData(
  object = HSC_MPP_CMP,
  metadata = hsc_mpp.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Lineage"
)
desired_levels_hsc_mpp <- c('HSC/MPP 1', 'HSC/MPP 2', 'HSC/MPP 3', 'HSC/MPP 4', 
                            'CMP')
Idents(HSC_MPP_CMP) <- factor(x = Idents(HSC_MPP_CMP), levels = desired_levels_hsc_mpp)
HSC_MPP_GMP_Dimplot <- DimPlot(HSC_MPP_CMP) + NoAxes()
LineageUMAP <- FeaturePlot(HSC_MPP_CMP, c("Lineage"), pt.size = 0.1) & scale_color_viridis_c()
HSC_MPP_GMP_Dimplot + LineageUMAP



