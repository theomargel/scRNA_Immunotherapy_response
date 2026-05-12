library(limma)
library(openxlsx)
library(Seurat)
library(data.table)
library(Matrix)

#The tpm counts txt file was downloaded from GEO and saved as an excel TPM.csv (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575)
GSE120575_counts <- fread("C:\\Users\\Sade-Feldman M, 2018\\TPM.csv")
#Keep only Pre-Treatment samples
GSE120575_counts_subset <- GSE120575_counts_subset[,GSE120575_counts_subset[,!grepl("^Post",GSE120575_counts_subset[2,])]]
cols2 <- !grepl("^Post", GSE120575_counts[2,])
table(cols2)
GSE120575_counts_subset <- GSE120575_counts[,..cols2]
#load metadata (this is a table derived from GSE120575_patient_ID_single_cells.txt.gz in the GEO database after removeing the first rows)
GSE120575_metadata <- fread("C:\\Users\\Sade-Feldman M, 2018\\metadata.csv")
GSE120575_metadata <- GSE120575_metadata[,-1]
colnames(GSE120575_metadata) <- as.character(GSE120575_metadata[1,])
GSE120575_metadata <- GSE120575_metadata[-1,]
GSE120575_counts_subset <- GSE120575_counts_subset[-2,]
colnames(GSE120575_counts_subset) <- as.character(GSE120575_counts_subset[1,])
GSE120575_counts_subset <- GSE120575_counts_subset[-1,]
GSE120575_counts_subset <- as.data.frame(GSE120575_counts_subset)
rownames(GSE120575_counts_subset) <- as.character(GSE120575_counts_subset[,1])
GSE120575_counts_subset <- GSE120575_counts_subset[,-1]

length(GSE120575_counts_subset1[,1])
nrow(GSE120575_counts_subset1)


mat_GSE <- as.matrix(GSE120575_counts_subset)
rownames(mat_GSE) <- rownames(GSE120575_counts_subset)
mode(mat_GSE) <- "numeric"
sparse_mat <- Matrix(mat_GSE, sparse = TRUE)
GSE120575_metadata <- as.data.frame(GSE120575_metadata)

#keep only pre-treatment samples 
GSE120575_metadata_rows <- GSE120575_metadata[grepl("^Pre_",GSE120575_metadata$Patient),]
GSE120575_metadata_rows <- GSE120575_metadata[!is.na(GSE120575_metadata$Patient),]
GSE120575_metadata_rows <- as.data.frame(GSE120575_metadata_rows)

row.names(GSE120575_metadata_rows) <- GSE120575_metadata_rows$title
row.names(GSE120575_metadata_rows) <- gsub("-",".",row.names(GSE120575_metadata_rows))
identical(row.names(GSE120575_metadata_rows), colnames(GSE120575_counts_subset))
setdiff(colnames(GSE120575_counts_subset), row.names(GSE120575_metadata_rows))
setdiff(row.names(GSE120575_metadata_rows), colnames(GSE120575_counts_subset))
 
row.names(GSE120575_metadata_rows) <- ifelse(
  row.names(GSE120575_metadata_rows) %in% names(name_corrections),
  name_corrections[row.names(GSE120575_metadata_rows)],
  row.names(GSE120575_metadata_rows)
)

GSE120575_metadata_rows$title = NULL
GSE120575_metadata_rows$`Sample name` = NULL
GSE120575_metadata_rows$`source name` = NULL
#Create seurat object 
seurat_object <- CreateSeuratObject( sparse_mat, assay = "RNA", meta.data = GSE120575_metadata_rows, project = "GSE120575" )

# The [[ operator can add columns to object metadata.  
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
# Show QC metrics for the first 5 cells
head(seurat_object@meta.data, 5)
seurat_object <- subset(seurat_object, features = rownames(seurat_object)[rowSums(seurat_object@assays$RNA$counts) > 0])
# Visualize QC metrics as a violin plot
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Patient")

FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

save(seurat_object, file="C:\\UsersSade-Feldman M, 2018\\seurat_object.RData")

#Filtering genes with less than 500 genes  and more than  5000 features (possible doublets)  and less than 20% mt genes (dying cells)
fil_seurat_object <- subset(seurat_object, subset = nFeature_RNA > 500 & percent.mt < 20 &  nFeature_RNA < 5000)

VlnPlot(fil_seurat_object, features = "nCount_RNA", layer="counts", group.by="Response",raster=FALSE,alpha=0.2, fill.by = "Patient") 
#scale_fill_manual(values=seurat_object@meta.data$Response) 


#normalization
fil_seurat_object <- NormalizeData(fil_seurat_object)

fil_seurat_object <- FindVariableFeatures(fil_seurat_object, selection.method = "vst", nfeatures = 3000)

# plot variable features with and without labels
VariableFeaturePlot(fil_seurat_object)


# scaling so that the mean of each gene equals to 0 
fil_seurat_object <- ScaleData(fil_seurat_object, features = NULL)
#runing PCA to reduce dimensions, taking into account only Variable Features 
fil_seurat_object <- RunPCA(fil_seurat_object, features = VariableFeatures(object = fil_seurat_object))

#visualize cells
DimPlot(fil_seurat_object, reduction = "pca", group.by = "Patient") 



#visualize genes
VizDimLoadings(fil_seurat_object, reduction = "pca", dims = 1:2)
ElbowPlot(fil_seurat_object)

#data integration
fil_seurat_object <-harmony::RunHarmony(
  object = fil_seurat_object,
  group.by.vars = "Patient",
)

# KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods 
fil_seurat_object <- FindNeighbors(fil_seurat_object, dims = 1:20, reduction = "harmony")
fil_seurat_object <- FindClusters(fil_seurat_object, resolution = 0.4)

# Look at cluster IDs of the first 5 cells
head(Idents(fil_seurat_object), 5)

fil_seurat_object <- RunUMAP(fil_seurat_object, reduction = "harmony", reduction.name = "harmony.umap", dims = 1:20)



ref <- celldex::BlueprintEncodeData()

expr_matrix <- GetAssayData(fil_seurat_object, assay = "RNA", slot = "data")


singleR_results_blueprint_f <- SingleR::SingleR(test = expr_matrix, ref = ref, labels = ref$label.main )

fil_seurat_object$cell_type <- singleR_results_blueprint_f$labels


relabel_as_monocytes <- c("DC", "Macrophages", "Neutrophils")

# Reassign these to "Monocytes"
fil_seurat_object$cell_type[fil_seurat_object$cell_type %in% relabel_as_monocytes] <- "Monocytes"






cell_type_markers_tissue <- FindAllMarkers(fil_seurat_object, group.by = "cell_type", test.use = "wilcox",slot = "data" ,min.pct = 0.4, only.pos = TRUE)

library(openxlsx)

# Split the dataframe by cell type
celltype_list_tissue <- split(cell_type_markers_tissue, cell_type_markers_tissue$cluster)

 wb <- createWorkbook()

 for (ct in names(celltype_list_tissue)) {
  addWorksheet(wb, sheetName = substr(ct, 1, 31))  
  writeData(wb, sheet = ct, celltype_list_tissue[[ct]])
}

 saveWorkbook(wb, file = "C:/Users/cell_type_markersSADE.xlsx", overwrite = TRUE)


library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(viridis)


library(dplyr)
library(dplyr)

top_markers_all_tissue <- cell_type_markers_tissue %>%
  filter(avg_log2FC > 0) %>%
  arrange(cluster, p_val_adj, desc(avg_log2FC))

top_genes_list_tissue <- list()
used_genes_tissue <- c()

for (ct in unique(top_markers_all_tissue$cluster)) {
  ct_genes <- top_markers_all_tissue %>%
    filter(cluster == ct, !gene %in% used_genes_tissue) %>%
    slice_head(n = 8)
  top_genes_list_tissue[[ct]] <- ct_genes$gene
  used_genes_tissue <- union(used_genes_tissue, ct_genes$gene)
}

 top_genes_tissue <- unlist(top_genes_list_tissue)
Idents(fil_seurat_object) <- "cell_type"

set.seed(111)


cells_by_type_tissue <- lapply(names(top_genes_list_tissue), function(ct) {
  all_cells <- WhichCells(fil_seurat_object, idents = ct)
  sample(all_cells, size = min(700, length(all_cells)))
})
names(cells_by_type_tissue) <- names(top_genes_list_tissue)

cells_ordered_tissue <- unlist(cells_by_type_tissue)

expr_mat_tissue <- GetAssayData(fil_seurat_object, slot = "data")[top_genes_tissue, cells_ordered_tissue ]
expr_scaled_tissue <- t(scale(t(as.matrix(expr_mat_tissue))))

# Annotations
cell_types_ordered_tissue <- fil_seurat_object@meta.data[cells_ordered_tissue, "cell_type"]
cell_type_levels_tissue <- names(top_genes_list_tissue)
cell_types_ordered_tissue <- factor(cell_types_ordered_tissue, levels = cell_type_levels_tissue)
cell_type_colors_tissue <- setNames(rainbow(length(cell_type_levels_tissue)), cell_type_levels_tissue)

ha <- HeatmapAnnotation(
  CellType = cell_types_ordered_tissue,
  col = list(CellType = cell_type_colors_tissue),
  show_annotation_name = FALSE
)

# Row and column splits  
row_split <- rep(names(top_genes_list_tissue), times = sapply(top_genes_list_tissue, length))
column_split <- cell_types_ordered_tissue

bmp("C:/Users/Supplementary_Figure2b.bmp", width = 2200, height = 1600, res = 200)
col_fun = colorRamp2(
  breaks = c(-1, 0, 2),  # Explicit breaks at -1, 0, and 2
  colors = c("magenta3", "gray21", "yellow")  # Pink, black, then yellow
)
Heatmap(expr_scaled_tissue,
        name = "Expression",
        top_annotation = ha,
        col = col_fun,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_column_names = FALSE,
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 8),
        column_split = cell_types_ordered_tissue,
        column_title = NULL,
        row_title = NULL,
        border = TRUE)

dev.off()

 
ct_mel_tis <- DimPlot(fil_seurat_object, reduction = "harmony.umap", group.by = "cell_type", label = FALSE) + 
  ggtitle("Validation Melanoma Tissue - Cell types") +
  xlab("UMAP1") +
  ylab("UMAP2")
ggsave("C:/Users/Figure2B.bmp", plot = ct_mel_tis, width = 8, height = 6, dpi = 300)

# Create a unique patient-celltype label to get a pseudobulk object that has the relevant metadata info
fil_seurat_object$Patient_celltype_manual <- paste(fil_seurat_object$Patient, fil_seurat_object$cell_type, sep = "_")

#  Aggregate expression per patient-celltype (NOTE:: DATA already log2 TPM)
cell_type_patient_pseudob <- AggregateExpression(
  fil_seurat_object,
  assays = "RNA",
  return.seurat = TRUE,
  group.by = c("Patient_celltype_manual", "Response", "cell_type" )
)

 Idents(cell_type_patient_pseudob) <- "cell_type"
cell_types <- unique(fil_seurat_object$cell_type)

 tpm_markers_by_cell_type <- list()

#  Loop through each cell type
for (cell_type in cell_types) {
  cat("Processing:", cell_type, "\n")
  
  # Subset pseudobulk samples for this cell type
  cell_type_obj <- subset(cell_type_patient_pseudob, idents = cell_type)
  
  # Extract log2 TPM matrix !! 
  log_tpm <- as.matrix(GetAssayData(cell_type_obj, assay = "RNA", slot = "data"))
  
   metadata <- cell_type_obj@meta.data
  metadata <- metadata[colnames(log_tpm), c("Response")]
  metadata <- as.data.frame(metadata)
  metadata$Response <-metadata$metadata 
  metadata$metadata = NULL
  # Filter out cell types with too few R/NR samples
  r_samples <- sum(metadata$Response == "Responder")
  nr_samples <- sum(metadata$Response == "Non-responder")
  if (r_samples < 3 || nr_samples < 3) {
    cat("Skipping", cell_type, "- Not enough R (", r_samples, ") or NR (", nr_samples, ") samples.\n")
    next
  }
  
   design <- model.matrix(~  Response, data = metadata)
  
  # Compute array weights and fit linear model with limma, this procedure is the best for continuous data like TPM
  weights <- arrayWeights(log_tpm, design)
  
  # 
  fit <- lmFit(log_tpm, design, weights = weights)
  fit <- eBayes(fit, trend = TRUE)
  
  res <- topTable(fit, number = Inf, sort.by = "P")
  res$cell_type <- cell_type
  res$genes <- rownames(res)
  res$comparison <- "R_vs_NR"

  tpm_markers_by_cell_type[[paste0("cell_type_", cell_type)]] <- res
  
  # Clean up
  rm(cell_type_obj)
  gc()
}

wb <- createWorkbook()

for (name in names(tpm_markers_by_cell_type)) {
  truncated_name <- substr(name, 1, 31)
  addWorksheet(wb, truncated_name)
  writeData(wb, truncated_name, tpm_markers_by_cell_type[[name]])
}

# Save workbook
saveWorkbook(wb, file = "C:\\Users\\supplementary_table_2.xlsx", overwrite = TRUE)
