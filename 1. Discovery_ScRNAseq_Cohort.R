library(stringr)
library(data.table)
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(SingleR)
library(harmony)
library(loomR)
library(SeuratDisk)
library(openxlsx)
library(DESeq2)
library(GENIE3)
library(Seurat)

 
#Schlenker metadata from ArrayExpress (https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-13770?query=E-MTAB-13770)
ramona_metadata <- read.delim2("C:/E-MTAB-13770.tsv")

# loom_dir is the directory where TILs scRNA data from Schlenker is stored in the internal PC
loom_dir <- "C:/Users/tils"
loom_files <- list.files(loom_dir, pattern = "\\.loom$", full.names = TRUE)

 
seurat_list <- list()

# Loop through each loom file
for (file_path in loom_files) {
  
  # Extract Patient ID from filename: P123_sample.loom -> P123
  file_name <- basename(file_path)
  patient_id <- str_extract(file_name, "P\\d+")
 
  loom <- Connect(filename = file_path)
  
  # Convert to Seurat object 
  seu <- as.Seurat(loom, counts = "spliced")
  
  # Add patient ID as metadata
  seu$Patient <- patient_id
  
   loom$close_all()
  
  # Add to list
  seurat_list[[length(seurat_list) + 1]] <- seu
  
}

# Merge all Seurat objects of Schlenker et al. 
schlenker_object <- merge(seurat_list[[1]], y = seurat_list[-1])


ramona_metadata$Patient <- paste0("P",ramona_metadata$individual)
pts <- unique(schlenker_object$Patient)
#keep only TILs samples, not PBMCs and only patients who were included in the Seurat object
ramona_metadata1 <- ramona_metadata[ramona_metadata$fraction == "CD45+",]
ramona_metadata1 <- ramona_metadata1[ramona_metadata1$Patient %in% pts, ]
ramona_metadata2 <- ramona_metadata1[ramona_metadata1$Characteristics.Table.source_chars.E.MTAB.13770. != "source_22",]

seurat_meta <- schlenker_object@meta.data
seurat_meta$cell_id <- rownames(seurat_meta)

new_meta <- left_join(seurat_meta, ramona_metadata2, by = "Patient")

rownames(new_meta) <- new_meta$cell_id
new_meta$cell_id <- NULL  

# Add to Seurat object
schlenker_object <- AddMetaData(schlenker_object, metadata = new_meta)
#Add Tumor Free samples to Responders
schlenker_object$Response <- ifelse(
  grepl("^NR_", schlenker_object$response.to.treatment), "NR",
  ifelse(schlenker_object$response.to.treatment == "TF", "R", NA)
)
schlenker_object$Study <- "Schlenker"
schlenker_object$Patient <- paste0("Ramona_",schlenker_object$Patient)

#Pozniak et al. dataset
#Load the Entire_TME object downloaded from KU Leuven server (https://rdr.kuleuven.be/dataset.xhtml?persistentId=doi:10.48804/GSAXBN)
Entire_TME <- readRDS("C:/Users/Entire_TME.rds")
pozniak_object = Entire_TME
#Keep only pre treatment samples
pozniak_object <- subset(pozniak_object, subset = `BT/OT` == "BT")
#Match to binary response based on the Supplementary table 1 information included in the original publication
pozniak_object$Response <- ifelse(pozniak_object$orig.ident %in% c("sc5rCMA070","scrCMA068","scrCMA109","sc5rCMA141","sc5rCMA144","sc5rCMA192","scrCMA038","scrCMA093","scrCMA094","scrCMA119"
) , "R" ,"NR")
pozniak_object$Patients <- pozniak_object$orig.ident
gene_names <- rownames(pozniak_object)
pozniak_object$Study <- "Pozniak"

seurat_list <- list(pozniak_object,schlenker_object)


# Get the common features (genes) across both Seurat objects
common_genes <- Reduce(intersect, lapply(seurat_list, function(seu) rownames(seu)))

# Make sure each Seurat object has the same set of genes
rna_list <- lapply(seurat_list, function(seu) {
  GetAssayData(seu, assay = "RNA", slot = "counts")[common_genes, ]
})

# Combine the RNA count matrices with the common genes
combined_rna <- do.call(cbind, rna_list)
# Replace NA with 0 
combined_rna[is.na(combined_rna)] <- 0

sparse_mat <- Matrix(combined_rna, sparse = TRUE)
seurat_merged <- CreateSeuratObject(counts = sparse_mat, project = "scRNA_Disc")

#Get union of all metadata column names
all_columns <- Reduce(union, lapply(seurat_list, function(seu) colnames(seu@meta.data)))

# Merge metadata and add NAs if they do not exist in the different metadata 
metadata_list <- lapply(seurat_list, function(seu) {
  meta <- seu@meta.data
  meta$cell_id <- rownames(meta)   
  
  missing_cols <- setdiff(all_columns, colnames(meta))
  for (col in missing_cols) {
    meta[[col]] <- NA
  }
  
  meta <- meta[, c(all_columns, "cell_id")]
  return(meta)
})

# Combine metadata
combined_metadata <- do.call(rbind, metadata_list)

# Match metadata to seurat_merged columns (cells)
combined_metadata <- combined_metadata[match(colnames(seurat_merged), combined_metadata$cell_id), ]
#Assign metadata
combined_metadata$cell_id <- NULL
seurat_merged@meta.data <- combined_metadata

seurat_merged[["percent.mt"]] <- PercentageFeatureSet(seurat_merged, pattern = "^MT-")
all_counts <- Matrix::colSums(GetAssayData(seurat_merged, slot = "counts"))

# Calculate HB genes expression per cell
hb_genes <- grep("^HB[ABDEGMZ]", rownames(seurat_merged), value = TRUE)
hb_expr <- Matrix::colSums(GetAssayData(seurat_merged, slot = "counts")[hb_genes,]) 

total_expr <- Matrix::colSums(GetAssayData(seurat_merged, slot = "data"))

# Percent HB expression per cell
percent_hb <- 100 * hb_expr / all_counts

# Filter cells with >5% HB gene expression to avoid the presence of red blood cells
seurat_no_hb <- subset(seurat_merged, cells = names(percent_hb[percent_hb < 5]))

seurat_merged$percent.hb <- percent_hb 

VlnPlot(seurat_merged, features = "percent.hb")

seurat_merged[["nCount_RNA"]] <- Matrix::colSums(GetAssayData(seurat_merged, assay = "RNA", slot = "counts"))
seurat_merged[["nFeature_RNA"]] <- Matrix::colSums(GetAssayData(seurat_merged, assay = "RNA", slot = "counts") > 0)

VlnPlot(seurat_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Patient")
#filter to exclude doublets ( based on nFeatures), cells with reduced info, low quality cells and erythrocytes
fil_seurat_merged <- subset(seurat_merged, subset = nFeature_RNA > 500 & percent.mt < 20 &  percent_hb < 5 & nFeature_RNA < 5000)



# Normalize the combined RNA assay, FindVariables and Scale with default settings
fil_seurat_merged <- NormalizeData(fil_seurat_merged, assay = "RNA")
fil_seurat_merged <- FindVariableFeatures(fil_seurat_merged)
fil_seurat_merged <- ScaleData(fil_seurat_merged)


# Run PCA on the normalized RNA assay 
fil_seurat_merged <- RunPCA(fil_seurat_merged, npcs = 30, verbose = T)

# Run Harmony based on Study AND Patients' samples to harmonize the UMAP space so that we can do good cell type annotation
fil_seurat_merged <-harmony::RunHarmony(
  object = fil_seurat_merged,
  group.by.vars = c("Study","Patient"),
)

# UMAP after Harmony
fil_seurat_merged <- RunUMAP(fil_seurat_merged, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony_pat")

fil_seurat_merged <- FindNeighbors(fil_seurat_merged, reduction = "harmony", dims = 1:30)
fil_seurat_merged <- FindClusters(fil_seurat_merged, resolution = 0.4)

cluster_markers <- FindAllMarkers(fil_seurat_merged, only.pos = TRUE,min.pct = 0.6, slot = "data",verbose = TRUE)

# Load Tirosh markers from the responding publication (https://doi.org/10.1126/science.aad0501)
tirosh_markers <- read.xlsx("C:\\Users\\tirosh_markers.xlsx")
#Load AUCell to match Non-Immune cell types
library(AUCell)
s.data <- GetAssayData(object = fil_seurat_merged, layer = "counts")
cell_rankings <- AUCell_buildRankings(s.data)

melanoma_markers <- tirosh_markers$melanoma
melanoma_markers <- melanoma_markers[!is.na(melanoma_markers)]
#calculate AUC score for each gene 
cells_mel_AUC <- AUCell_calcAUC(melanoma_markers, cell_rankings)

#get AUC histogram to check if threshold is correctly assigned. In this case it is
cells_mel_assignment <- AUCell_exploreThresholds(cells_mel_AUC, plotHist = TRUE, assign = TRUE)

#so we are using this to see the cells that match to melanoma cells
mel.cells <- cells_mel_assignment$geneSet$assignment 
#We do the same with Markers of Cancer Associated Fibroblasts
caf_markers <- tirosh_markers$CAFs
caf_markers <- caf_markers[!is.na(caf_markers)]
#calculate AUC score for each gene 
cells_caf_AUC <- AUCell_calcAUC(caf_markers, cell_rankings)

#get AUC histogram to check if threshold is correctly inputed.  
cells_caf_assignment <- AUCell_exploreThresholds(cells_caf_AUC, plotHist = TRUE, assign = TRUE)

#so we are using this to see the cells that match to ca fibroblasts
caf.cells <- cells_caf_assignment$geneSet$assignment 

#We do the same with Markers of Endothelial cells

endoth_markers <- tirosh_markers$Endothelial.cells
endoth_markers <- endoth_markers[!is.na(endoth_markers)]
#calculate AUC score for each gene 
cells_endoth_AUC <- AUCell_calcAUC(endoth_markers, cell_rankings)

#get AUC histogram to check if threshold is correctly inputed.  
cells_endoth_assignment <- AUCell_exploreThresholds(cells_endoth_AUC, plotHist = TRUE, assign = TRUE)

#so we are using this to see the cells that match to endothelial cells
endoth.cells <- cells_endoth_assignment$geneSet$assignment 

#Add the cells we found, to the cells they are mapped and View them compared to the Seurat Clusters 
fil_seurat_merged$tirosh_type <- ifelse(colnames(fil_seurat_merged) %in% caf.cells, "CAFs", ifelse(colnames(fil_seurat_merged) %in% mel.cells, "Melanoma", ifelse(colnames(fil_seurat_merged) %in% endoth.cells, "Endothelial", "Other")))

DimPlot(object = fil_seurat_merged,reduction = "umap_harmony_pat", group.by = c("tirosh_type","seurat_clusters"), label = TRUE)

#The cells map specific clusters with great precision so we assign the clusters to these cell types:

fil_seurat_merged$cell.type1 <- ifelse(fil_seurat_merged$seurat_clusters %in% c(15), "Endothelial", 
                                       ifelse(fil_seurat_merged$seurat_clusters %in% c(10, 13,17), "CAFs", 
                                              ifelse(fil_seurat_merged$seurat_clusters == 4, "Melanoma", "Other")))


ptprc_feat <- FeaturePlot(fil_seurat_merged, features = "PTPRC", reduction = "umap_harmony_pat") +
  ggtitle("PTPRC expression") +
  xlab("UMAP1") +
  ylab("UMAP2")

ggsave("C:/Users/Supplementary_Figure_1b.bmp", plot = ptprc_feat, width = 8, height = 6, dpi = 300)

#Rest of the clusters are PTPRC (CD45+) positive meaning they are Immune cells.
#We will annotate them using BlueprintEndodeData as reference from the SingleR R package
ref <- celldex::BlueprintEncodeData()
expr_matrix <- LayerData(fil_seurat_merged, assay = "RNA", layer = "data")

singleR_results_blueprint_m <- SingleR::SingleR(test = expr_matrix, ref = ref, labels = ref$label.main )
#singleR_results_blueprint_f <- SingleR::SingleR(test = expr_matrix, ref = ref, labels = ref$label.fine )

fil_seurat_merged$SingleR_blueprint_m <- singleR_results_blueprint_m$labels
#fil_seurat_merged$SingleR_blueprint_f <- singleR_results_blueprint_f$labels
#Assign the cells to their SingleR label based on BluePrintEncode reference label MAIN not fine
fil_seurat_merged$cell.type_m1 <- ifelse(fil_seurat_merged$cell.type1 == "Other",fil_seurat_merged$SingleR_blueprint_m, fil_seurat_merged$cell.type1  ) 

DimPlot(object = fil_seurat_merged,reduction = "umap_harmony_pat", group.by = c("cell.type_m1","seurat_clusters"), label = TRUE)
fil_seurat_merged$cell.type_m1[fil_seurat_merged$cell.type_m1 == "CAFs"] <- "Fibroblasts"

table(fil_seurat_merged$cell.type_m1)
#Remove cell types with less than 200 cells to reduce unwanted noise in the differential analyses
cell_counts <- table(fil_seurat_merged$cell.type_m1)

fil_seurat_merged <- subset(
  fil_seurat_merged,
  subset = !(cell.type_m1 %in% names(cell_counts[cell_counts < 200]))
)



cell.type_m1_markers_tissue <- FindAllMarkers(fil_seurat_merged, group.by = "cell.type_m1", test.use = "wilcox",slot = "data" ,min.pct = 0.4, only.pos = TRUE)

library(openxlsx)

# Split the dataframe by cell type
celltype_list_tissue <- split(cell.type_m1_markers_tissue, cell.type_m1_markers_tissue$cluster)

 wb <- createWorkbook()

 for (ct in names(celltype_list_tissue)) {
  addWorksheet(wb, sheetName = substr(ct, 1, 31))  # Excel sheet names limited to 31 chars
  writeData(wb, sheet = ct, celltype_list_tissue[[ct]])
}

# Save the workbook showcasing the marekrs of each assigned cell type
saveWorkbook(wb, file = "C:/Users/celltypemarkers.xlsx", overwrite = TRUE)


library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(viridis)


library(dplyr)
library(dplyr)

top_markers_all_tissue <- cell.type_m1_markers_tissue %>%
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
Idents(fil_seurat_merged) <- "cell.type_m1"

set.seed(111)


cells_by_type_tissue <- lapply(names(top_genes_list_tissue), function(ct) {
  all_cells <- WhichCells(fil_seurat_merged, idents = ct)
  sample(all_cells, size = min(700, length(all_cells)))
})
names(cells_by_type_tissue) <- names(top_genes_list_tissue)

cells_ordered_tissue <- unlist(cells_by_type_tissue)

expr_mat_tissue <- GetAssayData(fil_seurat_merged, slot = "data")[top_genes_tissue, cells_ordered_tissue ]
expr_scaled_tissue <- t(scale(t(as.matrix(expr_mat_tissue))))

 cell.type_m1s_ordered_tissue <- fil_seurat_merged@meta.data[cells_ordered_tissue, "cell.type_m1"]
cell.type_m1_levels_tissue <- names(top_genes_list_tissue)
cell.type_m1s_ordered_tissue <- factor(cell.type_m1s_ordered_tissue, levels = cell.type_m1_levels_tissue)
cell.type_m1_colors_tissue <- setNames(rainbow(length(cell.type_m1_levels_tissue)), cell.type_m1_levels_tissue)

ha <- HeatmapAnnotation(
  CellType = cell.type_m1s_ordered_tissue,
  col = list(CellType = cell.type_m1_colors_tissue),
  show_annotation_name = FALSE
)

 row_split <- rep(names(top_genes_list_tissue), times = sapply(top_genes_list_tissue, length))
column_split <- cell.type_m1s_ordered_tissue

bmp("C:/Users/supplementary_figure2A.bmp", width = 2200, height = 1600, res = 200)
col_fun = colorRamp2(
  breaks = c(-1, 0, 2),  
  colors = c("magenta3", "gray21", "yellow")  
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
        column_split = cell.type_m1s_ordered_tissue,
        column_title = NULL,
        row_title = NULL,
        border = TRUE)

dev.off()








study_mel_tis <- DimPlot(fil_seurat_merged, reduction = "umap_harmony_pat", group.by = "Study", label = TRUE)
DimPlot(fil_seurat_merged, reduction = "umap_harmony_pat", group.by = "cell.type_m1", label = TRUE)
#match some Response data which were not parsed from metadata probably due to added " " in the excel cell 
fil_seurat_merged$Response <- ifelse(fil_seurat_merged$Patient == "Ramona_P72", "R", fil_seurat_merged$Response )
fil_seurat_merged$Response <- ifelse(fil_seurat_merged$Patient == "Ramona_P63", "R", fil_seurat_merged$Response )
fil_seurat_merged$Response <- ifelse(fil_seurat_merged$Patient == "Ramona_P40", "R", fil_seurat_merged$Response )

#FINAL UMAPs
ct_mel_tis <- DimPlot(fil_seurat_merged, reduction = "umap_harmony_pat", group.by = "cell.type_m1", label = FALSE) + 
  ggtitle("Tissue Melanoma - Cell types") +
  xlab("UMAP1") +
  ylab("UMAP2")
ggsave("C:/Users/Figure2A.bmp", plot = ct_mel_tis, width = 8, height = 6, dpi = 300)

study_mel_tis <- DimPlot(fil_seurat_merged, reduction = "umap_harmony_pat", group.by = "Study") + 
  ggtitle("Tissue Melanoma - Studies") +
  xlab("UMAP1") +
  ylab("UMAP2")
ggsave("C:/Users/Supplementary_Figure1A.bmp", plot = study_mel_tis, width = 8, height = 6, dpi = 300)
#Create the PSEUDOBULK objects
cell_type_patient_pseudob <- AggregateExpression(fil_seurat_merged, assays = "RNA",return.seurat = T, group.by = c("Patient_celltype1","Response","cell.type_m1",
                                                                                                                   "Study"))
Idents(cell_type_patient_pseudob) <- "cell.type_m1"
cell_types <- unique(fil_seurat_merged$cell.type_m1)

# Initialize list to store DESeq2 results
nde_markers_by_cell_type <- list()

# Loop through each cell type to find markers comparing R vs NR using DESeq2
for (cell_type in cell_types) {
  # Subset cells in the current cell_type
  cell_type_cells <- WhichCells(cell_type_patient_pseudob, id = cell_type)
  
  cell_type_obj <- subset(cell_type_patient_pseudob, cells = cell_type_cells)
  
   Idents(cell_type_obj) <- "Response"
   
  # We make sure counts matrix is in integer mode
  counts_matrix <- as.matrix(GetAssayData(cell_type_obj, assay = "RNA", slot = "counts"))
  counts_matrix <- round(counts_matrix)
  
  metadata <- as.data.frame(cell_type_obj@meta.data)
  metadata$Response <- factor(metadata$Response)
  metadata$Study <- factor(metadata$Study)
  
  # Ensure the order of metadata matches the columns of the counts matrix
  metadata <- metadata[match(colnames(counts_matrix), rownames(metadata)), ]

    # Create DESeq2 design formula, including Study as a covariate to exclude the technical interstudy batch efffects from the results 
  dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                colData = metadata,
                                design = ~ Study + Response)  # Include Study as covariate
  
  # Run DESeq2 analysis
  dds <- DESeq(dds)
  
  #  Get DESeq2 results (R vs NR comparison)
  res <- results(dds, contrast = c("Response", "R", "NR"))
  
  res$cell_type <- cell_type
  res$genes <- rownames(res)
  res$comparison <- "R_vs_NR"
  nde_markers_by_cell_type[[paste0("cell_type_", cell_type)]] <- as.data.frame(res)
  
  
  # Clean up  
  rm(cell_type_obj)
  gc()
}

wb <- createWorkbook()



 for (name in names(nde_markers_by_cell_type)) {
  truncated_name <- gsub("cell_type_","", name)
  addWorksheet(wb, truncated_name) 
  writeData(wb, truncated_name, nde_markers_by_cell_type[[name]])  
}

saveWorkbook(wb, file =
               "C:/Users/supplementary_table_1.xlsx", overwrite = T)


Tcells <-subset(fil_seurat_merged, subset = (cell.type_m1 %in% c( "CD4+ T-cells"  , "CD8+ T-cells")))

library(ProjecTILs)
list.reference.maps()
list_of_ref.maps <- get.reference.maps(collection = "human")

list_of_ref_maps <- get.reference.maps(
  reference = c("CD4", "CD8"))

Tcells <- NormalizeData(Tcells)
Tcells <- FindVariableFeatures(Tcells, selection.method = "vst", nfeatures = 3000)
Tcells <- ScaleData(Tcells)


Tcells <- RunPCA(Tcells, npcs = 30, verbose = T)
Tcells <-harmony::RunHarmony(
  object = Tcells,
  group.by.vars = c("Study","Patient"))

Tcells <- RunUMAP(Tcells, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony_pat")

CD4 <- subset(Tcells, subset = cell.type_m1 == "CD4+ T-cells")

CD4 <- NormalizeData(CD4)
CD4 <- FindVariableFeatures(CD4, selection.method = "vst", nfeatures = 3000)
CD4 <- ScaleData(CD4)


CD4 <- RunPCA(CD4, npcs = 30, verbose = T)
CD4 <-harmony::RunHarmony(
  object = CD4,
  group.by.vars = c("Study","Patient"))

CD4 <- RunUMAP(CD4, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony_pat")


CD8 <- subset(Tcells, subset = cell.type_m1 == "CD8+ T-cells")
CD8 <- NormalizeData(CD8)
CD8 <- FindVariableFeatures(CD8, selection.method = "vst", nfeatures = 3000)
CD8 <- ScaleData(CD8)


CD8 <- RunPCA(CD8, npcs = 30, verbose = T)
CD8 <-harmony::RunHarmony(
  object = CD8,
  group.by.vars = c("Study","Patient"))

CD8 <- RunUMAP(CD8, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony_pat")

DimPlot(Tcells, group.by = "cell.type_m1", reduction = "umap_harmony_pat", label = TRUE)

# Extract CD4 and CD8 references from the list
ref_cd4 <- list_of_ref_maps$human$CD4
ref_cd8 <- list_of_ref_maps$human$CD8

CD4_Proj <- Run.ProjecTILs(query = CD4, ref = ref_cd4, filter.cells = F)
DimPlot(CD4_Proj, group.by = "functional.cluster", label = TRUE, repel = TRUE)
plot.projection(ref_cd4, CD4_Proj, linesize = 0.5, pointsize = 0.5)
 

CD8_Proj <- Run.ProjecTILs(query = CD8, ref = ref_cd8, filter.cells = F)
DimPlot(CD8_Proj, group.by = "functional.cluster", label = TRUE, repel = TRUE)
plot.projection(ref_cd8, CD8_Proj, linesize = 0.5, pointsize = 0.5)

CD4meta <- CD4_Proj@meta.data
CD8meta <- CD8_Proj@meta.data
tcm <- rbind(CD4meta, CD8meta)

TcellsMeta <- Tcells@meta.data
tcmr= tcm[rownames(TcellsMeta),]

Tcells <- AddMetaData(object = Tcells, metadata = tcmr)




DimPlot(Tcells, group.by = "functional.cluster",reduction = "umap_harmony_pat", label = TRUE, repel = TRUE)

unwanted_types <- c(
  "Class-switched memory B-cells",
  "CLP",
  "DC",
  "Macrophages",
  "Memory B-cells",
  "GMP",
  "Monocytes",
  "naive B-cells",
  "NK cells", 'Plasma cells'
)

# Remove these from the dataset
Tcells <- subset(
  Tcells, 
  subset = !(SingleR_blueprint_f %in% unwanted_types)
)

# Create a named vector mapping old names to new names
cell_type_map <- c(
  "CD4.CTL_EOMES"   = "CD4+ Cytotoxic EOMES+",
  "CD4.CTL_Exh"     = "CD4+ Cytotoxic Exhausted",
  "CD4.CTL_GNLY"    = "CD4+ Cytotoxic GNLY+",
  "CD4.Memory"      = "CD4+ Memory",
  "CD4.NaiveLike"   = "CD4+ Naive-like",
  "CD4.Tfh"         = "CD4+ Tfh",
  "CD4.Th17"        = "CD4+ Th17",
  "CD4.Treg"        = "CD4+ Treg",
  "CD8.CM"          = "CD8+ Central Memory",
  "CD8.EM"          = "CD8+ Effector Memory",
  "CD8.MAIT"        = "CD8+ MAIT Cell",
  "CD8.NaiveLike"   = "CD8+ Naive-like",
  "CD8.TEMRA"       = "CD8+ TEMRA Cell",
  "CD8.TEX"         = "CD8+ Exhausted",
  "CD8.TPEX"        = "CD8+ Progenitor Exhausted"
)


for (fname in names(cell_type_map)) {
  Tcells$cell.type[Tcells$cell.type == fname] <- cell_type_map[fname]
}
cell_type_plot <- DimPlot(Tcells, reduction = "umap_harmony_pat", group.by = "cell.type") +
  ggtitle("T-cell Subsets - Discovery Cohort") +
  xlab("UMAP1") +
  ylab("UMAP2")

ggsave("C:/Users/Figure_3A.bmp", plot = cell_type_plot, width = 8, height = 6, dpi = 300)

SaveSeuratRds(Tcells, "C:/Users/TcellsObj.Rdata")

library(Seurat)
library(openxlsx)

# genes of interest
genes <- unique(c(
  "LGALS1","IL12RB2","TNFRSF18","TNFRSF9","GK","UHRF2",
  "TP63","VCAM1","DUSP4","TNFRSF4","LINGO1","SELL","CD96","EPB41"
))


 Idents(Tcells) <- "cell.type"

cell_types <- levels(Tcells)

wb <- createWorkbook()

for (ct in cell_types) {
  
  # differential expression: ct vs all others
  de <- FindMarkers(
    object = Tcells,
    ident.1 = ct,
    logfc.threshold = 0,
    min.pct = 0.1,
    test.use = "wilcox"
  )
  
  # keep only genes of interest that exist
  de_sub <- de[intersect(genes, rownames(de)), , drop = FALSE]
  
  if (nrow(de_sub) == 0) next
  
   df <- data.frame(
    gene = rownames(de_sub),
    avg_log2FC = de_sub$avg_log2FC,
    p_val = de_sub$p_val,
    p_val_adj = de_sub$p_val_adj,
    pct_in_cluster = de_sub$pct.1,
    pct_other = de_sub$pct.2
  )
  
  df <- df[order(df$avg_log2FC, decreasing = TRUE), ]
  
  sheet_name <- substr(gsub("[^[:alnum:] ]", "", ct), 1, 31)
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, df)
}

saveWorkbook(wb, "C:/Users//supplementary_table4.xlsx", overwrite = TRUE)

library(Seurat)
library(DESeq2)
library(dplyr)

genes_focus <- c("TNFRSF18","TNFRSF9","TNFRSF4","SDC4",
                 "DUSP4","VCAM1","BATF", "IL12RB2")

celltypes_focus <- c(
  "CD4+ Cytotoxic Exhausted",
  "CD4+ Memory",
  "CD4+ Naive-like",
  "CD4+ Treg",
  "CD8+ Central Memory",
  "CD8+ Effector Memory",
  "CD8+ Exhausted",
  "CD8+ MAIT Cell",
  "CD8+ Naive-like",
  "CD8+ Progenitor Exhausted",
  "CD8+ TEMRA Cell"
)

Tcells$Patient_celltype <- paste(Tcells$Patient, Tcells$cell.type, sep="_")

pb_counts_obj <- AggregateExpression(
  Tcells,
  assays = "RNA",
  group.by = c("Patient_celltype","cell.type","Patient", "Response"),
  slot = "counts",
  return.seurat = TRUE
)


all_counts <- GetAssayData(pb_counts_obj, slot="counts")  
meta <- pb_counts_obj@meta.data

meta <- meta[meta$cell.type %in% celltypes_focus, ]
all_counts <- all_counts[, rownames(meta)]


dds <- DESeqDataSetFromMatrix(
  countData = all_counts,
  colData = meta,
  design = ~ 1  # no design, just normalize
)

dds <- estimateSizeFactors(dds)

norm_counts <- counts(dds, normalized=TRUE)


norm_counts_focus <- norm_counts[rownames(norm_counts) %in% genes_focus, ]



library(ggplot2)
library(dplyr)
library(tidyr)


colnames(norm_counts_focus) <- as.character(colnames(norm_counts_focus))

meta$Patient_celltype <- as.character(meta$Patient_celltype)
library(tibble)
library(stringr)

 df_long <- as.data.frame(norm_counts_focus) %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "Patient_celltype",
    values_to = "expression"
  )


df_long <- df_long %>%
  mutate(
    Patient = str_extract(Patient_celltype, "^[^_]+"),         
    cell.type = str_replace(Patient_celltype, "^[^_]+_", "")   
  )

# Order cell types from CD4 to CD8

celltypes_focus <- c(
  "CD4+ Cytotoxic Exhausted",
  "CD4+ Memory",
  "CD4+ Naive-like",
  "CD4+ Treg",
  "CD8+ Central Memory",
  "CD8+ Effector Memory",
  "CD8+ Exhausted",
  "CD8+ MAIT Cell",
  "CD8+ Naive-like",
  "CD8+ Progenitor Exhausted",
  "CD8+ TEMRA Cell"
)
df_long$cell.type <- NULL

meta$Patient = meta$Patient_celltype

df_long <- df_long %>%
  left_join(meta[, c("Patient","cell.type","Response")], by = "Patient")

# Plot: violins split by Response within each cell type
TCELGOF_Response <- ggplot(df_long, aes(x = cell.type, y = expression, fill = Response)) +
  geom_violin(position = position_dodge(width = 0.9), scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), color = "black", outlier.shape = NA) +
  facet_wrap(~gene, scales = "free_y", ncol = 4) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1),
    strip.text = element_text(size = 9)
  ) +
  labs(
    x = "",
    y = "DESeq2-normalized counts",
    fill = "Response"
  ) 

 


TCELGOF = ggplot(df_long, aes(x = cell.type, y = expression, fill = cell.type)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.08, outlier.shape = NA, color = "black") +
  facet_wrap(~gene, scales = "free_y", ncol = 4) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1),
    legend.position = "none",
    strip.text = element_text(size = 9)
  ) +
  labs(
    x = " ",
    y = "DESeq2-normalized counts"
  )
ggsave("C:/Users/Figure3B.bmp", plot = TCELGOF, width = 8, height = 6, dpi = 300)


