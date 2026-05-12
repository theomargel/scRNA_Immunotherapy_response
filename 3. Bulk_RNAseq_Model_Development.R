library(verification)
library(readODS)
library(readxl)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(sva)
library(pROC)
library(caret)
library(psych)
library(openxlsx)
library(readxl)
library(dplyr)

#SET PATHS  data where downloaded by the servers mentioned in the Data Availability section
input_dir <- "C:/dummy_path/input/"
output_dir <- "C:/dummy_path/output/"

#  Trilla-Fuertes Data Preparation
TrillaCounts <- read.xlsx(paste0(input_dir, "Trilla/RNA_raw_clean_counts.xlsx"))
TrillaMetadata <- read.xlsx(paste0(input_dir, "Trilla/metadata_Trilla_Fuertes.xlsx"))

rownames(TrillaCounts) <- TrillaCounts$Symbol
TrillaCounts$Symbol <- NULL
TrillaMetadata$Study <- "Trilla"
TrillaMetadata$Response <- ifelse(TrillaMetadata$`Better.response.to.anti-PD-1.treatment.1:.CR.2:.PR.3:.SS.4:.PD.9996:.NE` %in% c(1,2), "R", 
                                  ifelse(TrillaMetadata$`Better.response.to.anti-PD-1.treatment.1:.CR.2:.PR.3:.SS.4:.PD.9996:.NE` %in% c(3,4), "NR", NA))

TrillaMetadata <- TrillaMetadata[!is.na(TrillaMetadata$Response),]
ptsf_trilla <- intersect(TrillaMetadata$Unique.Sample.ID, colnames(TrillaCounts))
TrillaCounts <- TrillaCounts[, ptsf_trilla]
TrillaMetadata <- TrillaMetadata[TrillaMetadata$Unique.Sample.ID %in% ptsf_trilla,]

#  Hugo Data Preparation
HugoMetadata <- read.xlsx(paste0(input_dir, "Hugo/Metadata.xlsx"))
HugoCounts <- read.xlsx(paste0(input_dir, "Hugo/clean/GSE78220_raw_counts_GRCh38NCBI_cleaned.xlsx"))

rownames(HugoCounts) <- HugoCounts$Symbol
HugoCounts$Symbol <- NULL
HugoMetadata$Study <- "Hugo"
HugoMetadata$Response <- HugoMetadata$Binary.Response
HugoMetadata$Unique.Sample.ID <- HugoMetadata$Run 

# Hossain Data Preparation
HossainCounts <- read.xlsx(paste0(input_dir, "Hossain/Raw_Clean_Counts.xlsx"))
HossainMetadata <- read.xlsx(paste0(input_dir, "Hossain/Metadata_Hossain.xlsx"))

rownames(HossainCounts) <- HossainCounts$Symbol
HossainCounts$Symbol <- NULL
HossainMetadata$Study <- "Hossain"
HossainMetadata$Response <- HossainMetadata$Response.based.on.RECIST
HossainMetadata$Unique.Sample.ID <- HossainMetadata$Sample.Unique.ID

# Riaz Data Preparation
RiazMetadata <- read.xlsx(paste0(input_dir, "Riaz/GSE91061_metadata.xlsx"))
RiazCounts <- read.xlsx(paste0(input_dir, "Riaz/clean/GSE91061_raw_counts_GRCh38NCBI_cleaned.xlsx"))

rownames(RiazCounts) <- RiazCounts$Symbol
RiazCounts$Symbol <- NULL
RiazMetadata$Study <- "Riaz"
RiazMetadata <- RiazMetadata[RiazMetadata$`Binary.Response.(CRPR/SDPD)` != "NE" & RiazMetadata$Sample.collection == "PRE",]
RiazMetadata$Response <- RiazMetadata$`Binary.Response.(CRPR/SDPD)`
RiazMetadata$Unique.Sample.ID <- RiazMetadata$UniqueSampleID

ptsf_riaz <- intersect(RiazMetadata$Unique.Sample.ID, colnames(RiazCounts))
RiazCounts <- RiazCounts[, ptsf_riaz]
RiazMetadata <- RiazMetadata[RiazMetadata$Unique.Sample.ID %in% ptsf_riaz,]



common_genes <- Reduce(intersect, list(rownames(RiazCounts), rownames(HugoCounts), rownames(HossainCounts), rownames(TrillaCounts)))

metan <- rbind(RiazMetadata[, c("Unique.Sample.ID","Response","Study")], 
               HossainMetadata[, c("Unique.Sample.ID","Response","Study")],
               HugoMetadata[, c("Unique.Sample.ID","Response","Study")],
               TrillaMetadata[, c("Unique.Sample.ID","Response","Study")])

metan$Response <- ifelse(metan$Response == "N", "NR", metan$Response)

metan_counts <- cbind(RiazCounts[common_genes,], 
                      HossainCounts[common_genes,],
                      HugoCounts[common_genes,], 
                      TrillaCounts[common_genes,])

metan_counts$gene <- rownames(metan_counts)

# Create the Excel files 
write.xlsx(metan_counts, paste0(output_dir, "4datasetsallgenesNEW.xlsx"))
write.xlsx(metan, paste0(output_dir, "4datasetsmetadataNEW1.xlsx"))

data <- read_excel(paste0(output_dir, "4datasetsallgenesNEW.xlsx"))
data <- as.data.frame(data)
rownames(data) <- data$gene
data$gene <- NULL

metadata <- read_excel(paste0(output_dir, "4datasetsmetadataNEW1.xlsx"))
metadata <- metadata %>%
  dplyr::select(Unique.Sample.ID, Response, Study)

metadata$Study <- as.factor(metadata$Study)
metadata$Response <- as.factor(metadata$Response)
rownames(metadata) <- metadata$Unique.Sample.ID

######TRAINING DATA PREPARATION AND BIOMARKER DISCOVERY  (MELANOMA BULK RNA-SEQ)

 

#arrange metadata in same sample order as data
metadata <- metadata[colnames(data), ]
rownames(metadata) <- metadata$Unique.Sample.ID

#PCA before batch correction and deseq2 normalisation
data_t <- t(data)
pca <- prcomp(data_t, center = TRUE, scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  study = metadata$Study
)
ggplot(pca_df, aes(x = PC1, y = PC2, color = study)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_classic() +
  labs(
    title = "PCA of Gene Expression Before Normalisation",
    x = paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca)$importance[2, 2], 1), "%)")
  )

#remove batch effects
batch <- metadata$Study
group <- metadata$Response
levels(metadata$Response)
data_combat <- ComBat_seq(as.matrix(data), batch = batch, group=group)
data <- data_combat

#deseq2 normalisation and DE analysis
stopifnot(identical(colnames(data), rownames(metadata)))
levels(metadata$Response)
dds <- DESeqDataSetFromMatrix(
  countData = data,
  colData   = metadata,
  design    = ~ Response
)
dds <- DESeq(dds)

data_norm <- counts(dds, normalized = TRUE)

de_deseq2 <- as.data.frame(results(dds))
colnames(de_deseq2) <- paste0("bulk_melanoma_deseq2_", colnames(de_deseq2))
de_deseq2$genes <- row.names(de_deseq2)
write.xlsx(de_deseq2, "supplementary_table_5.xlsx")
#PCA after correction/normalisation
data_norm_t <- t(data_norm)
pca_norm <- prcomp(data_norm_t, center = TRUE, scale. = TRUE)
pca_norm_df <- data.frame(
  PC1 = pca_norm$x[, 1],
  PC2 = pca_norm$x[, 2],
  study = metadata$Study
)
ggplot(pca_norm_df, aes(x = PC1, y = PC2, color = study)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_classic() +
  labs(
    title = "PCA of Gene Expression After Batch Effect Removal and Normalisation",
    x = paste0("PC1 (", round(100 * summary(pca_norm)$importance[2, 1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca_norm)$importance[2, 2], 1), "%)")
  )

###selection of top genes based on scRNA analysis
#load validated genes from melanoma scRNA-seq analysis
scrna_path <- "C:/Users/supplementary_table_3.ods"
scrna_sheet_names <- ods_sheets(scrna_path)
scrna_list <- lapply(scrna_sheet_names, function(s) read_ods(scrna_path, sheet = s))
names(scrna_list) <- scrna_sheet_names
scrna_list <- scrna_list[sapply(scrna_list, function(df) is.data.frame(df) && nrow(df) > 0)]
scrna_df <- do.call(
  rbind,
  lapply(names(scrna_list), function(nm) {
    df <- scrna_list[[nm]]
    df$single_cell_source <- nm
    df
  })
)

#keep melanoma scRNA-seq validated genes that have p < 0.05 and agree in FC in melanoma bulk rna-seq analysis
top_genes_df <- inner_join(scrna_df, de_deseq2, by = c("Genes" = "genes")) %>%
  filter(bulk_melanoma_deseq2_pvalue < 0.05) %>%
  filter(
    (`Discovery: log2FoldChange` > 0 & bulk_melanoma_deseq2_log2FoldChange > 0) |
      (`Discovery: log2FoldChange` < 0 & bulk_melanoma_deseq2_log2FoldChange < 0)
  )

top_genes <- unique(top_genes_df$Genes)



######MODEL TRAINING (LOGISTIC REGRESSION-CROSSVALIDATION) (MELANOMA BULK RNA-SEQ)

#prepare caret input
df_model_input <- as.data.frame(data_norm_t)
df_model_input$Unique.Sample.ID <- rownames(df_model_input)
df_model_input <- inner_join(df_model_input, metadata)
rownames(df_model_input) <- df_model_input$Unique.Sample.ID
df_model_input <- df_model_input %>%
  dplyr::select(c(any_of(top_genes), Response, Unique.Sample.ID))
df_model_input$Response <- factor(
  df_model_input$Response,
  levels = c("NR", "R")
)
levels(df_model_input$Response)


#training of the model using cross-validation
k <- 10
ctrl <- trainControl(
  method = "cv",
  number = k,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

set.seed(123)
fit <- train(
  x = df_model_input %>% dplyr::select(-Response, -Unique.Sample.ID),
  y = df_model_input$Response,
  method = "glm",
  family = binomial(),
  trControl = ctrl,
  metric = "ROC"
)

#score predictions in total cross-validation
pred_prob_cv <- fit$pred
pred_prob_cv$Unique.Sample.ID <- df_model_input$Unique.Sample.ID[pred_prob_cv$rowIndex]

#model performance in total cross-validation
roc_obj_cv <- roc(
  response  = pred_prob_cv$obs,
  predictor = pred_prob_cv$R,
  levels    = c("NR", "R")
)
auc(roc_obj_cv)
ci.auc(roc_obj_cv)
scores_responders_cv <- pred_prob_cv%>%
  filter(obs == "R")
scores_non_responders_cv <- pred_prob_cv%>%
  filter(obs == "NR")
plot(
  roc_obj_cv,
  col = "blue",
  lwd = 3
)
pred_prob_cv$obs_01 <- ifelse(pred_prob_cv$obs == "R", 1, 0)
roc.area(pred_prob_cv$obs_01, pred_prob_cv$R)
#write_xlsx(pred_prob_cv, "probscores_crossvalidation_v2.xlsx")

#default threshold = 0.5 in cross-validation
pred_prob_cv$pred_class_0.5 <- ifelse(pred_prob_cv$R >= 0.5, "R", "NR")
pred_prob_cv$pred_class_0.5 <- factor(
  pred_prob_cv$pred_class_0.5,
  levels = c("NR", "R")
)
cm_cv_thres_0.5 <- confusionMatrix(
  data = pred_prob_cv$pred_class_0.5,
  reference = pred_prob_cv$obs,
  positive = "R"
)
print(cm_cv_thres_0.5)

#threshold for 90% sensitivity in cross-validation
roc_df_cv <- data.frame(
  threshold   = roc_obj_cv$thresholds,
  sensitivity = roc_obj_cv$sensitivities,
  specificity = roc_obj_cv$specificities
)
#after evaluation of roc_df_cv: the threshold closer to 0.9 sensitivity at cv is  0.190189029
pred_prob_cv$pred_class_cv_sens_90 <- ifelse(pred_prob_cv$R >= 0.190189029, "R", "NR")
pred_prob_cv$pred_class_cv_sens_90 <- factor(
  pred_prob_cv$pred_class_cv_sens_90,
  levels = c("NR", "R")
)
cm_cv_thres_cv_sesns_90 <- confusionMatrix(
  data = pred_prob_cv$pred_class_cv_sens_90,
  reference = pred_prob_cv$obs,
  positive = "R"
)
print(cm_cv_thres_cv_sesns_90)



#####MODEL VALIDATION (BLADDER CANCER BULK RNA-SEQ Mariathasan et al.)

#bladder cancer bulk RNA metadata
mariathasan_metadata <- read_excel("C:\\Users\\Mariathasan_pheno.xlsx")
mariathasan_metadata <- mariathasan_metadata%>%
  filter(!is.na(binaryResponse))
mariathasan_metadata <- mariathasan_metadata%>%  
  mutate(
    Response = case_when(
      binaryResponse == "CR/PR" ~ "R",
      binaryResponse == "SD/PD" ~ "NR",
    ))
mariathasan_metadata <- mariathasan_metadata%>%
  select(SampleID, Response)
mariathasan_metadata$Response <- factor(mariathasan_metadata$Response, levels = c("NR", "R"))
levels(mariathasan_metadata$Response)
row.names(mariathasan_metadata) <- mariathasan_metadata$SampleID
mariathasan_R <- mariathasan_metadata%>%
  filter(Response == "R")
mariathasan_NR <- mariathasan_metadata%>%
  filter(Response == "NR")

#bladder cancer bulk RNA data
mariathasan_data <- read_excel("C:\\Users\\Mariathasan\\clean\\Mariathasan_cleaned.xlsx")
mariathasan_data <- as.data.frame(mariathasan_data) %>%
  dplyr::select(all_of(c("Symbol", mariathasan_metadata$SampleID)))
row.names(mariathasan_data) <- mariathasan_data$Symbol
mariathasan_data$Symbol <- NULL

#align metadata order with data
mariathasan_metadata <- mariathasan_metadata[colnames(mariathasan_data), ]
row.names(mariathasan_metadata) <- mariathasan_metadata$SampleID

#deseq2 normalisation
stopifnot(identical(colnames(mariathasan_data), rownames(mariathasan_metadata)))
mariathasan_dds <- DESeqDataSetFromMatrix(
  countData = mariathasan_data,
  colData   = mariathasan_metadata,
  design = ~ 1
)
mariathasan_dds <- DESeq(mariathasan_dds)
mariathasan_data_norm <- counts(mariathasan_dds, normalized = TRUE)

#compare median gene expression in melanoma vs bladder dataset (for common genes)
row_medians_melanoma <- as.data.frame(apply(data_norm, 1, median, na.rm = TRUE))
colnames(row_medians_melanoma) <- "melanoma"
row_medians_melanoma$gene <- row.names(row_medians_melanoma)

row_medians_bladder <- as.data.frame(apply(mariathasan_data_norm, 1, median, na.rm = TRUE))
colnames(row_medians_bladder) <- "bladder"
row_medians_bladder$gene <- row.names(row_medians_bladder)

merged <- inner_join(row_medians_melanoma, row_medians_bladder)
melanoma_median <- median(merged$melanoma)
bladder_median <- median(merged$bladder)
norm_factor <- round(bladder_median/melanoma_median, 1)

#plot median gene expression distribution of melanoma and bladder cancer datasets 
#before normalisation of bladder cancer dataset by dividing it with norm`factor`
breaks <- seq(min(c(log2(merged$melanoma +1),log2(merged$bladder +1))), max(c(log2(merged$melanoma +1),log2(merged$bladder +1))), length.out = 30)
hist(log2(merged$melanoma +1),
     breaks = breaks,
     col = rgb(0,0,1,0.5),
     xlim = range(breaks),
     main = "Overlayed Histograms",
     xlab = "log2(gene_medians)")
hist(log2(merged$bladder +1),
     breaks = breaks,
     col = rgb(1,0,0,0.5),
     add = TRUE)
legend("topright",
       legend = c("melanoma","bladder"),
       fill = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)),
       cex=0.5)

#normalisation of bladder dataset based on the melanoma dataset
mariathasan_data_norm <- mariathasan_data_norm/norm_factor

#plot median gene expression distribution of melanoma and bladder cancer datasets 
#after normalisation of bladder cancer dataset by dividing it with norm`factor`
row_medians_bladder <- as.data.frame(apply(mariathasan_data_norm, 1, median, na.rm = TRUE))
colnames(row_medians_bladder) <- "bladder"
row_medians_bladder$gene <- row.names(row_medians_bladder)
merged <- inner_join(row_medians_melanoma, row_medians_bladder)
breaks <- seq(min(c(log2(merged$melanoma +1),log2(merged$bladder +1))), max(c(log2(merged$melanoma +1),log2(merged$bladder +1))), length.out = 30)
hist(log2(merged$melanoma +1),
     breaks = breaks,
     col = rgb(0,0,1,0.5),
     xlim = range(breaks),
     main = "Overlayed Histograms",
     xlab = "log2(gene_medians)")
hist(log2(merged$bladder +1),
     breaks = breaks,
     col = rgb(1,0,0,0.5),
     add = TRUE)
legend("topright",
       legend = c("melanoma","bladder"),
       fill = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)),
       cex=0.5)

#prepare as model validation cohort
df_model_input_val <- as.data.frame(t(mariathasan_data_norm))
df_model_input_val$SampleID <- row.names(df_model_input_val)
df_model_input_val <- inner_join(df_model_input_val, mariathasan_metadata, by = "SampleID") %>%
  dplyr::select(c(any_of(top_genes), "Response", "SampleID"))
df_model_input_val$Response <- factor(df_model_input_val$Response, levels = c("NR","R"))
rownames(df_model_input_val) <- df_model_input_val$SampleID

#validate model in bladder cancer RNA-seq dataset and add response info to the model predictions
pred_prob_val <- predict(fit, newdata = df_model_input_val, type = "prob")
pred_prob_val$SampleID <- rownames(pred_prob_val)
pred_prob_val <- inner_join(mariathasan_metadata, pred_prob_val)

#model performance in validation
roc_obj_val <- roc(
  response  = pred_prob_val$Response,
  predictor = pred_prob_val$R,
  levels    = c("NR", "R")
)
roc_obj_val$direction
auc(roc_obj_val)
ci.auc(roc_obj_val)
scores_responders_val <- pred_prob_val%>%
  filter(Response == "R")
scores_non_responders_val <- pred_prob_val%>%
  filter(Response == "NR")
plot(
  roc_obj_val,
  col = "blue",
  lwd = 3
)
pred_prob_val$obs_01 <- ifelse(pred_prob_val$Response == "R", 1, 0)
roc.area(pred_prob_val$obs_01, pred_prob_val$R)

#apply threshold predifined for 90% sensitivity in cross-validation
pred_prob_val$pred_class_cv_sens_90 <- ifelse(pred_prob_val$R >= 0.190189029, "R", "NR")
pred_prob_val$pred_class_cv_sens_90 <- factor(
  pred_prob_val$pred_class_cv_sens_90,
  levels = c("NR", "R")
)
cm_val_thres_cv_sens_90 <- confusionMatrix(
  data = pred_prob_val$pred_class_cv_sens_90,
  reference = pred_prob_val$Response,
  positive = "R"
)
print(cm_val_thres_cv_sens_90)



###model performance in low-TMB subgroup
#bladder cancer bulk RNA TMB metadata
mariathasan_metadata <- read_excel("C:\\Users\\Mariathasan\\Mariathasan_pheno.xlsx")
mariathasan_metadata <- mariathasan_metadata%>%
  select(SampleID, FMOne.mutation.burden.per.MB)
#merge TMB metadata with model scores
df <- inner_join(mariathasan_metadata, pred_prob_val, by="SampleID")
df <- df%>%
  filter(!is.na(FMOne.mutation.burden.per.MB))
#select low TMB-subgroup
df <- df%>%
 filter(FMOne.mutation.burden.per.MB <10)
#metrics of performance in the low-TMB subgroup
df$Response_01 <- ifelse(df$Response == "R", 1, 0)
roc_score <- roc(df$Response_01, df$R)
plot(roc_score, col="blue")
auc(roc_score)
ci.auc(roc_score)
roc.area(df$Response_01, df$R)
confusionMatrix(
  factor(df$pred_class_cv_sens_90),
  factor(df$Response),
  positive = "R"
)

###Correlation analysis 


# Target list (13BM)
genes13BM <- c("LGALS1", "UHRF2", "LINGO1", "SELL", "CD96", "EPB41", 
               "MPRIP", "IKZF1", "TLR10", "ST6GAL1", "ALDH5A1", "PDGFRB", "PLEC")



# Define the weights from your model output
weights <- c(
  CD96 = 5.507e-04, SELL = 1.489e-04, EPB41 = -4.547e-06, 
  LINGO1 = -2.103e-03, UHRF2 = -2.076e-03, LGALS1 = 2.069e-05, 
  PDGFRB = -1.243e-03, PLEC = -4.671e-05, ST6GAL1 = 3.589e-04, 
  TLR10 = -1.511e-04, IKZF1 = -1.031e-03, MPRIP = -1.223e-03, 
  ALDH5A1 = 3.846e-03
)
intercept <- 2.086e+00

# 1. FOR MELANOMA: Subset data_norm to just these genes and transpose
signature_genes <- intersect(names(weights), rownames(data_norm))
exp_matrix <- t(data_norm[signature_genes, setdiff(colnames(data_norm), "Genes")])

# 2. Calculate the score (Weighted Sum + Intercept)

sample_scores <- (exp_matrix %*% weights[signature_genes]) + intercept

# Prepare the full data for correlation (Transpose so genes are columns)
all_genes_mat <- t(data_norm[, setdiff(colnames(data_norm), "Genes")])

# Calculate correlation for all genes against the score


# Calculate correlations and p-values in one go
ctest <- corr.test(all_genes_mat, sample_scores, method = "spearman", adjust = "none")

cor_results_MEL <- data.frame(
  Gene = rownames(ctest$r),
  Correlation = as.numeric(ctest$r),
  p_value = as.numeric(ctest$p)
) %>%
  mutate(padj = p.adjust(p_value, method = "fdr")) %>%
  arrange(desc(Correlation))


# 1. FOR BLADDER Subset mariathasan_data_norm  to just these genes and transpose
signature_genes <- intersect(names(weights), rownames(mariathasan_data_norm ))
exp_matrix <- t(mariathasan_data_norm [signature_genes, setdiff(colnames(mariathasan_data_norm ), "Genes")])

# 2. Calculate the score (Weighted Sum + Intercept)
# This uses matrix multiplication (%*%) for speed
sample_scores <- (exp_matrix %*% weights[signature_genes]) + intercept

# Now sample_scores is a vector where each value is the "Model Score" for that sample

# Prepare the full data for correlation (Transpose so genes are columns)
all_genes_mat <- t(mariathasan_data_norm[, setdiff(colnames(mariathasan_data_norm ), "Genes")])

# Calculate correlation for all genes against the score
# This returns a matrix of R values
# Calculate correlations and p-values in one go
ctest <- corr.test(all_genes_mat, sample_scores, method = "spearman", adjust = "none")

cor_results_BLAD <- data.frame(
  Gene = rownames(ctest$r),
  Correlation = as.numeric(ctest$r),
  p_value = as.numeric(ctest$p)
) %>%
  mutate(padj = p.adjust(p_value, method = "fdr")) %>%
  arrange(desc(Correlation))

 

genes_posbl <- cor_results_BLAD %>%
  filter(Correlation > 0.3) %>%
  pull(Gene)

genes_negmel <- cor_results_MEL %>%
  filter(Correlation < -0.3) %>%
  pull(Gene)

genes_posmel <- cor_results_MEL %>%
  filter(Correlation > 0.3) %>%
  pull(Gene)

genes_negbl <- cor_results_BLAD %>%
  filter(Correlation < -0.3) %>%
  pull(Gene)

library(clusterProfiler)
library(org.Hs.eg.db)

pos_entrezbl <- bitr(genes_posbl,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)


neg_entrezmel <- bitr(genes_negmel,
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)

pos_entrezmel <- bitr(genes_posmel,
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)


neg_entrezbl <- bitr(genes_negbl,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)

ora_pos_blad <- enrichGO(
  gene = pos_entrezbl$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

ora_neg_blad <- enrichGO(
  gene = neg_entrezbl$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

ora_pos_mel <- enrichGO(
  gene = pos_entrezmel$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

ora_neg_mel <- enrichGO(
  gene = neg_entrezmel$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

library(writexl)

sheets <- list(
  "Bladder_Positive" = as.data.frame(ora_pos_blad),
  "Bladder_Negative" = as.data.frame(ora_neg_blad),
  "Melanoma_Positive" = as.data.frame(ora_pos_mel),
  "Melanoma_Negative" = as.data.frame(ora_neg_mel)
)

write_xlsx(sheets, "supplementary_table_6.xlsx")
 
