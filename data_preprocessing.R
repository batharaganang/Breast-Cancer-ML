########## Filtering for breast cancer samples

# Load target data
target_data <- read.table('TCGA_BRCA_Target.txt', sep='\t', header=TRUE)

# Converting target data to csv for future use
write.csv(target_data, "TCGA_BRCA_Target.csv", row.names=FALSE)

# Loading attribute data
data <- read.table('TCGA_BRCA_Attribute.txt', sep='\t', header=TRUE)

# Replacing "-" with "." in target data sampleID column
target_var <- gsub("-", ".", target_data$sampleID)

# Filtering breast cancer samples from the TPM normalized data
selected_data <- data[, colnames(data) %in% target_var]

# Saving the filtered result to csv file
write.csv(selected_data, "TCGA_BRCA_attribute_filtered.csv", row.names=FALSE)



########## Transcript to Genes

# Loading filtered attribute data
data <- read.csv('TCGA_BRCA_attribute_filtered.csv', header=TRUE)

# Converting log2 values back to its normal value
library(dplyr)
data_new <- data %>%
  mutate(across(-sample, ~ 2^(.) - 0.001))

# Changing minus values to 0 because transcripts cannot have negative value
data_new <- data_new %>%
  mutate(across(-sample, ~ ifelse(. < 0, 0, .)))

# Sanity Check
# col_sums <- colSums(data_new[, !colnames(data_new) %in% "sample"])

# Remove everything after the period in the ID column for convertion
data_new$sample <- sub("\\..*$", "", data_new$sample)

# Convert transcript ID to gene ID using biomart
library(biomaRt)
gene_ids <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters = "ensembl_gene_id",
  values = data_new$sample,
  mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
)

# Perform inner join on data_new and gene_ids
result <- merge(data_new, gene_ids, by.x = "sample", by.y = "ensembl_gene_id")

# Get rid of sample columns because we already have external_gene_name
result <- result[, -which(names(result) == "sample")]

# Get rid of empty external gene
result <- result[result$external_gene_name != "", ]

# Sum data with external_gene_name as the groupby
result_aggr <- aggregate(. ~ external_gene_name, data = result, FUN = sum)

row_sums <- rowSums(result_aggr[, !colnames(result_aggr) %in% "external_gene_name"])
#zero_count <- sum(row_sums == 0)
#print(zero_count)

# Filter data that have 0 sum of genes across sample (no such genes in breast cancer)
data_new_filtered <- subset(result_aggr, row_sums != 0)

# Save Data
write.csv(data_new_filtered, "TCGA_BRCA_attribute_genes.csv", row.names=FALSE)



########## Stratified Split and applying VST on train data

# Load Data
data <- read.csv('TCGA_BRCA_attribute_genes.csv', header=TRUE)
target_data <- read.table('TCGA_BRCA_Target.txt', sep='\t', header=TRUE)

# Setting up data for merging
target_data$sampleID <- gsub("-", ".", target_data$sampleID)
target_data <- target_data[,c("sampleID","PAM50Call_RNAseq")]
rownames(data) <- data$external_gene_name
data <- data[, -1]
data <- t(data)
data <- as.data.frame(data)

# Get rid of samples with '11' on the data
data <- data[!grepl("11", row.names(data)), ]

# Merge target and attribute data
data$sampleID <- row.names(data)
merged_data <- merge(data, target_data)

# Filter out where pam50 target is null
filtered_data <- merged_data[complete.cases(merged_data$PAM50Call_RNAseq), ]
filtered_data <- subset(merged_data, PAM50Call_RNAseq != "")
rownames(filtered_data) <- NULL  

library(caret)
# Perform stratified split
train_indices <- createDataPartition(filtered_data$PAM50Call_RNAseq, p = 0.8, list = FALSE)
# Create the training and testing datasets
train_data <- filtered_data[train_indices, ]
test_data <- filtered_data[-train_indices, ]

# Save test data
write.csv(test_data, "TCGA_BRCA_Test.csv", row.names=FALSE)

# Applying VST to train data
# Transforming train data back to original
rownames(train_data) <- train_data$sampleID
train_data <- train_data[, -1]
pam50 <- train_data$PAM50Call_RNAseq
train_data <- train_data[, -which(names(train_data) == "PAM50Call_RNAseq")]
train_data <- t(train_data)
train_data <- as.data.frame(train_data)

# Implement VST with seurat package
library(Seurat)

result <- FindVariableFeatures(
  train_data,
  selection.method = "vst",
  loess.span = 0.3,
  clip.max = "auto",
  mean.function = FastExpMean,
  dispersion.function = FastLogVMR,
  num.bin = 20,
  binning.method = "equal_width",
  nfeatures = 2000,
  verbose = TRUE)

# Get VST score for train data
train_data$vst_score <- result$vst.variance.standardized

# Sort train data by VST score
sorted_df <- train_data[order(-train_data$vst_score), ]
sorted_df$external_gene_name <- row.names(sorted_df)
rownames(sorted_df) <- NULL  
# sorted_df_test <- sorted_df[,c("external_gene_name","vst_score")]

# Save train data
write.csv(sorted_df, "TCGA_BRCA_Train.csv", row.names=FALSE)