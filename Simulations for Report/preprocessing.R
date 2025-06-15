BiocManager::install("barbieQ")
library(tidyverse)
data("monkeyHSPC", package = 'barbieQ')

# Formatting monkeyHSPC barbieQ object into a SummarizedExperiment object formatted for strainspy use
# Code inspired by read_sylph function from strainspy

###### Pre-processing inspired by Liyang's case study ########################

# Limit to months between 6 and 55
months <- monkeyHSPC$sampleMetadata$Months
monkeyHSPC <- monkeyHSPC[,6 < months & months < 55]
# Access barcodes that represent majority of abundance 
# Also must be present in at least 6 biological replicates
monkeyHSPC <- barbieQ::tagTopBarcodes(barbieQ = monkeyHSPC, nSampleThreshold = 6)

# Stacked barplot representing proportion of barcodes filtered out
barbieQ::plotBarcodeSankey(barbieQ = monkeyHSPC)
# Filtering out barcodes which don't contribute heavily
monkeyHSPC <- monkeyHSPC[SummarizedExperiment::rowData(monkeyHSPC)$isTopBarcode$isTop,]

############ Creating a Summarized Experiment formatted for StrainSpy ############

# Grab proportions of barcodes 
barcode_proportions <- as.data.frame(monkeyHSPC@assays@data@listData[["proportion"]])
temp <- barcode_proportions

# Check if each row is a unique barcode
# barcode_proportions |>
#   group_by(barcode) |>
#   summarise(n = n()) |>
#   summarise(max = max(n))

# Access the sample names
colnames(barcode_proportions) <- c(1:ncol(barcode_proportions))
barcode_proportions$barcode <- row.names(barcode_proportions)
# Remove barcode
sample <- as.character(c(1:(ncol(barcode_proportions)-1)))

# Create row_indices (each row is a unique barcode)
barcode_proportions <- tibble::rowid_to_column(barcode_proportions, "row_indices")
# Reformat the data.frame such that each row is a unique proportion
barcode_proportions  <- barcode_proportions |>
  pivot_longer(cols = all_of(sample), names_to = 'Sample', values_to = 'Proportion') |>
  mutate(col_indices = as.integer(Sample)) |>
  select(-Sample)
# barcode_proportions$Sample <- gsub('Gr_.','Gr',barcode_proportions $Sample)

# Access the metadata from monkeyHSPC
col_data <- monkeyHSPC@colData@listData[["sampleMetadata"]]
#col_data <- col_data[col_data@rownames %in% sample,]
col_data$Celltype <- gsub('Gr_.','Gr',col_data$Celltype)
col_data@rownames <- gsub('Gr_.','Gr',col_data@rownames)

col_data_2_groups <- col_data
col_data_2_groups$feature <- ifelse(col_data$Celltype == 'NK_CD56n_CD16p',col_data$Celltype,'Other')

n_rows <- max(barcode_proportions$row_indices)
n_cols <- max(barcode_proportions$col_indices)

# Creating a proportions summarized experiment that is asin transformed
se_HSPC <- SummarizedExperiment::SummarizedExperiment(
  assays = list(Matrix::sparseMatrix(
    i = barcode_proportions$row_indices,
    j = barcode_proportions$col_indices,
    x = barcode_proportions $Proportion,
    dims = c(n_rows, n_cols),
    repr = "R"
  )),
  colData = col_data
)

# A 2 group proportions SE
se_HSPC_2grp <- SummarizedExperiment::SummarizedExperiment(
  assays = list(Matrix::sparseMatrix(
    i = barcode_proportions$row_indices,
    j = barcode_proportions$col_indices,
    x = barcode_proportions$Proportion,
    dims = c(n_rows, n_cols),
    repr = "R"
  )),
  colData = col_data_2_groups
)

############################## Real Data #######################################

# Standard design across models
m_design <- as.formula("~ feature + Months")
#nom_design <- as.formula("~ Celltype")

# Zero-inflated beta fit 
zib_2grp <- strainspy::glmZiBFit(se = se_HSPC_2grp, design = design)

barcodes <- row.names(temp)

# Grab the minimum of the two p-values for each barcode
p.val.zib <- c()
beta.p <- beta_2grp@p_values$featureOther
zero.p <- beta_2grp@zi_p_values$featureOther
for(i in c(1:dim(se_HSPC_2grp)[[1]])) {
  p.val.zib[i] <- min(beta.p[i],zero.p[i])
}

# BH correction for the p-values
p.val.zib <- p.adjust(p.val.zib, method = 'BH')
p.val.zib[is.na(p.val.zib)] <- 1

# Find significant barcodes for ZiB
zib.sig <- barcodes[p.val.zib < 0.05]

# limma
design_limma <- model.matrix(~ se_HSPC_2grp$feature + se_HSPC_2grp$Months)
limma_dat <- as.matrix(se_HSPC_2grp@assays@data@listData[[1]])
limma_dat <- asin(sqrt(limma_dat))

limma_fit <- limma::lmFit(limma_dat, design_limma)
limma_fit <- limma::eBayes(limma_fit)

# Access p-values and apply BH correction
p.val_limma <- limma_fit$p.value[,2]
p.val_limma <- p.adjust(p.val_limma, method = 'BH')
limma.sig <- barcodes[p.val_limma < 0.05]

# Create venn diagram
VennDiagram::venn.diagram(
  x = list(limma.sig,zib.sig),
  category.names = c("limma" , "ZiB"),
  filename = 'venn.png',
  output=TRUE,
  fill = c("#7570b3", "#1b9e77"),  # Green and orange
  alpha = 0.6,
  cex = 1.25,
)

