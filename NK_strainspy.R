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
barcode_proportions$barcode <- row.names(barcode_proportions)

# Check if each row is a unique barcode
# barcode_proportions |>
#   group_by(barcode) |>
#   summarise(n = n()) |>
#   summarise(max = max(n))

# Access the sample names
sample <- colnames(barcode_proportions)
# Remove barcode
sample <- sample[-length(sample)]

# Create row_indices (each row is a unique barcode)
barcode_proportions <- tibble::rowid_to_column(barcode_proportions, "row_indices")
# Reformat the data.frame such that each row is a unique proportion
barcode_proportions  <- barcode_proportions |>
  select(all_of(sample),row_indices,barcode) |>
  pivot_longer(cols = all_of(sample), names_to = 'Sample', values_to = 'Proportion')
# Create col_indices according to grouping
barcode_proportions <- barcode_proportions |>
  group_by(Sample) |>
  mutate(col_indices = cur_group_id())
barcode_proportions$Sample <- gsub('Gr_.','Gr',barcode_proportions $Sample)

# Access the metadata from monkeyHSPC
col_data <- monkeyHSPC@colData@listData[["sampleMetadata"]]
col_data <- col_data[col_data@rownames %in% sample,]
col_data$Celltype <- gsub('Gr_.','Gr',col_data$Celltype)
col_data@rownames <- gsub('Gr_.','Gr',col_data@rownames)

col_data_2_groups <- col_data
col_data_2_groups$Celltype <- ifelse(col_data$Celltype == 'NK_CD56n_CD16p',col_data$Celltype,'Other')

n_rows <- max(barcode_proportions$row_indices)
n_cols <- max(barcode_proportions$col_indices)

# Creating a proportions summarized experiment that is asin transformed
se_HSPC_asin <- SummarizedExperiment::SummarizedExperiment(
  assays = list(Matrix::sparseMatrix(
    i = barcode_proportions$row_indices,
    j = barcode_proportions$col_indices,
    x = asin(sqrt(barcode_proportions $Proportion)),
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
    x = asin(sqrt(barcode_proportions $Proportion)),
    dims = c(n_rows, n_cols),
    repr = "R"
  )),
  colData = col_data_2_groups
)

############################## Fitting Models #######################################

# Standard design across models
design <- as.formula("~ Celltype + Months")

# Zero-inflated beta fit 
zib_2grp <- strainspy::glmZiBFit(se = se_HSPC_2grp, design = design)

# Ordinal beta fit
beta_2grp <- strainspy::glmFit(se = se_HSPC_2grp, design = design)


# Simulations via implantation
x <- as.vector(SummarizedExperiment::assay(se_HSPC_asin))
rescaled <- strainspy::rescale_beta(x)
strainspy::getZICoefficients(fit_asin)

# Comparison via AUROC (similar to Wirbel et al.)

