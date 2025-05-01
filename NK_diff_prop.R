BiocManager::install("barbieQ")
data("monkeyHSPC", package = 'barbieQ')
monkeyHSPC_clust <- barbieQ::clusterCorrelatingBarcodes(monkeyHSPC)

# Creating a SummarizedExperiment for strainspy usage
barcode_proportions <- as.data.frame(monkeyHSPC@assays@data@listData[["proportion"]])
barcode_proportions$barcode <- row.names(barcode_proportions)
sample <- colnames(barcode_proportions)
sample <- sample[-length(sample)]

barcode_proportions <- tibble::rowid_to_column(barcode_proportions, "row_indices")
barcode_proportions <- barcode_proportions |>
  pivot_longer(cols = all_of(sample), names_to = 'Sample', values_to = 'Proportion') |>
  group_by(Sample) |>
  mutate(col_indices = cur_group_id())

col_data <- monkeyHSPC@colData@listData[["sampleMetadata"]]
n_rows <- max(barcode_proportions$row_indices)
n_cols <- max(barcode_proportions$col_indices)



se_HSPC <- SummarizedExperiment::SummarizedExperiment(
  assays = list(Matrix::sparseMatrix(
    i = barcode_proportions$row_indices,
    j = barcode_proportions$col_indices,
    x = barcode_proportions$Proportion,
    dims = c(n_rows, n_cols),
    repr = "R"
  )),
  colData = col_data
)

se_HSPC_t <- strainspy::filter_by_presence(se_HSPC)

design <- as.formula(" ~ Celltype")

fit_p <- glmZiBFit(se = se_HSPC_t, design = design)

strainspy::summary(fit_p)
strainspy::top_hits(fit_p)

