BiocManager::install("barbieQ")
data("monkeyHSPC", package = 'barbieQ')
monkeyHSPC_clust <- barbieQ::clusterCorrelatingBarcodes(monkeyHSPC)

# Formatting monkeyHSPC barbieQ object into a SummarizedExperiment object formatted for strainspy use
# Code inspired by read_sylph function from strainspy

# Creating a data.frame with proportion data
barcode_proportions <- as.data.frame(monkeyHSPC@assays@data@listData[["proportion"]])
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
# Limit sample to just these months (featured from Liyang's case study)
months <- c('6.5m','12m','17m','27m','36m','42m')
sample <- sample[grepl(paste(months,collapse='|'), sample)]
m42_sample <- sample[grepl('42m', sample)]
sample <- sample[-c(length(sample),length(sample)-1)]

barcode_proportions <- tibble::rowid_to_column(barcode_proportions, "row_indices")
barcode_proportions_all <- barcode_proportions |>
  select(all_of(sample),row_indices,barcode) |>
  pivot_longer(cols = all_of(sample), names_to = 'Sample', values_to = 'Proportion')
barcode_proportions_42m <- barcode_proportions |>
  select(all_of(m42_sample),row_indices,barcode) |>
  pivot_longer(cols = all_of(m42_sample), names_to = 'Sample', values_to = 'Proportion')
barcode_proportions_all$Sample <- gsub('Gr_1','Gr',barcode_proportions_all$Sample)
barcode_proportions_all <- barcode_proportions_all |>
  group_by(Sample) |>
  mutate(col_indices = cur_group_id())
barcode_proportions_42m <- barcode_proportions_42m |>
  group_by(Sample) |>
  mutate(col_indices = cur_group_id())

# Access the metadata from monkeyHSPC
col_data <- monkeyHSPC@colData@listData[["sampleMetadata"]]
col_data <- col_data[col_data@rownames %in% sample,]
col_data$Celltype <- gsub('Gr_1','Gr',col_data$Celltype)
n_rows <- max(barcode_proportions_all$row_indices)
n_cols <- max(barcode_proportions_all$col_indices)

col_data_42m <- col_data[grepl('42m',col_data@rownames),]
n_rows_42 <- max(barcode_proportions_42m$row_indices)
n_cols_42 <- max(barcode_proportions_42m$col_indices)

# Creating a proportions summarized experiment
se_HSPC_n <- SummarizedExperiment::SummarizedExperiment(
  assays = list(Matrix::sparseMatrix(
    i = barcode_proportions_all$row_indices,
    j = barcode_proportions_all$col_indices,
    x = barcode_proportions_all$Proportion,
    dims = c(n_rows, n_cols),
    repr = "R"
  )),
  colData = col_data
)

se_HSPC_asin <- SummarizedExperiment::SummarizedExperiment(
  assays = list(Matrix::sparseMatrix(
    i = barcode_proportions_all$row_indices,
    j = barcode_proportions_all$col_indices,
    x = asin(sqrt(barcode_proportions_all$Proportion)),
    dims = c(n_rows, n_cols),
    repr = "R"
  )),
  colData = col_data
) 

se_HSPC_42m <- SummarizedExperiment::SummarizedExperiment(
  assays = list(Matrix::sparseMatrix(
    i = barcode_proportions_42m$row_indices,
    j = barcode_proportions_42m$col_indices,
    x = barcode_proportions_42m$Proportion,
    dims = c(n_rows_42, n_cols_42),
    repr = "R"
  )),
  colData = col_data_42m
)

se_HSPC_n <- strainspy::filter_by_presence(se_HSPC_n, rescale_abundance = T)
se_HSPC_asin <- strainspy::filter_by_presence(se_HSPC_asin, rescale_abundance = T)
se_HSPC_42m <- strainspy::filter_by_presence(se_HSPC_42m, rescale_abundance = T)

design <- as.formula(" ~ Celltype + Months")

fit_n <- strainspy::glmZiBFit(se = se_HSPC_n, design = design)
fit_asin <- strainspy::glmZiBFit(se = se_HSPC_asin, design = design)
fit_42m <- strainspy::glmZiBFit(se = se_HSPC_42m, design = as.formula(' ~ Celltype'))

strainspy::summary(fit_n)
strainspy::top_hits(fit_n, alpha = 0.1)

# FDR
strainspy::filter_by_presence()

# Benchmarking via simulation

# Volcano plots

# Ion know


