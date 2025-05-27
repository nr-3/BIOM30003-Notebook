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
col_data_2_groups$Celltype <- ifelse(col_data$Celltype == 'NK_CD56n_CD16p',col_data$Celltype,'Other')

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

############################## Fitting Models #######################################

# Standard design across models
m_design <- as.formula("~ Celltype + Months")
nom_design <- as.formula("~ Celltype")

# Zero-inflated beta fit 
zib_2grp <- strainspy::glmZiBFit(se = se_HSPC_2grp, design = design)

# Ordinal beta fit
beta_2grp <- strainspy::glmFit(se = se_HSPC_2grp, design = design)

nrow(temp)

######################## Simulation via implantation ################################

# Data.frame to hold the AUROC + FDR values for each simulation
sim_res <- data.frame()

# Sample size
n <- n_cols

# Which repetition
j <- 1

# Implantation parameters
beta_diff <- runif(1,0.5,0.9)
zero_effect <- 0.05

# Split samples into two random groups
sample <- c(1:n)
grp1 <- sort(sample(sample, n/2))
grp2 <- setdiff(sample, grp1)
simulation_se <- se_HSPC_2grp
simulation_se$sim <- NA
simulation_se[,grp1]$sim <- 'a'
simulation_se[,grp2]$sim <- 'b'

# 10% of the barcodes to be differentially abundant
implant <- round(0.1*n_rows)
implant <- sample(c(1:n_rows),implant)

true_signal <- rep(0, n_rows)
true_signal[implant] <- 1
grp1 <- simulation_se$sim=='a'

# Implant the signal in the barcodes randomly selected
for (i in implant) {
  s <- as.vector(simulation_se@assays@data@listData[[1]][i,])
  s <- rescale_beta(s[grp1], beta = beta_diff, zi = zero_effect)
  simulation_se@assays@data@listData[[1]][i,grp1] <- s$rescaled
}

sim_zib <- strainspy::glmZiBFit(se = simulation_se, design = as.formula('~sim'))
sim_ordb <- strainspy::glmFit(se = simulation_se, design = as.formula('~sim'))

# AUC Calculations (similar to Wirbel et al.'s SIMBA package)
p.val.zib <- sim_zib@zi_p_values$simb
p.val.ord <- sim_ordb@p_values$simb

auroc <- pROC::roc(predictor = -log10(p.val.ord + 1e-50),
             response = true_signal, levels = c(0,1), direction = '<')


sim_res <- rbind(sim_res, list(j, 'ordinal', get_AUC(p.val.ord, true_signal), 0))
sim_res <- rbind(sim_res, list(j, 'zib', get_AUC(p.val.zib, true_signal), 0))

ggplot(sim_res, aes(x = model, y = AUC)) +
  geom_col()
