############################## barbieQ #######################################

monkeyHSPC$sampleMetadata$Celltype <- ifelse(monkeyHSPC$sampleMetadata$Celltype == 'NK_CD56n_CD16p',
                                             monkeyHSPC$sampleMetadata$Celltype,
                                             'Other')

# Remove the monkey column
monkeyHSPC$sampleMetadata <- monkeyHSPC$sampleMetadata[,-1]
monkeyHSPC$sampleMetadata <- monkeyHSPC$sampleMetadata[,c(2,1)]

des_mat <- model.matrix(~ Celltype + Months, data = monkeyHSPC$sampleMetadata)

barbieQ::testBarcodeSignif(barbieQ = monkeyHSPC,
                           designMatrix = des_mat,
                           method = 'diffProp',
                           transformation = 'asin-sqrt')

# Having an issue with the "designFormula" part... seems to default to a contrast when I don't want it to

############################## fastANCOM #######################################

# Access proportion data
fancom_data <- as.matrix(monkeyHSPC@assays@data$proportion)
# 2 groups: NK_CD56n_CD16p and Other
fancom_grp <- col_data_2_groups$Celltype
# Months as additional variable (covariate)
fancom_conf <- col_data_2_groups$Months
# fastANCOM formats data with samples as rows and
# microbes (or in this case barcodes) as columns; so transpose is needed
fancom_data <- t(fancom_data)
fastANCOM_fit <- fastANCOM::fastANCOM(Y = fancom_data, x = fancom_grp, z = fancom_conf)$results$final

######################### General Linear Model #################################

glm_dat <- as.matrix(monkeyHSPC@assays@data$proportion)
glm_dat <- as.data.frame(t(glm_dat))
barcodes <- colnames(glm_dat)
glm_dat$Celltype <- monkeyHSPC$sampleMetadata$Celltype
glm_dat$Months <- monkeyHSPC$sampleMetadata$Months

glm_dat <- glm_dat |>
  pivot_longer(cols = all_of(barcodes), names_to = 'barcode', values_to = 'proportion')

lmerTest::lmer(proportion ~ Celltype + (1|Months), data = )

############################## Wilcoxon #########################################

wilcoxon_dat <- glm_dat
wilcoxon_fit <- coin::wilcox_test(proportion ~ as.factor(Celltype)|as.factor(Months), data = wilcoxon_dat)