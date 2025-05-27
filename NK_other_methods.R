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

######################### General Linear Model #################################

# Linear model fit using lm for each barcode

glm_dat <- as.matrix(simulation_se@assays@data@listData[[1]])
grp <- simulation_se$sim
num_barcodes <- dim(simulation_se)[1]

p.val_lm <- c()

for(i in c(1:num_barcodes)) {
  b_i <- glm_dat[i,]
  b_i <- data.frame(proportion = asin(sqrt(b_i)),
                      sim = grp)
  lm_sim <- lm(proportion ~ sim, data = b_i)
  aov_sim <- anova(lm_sim)
  p.val_lm[i] <- aov_sim$`Pr(>F)`
}

sim_res <- rbind(sim_res, list(j, 'LM', get_AUC(p.val_lm, true_signal), 0))

############################## Wilcoxon & T-test #########################################

# Used base R wilcox.test & t.test on each barcode

wil_dat <- as.matrix(simulation_se@assays@data@listData[[1]])
grp <- simulation_se$sim == 'a'
num_barcodes <- dim(simulation_se)[1]

p.val_wil <- c()
p.val_t <- c()

for(i in c(1:num_barcodes)) {
  b_i <- wil_dat[i,]
  b_i_grp1 <- asin(sqrt(b_i[grp]))
  b_i_grp2 <- asin(sqrt(b_i[!grp]))
  wil_sim <- wilcox.test(b_i_grp1,b_i_grp2,alternative = 'two.sided')
  t_sim <- t.test(b_i_grp1,b_i_grp2,alternative = 'two.sided')
  p.val_wil[i] <- wil_sim$p.value
  p.val_t[i] <- t_sim$p.value
}


sim_res <- rbind(sim_res, list(j, 'wilcoxon', get_AUC(p.val_wil, true_signal), 0))
sim_res <- rbind(sim_res, list(j, 't-test', get_AUC(p.val_t, true_signal), 0))


