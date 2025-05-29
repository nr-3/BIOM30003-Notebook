######################## Simulation via implantation ################################

# Data.frame to hold the AUROC + FDR values for each simulation
#sim_res_adj <- data.frame()

# Sample size
n <- n_cols
zero_effect <- 0.03

for(beta_diff in c(0.75,0.8,0.85,0.9,0.95)) {
  for(j in c(1:3)) {
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
    
    FDR_ord <- get_FDR(p.adjust(p.val.ord,method='BH'),true_signal)
    FDR_zib <- get_FDR(p.adjust(p.val.zib,method='BH'),true_signal)
    
    sim_res_adj <- rbind(sim_res_adj, list(j, 'ordinal', get_AUC(p.val.ord, true_signal), FDR_ord[1], FDR_ord[2], FDR_ord[3], beta_diff))
    sim_res_adj <- rbind(sim_res_adj, list(j, 'zib', get_AUC(p.val.zib, true_signal), FDR_zib[1], FDR_zib[2], FDR_zib[3], beta_diff))
    
    ############################## limma #######################################
    
    # Add random groups
    #monkeyHSPC$sampleMetadata$sim <- as.factor(simulation_se$sim)
    
    design_limma <- model.matrix(~ simulation_se$sim)
    limma_dat <- as.matrix(simulation_se@assays@data@listData[[1]])
    limma_dat <- asin(sqrt(limma_dat))
    
    limma_fit <- lmFit(limma_dat, design_limma)
    limma_fit <- eBayes(limma_fit)
    
    p.val_limma <- limma_fit$p.value[,2]
    
    adj_limma <- p.adjust(p.val_limma, method = 'BH')
    
    FDR_limma <- get_FDR(adj_limma,true_signal)
    
    sim_res_adj <- rbind(sim_res_adj, list(j, 'limma', get_AUC(p.val_limma, true_signal), FDR_limma[1], FDR_limma[2], FDR_limma[3], beta_diff))
    
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
    
    FDR_glm <- get_FDR(p.adjust(p.val_lm, method='BH'),true_signal)
    
    sim_res_adj <- rbind(sim_res_adj, list(j, 'LM', get_AUC(p.val_lm, true_signal), FDR_glm[1], FDR_glm[2], FDR_glm[3] , beta_diff))
    
    
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
    
    FDR_wil <- get_FDR(p.adjust(p.val_wil, method = 'BH'),true_signal)
    FDR_t <- get_FDR(p.adjust(p.val_t, method = 'BH'),true_signal)
    
    sim_res_adj <- rbind(sim_res_adj, list(j, 'wilcoxon', get_AUC(p.val_wil, true_signal), FDR_wil[1], FDR_wil[2], FDR_wil[3], beta_diff))
    sim_res_adj <- rbind(sim_res_adj, list(j, 't-test', get_AUC(p.val_t, true_signal), FDR_t[1], FDR_t[2], FDR_t[3], beta_diff))
    
    write.csv(sim_res_adj, file = 'C:/Users/nraja/Downloads/sim_res_adj.csv', row.names = F)
    View(sim_res_adj)
  }
}

sim_res_conf <- data.frame()

for(beta_diff in c(0.75)) {
  for(j in c(3)) {
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
    sim_res_conf <- lm_sim(sim_res = sim_res_conf, j)
    sim_res_conf <- limma_sim(sim_res = sim_res_conf, j)
    View(sim_res_conf)
    sim_res_conf <- strainspy_sim(sim_res = sim_res_conf, j)
  }
}

