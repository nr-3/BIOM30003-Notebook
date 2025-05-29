# Simulation helper functions

# Modified function from strainspy to account for scaling
rescale_beta <- function(x, beta = 0.95, zi=0.1) {
  
  # Ensure x is between 0 and 1
  if (any(x < 0 | x > 1)) stop("Values must be between 0 and 1")
  # Ensure beta is between 0 and 1
  if (any(beta < 0 | beta > 1)) stop("beta must be between 0 and 1")
  
  # Apply power transformation (Beta regression-style scaling)
  x_rescaled <- x*beta
  
  p1 <- sum(x_rescaled==0)/length(x_rescaled)
  p2 <- p1 + zi
  
  if (p2 > 1) {
    stop("ZI is too large. The total proportion of
                   zeros is greater than 1 after adjustment")
  }
  
  if(p2 > p1){ # when zi > 0
    # We cannot guarantee this won't fail for very small zi values, if that happens, return input
    sz <- rbinom(1, length(x_rescaled), p2) - sum(x_rescaled==0)
    nz = which(x_rescaled!=0)
    if(sz < length(nz) & sz > 0) {
      selected <- sample(which(x_rescaled!=0), sz, replace=FALSE)
      x_rescaled[selected] <- 0
    } 
  }
  
  # Added line to scale everything to sum up to one (proportion)
  x_rescaled <- x_rescaled/sum(x_rescaled)
  
  return(list(
    rescaled=x_rescaled,
    expected_beta=log(beta),
    expected_zi=log(p2/(1-p2)/(p1/(1-p1)))
  ))
}

# Get AUROC

get_AUC <- function(p.val, true_signal) {
  auroc <- pROC::roc(predictor = -log10(p.val + 1e-50),
                     response = true_signal, levels = c(0,1), direction = '<')
  return(auroc$auc)
}
library(ggrepel)
# Get FDR

# FDR Calculations (similar to Wirbel et al.'s SIMBA package)
get_FDR <- function(p.val, true_signal, alpha = 0.05) {
  
  # Remove any NA p-values
  if(is.na(sum(p.val))) {
    na_val <- is.na(p.val)
    true_signal <- true_signal[!na_val]
    p.val <- p.val[!na_val]
  }
  
  # True Positives - Is significant and signal was implanted
  TP <- sum(true_signal[p.val < alpha] == 1)
  # True Negatives - Is not significant and no signal implanted
  TN <- sum(true_signal[p.val > alpha] == 0)
  # Type1 Error - Is signficant, but no signal was implanted
  FP <- sum(true_signal[p.val < alpha] == 0)
  # Type2 Error - Is not significant, but a signal was implanted
  FN <- sum(true_signal[p.val > alpha] == 1)
  # FDR
  FDR <- FP/(FP+TP)
  FDR <- ifelse(is.na(FDR), 0, FDR)
  # Precision
  PR <- TP/(FP+TP)
  PR <- ifelse(is.na(PR), 0, PR)
  # recall
  R <- TP/(TP+FN)
  R <- ifelse(is.na(R), 0, R)
  return(c(FDR,PR,R))
}

plot_AUC(sim_res)

# Plot AUC
plot_AUC <- function(sim_res) {
  AUC_dat <- sim_res |>
    group_by(model, beta) |>
    summarise(AUC = mean(AUC))
  
  plot <- AUC_dat |>
    ggplot(aes(x=beta,y=AUC,color=model)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = scales::hue_pal()(length(unique(AUC_dat$model)))) +
    theme_minimal(base_size = 14) +
    labs(x = 'Strength of Signal',y='Mean AUROC') +
    theme(
      legend.position = 'right',
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
    
    
  return(plot)
  
}

plot_FDR <- function(sim_res) {
  FDR_dat <- sim_res |>
    group_by(model) |>
    mutate(model_id = cur_group_id()) |>
    ungroup() |>
    group_by(model_id, model, beta) |>
    summarise(FDR = mean(FDR))
  
  plot <- FDR_dat |>
    ggplot(aes(x=beta,y=1-FDR,group=model_id,color=factor(model))) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = scales::hue_pal()(length(unique(AUC_dat$model_id)))) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  
  return(plot)
  
}

plot_FDR(sim_res_conf)
plot_AUC(sim_res_conf)


lm_sim <- function(sim_res, j) {
  glm_dat <- as.matrix(simulation_se@assays@data@listData[[1]])
  grp <- simulation_se$sim
  num_barcodes <- dim(simulation_se)[1]
  
  p.val_lm <- c()
  
  for(i in c(1:num_barcodes)) {
    b_i <- glm_dat[i,]
    b_i <- data.frame(proportion = asin(sqrt(b_i)),
                      sim = grp, Months = simulation_se$Months)
    lm_sim <- lm(proportion ~ sim + Months, data = b_i)
    aov_sim <- anova(lm_sim)
    p.val_lm[i] <- aov_sim$`Pr(>F)`
  }
  
  FDR_glm <- get_FDR(p.adjust(p.val_lm, method='BH'),true_signal)
  
  sim_res <- rbind(sim_res, list(j, 'LM', get_AUC(p.val_lm, true_signal), FDR_glm[1], FDR_glm[2], FDR_glm[3] , beta_diff))
  
  return(sim_res)
}

limma_sim <- function(sim_res, j) {
  design_limma <- model.matrix(~ simulation_se$sim + simulation_se$Months)
  limma_dat <- as.matrix(simulation_se@assays@data@listData[[1]])
  limma_dat <- asin(sqrt(limma_dat))
  
  limma_fit <- lmFit(limma_dat, design_limma)
  limma_fit <- eBayes(limma_fit)
  
  p.val_limma <- limma_fit$p.value[,2]
  
  adj_limma <- p.adjust(p.val_limma, method = 'BH')
  
  FDR_limma <- get_FDR(adj_limma,true_signal)
  
  sim_res <- rbind(sim_res, list(j, 'limma', get_AUC(p.val_limma, true_signal), FDR_limma[1], FDR_limma[2], FDR_limma[3], beta_diff))
  
  return(sim_res)
}

strainspy_sim <- function(sim_res, j) {
  sim_zib <- strainspy::glmZiBFit(se = simulation_se, design = as.formula('~sim + Months'))
  sim_ordb <- strainspy::glmFit(se = simulation_se, design = as.formula('~sim + Months'))
  
  # AUC Calculations (similar to Wirbel et al.'s SIMBA package)
  p.val.zib <- sim_zib@zi_p_values$simb
  p.val.ord <- sim_ordb@p_values$simb
  
  FDR_ord <- get_FDR(p.adjust(p.val.ord,method='BH'),true_signal)
  FDR_zib <- get_FDR(p.adjust(p.val.zib,method='BH'),true_signal)
  
  sim_res <- rbind(sim_res, list(j, 'ordinal', get_AUC(p.val.ord, true_signal), FDR_ord[1], FDR_ord[2], FDR_ord[3], beta_diff))
  sim_res <- rbind(sim_res, list(j, 'zib', get_AUC(p.val.zib, true_signal), FDR_zib[1], FDR_zib[2], FDR_zib[3], beta_diff))
  
  return(sim_res)
}
