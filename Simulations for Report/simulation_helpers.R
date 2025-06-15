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

# AUC Calculations
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

### Plotting functions
model_colors <- c("ordinal" = "#1b9e77", 
               "zib"     = "#d95f02", 
               "limma"   = "#7570b3", 
               "LM"      = "#e7298a", 
               "wilcoxon"= "#66a61e", 
               "t-test"  = "#e6ab02",
               "Logistic" = "#1f78b4")

# Plot AUC
plot_AUC <- function(sim_res) {
  
  # Calculate the mean for each model and scalar
  AUC_dat <- sim_res |>
    group_by(model, beta) |>
    summarise(
      mean = mean(AUC),
    )
  
  # Plot lines for each model's mean AUC across scalars
  plot <- AUC_dat |>
    ggplot(aes(x = factor(beta), y = mean, group = model, color = model)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.2) +
    scale_color_manual(values = model_colors) +
    ggthemes::theme_clean(base_size = 16) +
    theme(
      plot.background = element_blank(),
      legend.background = element_blank(),
      legend.title = element_blank()
    ) +
    labs(x = 'Scalar Factor', y = 'Mean AUROC')
    
    
  return(plot)
}

plot_FDR <- function(sim_res) {
  
  # Calculate the mean for each model and scalar
  FDR_dat <- sim_res |>
    group_by(model, beta) |>
    summarise(FDR = mean(FDR))
  
  # Plot lines for each model's mean precision across scalars
  plot <- FDR_dat |>
    ggplot(aes(x=factor(beta),y=1-FDR,group = model,color=model)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.2) +
    scale_color_manual(values = model_colors) +
    ggthemes::theme_clean(base_size = 16) +
    theme(plot.background = element_blank(),
          legend.background = element_blank(), legend.title = element_blank()) +
    labs(x = 'Scalar Factor', y = 'Mean Precision (1-FDR)')
  
  
  return(plot)
  
}

plot_recall <- function(sim_res) {
  
  # Calculate the mean for each model and scalar
  recall_dat <- sim_res |>
    group_by(model, beta) |>
    summarise(recall = mean(recall))
  
  # Plot lines for each model's mean recall across scalars
  plot <- recall_dat |>
    ggplot(aes(x=factor(beta),y=recall,group=model,color=model)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.2) +
    scale_color_manual(values = model_colors) +
    ggthemes::theme_clean(base_size = 16) +
    theme(plot.background = element_blank(),
          legend.background = element_blank(), legend.title = element_blank()) +
    labs(x = 'Scalar Factor', y = 'Mean Recall')
  
  
  return(plot)
  
}

# Model fitting functions

# Fits the linear model and logistic regression for each barcode
lm_sim <- function(simulation_se, sim_res, j, true_signal) {
  
  # Accesses implanted data
  glm_dat <- as.matrix(SummarizedExperiment::assays(simulation_se)[[1]])
  grp <- SummarizedExperiment::colData(simulation_se)$sim
  # Creates binaries for logistic regression
  binary_grp <- ifelse(grp == 'a', 1, 0)
  num_barcodes <- dim(simulation_se)[1]
  
  # Storage vectors for p values
  p.val_lm <- c()
  p.val_log <- c()
  
  # Loops through each barcode
  for(i in c(1:num_barcodes)) {
    # Accesses the data for each barcode
    b_i <- glm_dat[i,]
    # Creates a data frame with normalized proportions via asin
    b_i <- data.frame(asin_proportion = asin(sqrt(b_i)), proportion = b_i, 
                      sim = grp, sim_bi = binary_grp, Months = simulation_se$Months)
    # Fits the single barcode data for the two models
    lm_sim <- lm(asin_proportion ~ sim + Months, data = b_i)
    log_sim <- glm(sim_bi ~ proportion + Months, data = b_i)
    
    # Accesses the p value for each barcode and stores in vector
    aov_sim <- anova(lm_sim)
    p.val_lm[i] <- aov_sim$`Pr(>F)`
    p.val_log[i] <- summary(log_sim)$coefficients[2,4]
  }
  
  # Gets the FDR, precision, and recall for each model based on the implanted signal
  FDR_lm <- get_FDR(p.adjust(p.val_lm, method='BH'),true_signal)
  FDR_glm <- get_FDR(p.adjust(p.val_log, method = 'BH'), true_signal)
  
  # Updates the results data.frame
  sim_res <- rbind(sim_res, list(j, 'LM', get_AUC(p.val_lm, true_signal), FDR_lm[1], FDR_lm[2], FDR_lm[3] , beta_diff))
  sim_res <- rbind(sim_res, list(j, 'Logistic', get_AUC(p.val_log, true_signal), FDR_glm[1], FDR_glm[2], FDR_glm[3] , beta_diff))
  
  return(sim_res)
}

limma_sim <- function(simulation_se, sim_res, j, true_signal) {
  # Creates design matrix
  design_limma <- model.matrix(~ SummarizedExperiment::colData(simulation_se)$sim + SummarizedExperiment::colData(simulation_se)$Months)
  
  # Accesses and normalizes the data
  limma_dat <- as.matrix(SummarizedExperiment::assays(simulation_se)[[1]])
  limma_dat <- asin(sqrt(limma_dat))
  
  # Applies lmFit according to design matrix and then analysis via empirical Bayes
  limma_fit <- limma::lmFit(limma_dat, design_limma)
  limma_fit <- limma::eBayes(limma_fit)
  
  # Accesses and applies BH correction to the p values
  p.val_limma <- limma_fit$p.value[,2]
  adj_limma <- p.adjust(p.val_limma, method = 'BH')
  
  # Gets the FDR, precision, and recall
  FDR_limma <- get_FDR(adj_limma,true_signal)
  
  # Updates the simulation results
  sim_res <- rbind(sim_res, list(j, 'limma', get_AUC(p.val_limma, true_signal), FDR_limma[1], FDR_limma[2], FDR_limma[3], beta_diff))
  
  return(sim_res)
}

# Similar process as "lm_sim" except applied to Wilcoxon and t
wilcoxon_t_sim <- function(simulation_se, sim_res, j, true_signal) {
  wil_dat <- as.matrix(SummarizedExperiment::assays(simulation_se)[[1]])
  grp <- SummarizedExperiment::colData(simulation_se)$sim == 'a'
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
  
  sim_res <- rbind(sim_res, list(j, 'wilcoxon', get_AUC(p.val_wil, true_signal), FDR_wil[1], FDR_wil[2], FDR_wil[3], beta_diff))
  sim_res <- rbind(sim_res, list(j, 't-test', get_AUC(p.val_t, true_signal), FDR_t[1], FDR_t[2], FDR_t[3], beta_diff))
}


# Similar to limma_sim except applied to the StrainSpy methods
strainspy_sim <- function(simulation_se, sim_res, j, true_signal) {
  sim_zib <- strainspy::glmZiBFit(se = simulation_se, design = as.formula('~sim + Months'))
  sim_ordb <- strainspy::glmFit(se = simulation_se, design = as.formula('~sim + Months'))
  
  # AUC Calculations (similar to Wirbel et al.'s SIMBA package)
  p.val.zib <- c()
  beta.p <- sim_zib@p_values$simb
  zero.p <- sim_zib@zi_p_values$simb
  for(i in c(1:dim(simulation_se)[[1]])) {
    p.val.zib[i] <- min(beta.p[i],zero.p[i])
  }
  p.val.ord <- sim_ordb@p_values$simb
  
  FDR_ord <- get_FDR(p.adjust(p.val.ord,method='BH'),true_signal)
  FDR_zib <- get_FDR(p.adjust(p.val.zib,method='BH'),true_signal)
  
  sim_res <- rbind(sim_res, list(j, 'ordinal', get_AUC(p.val.ord, true_signal), FDR_ord[1], FDR_ord[2], FDR_ord[3], beta_diff))
  sim_res <- rbind(sim_res, list(j, 'zib', get_AUC(p.val.zib, true_signal), FDR_zib[1], FDR_zib[2], FDR_zib[3], beta_diff))
  
  return(sim_res)
}

# Loop which runs all sims
all_sims <- function(j, beta_diff, zero_effect) {
  sim2_res <- data.frame()
  
  n <- dim(se_HSPC_2grp)[[2]]
  n_rows <- dim(se_HSPC_2grp)[[1]]
  
  # Split samples into two random groups
  sample <- c(1:n)
  grp1 <- sort(sample(sample, n/2))
  grp2 <- setdiff(sample, grp1)
  simulation_se <- se_HSPC_2grp
  SummarizedExperiment::colData(simulation_se)$sim <- NA
  SummarizedExperiment::colData(simulation_se)$sim[grp1] <- 'a'
  SummarizedExperiment::colData(simulation_se)$sim[grp2] <- 'b'
  
  # 10% of the barcodes to be differentially abundant
  implant <- round(0.1*n_rows)
  implant <- sample(c(1:n_rows),implant)
  
  true_signal <- rep(0, n_rows)
  true_signal[implant] <- 1
  grp1 <- SummarizedExperiment::colData(simulation_se)$sim=='a'
  
  # Implant the signal in the barcodes randomly selected
  for (i in implant) {
    s <- as.vector(SummarizedExperiment::assays(simulation_se)[[1]][i,])
    s <- rescale_beta(s[grp1], beta = beta_diff, zi = zero_effect)
    SummarizedExperiment::assays(simulation_se)[[1]][i,grp1] <- s$rescaled
  }
  sim2_res <- lm_sim(simulation_se, sim_res = sim2_res, j, true_signal)
  sim2_res <- limma_sim(simulation_se, sim_res = sim2_res, j, true_signal)
  sim2_res <- wilcoxon_t_sim(simulation_se, sim_res = sim2_res,j, true_signal)
  sim2_res <- strainspy_sim(simulation_se, sim_res = sim2_res, j, true_signal)

  
  colnames(sim2_res) <- colnames
  
  return(sim2_res)
}
