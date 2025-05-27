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

# Get FDR

# FDR Calculations (similar to Wirbel et al.'s SIMBA package)
get_FDR <- function(p.val, true_signal, alpha = 0.05) {
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
}

