#' Generate a Volcano Plot from a betaGLM Object
#'
#' This function creates a volcano plot from the results of a betaGLM analysis,
#' highlighting the relationship between effect size (coefficients) and statistical significance (-log10 p-values).
#'
#' @param fit A `betaGLM` object containing model results.
#' @param coef Integer, specifying the coefficient index of the to extract top hits for plotting. Default is 2.
#' @param plot Logical, whether to return a ggplot object (default: TRUE). If FALSE, returns the processed data.
#'
#' @return A ggplot2 object if `plot = TRUE`, otherwise a dataframe containing the relevant statistics.

volcano_plot <- function(fit, coef = 2, plot = T) {
  
  # Validate input
  if(!inherits(fit, 'betaGLM')) {
    stop('Input must be a betaGLM object.')
  }
  
  # Access p values and effect sizes
  hits <- top_hits(fit, coef = coef, alpha = 1)
  
  # Checks if multiple models exist in the fit
  if ('zi_p_adjust' %in% colnames(hits)) {
    
    # Creating a tibble that can be facetted by Model
    p <- hits |>
      tidyr::pivot_longer(cols = c(p_value, zi_p_value), names_to = 'Model', values_to = 'P_value') |>
      dplyr::select(Genome_file, Model, P_value) |>
      dplyr::mutate(Model = ifelse(stringr::str_detect(Model, 'zi'), 'ZiB', 'Beta'))
    hits <- hits |>
      tidyr::pivot_longer(cols = c(coefficient,zi_coefficient), 
                          names_to = 'Model', values_to = 'Coefficient') |>
      dplyr::mutate(Model = ifelse(stringr::str_detect(Model, 'zi'), 'ZiB', 'Beta')) |>
      dplyr::inner_join(p, by = c('Genome_file' = 'Genome_file', 'Model' = 'Model'))
    
    # Return the data if requested
    if(!plot) {
      return(hits)
      
    }
    
    # Plot facetted volcano plots
    plot <- hits |> ggplot2::ggplot(aes(x = Coefficient, y = -log10(P_value), colour = p_adjust)) +
      ggplot2::geom_point() +
      ggplot2::facet_grid(~ Model, scale = 'free_x') +
      ggthemes::theme_few() +
      ggplot2::theme(panel.spacing = unit(1, "cm"), axis.title.x = NULL) 
  } 
  
  # Plot one-model volcano plot
  else {
    hits <- hits |>
      dplyr::mutate(log10p = log10(p_value))
    
    if(!plot) {
      return(hits)
    }
    
    plot <- hits |> ggplot2::ggplot(ggplot2::aes(x = coefficient, y = -log10padj)) +
      ggplot2::geom_point()
  }
  
  return(plot)
}

############################################################################################

volcano_plot(fit_p,plot = T)
  


# NOT FOR PACKAGE
# Initial volcano plot which manually calculates fold change (didn't work)

volcano_plot_initial <- function(se, fit, predictor, plot = T) {
  
  # Creating tidy tibble for plot
  idx <- which(colnames(tibble::as_tibble(se@colData@listData)) == predictor)
  case <- se@colData@listData[[idx]]
  
  # Removing NA values
  se <- se[,!is.na(case)]
  asy <- SummarizedExperiment::assay(se)
  case <- se@colData@listData[[idx]]
  unique_case <- unique(case)
  
  # Creating the fold change
  plot_data <-as.data.frame(as.matrix(asy))
  plot_data$case1mean <- rowMeans(plot_data[,case==unique_case[1]], na.rm = T)
  plot_data$case2mean <- rowMeans(plot_data[,case==unique_case[2]], na.rm = T)
  plot_data$fold_change <- plot_data$case1mean/plot_data$case2mean
  plot_data <- tibble::rownames_to_column(plot_data, 'Contig_name')
  
  # Including the p-values in the plot data
  hits <- top_hits(fit, alpha = 1)
  plot_data <- plot_data |>
    dplyr::select(Contig_name, fold_change) |>
    dplyr::inner_join(hits) |>
    dplyr::mutate(log10padj = log10(p_adjust))
  
  # Create a plot
  plot <- plot_data |> ggplot2::ggplot(ggplot2::aes(x = log2(fold_change),y = -log10padj)) +
    ggplot2::geom_point()
  
  return(plot)
}


