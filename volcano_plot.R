#' Generate a Volcano Plot from a betaGLM Object
#'
#' This function creates a volcano plot from the results of a betaGLM analysis,
#' highlighting the relationship between effect size (coefficients) and statistical significance (-log10 p-values).
#'
#' @param fit A `betaGLM` object containing model results.
#' @param coef Integer, specifying the phenotype index of the betaGLM to extract p values and effect size for plotting. Default is 2.
#' @param alpha Float indicating signficance level
#' @param palette A vector of user defined palette to colour plot. Default is the viridis palette.
#' @param label Logical, whether to label most signifcant points of the plot with Species name.
#' @param regex String, a regular expression which can be used to format the labels
#' @param plot Logical, whether to return a ggplot object (default: TRUE). If FALSE, returns the processed data.
#'
#'
#' @return A ggplot2 object if `plot = TRUE`, otherwise a dataframe containing the relevant statistics.

volcano_plot <- function(fit, coef = 2, alpha = 0.05, 
                         palette = c("#DAA520", "#3CB371", "#483D8B"), 
                         label = F, 
                         regex = '\\b[A-Z][a-z]+ [a-z]+(?: [a-z]+)?(?: [A-Z]{2,}[^\\s,]*)*', 
                         plot = T) {
  
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
      tidyr::pivot_longer(cols = c(p_value, zi_p_value), names_to = 'Model', values_to = 'p') |>
      dplyr::select(Genome_file, Model, p) |>
      dplyr::mutate(Model = ifelse(stringr::str_detect(Model, 'zi'), 'ZiB', 'Beta'))
    p_adj <- hits |>
      tidyr::pivot_longer(cols = c(p_adjust, zi_p_adjust), names_to = 'Model', values_to = 'p_adj') |>
      dplyr::select(Genome_file, Model, p_adj) |>
      dplyr::mutate(Model = ifelse(stringr::str_detect(Model, 'zi'), 'ZiB', 'Beta'))
    hits <- hits |>
      tidyr::pivot_longer(cols = c(coefficient,zi_coefficient), 
                          names_to = 'Model', values_to = 'Coefficient') |>
      dplyr::mutate(Model = ifelse(stringr::str_detect(Model, 'zi'), 'ZiB', 'Beta')) |>
      dplyr::inner_join(p, by = c('Genome_file', 'Model')) |>
      dplyr::inner_join(p_adj, by = c('Genome_file', 'Model'))
    
    # Return the data if requested
    if(!plot) {
      return(hits)
      
    }
    
    # Plot facetted volcano plots
    plot <- hits |> ggplot2::ggplot(aes(x = Coefficient, y = -log10(p), colour = p_adj)) +
      ggplot2::geom_point() +
      ggplot2::facet_grid(~ Model, scale = 'free_x') +
      ggthemes::theme_few() +
      ggplot2::scale_color_gradientn(colors = palette) + 
      ggplot2::theme(panel.spacing = unit(1, 'cm')) +
      ggplot2::xlab('Effect Size') +
      ggplot2::ylab('-log10(p)')
    
    if(label) {
      plot <- plot + ggrepel::geom_text_repel(data = hits[hits$p_adj < alpha,], 
                                              mapping = aes(label = str_extract(Contig_name, regex)),
                                              colour = 'black',
                                              size = 2.5,
                                              alpha = 0.8,
                                              max.overlaps = Inf)
    }
  }
  
  # Plot one-model volcano plot
  else {
    hits <- hits |>
      dplyr::mutate(log10p = log10(p_value)) |>
      dplyr::mutate(sig = ifelse(p_adjust < alpha, T, F))
    
    if(!plot) {
      return(hits)
    }
    
    plot <- hits |> ggplot2::ggplot(ggplot2::aes(x = coefficient, y = -log10p, color = p_adjust)) +
      ggplot2::geom_point() +
      ggthemes::theme_few() +
      ggplot2::scale_color_gradientn(colors = palette) +
      ggplot2::xlab('Effect Size') +
      ggplot2::ylab('-log10(p)')
  }
  
  # If labels are requested, then add.
  if(label) {
    plot <- plot + ggrepel::geom_text_repel(data = hits[hits$p_adjust < alpha,], 
                                            mapping = aes(label = stringr::str_extract(Contig_name, regex)),
                                            colour = 'black',
                                            size = 2.5,
                                            alpha = 0.8,
                                            max.overlaps = Inf)
  }
  
  return(plot)
}

############################################################################################

volcano_plot(fit_p, alpha = 0.20, label = T)


plot + ggrepel::geom_label_repel(data = plot_d[plot_d$signif,], mapping = aes(label = Genome_file))

plot_d |>
  dplyr::left_join(taxonomy, by = dplyr::join_by(Genome_file==Genome)) |>
  dplyr::filter(p < 0.001) |>
  dplyr::select(Genome_file, Genus, Species)


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


