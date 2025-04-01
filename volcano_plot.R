# Potential volcano plot??

# use geom jitter?

# Should log2fold change be in terms of difference in ANI?

volcano_plot <- function(se, fit, predictor, plot = T) {

  # Creating tidy tibble for plot
  asy <- SummarizedExperiment::assay(se)

  case <- which(colnames(tibble::as_tibble(se@colData@listData)) == predictor)
  case <- se@colData@listData[[case]]
  unique_case <- unique(case)

  plot_data <-as.data.frame(as.matrix(asy))
  plot_data$case1mean <- rowMeans(plot_data[,case==unique_case[1]], na.rm = T)
  plot_data$case2mean <- rowMeans(plot_data[,case==unique_case[2]], na.rm = T)
  plot_data$fold_change <- plot_data$case1mean/plot_data$case2mean
  plot_data <- tibble::rownames_to_column(plot_data, 'Contig_name')
  
  hits <- top_hits(fit, alpha = 1)
  plot_data <- plot_data |>
    dplyr::select(Contig_name, fold_change) |>
    dplyr::inner_join(hits) |>
    dplyr::mutate(log10padj = log10(p_adjust))
  
  plot_data$coef <- fit@coefficients@listData$Case_statusPD
  
  
  plot <- plot_data |> ggplot2::ggplot(ggplot2::aes(x = coef,y = -log10padj)) +
    ggplot2::geom_point()

  return(plot)
}