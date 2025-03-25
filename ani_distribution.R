# Boxplot of ANI distribution based on a specified predictor
# Parameters:
# se - SummarizedExperiment
# genome_file - String containing specified organism
# predictor - String containing the distinguishing predictor
# plot - boolean indicating whether a plot should be produced

ani_boxplot <- function(se, genome_file, predictor, plot = T) {
  
  # Validate input
  if(!inherits(se, 'SummarizedExperiment')) {
    stop('Input must be a SummarizedExperiment object.')
  }

  # Access assay data
  asy <- SummarizedExperiment::assay(se)
  
  # Access indices of inputted contig and predictor
  contig <- which(SummarizedExperiment::rowData(se)$Genome_file == genome_file)
  case <- which(colnames(tibble::as_tibble(se@colData@listData)) == predictor)
  
  # Validate that the contig/predictor exists
  if(length(contig) == 0 | length(case) == 0) {
    stop('Provide a valid predictor or contig (ex. \'GCF_025147485.1\')')
  }
  
  # Access the data for the predictor
  case <- se@colData@listData[[case]]
  
  # Create data.frame for plotting
  plot_data <- data.frame(sample = names(asy[contig,]), ani = unname(asy[contig,]), predictor = case)
  
  # Return data if requested
  if(!plot) {
    return(plot_data)
  }
  
  custom_colors <- c('#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5')
  
  custom_colors <- custom_colors[1:length(unique(plot_data$predictor))]
  
  
  plot <- ggplot2::ggplot(plot_data, aes(x=predictor,y=ani, fill = predictor)) +
    geom_boxplot(outliers=F) +
    ggbeeswarm::geom_quasirandom() +
    ggplot2::xlab(label = predictor) +
    ggplot2::scale_fill_manual(values = custom_colors) +
    ggplot2::scale_y_continuous(name = 'ANI', limits = c(97,100)) +
    theme_classic() +
    theme(legend.position='none')
  
  return(plot)

}

ani_boxplot(se, 'GCF_025147485.1', 'Case_status', plot=T)
