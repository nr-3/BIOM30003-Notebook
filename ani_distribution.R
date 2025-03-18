# Boxplot of ANI distribution based on a specified predictor
ani_distribution <- function(se, genome_file, predictor, plot = T) {
  
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
    stop('Provide a valid predictor or genome file')
  }
  
  # Access the data for the predictor
  case <- se@colData@listData[[case]]
  
  # Create data.frame for plotting
  plot_data <- data.frame(sample = names(asy[contig,]), ani = unname(asy[contig,]), predictor = case)
  
  # Return data if requested
  if(!plot) {
    return(plot_data)
  }
  
  custom_colors <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#b15928',
                     '#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5')
  
  custom_colors <- custom_colors[1:length(unique(plot_data$predictor))]
  
  
  plot <- ggplot2::ggplot(plot_data, aes(x=predictor,y=ani, fill = predictor)) +
    geom_boxplot() +
    ggplot2::xlab(label = predictor) +
    ggplot2::scale_fill_manual(values = custom_colors) +
    ggplot2::scale_y_continuous(name = 'ANI', limits = c(97,100)) +
    theme_classic() +
    theme(legend.position='none')
  
  return(plot)

}

ani_distribution(se, 'GCF_025147485.1', 'Case_status', plot=T)
