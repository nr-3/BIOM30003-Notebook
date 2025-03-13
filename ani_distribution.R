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
  if(contig == 0 | case == 0) {
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
  
  
  plot <- ggplot2::ggplot(plot_data, aes(x=predictor,y=ani)) +
    geom_boxplot() +
    ggplot2::xlab(label = predictor) +
    ggplot2::scale_y_continuous(name = 'ANI', limits = c(97,100))
  
  return(plot)

}

ANI_distribution(se, 'GCF_025147485.1', 'Case_status', plot=T)
