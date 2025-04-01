#' Perform PCA and Generate a Biplot
#'
#' This function performs Principal Component Analysis (PCA) on a SummarizedExperiment object
#' and generates a PCA biplot colored by a specified phenotype.
#'
#' @param se A `SummarizedExperiment` object containing expression data.
#' @param predictor A character string specifying the phenotype by which the plot is coloured
#' @param plot A logical value indicating whether to return the PCA plot (`TRUE`) or PCA data (`FALSE`). Default is `TRUE`.
#'
#' @return If `plot = TRUE`, returns a `ggplot2` PCA plot. If `plot = FALSE`, returns a tibble with PCA results.

pca_biplot <- function(se, predictor, plot = T) {

  asy<-SummarizedExperiment::assay(se)
  
  # PCA
  PC <- prcomp(t(as.matrix(asy)), scale. = T)
  
  # Calculate the percent variance from the first two PC
  pve <- PC$sdev^2/sum(PC$sdev^2)
  xLab <- paste0("PC1 (",round(pve[1]*100,2),"%)")
  yLab <- paste0("PC2 (",round(pve[2]*100,2),"%)")
  
  #Format PC data
  if(is.null(rownames(PC$x))){
    rownames(PC$x)<-1:nrow(PC$x)
  }
  dat <- tibble::as_tibble(PC$x)
  
  # Access the phenotype by which the plot is coloured
  col <- which(colnames(tibble::as_tibble(se@colData@listData)) == predictor)
  col <- se@colData@listData[[col]]
  dat$pred <- col
  
  # Return plot data if requested
  if(!plot) {
    return(dat)
  }
  
  # Create a PCA plot
  plot <- ggplot2::ggplot(
    data = dat, 
    mapping = aes(x=PC1, y=PC2, color = pred)
    ) +
    ggplot2::geom_point()+
    ggplot2::theme_bw() +
    theme_classic()+
    labs(x=xLab,y=yLab)
  
  # Add a gradient colour if phenotype is continuous
  if(is.numeric(col)) {
    plot <- plot + scale_color_viridis_c(option = 'plasma')
  }

  return(plot)
}

pca_biplot(se_lee, predictor = 'antibiotics_current_use')

PC <- prcomp(t(as.matrix(asy)), scale. = T)
class(PC)
View(head(PC$rotation))

autoplot(PC,loadings=T,loadings.label=T)

