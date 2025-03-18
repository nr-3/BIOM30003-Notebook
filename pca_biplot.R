# PCA plot

pca_biplot <- function(se, predictor) {

  asy<-SummarizedExperiment::assay(se)
  PC <- prcomp(t(as.matrix(asy)), scale. = T)
  
  pve <- PC$sdev^2/sum(PC$sdev^2)
  xLab <- paste0("PC1 (",round(pve[1]*100,2),"%)")
  yLab <- paste0("PC2 (",round(pve[2]*100,2),"%)")
  
  #Format input data
  if(is.null(rownames(PC$x))){
    rownames(PC$x)<-1:nrow(PC$x)
  }
  
  dat <- tibble::as_tibble(PC$x)
  
  col <- which(colnames(tibble::as_tibble(se@colData@listData)) == predictor)
  col <- se@colData@listData[[col]]
  
  dat$pred <- col
  
  # Create a PCA plot
  plot <- ggplot(
    data = dat, 
    mapping = aes(x=PC1, y=PC2, color = pred)
    ) +
    geom_point()+
    theme_bw()+
    coord_fixed()+
    labs(x=xLab,y=yLab)

  return(plot)
}
