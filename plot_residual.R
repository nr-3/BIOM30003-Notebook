#plot_residual <- function(object) {

# Access the residual data (need to double check if this is truly the data)
res_data <- as.tibble(fit@residuals@listData)

# Manipulate the data such that the residuals exist in one long column
test <- as.tibble(res_data) |>
  pivot_longer(cols = names(res_data)[1:ncol(res_data)], names_to = 'sample', values_to = 'residuals') |>
  # This is a temporary addition: the residuals are averaged out for each sample
  dplyr::group_by(sample) |>
  dplyr::summarise(residuals=mean(residuals))

# Plot the residuals with a horizontal line at zero
res_plot <- ggplot2::ggplot(test) +
  ggplot2::geom_point(aes(x = sample, y = residuals)) +
  ggplot2::geom_hline(aes(yintercept = 0), col='red', linetype = 'dotted') +
  ggplot2::theme(
    axis.text.x=ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    )
  
 # return(res_plot)

#}