# Manhattan plot
library(strainspy)

# Default values for the function (for now)
coef <- 2
taxonomy <- read_taxonomy(example_taxonomy_path)
method <- 'HMP'
alpha <- 0.05
levels <- c("Phylum", "Genus", "Species", "Strain")
colour_level <- 'Phylum'
plot <- T

hier_p <- hadjust(fit, taxonomy = taxonomy)

plot_data <- plot_manhattan(fit, method = 'holm', plot = F)


plot_data_tax <- plot_manhattan(fit, taxonomy = taxonomy, method = method, plot = F)
plot_data_tax <- plot_data_tax |>
  dplyr::arrange(Phylum)
  
plot_manhattan(fit, taxonomy = taxonomy, method = method, plot = T)

View(plot_data_tax)

# Revised manhattan with taxonomical data
gg <- ggplot2::ggplot(plot_data_tax, aes(x = Name, y = log_p_adjust, colour = Phylum)) +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_grid(Level ~ Model, scales = 'free_x')

gg

# Revised manhattan without taxonomical data (not much I can do)
gg <- ggplot2::ggplot(plot_data, aes)


fit
