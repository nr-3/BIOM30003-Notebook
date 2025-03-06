library(strainspy)
library(tidyverse)

se <- read_sylph('zeevi_query.tsv')
se <- filter_by_presence(se, min_nonzero = 30)

taxonomy <- read_taxonomy('sylph_taxonomy.tsv')

temp <- taxonomy |>
  select(Genome, Species)

rowData(se) <- inner_join(temp, as.data.frame(rowData(se)), by = join_by(Genome == Genome_file))
se@assays@data@listData[[1]]@x <- se@assays@data@listData[[1]]@x/100

#Each unique species in se
zeevi_species <- unique(rowData(se)$Species)

#Limit the number of species
zeevi_species <- sample(zeevi_species, size = 3)

#sub_se
sub_se <- se[rowData(se)$Species %in% zeevi_species,]

#@assays@data@listData[[1]]@x

rowData(sub_se)$test <- F

#for(sp in zeevi_species) {
#Subset for one species
sp<-sample(zeevi_species, size = 1)

sp_se <- sub_se[rowData(sub_se)$Species==sp,]
sp_ani <- sp_se@assays@data@listData[[1]]@x
#Sample ANIs for one group
control_ani <- sample(sp_ani, size = 50)
#Another sample for the other
exp_ani <- sample(sp_ani[!sp_ani %in% control_ani], size = 50)


colData(sp_se)$test <- F
colData(sp_se)$test[sample((1:nrow(colData(sp_se))),99)] <- T

design <- as.formula(' ~ test')
fit <- glmZiBFit(sp_se, design, nthreads = 4)
top_hits(fit, alpha=0.5)
plot_manhattan(fit,method = 'holm')