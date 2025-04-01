# Lee data
meta_data <- readr::read_csv('C:/Users/nraja/BIOM30003-Notebook/Lee/lee_PRJEB43119_metadata.csv')
lee_data <- readr::read_tsv('C:/Users/nraja/BIOM30003-Notebook/Lee/lee_sylph_profile_gtdb-r220-c200-dbv1.tsv/lee_sylph_profile_gtdb-r220-c200-dbv1.tsv')

join_table <- readr::read_tsv("C:/Users/nraja/Downloads/filereport_read_run_PRJEB43119.tsv")

join_table$fastq_ftp <- gsub('^(.*?/){6}','',join_table$fastq_ftp)
join_table$fastq_ftp <- gsub(';.*','',join_table$fastq_ftp)

join_table$fastq_ftp <- paste0('./',join_table$fastq_ftp)
join_table$sample_title <- gsub(' wgs','',join_table$sample_title)

meta_data <- join_table |>
  select(fastq_ftp,sample_title) |>
  inner_join(meta_data, by = join_by(sample_title == sample_id))

se_lee <- read_sylph('C:/Users/nraja/BIOM30003-Notebook/Lee/lee_sylph_profile_gtdb-r220-c200-dbv1.tsv/lee_sylph_profile_gtdb-r220-c200-dbv1.tsv', meta_data)
se_lee <- filter_by_presence(se_lee, min_nonzero = 30)

pca_biplot(se_lee, 'antibiotics_current_use')

se_lee_q <- read_sylph('C:/Users/nraja/BIOM30003-Notebook/Lee/lee_sylph_query_gtdb_220_id99.tsv\lee_sylph_query_gtdb_220_id99.tsv')

design <- as.formula(" ~ antibiotics_current_use")

fit <- glmZiBFit(se_lee, design)

taxonomy <- read_taxonomy(example_taxonomy_path)

volcano_plot(se = se_lee, fit = fit, predictor = 'antibiotics_current_use')
head(top_hits(fit, alpha = 1))

case <- which(colnames(tibble::as_tibble(se_lee@colData@listData)) == 'antibiotics_current_use')
case <- se_lee@colData@listData[[case]]
case
unique_case <- unique(case)
dat <- as.data.frame(as.matrix(SummarizedExperiment::assay(se_lee)))
dat[,case==unique_case[1]]
dim(dat)
is.na(case)