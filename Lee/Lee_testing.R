# Lee data
meta_data <- readr::read_csv('C:/Users/nraja/BIOM30003-Notebook/Lee/lee_PRJEB43119_metadata.csv')
lee_data <- readr::read_tsv('C:/Users/nraja/BIOM30003-Notebook/Lee/lee_sylph_profile_gtdb-r220-c200-dbv1.tsv/lee_sylph_profile_gtdb-r220-c200-dbv1.tsv')

# File names don't match between meta data and sylph file
# Join table created to match up file names
join_table <- readr::read_tsv("C:/Users/nraja/Downloads/filereport_read_run_PRJEB43119.tsv")

# Removing everything before the 6th `/` and after the `;`
join_table$fastq_ftp <- gsub('^(.*?/){6}','',join_table$fastq_ftp)
join_table$fastq_ftp <- gsub(';.*','',join_table$fastq_ftp)

# Format as it is in the sylph files
join_table$fastq_ftp <- paste0('./',join_table$fastq_ftp)
join_table$sample_title <- gsub(' wgs','',join_table$sample_title)

# Add the formatted file names back into the meta data
meta_data <- join_table |>
  select(fastq_ftp,sample_title) |>
  inner_join(meta_data, by = join_by(sample_title == sample_id))

# Read profile sylph file
se_lee <- read_sylph('C:/Users/nraja/BIOM30003-Notebook/Lee/lee_sylph_profile_gtdb-r220-c200-dbv1.tsv/lee_sylph_profile_gtdb-r220-c200-dbv1.tsv', meta_data)
se_lee <- filter_by_presence(se_lee, min_nonzero = 30)

# PCA plot check
pca_biplot(se_lee, "treatment")

# Read query sylph file
se_lee_q <- read_sylph('C:/Users/nraja/BIOM30003-Notebook/Lee/lee_sylph_query_gtdb_220_id99.tsv/lee_sylph_query_gtdb_220_id99.tsv', meta_data)
se_lee_q <- filter_by_presence(se_lee, min_nonzero = 30)

design <- as.formula(" ~ treatment")

# Fit using zero inflated beta
fit_p <- glmFit(se_lee, design)
fit_q <- glmZiBFit(se_lee_q, design)

taxonomy <- read_taxonomy(example_taxonomy_path)

# Volcano plot check
volcano_plot(fit_p)
fit_q@p_values

head(top_hits(fit, alpha = 1))