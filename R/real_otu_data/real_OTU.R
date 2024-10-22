#___________________________________
# Pre-processing of real OTU dataset
#___________________________________

# Load data ---------------------------------------------------------------
load("data/real/OTU/ageset.Rdata")
OTU_original <- data.obj$otu.tab
OTU_original <- t(OTU_original) # samples by OTUs
dim(OTU_original)

#  Filter OTU to remove less informative OTUs and reduce dimension --------
# We imposed two filters:
# (1) OTU prevalence < 10%
prevalence_threshold <- 0.1
OTU_prevalence <- apply(OTU_original != 0, 2, mean)
OTU_keep <- OTU_prevalence >= prevalence_threshold
# (2) Median non-zero counts < 10.
median_count_threshold <- 10
median_counts <- apply(OTU_original, 2, function(x) median(x[x != 0], na.rm = TRUE))
OTU_keep <- OTU_keep & median_counts >= median_count_threshold
# Perform filtering
OTU_filtered <- OTU_original[, OTU_keep]
dim(OTU_filtered)


# Add pseudo counts and normalize with centered log ratio transfor --------
OTU_filtered_pc <- OTU_filtered
OTU_filtered_pc[OTU_filtered_pc == 0] <- 0.5
OTU_filtered_n <- sweep(log(OTU_filtered_pc),
                        MARGIN=1,
                        STATS=rowMeans(log(OTU_filtered_pc)) )


# Form dataframe ----------------------------------------------------------
meta_data <- data.obj$meta.dat
country = meta_data$geo_loc_name
# Square-root transformation of the age
sqrt_age <- sqrt(meta_data$age)
df <- data.frame(country, OTU_filtered_n, sqrt_age) |>
  tibble::rownames_to_column("sample_id")
# Drop observations with missing age
df_complete <- tidyr::drop_na(df, sqrt_age)
dim(df_complete) # dropped 40 observations

# Write processed data and OTU phylo info matrix --------------------------
write.csv(df_complete, "data/real/OTU/normalized_data.csv", row.names = F)
OTU_kept_ids = colnames(OTU_filtered)
OTU_phylo = data.obj$otu.name
OTU_kept_phylo = OTU_phylo[OTU_kept_ids,]
write.csv(OTU_kept_phylo, "data/real/OTU/final_taxonomies.csv", row.names = T)

