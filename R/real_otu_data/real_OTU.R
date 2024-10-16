#rm(list=ls())
#library(edgeR)
#library(scran)
library(dplyr)
load("data/real/OTU/ageset.Rdata")

OTU_original <- data.obj$otu.tab
OTU_original <- t(OTU_original)
#filter OTU

# Filter based on OTU prevalence
prevalence_threshold <- 0.1
OTU_prevalence <- apply(OTU_original != 0, 2, mean)
OTU_keep <- OTU_prevalence >= prevalence_threshold

# Filter based on median non-zero counts
median_count_threshold <- 10
median_counts <- apply(OTU_original, 2, function(x) median(x[x != 0], na.rm = TRUE))
OTU_keep <- OTU_keep & median_counts >= median_count_threshold

# Subset the OTU table based on the filtered columns
OTU_filtered <- OTU_original[, OTU_keep]
dim(OTU_filtered)
# calculate 97% quantile for each taxon
quantiles <- apply(OTU_filtered, 2, quantile, probs = 0.97)

# Replace values greater than the 97th percentile with the 97th percentile value
for (i in 1:ncol(OTU_filtered)) {
  OTU_filtered[OTU_filtered[,i] > quantiles[i], i] <- quantiles[i]
}
dim(OTU_filtered)

OTU_filtered = sqrt(OTU_filtered)

OTU_filtered_N<- t(apply(OTU_filtered, 1, function(OTU_filtered) OTU_filtered/sum(OTU_filtered)))


OTU = OTU_filtered_N


meta_data <- data.obj$meta.dat
age <- meta_data['age']
age <-sqrt(age)
age <- data.frame(data_name = rownames(age), age = age$age)
OTU<-data.frame(data_name = rownames(OTU),OTU)
# Merge the data frames
merged_data <- merge(OTU, age, by = "data_name")



country = meta_data['geo_loc_name']
# Create a new column that maps each country name to an integer value
country$geo_loc_name <- ifelse(country$geo_loc_name == "USA", 1, 
                               ifelse(country$geo_loc_name == "Malawi", 2,
                                      ifelse(country$geo_loc_name == "Venezuela", 3, NA)))

country$geo_loc_name = as.factor(country$geo_loc_name)
country <- data.frame(data_name = rownames(country), country = country$geo_loc_name)
merged_data2 <- merge(country,merged_data, by = "data_name")


#merged_data2 <- merged_data2 %>% select(-data_name)
# Delete the column named "data_name"
merged_data2  <- merged_data2 [, !colnames(merged_data2 ) %in% "data_name"]


dim(merged_data2)

merged_data2 <- merged_data2[complete.cases(merged_data2$age), ]
dim(merged_data2)


data = merged_data2

write.csv(data, "data/real/OTU/tss_normalized_data.csv", row.names = F)

# get taxonomies of OTUs included in this final data set
OTU_ids = colnames(data)
OTU_ids = OTU_ids[2:(length(OTU_ids)-1)]
OTU_taxonomies = data.obj$otu.name
final_OTU_taxonomies = OTU_taxonomies[str_remove(OTU_ids, "X"),]
write.csv(final_OTU_taxonomies, "data/real/OTU/taxonomies.csv", row.names = T)

