####################################################################################
#load required packages

require(openxlsx)
require(phyloseq)
require(readxl)
require(vegan)
require(ape)
require(ggplot2)
require(RColorBrewer)
require(ggpubr)
require(cowplot)
library(scales)


####################################################################################
#set working directory

setwd("~/Dropbox/Mayo_RS/R/202205 SI transplant/2nd_round_code/github code/")
#setwd("./202205 SI transplant/") #change your working directory

file_name <- "./small_bowel_humanization_data_202210.xlsx"
meta_file_name <- "./metadata_10_26_22.xlsx"


data_counts <- as.data.frame(read_excel(file_name, sheet = 1), stringsAsFactors=F)

taxa_names_full <- data_counts$ASV_name #these are the species names
rownames(data_counts) <- taxa_names_full

data_counts <- data_counts[,2:ncol(data_counts)]


####################################################################################

metadata <- as.data.frame(read_excel(meta_file_name, sheet = 1), stringsAsFactors=F)

str(metadata)

#lets count the number of samples
table(metadata$origin) #mouse and human
table(metadata$sample_class)
#negative_controls        phenotypes         profiling 
#9                59                76 

head(metadata, 10)


####################################################################################
#subset samples for which we have metadata 
# & put samples and metadata in same order

data_counts <- data_counts[, colnames(data_counts) %in% metadata$sample_name]
data_counts <- data_counts[, order(colnames(data_counts))]
metadata <- metadata[order(metadata$sample_name),]

#inspect if they are in the same order
cbind(metadata$sample_name, colnames(data_counts))


####################################################################################
#check if there are samples we have to remove because they have very few sequencing 
#reads, especially check the negative controls, they should be low and removed

#one column is one sequencing file / stool sample.
#rows=species #cols=stoolsample
#the total number of sequencing reads is given by:
hist(colSums(data_counts))

#negative control rows
neg_control_names <- metadata$sample_name[which(metadata$sample_class == "negative_controls")]

hist(colSums(data_counts), las=1)
hist(colSums(data_counts)[neg_control_names], col="red", add=T)

colSums(data_counts)[neg_control_names]


metadata$nr_reads <- colSums(data_counts)


metadata[metadata$origin == "human" & metadata$sample_class == "phenotypes",]
metadata[metadata$origin == "human" & metadata$sample_class == "profiling",]


#manually inspecting read depth
metadata[metadata$human_source == "11.1",]
metadata[metadata$human_source == "11.2",]
metadata[metadata$human_source == "15.2",]

metadata[metadata$human_source == "SIBO201",] 
metadata[metadata$human_source == "DIET2004Z",] 
metadata[metadata$human_source == "DIET2013M",] 


####################################################################################
#see which taxa are found in negative controls

neg_control_taxa_names <- names(which(rowSums(data_counts[,neg_control_names] > 0) > 2))


####################################################################################
#remove samples with fewer reads than median of negative controls
#also remove negative controls

#remove negative controls
read_cutoff <- median(colSums(data_counts)[neg_control_names])

too_low_samples <- names(data_counts)[which(colSums(data_counts) <= read_cutoff)]
to_remove_names <- unique(c(neg_control_names, too_low_samples))

metadata <- metadata[-which(metadata$sample_name %in% to_remove_names),]
cols_to_rm <- which(colnames(data_counts) %in% to_remove_names)
data_counts <- data_counts[,-cols_to_rm]

dim(metadata); dim(data_counts)
#head(metadata); colnames(data_counts)


####################################################################################
#some inspection of what is annotated and present in a reasonable number of samples

#cutoff in >10% of samples
data_counts_sub <- data_counts[rowSums(data_counts > 1) > 0.1*ncol(data_counts),]
dim(data_counts_sub)
head(data_counts_sub, 2)

data_RA_sub <- sweep(data_counts_sub, MARGIN = 2, colSums(data_counts_sub), '/')
dim(data_RA_sub)
#350 119


####################################################################################
#get true genus level (collapse ASVs) and make relative abundance

genus_level <- as.character(sapply(taxa_names_full, function(x) strsplit(x, "\\.")[[1]][5]))

#remove all instances with _NUMBER, as this is another way of indicating ASVs
genus_level <- gsub("_\\d", "", genus_level)

agg <- aggregate(data_counts, by=list(genus_level),sum)

#remove first 4 and unclassified
agg_sub <- agg[5:nrow(agg),]
to_rm <- which(agg_sub$Group.1 == "unclassified")
agg_sub <- agg_sub[-to_rm,]

data_counts_genus <- sapply(agg_sub[,-1], function(y) {as.numeric(as.character(y))})
rownames(data_counts_genus) <- agg_sub$Group.1


#removing rows that miss >90%
data_counts_genus_sub <- data_counts_genus[rowSums(data_counts_genus > 1) > 0.1*ncol(data_counts_genus),]
data_RA_genus_sub <- sweep(data_counts_genus_sub, 2, colSums(data_counts_genus_sub), '/')

dim(data_RA_genus_sub)
#132 119


####################################################################################
#inspecting first 2 lines of the data
head(data_RA_sub,2)

#makes a new column for metadata combining source and sample type
metadata$combined <- paste(metadata$transplant_source, metadata$sample_type, sep=".")


####################################################################################
#subset by experiment type

head(metadata)

table(metadata$sample_class)

profiling_rows <- which(metadata$sample_class == "profiling")
phenotypes_rows <- which(metadata$sample_class == "phenotypes")

#dim(data_RA_sub)
#dim(data_counts_sub)
#dim(data_counts_genus_sub)
#dim(data_RA_genus_sub)

metadata_full <- metadata
metadata_phenotyping <- metadata[phenotypes_rows,]
metadata <- metadata[profiling_rows,]

dim(metadata_full); dim(metadata); dim(metadata_phenotyping)


data_RA_sub_full <- data_RA_sub
data_counts_sub_full <- data_counts_sub
data_counts_genus_sub_full <- data_counts_genus_sub
data_RA_genus_sub_full <- data_RA_genus_sub

dim(data_RA_sub_full); dim(data_counts_sub_full); dim(data_counts_genus_sub_full); dim(data_RA_genus_sub_full)


data_RA_sub <- data_RA_sub[,profiling_rows]
data_counts_sub <- data_counts_sub[,profiling_rows]
data_counts_genus <- data_counts_genus_sub[,profiling_rows]
data_RA_genus_sub <- data_RA_genus_sub[,profiling_rows]

dim(data_RA_sub); dim(data_counts_sub); dim(data_counts_genus); dim(data_RA_genus_sub)


data_RA_sub_phenotyping <- data_RA_sub_full[,phenotypes_rows]
data_counts_sub_phenotyping <- data_counts_sub_full[,phenotypes_rows]
data_counts_genus_sub_phenotyping <- data_counts_genus_sub_full[,phenotypes_rows]
data_RA_genus_sub_phenotyping <- data_RA_genus_sub_full[,phenotypes_rows]

dim(data_RA_sub_phenotyping); dim(data_counts_sub_phenotyping); dim(data_counts_genus_sub_phenotyping); dim(data_RA_genus_sub_phenotyping)



