################################################################################
#now look for engraftment


#----------------------------
#rarify and transform to relative abundance again
#for engraftment do less filtering

#ASV level
data_counts_asv <- data_counts[,metadata$sample_name]
data_counts_asv_sub <- data_counts_asv[rowSums(data_counts_asv > 1) > 2,]

raremax <- min(colSums(data_counts_asv_sub))
data_counts_asv_sub_rar <- t(rrarefy(t(data_counts_asv_sub), raremax)) 
data_mat_rar_ra <- sweep(data_counts_asv_sub_rar, MARGIN = 2, colSums(data_counts_asv_sub_rar), '/')


#genus level
data_counts_genus <- data_counts_genus_copy[,metadata$sample_name]
data_counts_genus_sub <- data_counts_genus[rowSums(data_counts_genus > 1) > 2,]

raremax <- min(colSums(data_counts_genus_sub))
data_counts_genus_sub_rar <- t(rrarefy(t(data_counts_genus_sub), raremax)) 
data_RA_sub_genus <- sweep(data_counts_genus_sub_rar, MARGIN = 2, colSums(data_counts_genus_sub_rar), '/')


#-------------------------------------------------------------------------------
#inspecting the relative abundance distribution (ylim set because of heavily zero inflated as expected)
hist(data_mat_rar_ra, breaks=50, ylim=c(0,200))

donors <- unique(metadata$human_source)


#using SI as input sample look for overlap per mouse of all detected taxa in the human sample
engraft_list_ASV <- list()
input_n_list_ASV <- list()
for (i in 1:length(donors)) {
  human_row <- intersect(grep("human.SI.original", metadata$combined), which(metadata$human_source == donors[i]))
  taxa_above_cutoff <- names(which(data_mat_rar_ra[,human_row] > 0))
  temp_rows <- intersect(grep("mouse", metadata$sample_type), which(metadata$human_source == donors[i]))
  
  engraft_list_ASV[[i]] <- colSums(data_mat_rar_ra[taxa_above_cutoff, temp_rows] > 0)
  input_n_list_ASV[[i]] <- length(taxa_above_cutoff)
}
names(engraft_list_ASV) <- donors
#$`11.2`
#SB10 SB11 SB12 SB13 SB14 SB15 SB16  SB5  SB6 SB64 SB65 SB66 SB67 SB68 SB69  SB7 SB70 SB71 SB72 SB73 SB74 SB75  SB8  SB9 
#   9    5    8    9    8    8    7    9   11    8   10    6    6    8    8    9    7    5    5    6    5    5   11    6 

#$`15.2`
#SB109 SB110 SB111 SB112 SB113 SB114 SB115 SB116 SB117 SB118 SB119  SB83  SB84  SB85  SB86  SB87  SB88  SB90  SB91  SB92  SB93 
# 10    17    13    34    34    33    40    32    28    30    32     4     5     5    25    23    19    42    29    28    30 

#$`11.1`
#SB19 SB20 SB21 SB22 SB23 SB24 SB25 SB26 SB27 SB28 SB29 SB30 SB41 SB42 SB43 SB44 SB45 SB46 SB47 SB48 SB49 SB50 SB51 SB52 
# 31   43   32   45   31   41   41   41   41   55   49   49   13   10   15   14   24   26   22   21   22   25   24   24 

names(input_n_list_ASV) <- donors
#$`11.2`
#[1] 75

#$`15.2`
#[1] 147

#$`11.1`
#[1] 124


agg_list_ASVs <- list() #aggregate this by combined
for (i in 1:length(donors)) {
  percentages <- engraft_list_ASV[[i]] / input_n_list_ASV[[i]]
  agg_list_ASVs[[i]] <- aggregate(percentages, list(metadata[names(percentages), "combined"]), mean)
}
names(agg_list_ASVs) <- donors

#$`11.2`
#Group.1          x
#1      human.Mix.mouse.SI 0.10666667
#2   human.Mix.mouse.Stool 0.07000000
#3       human.SI.mouse.SI 0.13333333
#4    human.SI.mouse.Stool 0.10000000
#5    human.Stool.mouse.SI 0.09333333
#6 human.Stool.mouse.Stool 0.09333333

#$`15.2`
#Group.1          x
#1      human.Mix.mouse.SI 0.23979592
#2   human.Mix.mouse.Stool 0.15192744
#3       human.SI.mouse.SI 0.09070295
#4    human.SI.mouse.Stool 0.03174603
#5    human.Stool.mouse.SI 0.20748299
#6 human.Stool.mouse.Stool 0.21938776

#$`11.1`
#Group.1         x
#1      human.Mix.mouse.SI 0.3911290
#2   human.Mix.mouse.Stool 0.1915323
#3       human.SI.mouse.SI 0.3104839
#4    human.SI.mouse.Stool 0.1048387
#5    human.Stool.mouse.SI 0.3044355
#6 human.Stool.mouse.Stool 0.1875000



#engraftment for genus level
engraft_list_genus <- list()
input_n_list_genus <- list()
for (i in 1:length(donors)) {
  human_row <- intersect(grep("human.SI.original", metadata$combined), which(metadata$human_source == donors[i]))
  taxa_above_cutoff <- names(which(data_RA_sub_genus[,human_row] >  0))
  temp_rows <- intersect(grep("mouse", metadata$sample_type), which(metadata$human_source == donors[i]))
  
  engraft_list_genus[[i]] <- colSums(data_RA_sub_genus[taxa_above_cutoff, temp_rows] > 0)
  input_n_list_genus[[i]] <- length(taxa_above_cutoff)
}
names(engraft_list_genus) <- donors
#$`11.2`
#SB10 SB11 SB12 SB13 SB14 SB15 SB16  SB5  SB6 SB64 SB65 SB66 SB67 SB68 SB69  SB7 SB70 SB71 SB72 SB73 SB74 SB75  SB8  SB9 
# 12   10   14   13   11   11   15   11   13    8    8    4   11    9    9   13   10    8    7   10    7    7   12   10 

#$`15.2`
#SB109 SB110 SB111 SB112 SB113 SB114 SB115 SB116 SB117 SB118 SB119  SB83  SB84  SB85  SB86  SB87  SB88  SB90  SB91  SB92  SB93 
# 10    23    16    34    32    31    32    29    25    31    27     5     7     6    25    24    22    37    29    24    29 

#$`11.1`
#SB19 SB20 SB21 SB22 SB23 SB24 SB25 SB26 SB27 SB28 SB29 SB30 SB41 SB42 SB43 SB44 SB45 SB46 SB47 SB48 SB49 SB50 SB51 SB52 
# 28   36   33   36   28   40   41   37   35   42   41   43   13   12   13   13   22   24   23   22   24   22   23   22 

names(input_n_list_genus) <- donors
#$`11.2`
#[1] 45

#$`15.2`
#[1] 80

#$`11.1`
#[1] 76


agg_list_genus <- list() #aggregate this by combined
for (i in 1:length(donors)) {
  percentages <- engraft_list_genus[[i]] / input_n_list_genus[[i]]
  agg_list_genus[[i]] <- aggregate(percentages, list(metadata[names(percentages), "combined"]), mean)
}
names(agg_list_genus) <- donors
#$`11.2`
#Group.1         x
#1      human.Mix.mouse.SI 0.2777778
#2   human.Mix.mouse.Stool 0.1722222
#3       human.SI.mouse.SI 0.2722222
#4    human.SI.mouse.Stool 0.1722222
#5    human.Stool.mouse.SI 0.2555556
#6 human.Stool.mouse.Stool 0.2000000

#$`15.2`
#Group.1         x
#1      human.Mix.mouse.SI 0.4031250
#2   human.Mix.mouse.Stool 0.2958333
#3       human.SI.mouse.SI 0.2041667
#4    human.SI.mouse.Stool 0.0750000
#5    human.Stool.mouse.SI 0.3500000
#6 human.Stool.mouse.Stool 0.3718750

#$`11.1`
#Group.1         x
#1      human.Mix.mouse.SI 0.5296053
#2   human.Mix.mouse.Stool 0.2993421
#3       human.SI.mouse.SI 0.4802632
#4    human.SI.mouse.Stool 0.1677632
#5    human.Stool.mouse.SI 0.4375000
#6 human.Stool.mouse.Stool 0.2993421


#-------------------------------------------------------------------------------
#now Stool as input sample

#using SI as input sample
engraft_list_ASV <- list()
input_n_list_ASV <- list()
for (i in 1:length(donors)) {
  human_row <- intersect(grep("human.Stool.original", metadata$combined), which(metadata$human_source == donors[i]))
  #get names of ASVs above 0.1% and then count how many of these are also in the mouse samples
  taxa_above_cutoff <- names(which(data_mat_rar_ra[,human_row] > 0))
  temp_rows <- intersect(grep("mouse", metadata$sample_type), which(metadata$human_source == donors[i]))
  
  engraft_list_ASV[[i]] <- colSums(data_mat_rar_ra[taxa_above_cutoff, temp_rows] > 0)
  input_n_list_ASV[[i]] <- length(taxa_above_cutoff)
}
names(engraft_list_ASV) <- donors
names(input_n_list_ASV) <- donors

agg_list_ASVs <- list() #aggregate this by combined
for (i in 1:length(donors)) {
  percentages <- engraft_list_ASV[[i]] / input_n_list_ASV[[i]]
  agg_list_ASVs[[i]] <- aggregate(percentages, list(metadata[names(percentages), "combined"]), mean)
}
names(agg_list_ASVs) <- donors



#engraftment for genus level
engraft_list_genus <- list()
input_n_list_genus <- list()
for (i in 1:length(donors)) {
  human_row <- intersect(grep("human.Stool.original", metadata$combined), which(metadata$human_source == donors[i]))
  #get names of ASVs above 0.1% and then count how many of these are also in the mouse samples
  taxa_above_cutoff <- names(which(data_RA_sub_genus[,human_row] >  0))
  temp_rows <- intersect(grep("mouse", metadata$sample_type), which(metadata$human_source == donors[i]))
  
  engraft_list_genus[[i]] <- colSums(data_RA_sub_genus[taxa_above_cutoff, temp_rows] > 0)
  input_n_list_genus[[i]] <- length(taxa_above_cutoff)
}
names(engraft_list_genus) <- donors
names(input_n_list_genus) <- donors

agg_list_genus <- list() #aggregate this by combined
for (i in 1:length(donors)) {
  percentages <- engraft_list_genus[[i]] / input_n_list_genus[[i]]
  agg_list_genus[[i]] <- aggregate(percentages, list(metadata[names(percentages), "combined"]), mean)
}
names(agg_list_genus) <- donors

#stool >60% engraftment


#-------------------------------------------------------------------------------
#see if at least the most abundant taxa in human donors are engrafted; most abundant genus

#engraftment for genus level with SI as input sample
engraft_list_genus_max <- list()
engraft_list_genus_any_mouse <- list()

for (i in 1:length(donors)) {
  human_row <- intersect(grep("human.SI.original", metadata$combined), which(metadata$human_source == donors[i]))
  #get the highest abundant taxa
  taxa_max <- names(data_RA_sub_genus[,human_row])[which(data_RA_sub_genus[,human_row] == max(data_RA_sub_genus[,human_row]))]
  temp_rows <- intersect(grep("mouse", metadata$sample_type), which(metadata$human_source == donors[i]))
  temp_test <- data_RA_sub_genus[taxa_max, temp_rows] > 0
  temp_test[temp_test == TRUE] <- "detected"
  temp_test[temp_test == FALSE] <- "not_detected"
  engraft_list_genus_max[[i]] <- temp_test
  
  
  #sub up the mice by sample type
  #get names of ASVs above 0.1% and then count how many of these are also in the mouse samples
  taxa_above_cutoff <- names(data_RA_sub_genus[,human_row])[data_RA_sub_genus[,human_row] > 0.001]
  
  temp_test_SI <- table(rowSums(data_RA_sub_genus[taxa_above_cutoff, temp_rows[metadata$sample_type[temp_rows] == "mouse.SI"]]) > 0)
  names(temp_test_SI)[names(temp_test_SI) == TRUE] <- "detected"
  names(temp_test_SI)[names(temp_test_SI) == FALSE] <- "not_detected"
  
  temp_test_Stool <- table(rowSums(data_RA_sub_genus[taxa_above_cutoff, temp_rows[metadata$sample_type[temp_rows] == "mouse.Stool"]]) > 0)
  names(temp_test_Stool)[names(temp_test_Stool) == TRUE] <- "detected"
  names(temp_test_Stool)[names(temp_test_Stool) == FALSE] <- "not_detected"
  
  temp_test <- rbind(temp_test_SI, temp_test_Stool)
  engraft_list_genus_any_mouse[[i]] <- temp_test
  
}
names(engraft_list_genus_max) <- donors
names(engraft_list_genus_any_mouse) <- donors


agg_list_genus_max <- list() #aggregate this by combined
for (i in 1:length(donors)) {
  agg_list_genus_max[[i]] <- table(engraft_list_genus_max[[i]], metadata[names(engraft_list_genus_max[[i]]), "combined"])
  
}
names(agg_list_genus_max) <- donors


#proportion of mice that have contain the most abundant genus from the human donor sample in their SI (at any level)
#$`11.2`
#                   human.Mix.mouse.SI human.Mix.mouse.Stool human.SI.mouse.SI human.SI.mouse.Stool human.Stool.mouse.SI human.Stool.mouse.Stool
#detected                      4                     1                 4                    4                    4                       3
#not_detected                  0                     3                 0                    0                    0                       1

#$`15.2`
#                   human.Mix.mouse.SI human.Mix.mouse.Stool human.SI.mouse.SI human.SI.mouse.Stool human.Stool.mouse.SI human.Stool.mouse.Stool
#detected                      4                     0                 3                    1                    1                       1
#not_detected                  0                     3                 0                    2                    3                       3

#$`11.1`
#                   human.Mix.mouse.SI human.Mix.mouse.Stool human.SI.mouse.SI human.SI.mouse.Stool human.Stool.mouse.SI human.Stool.mouse.Stool
#detected                  4                     4                 4                    4                    4                       4


#engraft_list_genus_any_mouse genus above 0.1% in human samples and then count how many of these are also in ANY of the mouse samples
engraft_list_genus_any_mouse

#$`11.2`
#                   not_detected detected
#temp_test_SI              15       13    13/28 = 46.4%
#temp_test_Stool           22        6

#$`15.2`
#                   not_detected detected
#temp_test_SI              18       29    29/49 = 59.2%
#temp_test_Stool           25       22

#$`11.1`
#                   not_detected detected
#temp_test_SI               9       21    21/30 = 70%
#temp_test_Stool           21        9



#-------------------------------------------------------------------------------
#engraftment for genus level with Stool as input sample

engraft_list_genus_max <- list()
engraft_list_genus_any_mouse <- list()

for (i in 1:length(donors)) {
  human_row <- intersect(grep("human.Stool.original", metadata$combined), which(metadata$human_source == donors[i]))
  #get the highest abundant taxa
  taxa_max <- names(data_RA_sub_genus[,human_row])[which(data_RA_sub_genus[,human_row] == max(data_RA_sub_genus[,human_row]))]
  temp_rows <- intersect(grep("mouse", metadata$sample_type), which(metadata$human_source == donors[i]))
  temp_test <- data_RA_sub_genus[taxa_max, temp_rows] > 0
  temp_test[temp_test == TRUE] <- "detected"
  temp_test[temp_test == FALSE] <- "not_detected"
  engraft_list_genus_max[[i]] <- temp_test
  
  
  #sub up the mice by sample type
  #get names of ASVs above 0.1% and then count how many of these are also in the mouse samples
  taxa_above_cutoff <- names(data_RA_sub_genus[,human_row])[data_RA_sub_genus[,human_row] > 0.001]
  
  temp_test_SI <- table(rowSums(data_RA_sub_genus[taxa_above_cutoff, temp_rows[metadata$sample_type[temp_rows] == "mouse.SI"]]) > 0)
  names(temp_test_SI)[names(temp_test_SI) == TRUE] <- "detected"
  names(temp_test_SI)[names(temp_test_SI) == FALSE] <- "not_detected"
  
  temp_test_Stool <- table(rowSums(data_RA_sub_genus[taxa_above_cutoff, temp_rows[metadata$sample_type[temp_rows] == "mouse.Stool"]]) > 0)
  names(temp_test_Stool)[names(temp_test_Stool) == TRUE] <- "detected"
  names(temp_test_Stool)[names(temp_test_Stool) == FALSE] <- "not_detected"
  
  temp_test <- rbind(temp_test_SI, temp_test_Stool)
  engraft_list_genus_any_mouse[[i]] <- temp_test
  
}
names(engraft_list_genus_max) <- donors
names(engraft_list_genus_any_mouse) <- donors


agg_list_genus_max <- list() #aggregate this by combined
for (i in 1:length(donors)) {
  agg_list_genus_max[[i]] <- table(engraft_list_genus_max[[i]], metadata[names(engraft_list_genus_max[[i]]), "combined"])
  
}
names(agg_list_genus_max) <- donors
#$`11.2`
#                    human.Mix.mouse.SI human.Mix.mouse.Stool human.SI.mouse.SI human.SI.mouse.Stool human.Stool.mouse.SI human.Stool.mouse.Stool
#detected                      4                     4                 3                    1                    4                       4
#not_detected                  0                     0                 1                    3                    0                       0

#$`15.2`
#                    human.Mix.mouse.SI human.Mix.mouse.Stool human.SI.mouse.SI human.SI.mouse.Stool human.Stool.mouse.SI human.Stool.mouse.Stool
#detected                      4                     3                 0                    1                    3                       3
#not_detected                  0                     0                 3                    2                    1                       1

#$`11.1`
#                    human.Mix.mouse.SI human.Mix.mouse.Stool human.SI.mouse.SI human.SI.mouse.Stool human.Stool.mouse.SI human.Stool.mouse.Stool
#detected                      4                     4                 4                    3                    4                       4
#not_detected                  0                     0                 0                    1                    0                       0


#engraft_list_genus_any_mouse
engraft_list_genus_any_mouse
#$`11.2`
#                   not_detected detected
#temp_test_SI               2       37    
#temp_test_Stool            7       32    32/39 = 82.1%

#$`15.2`
#                   not_detected detected
#temp_test_SI               5       43    
#temp_test_Stool            8       40    40/48 = 83.3%

#$`11.1`
#                   not_detected detected
#temp_test_SI               2       39    
#temp_test_Stool            9       32    32/41 = 78.0%



