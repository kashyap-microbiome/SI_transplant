permanova_function <- function(tax_df, method, groups, meta) {
  
  require(vegan)
  require(ape)
  require(RColorBrewer)
  require(ggpubr)
  require(cowplot)
  require(reshape2)
  require(dplyr)
  
  
  beta_div <- as.matrix(vegdist(t(tax_df), method = method))
  
  #Adonis to check differences in centroid
  beta_dist = as.dist(beta_div)
  col_ind <- which(colnames(meta) == groups)

  set.seed(42)
  
  ad <- adonis(beta_div ~ meta[,col_ind], permutations=999)
  p_val <- ad$aov.tab[1,6]
  r_sq <- ad$aov.tab[1,5] 
  
  #Run Stats for diff. dispersion
  disp <- betadisper(beta_dist, meta[,col_ind])
  p_val_disp <- permutest(disp)$tab[1, 6]
  
  
  ######################################################################
  #beta diversity boxplots
  
  sel_colors <- brewer.pal(n = 8, "Dark2")[c(3,2,1,4:8)]
  
  wu.m = reshape2::melt(as.matrix(beta_dist))
  wu.m = wu.m %>%
    filter(as.character(Var1) != as.character(Var2)) %>%
    mutate_if(is.factor, as.character)
  
  #add metadata to it
  wu.m$meta1 <- sapply(wu.m$Var1, function(x) meta[,col_ind][which(colnames(tax_df) == x)])
  wu.m$meta2 <- sapply(wu.m$Var2, function(x) meta[,col_ind][which(colnames(tax_df) == x)])
  
  wu.m$comb <- paste(wu.m$meta1, wu.m$meta2, sep="_")
  
  #this will result in 4 groups with two variables, for example 0 vs 1 and 1 vs 0 are the same and one should be removed
  wu.m_sub <- wu.m[wu.m$comb %in% c("0_0", "1_1"),]
  wu.m_sub$comb <- as.factor(paste(wu.m_sub$meta1, wu.m_sub$meta2, sep="_"))
  
  beta_boxplot <- ggplot(wu.m_sub, aes(x = comb, y = value)) +
    geom_boxplot(fill=c(sel_colors[1], sel_colors[2]), outlier.color=NA, notch=TRUE) +
    labs(y = "Bray Curtis dissimilarity", x= "Groups") +
    guides(color="none") +
    stat_compare_means(method = "t.test", label = "p.signif", label.x.npc = "center") +
    theme_cowplot()
  
  
  ######################################################################
  #PCOA of beta div
  
  PCOA <- pcoa(beta_div)$vectors #from ape package
  variances_explained <- pcoa(beta_div)$values$Relative_eig *100
  
  PCOA <- t(apply(PCOA, 1, as.numeric))
  PCOA <- as.data.frame(PCOA)
  colnames(PCOA) <- paste("PC",1:ncol(PCOA), sep="")
  
  PCOA <- cbind(PCOA, meta[,col_ind])
  colnames(PCOA)[ncol(PCOA)] <- "groups"
  
  centroids <- aggregate(cbind(PCOA$PC1,PCOA$PC2) ~ PCOA$groups ,PCOA,mean)
  colnames(centroids) <- c('groups', "PC1", "PC2")
  
  ######################################################################
  #pairwise adonis for differences: all pairwise PERMANOVA comparisons
  groups <- as.character(PCOA$groups)
  n <- length(unique(groups))
  
  #very small groups are not reliable; e.g. do not look at data for 4.
  
  pair_ad_list <- list() #initiate an empty list to store the data in
  for(i in 1:(n-1)){
    for(u in (i+1):n) {
      t1 <- unique(groups)[[i]]
      t2 <- unique(groups)[[u]]
      sams <- rownames(PCOA[groups == t1 | groups == t2,])
      beta_dist = beta_div[sams,sams]
      ad_2 = adonis(beta_dist ~ PCOA[sams,'groups'], data=PCOA, permutations=999)
      p_val2 <- ad_2$aov.tab[1,6]
      r_sq2 <- ad_2$aov.tab[1,5]
      pair_ad_list[[paste(t1, "_vs_", t2, sep="")]] <- c("pval", p_val2, "r_sq", r_sq2)
    }
  }
  
  
  permanova_res_list <- list(permanova = ad, perm_p = p_val, perm_rsq = r_sq, 
                             dispersion = disp, disp_p = p_val_disp, pairwise_perm = pair_ad_list)
  
  
  ######################################################################
  #plot PCoA
  
  pc1_lab <- paste("PC1 ", round(variances_explained[1], digits=1), "%", sep="")
  pc2_lab <- paste("PC2 ", round(variances_explained[2], digits=1), "%", sep="")
  
  PCoA <- ggplot(PCOA) +
    geom_point(size = 2.5, aes_string(x = "PC1", y = "PC2", color = as.factor(PCOA$groups), alpha=0.25)) + 
    scale_color_manual(values=sel_colors) +
    guides(color=guide_legend(nrow=3), alpha="none") +
    theme(legend.title=element_blank()) +
    stat_chull(aes(x=PC1,y=PC2,color=as.factor(groups), fill=as.factor(groups)), alpha = 0.01, geom = "polygon") +
    scale_fill_manual(values=sel_colors) +
    labs(x= pc1_lab, y= pc2_lab) +
    theme_cowplot()
  
  
  ######################################################################
  #combining results
  
  res_list <- list(permanova_res_list, PCoA, beta_boxplot, wu.m)
  
  return(res_list)  

  
}