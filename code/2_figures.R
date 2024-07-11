#ALPHA DIVERSITY

metadata$group <- with(metadata, 
                       ifelse(combined=="human.SI.original", "Donor SI", 
                              ifelse(combined=="human.Stool.original", "Donor feces", 
                                     ifelse(combined=="human.Stool.mouse.Stool", "Mouse feces (input feces)",
                                            ifelse(combined=="human.SI.mouse.SI", "Mouse SI (input SI)",
                                                   ifelse(combined=="human.SI.mouse.Stool", "Mouse feces (input SI)",
                                                          ifelse(combined=="human.Mix.mouse.SI", "Mouse SI (input SI+feces)",
                                                                 ifelse(combined=="human.Mix.mouse.Stool", "Mouse feces (input SI+feces)",
                                                                        ifelse(combined=="human.Stool.mouse.SI", "Mouse SI (input feces)",NA)))))))))


samples_of_int <- metadata$sample_name[metadata$group %in% c("Donor SI", "Donor feces", "Mouse SI (input SI)", "Mouse SI (input feces)","Mouse SI (input SI+feces)")]

rownames(metadata) <- metadata$sample_name

data_RA_sub_diversity <- data_RA_sub[,samples_of_int]
data_counts_sub_diversity <- data_counts_sub[,samples_of_int]
metadata_diver <- metadata[samples_of_int,]

species <- specnumber(data_counts_sub_diversity, MARGIN=2) 
raremax <- min(colSums(data_counts_sub_diversity))
data_mat_rar <- t(rrarefy(t(data_counts_sub_diversity), raremax)) 
data_mat_rar_ra <- sweep(data_mat_rar, 2, colSums(data_mat_rar), '/')

alpha <-vegan::diversity(data_mat_rar_ra, MARGIN=2, index = "shannon")

metadata_diver$group <- factor(metadata_diver$group, levels= c(
  "Mouse SI (input SI+feces)",
  "Mouse SI (input feces)",
  "Donor feces",
  "Mouse SI (input SI)",
  "Donor SI"))

my_colors <- c(
  "Mouse SI (input SI+feces)" = "#1b9e77", 
  "Mouse SI (input feces)" = "#7570b3",
  "Donor feces" = "#a6761d","Mouse SI (input SI)" = "#d95f02", 
  "Donor SI" = "#3399FF")



figure_alpha <- ggplot(metadata_diver, aes(x = group, y = alpha, color = group)) +
  geom_boxplot(outlier.shape = NA, fill = "black", color="black", varwidth = T, alpha = 0.1, position = position_dodge(width = 0.8), size = 1, lwd=0.7) +
  stat_boxplot(geom ='errorbar',  color="black", width=0.3, size=1) + 
  geom_jitter(shape = 16, size = 3, position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8), alpha = 1) +
  scale_color_manual(values = my_colors) +
  labs(x = "", y = "Shannon diversity") +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 16, color="black"),
    axis.title.y = element_text(size = 14, color="black"), 
    axis.text.x = element_text(size = 16, color="black"), 
    axis.text.y = element_text(size = 16, color="black"),
    axis.line = element_line(size = 1.2),
    axis.line.x.top = element_line(color = "black", size = 1.2),
    axis.line.y.right = element_line(color = "black", size = 1.2),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),axis.ticks = element_line(color = "black", size = 1.2),
    axis.ticks.length = unit(0.2, "cm") )+ 
  scale_y_continuous(labels = scales::label_number())+
  geom_segment(aes(x = 5, xend = 4, y = 3.5, yend = 3.5), color = "black", size = 0.8, lineend = "square") +
  geom_segment(aes(x = 5, xend = 3, y = 3.9, yend = 3.9), color = "black", size = 0.8, lineend = "square") +
  geom_segment(aes(x =5, xend = 2, y = 4.3, yend = 4.3), color = "black", size = 0.8, lineend = "square") +
  geom_segment(aes(x =5, xend = 1, y = 4.7, yend =4.7), color = "black", size = 0.8, lineend = "square") +
  annotate("text", x = 4.5, y = 3.75, label = "ns", size = 6, angle=270) +
  annotate("text", x = 4, y = 4, label = "**", size = 7, angle=270) +
  annotate("text", x = 3.5, y = 4.4, label = "***", size = 7, angle=270) +
  annotate("text", x = 3, y = 4.8, label = "**", size = 7, angle=270) +
  geom_segment(aes(x = 5, xend = 5, y = 2.6 , yend = 4.7), color = "black", size = 0.8) +
  geom_segment(aes(x = 4, xend = 4, y = 3.3, yend = 3.5), color = "black", size = 0.8) +
  geom_segment(aes(x = 3, xend = 3, y = 3.9, yend = 3.9), color = "black", size = 0.8) +
  geom_segment(aes(x = 2, xend = 2, y = 4.2, yend = 4.3), color = "black", size = 0.8)+
  geom_segment(aes(x = 1, xend = 1, y =3.9 , yend = 4.7), color = "black", size = 0.8)


figure_alpha


fit <- aov(alpha ~ metadata_diver$group)
TukeyHSD(fit)





#BETA DIVERSITY

set.seed(42)

beta_div<- as.matrix(vegdist(t(data_mat_rar), method = "bray"))
#cbind(rownames(beta_div), rownames(metadata_diver))

PCOA <- as.data.frame(pcoa(beta_div)$vectors)

var_expl <- pcoa(beta_div)$values$Relative_eig *100 

colnames(PCOA) <- paste("PC",1:ncol(PCOA), sep="")

PCOA <- cbind(PCOA, combined = metadata_diver$combined)
PCOA <- cbind(PCOA, group = metadata_diver$group)

PCOA$PC1 <- as.numeric(PCOA$PC1)
PCOA$PC2 <- as.numeric(PCOA$PC2)

range(PCOA$PC1) 
range(PCOA$PC2)

centroids <- aggregate(cbind(PCOA$PC1,PCOA$PC2) ~ PCOA$combined,PCOA,mean)
colnames(centroids) <- c("combined", "PC1", "PC2")

pc1_lab <- paste("PC1 ", round(var_expl[1], digits=1), "%", sep="")
pc2_lab <- paste("PC2 ", round(var_expl[2], digits=1), "%", sep="")

table(PCOA$group)


figure_beta<-ggplot(PCOA, aes(x = PC1, y = PC2, color = group)) + 
  geom_point(size = 3, alpha=0.9, show.legend= TRUE) + 
  scale_color_manual(values=my_colors) + 
  scale_fill_manual(values=my_colors) + 
  #guides(color=guide_legend(nrow=8), alpha=TRUE) +
  stat_chull(aes(x=PC1, y=PC2, color=group, fill=group), alpha = 0.15, geom = "polygon",show.legend=TRUE, lwd=0.75) +
  labs(x= pc1_lab, y= pc2_lab) +
  theme_cowplot() +
  theme(legend.title=element_blank(), legend.direction = "vertical", legend.position = c(0.1,0.9))+
  theme_minimal()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.text=element_text(size=14),
        axis.line = element_line(color = "black"),
        axis.line.x = element_line(color = "black", size = 1.2),
        axis.line.y = element_line(color = "black", size = 1.2),
        axis.text = element_text(size = 14,color = "black"),
        axis.title = element_text(size = 16,color = "black")) + 
  theme(axis.ticks = element_line(linewidth = 1),axis.ticks.length = unit(0.2, "cm"))+ 
  scale_x_continuous(
    breaks = seq(-0.4,0.6, by = 0.2), 
    limits = c(-0.4, 0.6),
    labels = label_number()) + scale_y_continuous(
      breaks = seq(-0.4, 0.4, by = 0.2), 
      limits = c(-0.4,0.4),
      labels = label_number()
    ) + coord_fixed(ratio = 1)



figure_beta

