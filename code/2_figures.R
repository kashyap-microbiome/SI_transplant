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

#
data_RA_sub_diversity <- data_RA_sub[,samples_of_int]
data_counts_sub_diversity <- data_counts_sub[,samples_of_int]
metadata_diver <- metadata[samples_of_int,]


species <- specnumber(data_counts_sub_diversity, MARGIN=2) 
raremax <- min(colSums(data_counts_sub_diversity))
data_mat_rar <- t(rrarefy(t(data_counts_sub_diversity), raremax)) 
data_mat_rar_ra <- sweep(data_mat_rar, 2, colSums(data_mat_rar), '/')

alpha <- vegan::diversity(data_mat_rar_ra, MARGIN=2, index = "shannon")

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


metadata_diver$alpha <- as.numeric(alpha)


fit <- aov(alpha ~ metadata_diver$group)
TukeyHSD(fit)


# Create data frames for annotations
sig_bars <- data.frame(
  x = c(5, 5, 5, 5, 5, 4, 3, 2, 1),
  xend = c(4, 3, 2, 1, 5, 4, 3, 2, 1),
  y = c(3.5, 3.9, 4.3, 4.7, 2.6, 3.3, 3.9, 4.2, 3.9),
  yend = c(3.5, 3.9, 4.3, 4.7, 4.7, 3.5, 3.9, 4.3, 4.7)
)

star_labels <- data.frame(
  x = c(4.5, 4, 3.5, 3),
  y = c(3.69, 3.95, 4.35, 4.75),
  label = c("ns", "**", "***", "**")
)

# Improved plot code
figure_alpha <- ggplot(metadata_diver, aes(x = group, y = alpha)) +
  # Boxplot with consistent styling
  geom_boxplot(
    aes(fill = group),
    outlier.shape = NA,
    color = "black",
    alpha = 0.1,
    varwidth = TRUE,
    position = position_dodge(width = 0.8),
    size = 0.7
  ) +
  # Error bars
  stat_boxplot(
    geom = 'errorbar',
    color = "black",
    width = 0.3,
    size = 1
  ) + 
  # Jitter points with group coloring
  geom_jitter(
    aes(color = group, shape = human_source),  # Map shape to human_source
    size = 3,
    position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8),
    alpha = 1
  ) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  scale_shape_manual(values = c(16, 17, 15)) + 
  
  # Axis labels and theme
  labs(x = "", y = "Shannon diversity") +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16, face = "bold"),      
    axis.text = element_text(size = 14, face = "bold"), 
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_blank(), 
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.line = element_line(size = 1.2, color = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_line(color = "black", size = 1.2),
    axis.ticks.length = unit(0.2, "cm")
  ) +
  # Significance annotations
  geom_segment(
    data = sig_bars,
    aes(x = x, xend = xend, y = y, yend = yend),
    color = "black",
    size = 0.8,
    lineend = "square"
  ) +
  geom_text(
    data = star_labels,
    aes(x = x, y = y, label = label),
    size = 9,
    angle = 270
  )

figure_alpha


figure_alpha <- figure_alpha +
  scale_x_discrete(labels = c("Mouse SI (input SI+feces)" = "Mouse SI\n(input human SI + feces)", 
                              "Mouse SI (input feces)" = "Mouse SI\n(input human feces)",
                              "Donor feces" = "Human donor feces",
                              "Mouse SI (input SI)" = "Mouse SI\n(input human SI)",
                              "Donor SI" = "Human donor SI"))





###

#now make the alpha diversity plot based on number of detected species
observed_ASVs <- colSums(data_mat_rar_ra > 0.0005) #above 0.05

hist(observed_ASVs)

metadata_diver$observed_ASVs <- as.numeric(observed_ASVs)

fit <- aov(observed_ASVs ~ metadata_diver$group)
TukeyHSD(fit)




# Original data with updated x coordinates
sig_bars <- data.frame(
  x = c(5, 5, 5, 5, 1, 4, 3, 2, 5, 2),
  xend = c(4, 3, 2, 1, 1, 4, 3, 2, 5, 2),
  y = c(85, 100, 155, 165, 165, 85, 95, 110, 165, 155),
  yend = c(85, 100, 155, 165, 160, 80, 100, 100, 70, 150)
)


star_labels <- data.frame(
  x = c(4.5, 4, 3.5, 3),
  y = c(93, 102, 156, 167),
  label = c("ns", "**", "***", "***")
)


# Improved plot code
figure_observed <- ggplot(metadata_diver, aes(x = group, y = observed_ASVs)) +
  # Boxplot with consistent styling
  geom_boxplot(
    aes(fill = group),
    outlier.shape = NA,
    color = "black",
    alpha = 0.1,
    varwidth = TRUE,
    position = position_dodge(width = 0.8),
    size = 0.7
  ) +
  # Error bars
  stat_boxplot(
    geom = 'errorbar',
    color = "black",
    width = 0.3,
    size = 1
  ) + 
  # Jitter points with group coloring
  geom_jitter(
    aes(color = group, shape = human_source),  # Map shape to human_source
    size = 3,
    position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8),
    alpha = 1
  ) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  scale_shape_manual(values = c(16, 17, 15)) + 
  
  # Axis labels and theme
  labs(x = "", y = "Observed ASVs") +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16, face = "bold"),      
    axis.text = element_text(size = 14, face = "bold"), 
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_blank(), 
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.line = element_line(size = 1.2, color = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_line(color = "black", size = 1.2),
    axis.ticks.length = unit(0.2, "cm")
  ) +
  # Significance annotations
  geom_segment(
    data = sig_bars,
    aes(x = x, xend = xend, y = y, yend = yend),
    color = "black",
    size = 0.8,
    lineend = "square"
  ) +
  geom_text(
    data = star_labels,
    aes(x = x, y = y, label = label),
    size = 9,
    angle = 270
  )

figure_observed




figure_observed <- figure_observed +
  scale_x_discrete(labels = c("Mouse SI (input SI+feces)" = "", 
                              "Mouse SI (input feces)" = "",
                              "Donor feces" = "",
                              "Mouse SI (input SI)" = "",
                              "Donor SI" = ""))

# Combine the two plots side by side
combined_plot <- figure_alpha + figure_observed + plot_layout(ncol = 2)

plot_name <- "../figures/alpha combined.png"
ggsave(plot_name, combined_plot, width = 11, height = 5, units = "in", dpi = 300)





###

#BETA DIVERSITY

set.seed(42)

beta_div <- as.matrix(vegdist(t(data_mat_rar), method = "bray"))
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


PCOA$human_source <- metadata_diver$human_source



figure_beta <- ggplot(PCOA, aes(x = PC1, y = PC2, color = group, shape = human_source)) + 
  geom_point(size = 3, alpha = 0.9, show.legend = TRUE) + 
  scale_color_manual(values = my_colors) + 
  scale_fill_manual(values = my_colors) + 
  stat_chull(
    aes(x = PC1, y = PC2, group = group, color = group, fill = group), 
    alpha = 0.15, geom = "polygon", show.legend = TRUE, lwd = 0.75
  ) +
  labs(x = pc1_lab, y = pc2_lab) +
  theme_cowplot() +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.direction = "vertical",
    legend.position = "right",  # <-- This moves the legend to the right
    legend.text = element_text(size = 14),
    axis.line = element_line(color = "black"),
    axis.line.x = element_line(color = "black", size = 1.2),
    axis.line.y = element_line(color = "black", size = 1.2),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  ) +
  theme(axis.ticks = element_line(linewidth = 1), 
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 16, face = "bold"),      
        axis.text = element_text(size = 14, face = "bold"), 
        plot.title = element_text(size = 16, face = "bold")) +
  scale_x_continuous(
    breaks = seq(-0.4, 0.6, by = 0.2), 
    limits = c(-0.4, 0.6),
    labels = scales::label_number()
  ) +
  scale_y_continuous(
    breaks = seq(-0.4, 0.4, by = 0.2), 
    limits = c(-0.4, 0.4),
    labels = scales::label_number()
  ) +
  coord_fixed(ratio = 1)



figure_beta <- figure_beta #+ guides(shape = "none")


plot_name <- "../figures/beta.png"
ggsave(plot_name, figure_beta, width = 8, height = 4.5, units = "in", dpi = 300)



