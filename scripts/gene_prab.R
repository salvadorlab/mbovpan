# Author: Noah Legall
# Purpose: Creation of Pangenome curve figure for mbovpan 
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(ggdendro)
library(ggtree)
library(ggnewscale)

gene_pres_abs <- read.csv("../noahaus/mbovpan/auxilary/gene_presence_absence.csv", header = TRUE, stringsAsFactors = FALSE, row.names = "Gene")
accessory_genome <- gene_pres_abs[!(is.na(gene_pres_abs$Accessory.Fragment)),]
core_genome <- gene_pres_abs[is.na(gene_pres_abs$Accessory.Fragment),]
auxil <- gene_pres_abs %>% select(2:14)
accessory_pa <- accessory_genome %>% select(14:(ncol(accessory_genome)))

accessory_pa[!(accessory_pa=="")] <- 1
accessory_pa[accessory_pa==""] <- 0 

num_col <- ncol(accessory_pa)
task <- function(x){sum(as.numeric(as.character(x)))}
accessory_pa$pr <- apply(accessory_pa,1,task)
accessory_pa$perc_pr <- accessory_pa$pr/num_col
accessory_pa <- accessory_pa %>% filter(perc_pr >= 0.15 & perc_pr <= 0.99) %>% select(-c("pr","perc_pr"))

accessory_pa$gene_id = rownames(accessory_pa)
accessory_pa_long <- accessory_pa %>% gather(sample,prab,-gene_id)
accessory_matrix <- as.matrix(accessory_pa %>% select(-gene_id))

# perform clustering 
accessory_transpose <- t(accessory_matrix)


accessory_matrix <- as.matrix(accessory_pa %>% select(-gene_id))
rownames(accessory_matrix) <- gene_id
accessory_dendro <- as.dendrogram(hclust(d = dist(accessory_transpose, method = "binary"), method = "ward.D"))

dendro_plot <- ggdendrogram(data = accessory_dendro, rotate = TRUE) + 
  theme(axis.text.y = element_blank())

# Preview the plot
print(dendro_plot)

# reorder the rows 
accessory_order <- order.dendrogram(accessory_dendro)
accessory_pa_long$sample <- factor(x = accessory_pa_long$sample,
                               levels = rownames(accessory_transpose)[accessory_order], 
                               ordered = TRUE)

t1 <- gheatmap(ggtree(accessory_dendro),accessory_transpose, offset = 0.5, colnames = FALSE) +  
  scale_fill_manual(values = c("gray75","darkblue"), name = "Presence/Absence")




#heatmap_plot <- accessory_pa_long %>% ggplot(aes(x = gene_id,y = sample, fill = prab)) + 
#  geom_tile() +
#  scale_fill_manual(values = c("white","steelblue")) +
#  theme_minimal() + 
#  theme(axis.text = element_blank(),
#        axis.ticks = element_blank(),
#        legend.position = "top")
#grid.newpage()
#print(heatmap_plot, 
#      vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
#print(dendro_plot, 
#      vp = viewport(x = 0.90, y = 0.465, width = 0.2, height = 1.1))

# use ggplot to make the heatmap 
#accessory_gg <- accessory_pa_long %>% ggplot(aes(x = gene_id,y = sample, fill = prab)) + 
#  geom_tile() +
#  scale_fill_manual(values = c("white","steelblue")) +
#  theme_minimal() + 
#  theme(axis.text = element_blank(),
#        axis.ticks = element_blank())