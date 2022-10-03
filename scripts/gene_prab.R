# Author: Noah Legall
# Purpose: Creation of Pangenome curve figure for mbovpan 
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(ggdendro)
library(ggtree)
library(ggnewscale)

# load in the general metadata
args = commandArgs(trailingOnly=TRUE)




# load in the gene presence absence data, keep only accessory
gene_pres_abs <- read.csv("gene_presence_absence.csv", header = TRUE, stringsAsFactors = FALSE, row.names = "Gene")
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
rownames(accessory_transpose) <- gsub(".annot","",rownames(accessory_transpose))


accessory_matrix <- as.matrix(accessory_pa %>% select(-gene_id))
rownames(accessory_matrix) <- gene_id
accessory_dendro <- as.dendrogram(hclust(d = dist(accessory_transpose, method = "binary"), method = "ward.D"))

# reorder the rows 
accessory_order <- order.dendrogram(accessory_dendro)
accessory_pa_long$sample <- factor(x = accessory_pa_long$sample,
                               levels = rownames(accessory_transpose)[accessory_order], 
                               ordered = TRUE)

pdf("gene_prab_figures.pdf")

if(length(args[1]) != 0){
  isolate_dat <- read.csv(args[1], stringsAsFactors = FALSE)
  for(i in 1:length(colnames(isolate_dat))){
    if(i == "Name"){
      next
    }
    else{  
      ad_gg <- ggtree(accessory_dendro)
      ad_gg[["data"]]$label <- gsub(".annot","",ad_gg$data$label)

      ad_gg <- ad_gg %<+% isolate_dat + 
        geom_tippoint(aes(color = i))

      t1 <- gheatmap(ad_gg,accessory_transpose, offset = 0.5, colnames = FALSE) +  
        scale_fill_manual(values = c("gray75","darkblue"), name = "Presence/Absence")
      plot(t1)
    }
}
  
} else{
  ad_gg <- ggtree(accessory_dendro)
  ad_gg[["data"]]$label <- gsub(".annot","",ad_gg$data$label)
  
  t1 <- gheatmap(ad_gg,accessory_transpose, offset = 0.5, colnames = FALSE) +  
    scale_fill_manual(values = c("gray75","darkblue"), name = "Presence/Absence")
  plot(t1) 
}
dev.off()