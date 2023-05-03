# Author: Noah Legall
# Purpose: Output a presence absence matrix of M. bovis virulence genes 
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(ggdendro)
library(ggtree)
library(ggnewscale)

# load in the general metadata
args = commandArgs(trailingOnly=TRUE)
print("testing")
gene_pres_abs <- read.csv("mbov_virulent_prab.csv", header = TRUE, stringsAsFactors = FALSE, row.names = "Gene")

# load in the gene presence absence data, keep only accessory
# we should already have access to this in our directory
print(args[1])
test = toString(args[1])
gene_pres_abs <- read.csv(test, header = TRUE, stringsAsFactors = FALSE, row.names = "Gene")
head(gene_pres_abs)
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
rownames(accessory_transpose) <- gsub("_trimmed_R1.scaffold.annot","",rownames(accessory_transpose))


accessory_matrix <- as.matrix(accessory_pa %>% select(-gene_id))
rownames(accessory_matrix) <- accessory_pa$gene_id
accessory_dendro <- as.dendrogram(hclust(d = dist(accessory_transpose, method = "binary"), method = "ward.D"))

ad_gg <- ggtree(accessory_dendro)
ad_gg[["data"]]$label <- gsub(".annot","",ad_gg$data$label)
print("testing")
pdf("gene_prab_figures.pdf")

if(length(args[2]) != 0){
  print("this shouldn't work")
  isolate_dat <- read.csv(args[2], stringsAsFactors = FALSE, check.names = FALSE)
  
  for(i in 1:length(colnames(isolate_dat))){
    if(colnames(isolate_dat)[i] == "Name" || length(unique(isolate_dat[,i])) == 1 ){
      next
    }
    else{
      
      mbov_tree <- function(mydata,metadata){
      ad_gg <- ggtree(mydata)
      ad_gg[["data"]]$label <- gsub(".annot","",ad_gg$data$label)
      
      
      ad_gg <- ad_gg %<+% isolate_dat +
        ggtree::vexpand(.1, -1)
      row_id <- subset(ad_gg[["data"]], isTip == TRUE)$label
      ad_gg_onlytip <- as.data.frame(subset(ad_gg[["data"]], isTip == TRUE)[,metadata])
      rownames(ad_gg_onlytip) <- row_id
      print(head(ad_gg_onlytip))
      
    
      t1 <- gheatmap(ad_gg, ad_gg_onlytip, width = 0.3, colnames = FALSE) +
        scale_fill_manual(values = hcl.colors(length(unique(ad_gg_onlytip[,metadata])),palette = "Zissou 1"), name = metadata)
      
      t1_scaled <- t1 + new_scale_fill()
      t2 <- gheatmap(t1_scaled,accessory_transpose, offset = 3, colnames_angle=45,hjust = 1,colnames_offset_y = -2.5,font.size = 3) +  
        scale_fill_manual(values = c("gray75","darkblue"), name = "Presence/Absence")
      Sys.sleep(1) #gives program the time to make the figure
      t2
      
      }
      
      plot(mbov_tree(accessory_dendro,colnames(isolate_dat)[i]))
    
  
} 
   }
  } else{
  ad_gg <- ggtree(accessory_dendro)
  ad_gg[["data"]]$label <- gsub(".annot","",ad_gg$data$label)
  
  t1 <- gheatmap(ad_gg,accessory_transpose, offset = 0.5, colnames = FALSE) +  
    scale_fill_manual(values = c("gray75","darkblue"), name = "Presence/Absence")
  plot(t1) 
  
  mbov_tree <- function(mydata,metadata){
    ad_gg <- ggtree(mydata)
    ad_gg[["data"]]$label <- gsub(".annot","",ad_gg$data$label)
    
    
    ad_gg <- ad_gg %<+% isolate_dat
    
    t2 <- gheatmap(ad_gg,accessory_transpose, offset = 0.5, colnames_angle=-45) +  
      scale_fill_manual(values = c("gray75","darkblue"), name = "Presence/Absence")
    Sys.sleep(1) #gives program the time to make the figure
    t2
    
  }
}
dev.off()

### Finish Code and update mbovpan 