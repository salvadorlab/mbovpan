library(ggplot2)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

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

pa_transpose <- t(data.matrix(accessory_pa))
isolate_ids <- gsub("_trimmed_R1.scaffold.annot","",row.names(pa_transpose))
prab_pca <- prcomp(pa_transpose)
variance <- (prab_pca$sdev)^2
loadings <- prab_pca$rotation
scores <- prab_pca$x 

scores4 <- as.data.frame(scores[,1:4])
scores4$label <- isolate_ids
isolate_dat <- read.csv(args[1], stringsAsFactors = FALSE)
scores4 <- scores4 %>% left_join(isolate_ids,by = c("label" = "Name"))

# Add code to make a PCA for each pairwise comp. 

pdf("pca_figures.pdf")

if(length(args[1]) != 0){  
  for(k in 1:length(colnames(isolate_dat))){
    if(colnames(isolate_dat)[k] == "Name"){
      print("identified Name column")
      next
    }
    else{  
      for(i in 1:4){
        for(j in 1:4){
          if(i > j){

          scatterplot <- function(data_used, x.variable, y.variable, fill.variable) {
          ggplot(data_used,aes_(x=x.variable,y=y.variable)) +
            geom_point(aes_(color = as.factor(fill.variable) )) +
            theme_minimal() + 
            ggtitle("M. bovis pangenome (15% to 99% PRAB)") + 
            xlab(names(scores4)[i]) +
            ylab(names(scores4)[j]) + 
            labs(color = colnames(isolate_dat)[k])
          }
          
          plot(scatterplot(scores4, colnames(scores4)[i], colnames(scores4)[i], names(scores4[,colnames(isolate_dat)[k]])[1]))
          
          } else {
            next
          }
        }
      }
    }
  }
  
} 