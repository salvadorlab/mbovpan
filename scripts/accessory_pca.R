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
isolate_ids <- data.frame(Samples = gsub("_trimmed_R1.scaffold.annot","",row.names(pa_transpose)))
prab_pca <- prcomp(pa_transpose)
variance <- (prab_pca$sdev)^2
loadings <- prab_pca$rotation
rownames(loadings) <- colnames(pa_transpose)
scores <- prab_pca$x 

scores4 <- as.data.frame(scores[,1:4])
head(scores4)

# Add code to make a PCA for each pairwise comp. 

pdf("pca_figures.pdf")

if(length(args[1]) != 0){
  print(args[1])
  isolate_dat <- read.csv(args[1], stringsAsFactors = FALSE)
  print(colnames(isolate_dat))
  scores4$Name <- gsub("1.annot","",rownames(scores4))
  print(scores4$Name)
  scores4 <- scores4 %>% left_join(isolate_dat,by = "Name")
  print(scores4)
  
  test <- ggplot(scores4,aes(x=PC1,y=PC2)) +
    geom_point(aes(colour = factor(species))) +
    theme_minimal() + 
    ggtitle("M. bovis pangenome (15% to 99% PRAB)")
  
  plot(test)
  
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
          
          plot(scatterplot(scores4, scores4[,i], scores4[,j], scores4[,colnames(isolate_dat)[k]]))
          
          } else {
            next
          }
        }
      }
    }
  }
  
} else{
  for(i in 1:4){
    for(j in 1:4){
      if(i > j){
      a <- ggplot(as.data.frame(scores4),aes(x=colnames(scores4)[i],y=colnames(scores4)[j])) +
        geom_point() +
        theme_minimal() + 
        ggtitle("M. bovis pangenome (15% to 99% PRAB)")
      plot(a)
      } else {
        next 
      }
    }
  }
}