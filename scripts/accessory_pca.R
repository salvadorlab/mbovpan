library(ggplot2)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

gene_pres_abs <- read.csv(args[1], header = TRUE, stringsAsFactors = FALSE, row.names = "Gene")
accessory_pa <- gene_pres_abs %>% select(14:(ncol(gene_pres_abs)))

accessory_pa[!(accessory_pa=="")] <- 1
accessory_pa[accessory_pa==""] <- 0 

num_col <- ncol(accessory_pa)
task <- function(x){sum(as.numeric(as.character(x)))}
accessory_pa$pr <- apply(accessory_pa,1,task)
accessory_pa$perc_pr <- accessory_pa$pr/num_col
accessory_pa <- accessory_pa %>% filter(perc_pr >= 0.15 & perc_pr <= 0.99) %>% select(-c("pr","perc_pr"))

print("percentages are calculated")

pa_transpose <- t(data.matrix(accessory_pa))
isolate_ids <- row.names(pa_transpose)

#remove random '1' at the end of the string
isolate_ids <- substr(isolate_ids, 1, nchar(isolate_ids)-1)
print(isolate_ids)
prab_pca <- prcomp(pa_transpose)
variance <- (prab_pca$sdev)^2
loadings <- prab_pca$rotation
scores <- prab_pca$x 

scores4 <- as.data.frame(scores[,1:4])
scores4$label <- isolate_ids
isolate_dat <- read.csv(args[2], stringsAsFactors = FALSE, check.names = FALSE)
scores4 <- scores4 %>% left_join(isolate_dat,by = c("label" = "Name"))
head(scores4)

# Add code to make a PCA for each pairwise comp. 

pdf("pca_figures.pdf")

if(length(args[1]) != 0){  
  for(k in 1:length(colnames(isolate_dat))){
    if(colnames(isolate_dat)[k] == "Name" || colnames(isolate_dat)[k] == "Index"){
      print("identified Name column")
      next
    }
    else{  
      for(i in 1:4){
        for(j in 1:4){
          if(i > j){

          scatterplot <- function(data_used, x.variable, y.variable, fill.variable) {
          ggplot(data_used,aes(x=x.variable,y=y.variable)) +
            geom_point(aes(color = fill.variable)) +
            theme_minimal() + 
            xlab(paste(names(scores4)[i],"(",signif(variance[i],3),"% )")) +
            ylab(paste(names(scores4)[j],"(",signif(variance[j],3),"% )")) + 
            labs(color=colnames(isolate_dat)[k])  +
              theme(
                axis.text = element_text(size = 13),
                axis.title = element_text(size = 15),
                axis.line = element_line(color = "black")
              )
          }
          
          plot(scatterplot(scores4, scores4[,i], scores4[,j], scores4[,colnames(isolate_dat)[k]]))
          
          } else {
            next
          }
        }
      }
    }
  }
  
} 