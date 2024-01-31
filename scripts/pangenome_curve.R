#install.packages("ggplot2")
#install.packages("dplyr")

library(ggplot2)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

gene_pres_abs <- read.csv(args[1], header = TRUE, stringsAsFactors = FALSE, row.names = "Gene")
genes <- rownames(gene_pres_abs)
accessory_pa <- gene_pres_abs %>% select(3:(ncol(gene_pres_abs)))

accessory_pa[!(accessory_pa=="")] <- 1
accessory_pa[accessory_pa==""] <- 0 

task <- function(x){sum(as.numeric(as.character(x)))}

#compute number of conserved genes
genome_num <- c()
total_gene_num <- c()
conserved_gene_num <- c()

num_iter = 100

if(ncol(accessory_pa) < 100){
  num_iter = ncol(accessory_pa)
}



#outer loop is for changing the genome size
for(i in 2:100){
  genome_num <- append(genome_num,i)
  conserved_gene_avg = 0
  total_gene_avg = 0
  
  #inner loop is for going through ten iterations to ultimately average the results
  for(j in 1:5){
    rand_columns <- sample(colnames(accessory_pa))[1:i]
    prab_subset <- accessory_pa %>% select(all_of(rand_columns))
    prab_subset$presence_num <- apply(prab_subset,1,task)
    conserved_gene_avg = conserved_gene_avg + nrow(prab_subset %>% filter(presence_num == i))
    total_gene_avg = total_gene_avg + nrow(prab_subset %>% filter(presence_num > 0))
  }
  
  total_gene_num <- append(total_gene_num, total_gene_avg/5)
  conserved_gene_num <- append(conserved_gene_num, conserved_gene_avg/5)
  
}

plot_data <- data.frame(genome = genome_num, conserved = conserved_gene_num, total = total_gene_num)

type_vector <- append(rep("Conserved genes",99),rep("Total genes",99))
data_vector <- append(plot_data$conserved,plot_data$total)
genome_num_vector <- append(plot_data$genome,plot_data$genome)

plot_data <- data.frame(genome = genome_num_vector, group = type_vector, value = data_vector)

pdf("pangenome_curve.pdf")

ggplot(data=plot_data, aes(x=genome, y=value, group=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group)) +
  xlab("No. of genomes") +
  ylab("No. of genes") +
  theme_minimal() + 
  theme(
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15),
    axis.line = element_line(color = "black")
  ) +
  scale_color_discrete(name = "Key")

dev.off()
