# Author: Noah Legall
# Purpose: Creation of Pangenome curve figure for mbovpan 
library(ggplot2)

no_genes_pangen <- read.table("number_of_genes_in_pan_genome.Rtab")
boxplot(no_genes_pangen, data=no_genes_pangen, main="No. of genes in the pan-genome",
        xlab="No. of genomes", ylab="No. of genes",varwidth=TRUE, ylim=c(0,max(no_genes_pangen)), outline=FALSE)
no_conserved_genes <- read.table("number_of_conserved_genes.Rtab")
boxplot(no_conserved_genes, data=no_conserved_genes, main="Number of conserved genes",
        xlab="No. of genomes", ylab="No. of genes",varwidth=TRUE, ylim=c(0,max(no_conserved_genes)), outline=FALSE)
no_genes_pangen$Type <- "total genes"
no_conserved_genes$Type <- "core genes"
pangene_curve <- rbind(no_genes_pangen,no_conserved_genes)
pangene_curve$id <- rownames(pangene_curve)
pangene_curve <- reshape(data=pangene_curve, idvar="id",
                         varying = colnames(pangene_curve)[!(colnames(pangene_curve) %in% c("id","Type"))],
                         v.name=c("value"), times=colnames(pangene_curve)[!(colnames(pangene_curve) %in% c("id","Type"))], direction="long")


# White Background
# Lines on left and bottom border
# Change legend title
pan_gg <- ggplot(pangene_curve, aes(x=time, y=value,color=Type)) + 
  geom_boxplot(position=position_dodge(0),alpha = 0.8) +
  scale_x_discrete(breaks=NULL) +
  labs(x="Number of Genomes",
       y="Number of Genes") +
  theme(axis.title = element_text(size = 12),
        panel.background = element_rect(fill = 'white'),
        panel.border = element_blank(), 
        axis.line = element_line(size = 1.5),
        legend.title = element_text(size = 12),
        legend.background = element_blank())
ggsave("pangenome_curve.png",plot = pan_gg,device="png")# Author: Noah Legall
# Purpose: Creation of Pangenome curve figure for mbovpan 
library(ggplot2)

no_genes_pangen <- read.table("number_of_genes_in_pan_genome.Rtab")
boxplot(no_genes_pangen, data=no_genes_pangen, main="No. of genes in the pan-genome",
        xlab="No. of genomes", ylab="No. of genes",varwidth=TRUE, ylim=c(0,max(no_genes_pangen)), outline=FALSE)
no_conserved_genes <- read.table("number_of_conserved_genes.Rtab")
boxplot(no_conserved_genes, data=no_conserved_genes, main="Number of conserved genes",
        xlab="No. of genomes", ylab="No. of genes",varwidth=TRUE, ylim=c(0,max(no_conserved_genes)), outline=FALSE)
no_genes_pangen$Type <- "total genes"
no_conserved_genes$Type <- "core genes"
pangene_curve <- rbind(no_genes_pangen,no_conserved_genes)
pangene_curve$id <- rownames(pangene_curve)
pangene_curve <- reshape(data=pangene_curve, idvar="id",
                         varying = colnames(pangene_curve)[!(colnames(pangene_curve) %in% c("id","Type"))],
                         v.name=c("value"), times=colnames(pangene_curve)[!(colnames(pangene_curve) %in% c("id","Type"))], direction="long")


# White Background
# Lines on left and bottom border
# Change legend title
pan_gg <- ggplot(pangene_curve, aes(x=time, y=value,color=Type)) + 
  geom_boxplot(position=position_dodge(0),alpha = 0.8) +
  scale_x_discrete(breaks=NULL) +
  labs(x="Number of Genomes",
       y="Number of Genes") +
  theme(axis.title = element_text(size = 12),
        panel.background = element_rect(fill = 'white'),
        panel.border = element_blank(), 
        axis.line = element_line(size = 1.5),
        legend.title = element_text(size = 12),
        legend.background = element_blank())
ggsave("pangenome_curve.png",plot = pan_gg,device="png")
