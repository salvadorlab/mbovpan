# Return a filtered PRAB matrix
import pandas as pd

mbov_prab = pd.read_csv("..\\auxilary\gene_presence_absence.csv")
mbov_virulence = pd.read_csv("..\\auxilary\M_bovis_virulence_genes.csv")

# extract virulence genes as a list of unique elements
virulent_gene_list = list(set(mbov_virulence.loc[:,"Accronyme"].to_list()))

#return only rows that have the gene name matching in the annotation!
filtered_list = [ i for i in range(len(mbov_prab.index)) if str(mbov_prab.loc[i,"Gene"]).split("_")[0] in '\t'.join(virulent_gene_list)]

#filter the rows with the passing indicies
mbov_prab.iloc[filtered_list].to_csv("mbov_virulent_prab.csv",index=False)