# Return a filtered PRAB matrix
import pandas as pd
import sys

mbov_prab = pd.read_csv(sys.argv[1],low_memory=False)
mbov_virulence = pd.read_csv(sys.argv[2])

print("files read successfully")

# extract virulence genes as a list of unique elements
virulent_gene_list = list(set(mbov_virulence.loc[:,"Accronyme"].to_list()))

#return only rows that have the gene name matching in the annotation!
filtered_list = [] 
for i in range(len(mbov_prab.index)):
    #print(mbov_prab.loc[i,"Gene"])
    #print(mbov_prab.loc[i,"Gene"].split("_"))
    #print(mbov_prab.loc[i,"Gene"].split("_")[0])
    for gene in mbov_prab.loc[i,"Gene"].split("~~~"):
        if gene.split("_")[0] in virulent_gene_list:
            filtered_list.append(i)
            mbov_prab.loc[i,"Gene"] = gene.split("_")[0]
            print("part of the list")
#filter the rows with the passing indicies
mbov_prab.iloc[filtered_list].drop_duplicates().to_csv("mbov_virulent_prab.csv",index=False)
print(mbov_prab.iloc[filtered_list])
print("prab created successfully")