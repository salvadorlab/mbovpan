# Return a filtered PRAB matrix
import pandas as pd
import sys

mbov_prab = pd.read_csv(sys.argv[1],low_memory=False)
mbov_virulence = pd.read_csv(sys.argv[2])

print("files read successfully")

# extract virulence genes as a list of unique elements
virulent_gene_list = list(set(mbov_virulence.loc[:,"Accronyme"].to_list()))
print(virulent_gene_list)

#return only rows that have the gene name matching in the annotation!
filtered_list = [] 
for i in range(len(mbov_prab.index)):
    print(mbov_prab.loc[i,"Gene"])
    print(mbov_prab.loc[i,"Gene"].split("_"))
    print(mbov_prab.loc[i,"Gene"].split("_")[0])
    if str(mbov_prab.loc[i,"Gene"]).split("_")[0] in '\t'.join(virulent_gene_list):
        filtered_list.append(i)
        print("part of the list")
#filter the rows with the passing indicies
mbov_prab.iloc[filtered_list].to_csv("mbov_virulent_prab.csv",index=False)
print("prab created successfully")