import pandas as pd
import sys 

sys.argv[1]

df = pd.read_csv("mb.out.csv")
#print("will parse blast results")
#print(df)

# create a column of length % based on the subject (should be more than 75%)
df["percent_length"] = df["qlen"]/df["slen"]
filtered_df = df[df['percent_length'] >= 0.75]

#print(filtered_df)

# create a dictionary with key of subject and value of tuple (index,percent) ONLY if it is the max percent value 
# must go through list once, then filter the df by the dictionary values

index_to_keep = {}

for index, row in filtered_df.iterrows():
    if row['sseqid'] in index_to_keep:
        if index_to_keep[row['sseqid']][1] < row['percent_length']:
            index_to_keep[row['sseqid']] = (index,row['percent_length'])
        else:
            continue
    else:
        index_to_keep[row['sseqid']] = (index,row['percent_length'])
    
#print(index_to_keep)

# extract the indicies and filter dataframe
indicies = [x[0] for x in index_to_keep.values()]
#print(indicies)
#print(len(indicies))

final_df = filtered_df.loc[filtered_df.index.isin(indicies)]
#print(final_df)

# now we need to only keep the gene_presence_absence rows that show up in the filtered list and save the output as a new file to do pca on
prab = pd.read_csv(sys.argv[1])
#print(prab)

final_prab = prab[prab['Gene'].isin(final_df['qseqid'])]
#print(final_prab)
final_prab.to_csv('mbovis_filtered_cogs.csv', index=False)
