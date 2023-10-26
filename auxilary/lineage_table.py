# Return a filtered PRAB matrix
import pandas as pd
import sys
import json
import glob

#work with all the json files
jsons = glob.glob("*.results.json")

print("id,main lineage,sub lineage")
for j in jsons:
    print(j)
    f = open(j)
    data = json.load(f)

    lineage_str = "{},{},{}".format(data["id"],data["main_lin"],data["sublin"])
    print(lineage_str)
f.close()


#work with the spoligotyping files
spotyping_csv = pd.read_table(sys.argv[1],header=None)

for i in spotyping_csv.index:
    spotyping_csv[0][i] = spotyping_csv[0][i].split("&")


print("files read successfully")
#print(lineage_csv)
print(spotyping_csv)


