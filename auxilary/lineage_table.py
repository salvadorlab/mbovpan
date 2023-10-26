# Return a filtered PRAB matrix
import pandas as pd
import sys
import json
import glob

#work with all the json files
jsons = glob.glob("*.results.json")

for j in jsons:
    print(j)
    f = open(j)
    data = json.load(f)

    lineage_str = "{},{},{}".format(data["id"],data["main_lin"],data["sublin"])
    print(lineage_str)
f.close()


# lineage_csv = pd.read_json(sys.argv[1], lines=True, encoding='utf-8-sig')
spotyping_csv = pd.read_csv(sys.argv[1])

print("files read successfully")
#print(lineage_csv)
print(spotyping_csv)


