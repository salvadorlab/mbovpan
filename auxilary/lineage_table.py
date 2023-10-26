# Return a filtered PRAB matrix
import pandas as pd
import sys
import json

lineage_csv = open(sys.argv[1])
data = json.load(lineage_csv)

for i in data['id']:
    print(i)


# lineage_csv = pd.read_json(sys.argv[1], lines=True, encoding='utf-8-sig')
spotyping_csv = pd.read_csv(sys.argv[2])

print("files read successfully")
print(lineage_csv)
print(spotyping_csv)


