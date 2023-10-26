# Return a filtered PRAB matrix
import pandas as pd
import sys

lineage_csv = pd.read_json(sys.argv[1])
spotyping_csv = pd.read_csv(sys.argv[2])

print("files read successfully")
print(lineage_csv)
print(spotyping_csv)


