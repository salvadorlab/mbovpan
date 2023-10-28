# Return a filtered PRAB matrix
import pandas as pd
import json
import glob

#work with all the json files and log files
#sort to make sure they are operating on the same index
jsons = sorted(glob.glob("*.results.json"))
spotypes = sorted(glob.glob("*.log"))


#variables to print 
id = []
main_lin = []
sub_lin = []
spotype_binary = []
spotype_octo = []
warning = []


print("id,main lineage,sub lineage,spoligotype (binary),spoligotype (octal), note")
for i in range(len(jsons)):
    print(i) # test the value of the index
    json_file = open(jsons[i])
    spo_file = open(spotypes[i])

    # data of the lineages
    data = json.load(json_file)
    id.append(data["id"])
    main_lin.append(data["main_lin"])
    sub_lin.append(data["sublin"])

    #data of the spoligotype
    first_line = spo_file.readline().strip('\n').split()
    spotype_binary.append(first_line[1])
    spotype_octo.append(first_line[2])

    #good to notify the user if a spotyping pattern is not M. bovis
    #this code checks the final 5 spotype spacers and checks that they
    #are '00000', which all M. bovis possess
    if(first_line[1][len(first_line) - 4:len(first_line) - 1] == "00000"):
        warning.append(" ")
    else:
        warning.append("NOT M. BOVIS")
    
    json_file.close()
    spo_file.close()


for i in range(len(id)):
    print("{},{},{},{},{},{}".format(id[i],main_lin[i],sub_lin[i],spotype_binary[i],spotype_octo[i],warning[i]))



