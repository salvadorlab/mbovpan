# Return a filtered PRAB matrix
import json
import glob

#work with all the json files and log files
#sort to make sure they are operating on the same index
print("welcome!")
spotypes = sorted(glob.glob("*.log"))

for i in range(len(spotypes)):
    spo_file = open(spotypes[i])


    #data of the spoligotype
    first_line = spo_file.readline().strip('\n').split()

    #good to notify the user if a spotyping pattern is not M. bovis
    #this code checks the final 5 spotype spacers and checks that they
    #are '00000', which all M. bovis possess
    if(first_line[1][len(first_line[1]) - 5:len(first_line[1])] == "00000"):
        print("0")
    else:
        print("1")

    spo_file.close()



