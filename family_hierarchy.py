
import os
import os.path as ospath
import pandas as pd

#creates the csv files that associates each family to a list of the genuses inside it
file = "Plantae_Genus.txt"
with open(file, "r") as f:
    seq = f.read()

undesirables = ["_", "\n", "\r"]
for undesirable in undesirables:
    seq= seq.replace(undesirable, " ")

hierarchy = pd.DataFrame(columns = ["Family", "genuses"])
separated = seq.split()

current_family = "None"
current_genuses = []


for i in range(0,len(separated)):
    if separated[i] == "Family:" :
        hierarchy.loc[i] = current_family, current_genuses
        current_genuses = []
        current_family = separated[i+1]
    if separated[i] == "Genus:":
        current_genuses.append(separated[i+1])

hierarchy.index = range(0,len(hierarchy))
hierarchy.drop(0)

hierarchy.to_csv(os.getcwd()+ "/hierarchy_report.csv")





