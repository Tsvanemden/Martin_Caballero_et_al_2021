import os
import pandas as pd

source = "testA.txt"
dest = "testD.txt"

SraRunTable = pd.read_csv("SraRunTable.txt", sep = "\t")

source = SraRunTable["sraName"]
dest = SraRunTable["sampleName"]

for i, r in source.items():  
    oldname = "sraFiles/" + source[i] + ".fastq.gz"
    newname = "sraFiles/" + dest[i] + ".fastq.gz"
    os.rename(oldname, newname)
    print("renamed " + oldname + " to " + newname)