# Convert copy number ratios to format expected by circos


import os
import argparse
import pandas as pd
parser = argparse.ArgumentParser()
parser.add_argument('-i', type = str, help='Control-Freec output')
parser.add_argument('-o', type = str, help='Input for circos')
parser.add_argument('--cnv', type = str, help='CNVs identified by Control-FREEC.')
args = parser.parse_args()

#args.i = "/home/e840r/Documents/WGS/viz/circos/tumor01_C010-AML-16PB5058_merged.mdup.bam_ratio.txt"
#args.o = "/home/e840r/Documents/WGS/viz/circos/16PB5058cn.txt"

CNVs=[]
with open(args.cnv,"r") as infile:
    for line in infile:
        linesplit = line.rstrip("\n").split("\t")
        CNVs.append((linesplit[0],int(linesplit[1]),int(linesplit[2]),linesplit[4]))

def bin2CNV(chr,start):
    result="normal"
    for CNV in CNVs:
        if str(chr)==str(CNV[0]) and start>=CNV[1] and start<=CNV[2]:
            result = CNV[3]
    return result


df = pd.read_csv(args.i,sep="\t")
binsize = df.loc[1,"Start"] - df.loc[0,"Start"]


with open(args.o,"w") as outfile:
    for i in range(df.shape[0]):
        if df.loc[i,"Ratio"] >0:
            cn_state = bin2CNV(df.loc[i,"Chromosome"],int(df.loc[i,"Start"]))
            if cn_state=="gain":
                color = "red"
                stroke_color = "dred"
            elif cn_state=="loss":
                color = "blue"
                stroke_color = "dblue"
            else:
                color = "grey"
                stroke_color = "black"

            tmp= outfile.write(str(df.loc[i,"Chromosome"]) + "\t" +str(df.loc[i,"Start"]) + "\t"+str(df.loc[i,"Start"]+binsize) 
             +"\t" + str(2*df.loc[i,"Ratio"])+"\tcolor=" +color+",stroke_color=" +stroke_color+"\n")

