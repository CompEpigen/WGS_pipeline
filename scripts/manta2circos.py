# Convert VCF of structural variants to format expected by circos to draw arcs between SVs


import os
import argparse
import pandas as pd
import vcfpy
parser = argparse.ArgumentParser()
parser.add_argument('-i', type = str, help='VCF of structural variants')
parser.add_argument('-o', type = str, help='Link file for circos')
args = parser.parse_args()

#args.i = "/home/e840r/Documents/WGS/viz/WGS_data/filtered_VCF/C010-AML-15KM18875.vcf"
#args.o = "/home/e840r/Documents/WGS/viz/circos/tmp.txt"

reader = vcfpy.Reader.from_path(args.i)

color_mapping = {"BND":"green","TRA":"green","DEL":"blue","DUP":"red","INV":"black"}

with open(args.o,"w") as outfile:
    for record in reader:
        chr = record.CHROM
        pos = record.POS
        if record.INFO["SVTYPE"] in ["BND","TRA"]:
            chr2 = record.ALT[0].mate_chrom
            pos2 = record.ALT[0].mate_pos
        else:
            chr2 = chr
            pos2 = record.INFO["END"]
        
        low_mappability = "MAPPABILITY" in record.FILTER
        if record.INFO["SVTYPE"] == "DEL":
            color="DEL"
        elif record.INFO["SVTYPE"] == "DUP":
            color = "DUP"
        elif record.INFO["SVTYPE"] =="INV":
            color = "INV"
        else:  #BND
            if record.CHROM != record.ALT[0].mate_chrom:
                color = "TRA"
            else:
                if record.POS > record.ALT[0].mate_pos:
                    continue
                if record.ALT[0].orientation =="-" and record.ALT[0].mate_orientation == "+":
                    color = "DEL"
                elif record.ALT[0].orientation =="+" and record.ALT[0].mate_orientation == "-":
                    color = "DUP"
                else:
                    color = "INV"
        color = color_mapping[color]
        #if low_mappability:
        #    thickness = "1.0p"
        #    color = "vl"+color
        #else:
        thickness = "3p"
        tmp = outfile.write(chr+"\t"+str(pos)+"\t"+str(pos) + "\t" + chr2+"\t"+str(pos2)+"\t"+str(pos2)+"\t" + "color="+color+",thickness="+thickness+"\n")
            
        

