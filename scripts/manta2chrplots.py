# Takes as input a VCF of SVs (tested on manta output) and produces a tsv file containing the following columns:
# chr, start, end, length, type, chr2
# This is used for plotting the SVs.

import os
import argparse
import vcfpy

parser = argparse.ArgumentParser()
parser.add_argument('-i', type = str, help='Input VCF file')
parser.add_argument('-o', type = str, help='Output tsv file')
args = parser.parse_args()

#args.i="/home/e840r/Documents/WGS/viz/WGS_data/filtered_VCF/C010-AML-15KM12995.vcf"
#args.o = "tmp.tsv"

reader = vcfpy.Reader.from_path(args.i)


c=0
with open(args.o,"w") as outfile:
    tmp = outfile.write("\t".join(["chr","start","end","length","type","chr2","color"]) + "\n")
    for record in reader:
        low_mappability = "MAPPABILITY" in record.FILTER
        chr = record.CHROM
        start = record.POS
        #end = record.INFO["END"]
        type = record.INFO["SVTYPE"]
        if type in ["DEL","DUP","INV","INS"]:
            end = record.INFO["END"]
            chr2 = chr
            length = end-start
        else:
            type = "BND"
            chr2 = record.ALT[0].mate_chrom
            if chr==chr2:
                end = record.ALT[0].mate_pos
                if record.POS > record.ALT[0].mate_pos:
                    continue
                length = record.ALT[0].mate_pos - record.POS
                if record.ALT[0].orientation =="-" and record.ALT[0].mate_orientation == "+":
                    type = "DEL"
                elif record.ALT[0].orientation =="+" and record.ALT[0].mate_orientation == "-":
                    type = "DUP"
                else:
                    type = "INV"
            else:
                end = start
                length = 1000000
                #continue
        #if low_mappability:
        #    color = "grey"
        #else:
        if type=="DEL": color = "blue"
        elif type=="DUP": color = "red"
        elif type=="INV": color = "black"
        elif type =="ITX": color = "black"
        else: color = "green"
        tmp = outfile.write("\t".join([str(chr),str(start),str(end),str(length),str(type),str(chr2),color]) + "\n")