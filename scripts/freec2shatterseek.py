import os
import argparse
import vcfpy

parser = argparse.ArgumentParser()
parser.add_argument('-i', type = str, help='Input VCF file')
parser.add_argument('-o', type = str, help='Output tsv file')
args = parser.parse_args()


#args.i="/home/e840r/Documents/AML/out/CNV/FREEC/C010-AML-15PB19457/C010-AML-15PB19457.bam_ratio.txt"
#args.o = "/home/e840r/Documents/WGS/SVs/shatterseek/15PB19457_shatterseekCNV.tsv"


with open(args.o,"w") as outfile:
    tmp = outfile.write("\t".join(["chr","start","end","total_cn"]) + "\n")
    with open(args.i,"r") as infile:
        tmp=infile.readline()
        current_chr="1"
        current_start="1"
        current_end="10000"
        current_cn="2"
        for line in infile:
            (chr,start,ratio,medianRatio,cn) = line.rstrip("\n").split("\t")
            if chr!=current_chr or cn!=current_cn:
                if current_chr!="Y":
                    tmp = outfile.write("\t".join([current_chr,current_start,current_end,current_cn]) + "\n")
                current_chr = chr
                current_start = str(int(start))
                current_end = str(int(start)+9999)
                current_cn=cn
            else:
                current_end= str(int(start)+9999)
        #last line
        if current_chr!="Y":
            tmp = outfile.write("\t".join([current_chr,current_start,current_end,current_cn]) + "\n")
