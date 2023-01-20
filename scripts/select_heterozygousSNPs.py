import os
import argparse
import vcfpy

parser = argparse.ArgumentParser()
parser.add_argument('-i', type = str, help='Input VCF file')
parser.add_argument('--chr', type = str, help='Chromosome to select.')
parser.add_argument('-o', type = str, help='Output bed file')
args = parser.parse_args()

reader = vcfpy.Reader.from_path(args.i)

with open(args.o,"w") as outfile:
    for record in reader:
        chr = record.CHROM
        if args.chr is not None and chr!=args.chr: continue
        pos = record.POS
        #print(record)
        #print(record.calls[0])
        AD = record.calls[0].data["AD"]
        if ( (AD[0]+AD[1]) > 10) and ( AD[0] / (AD[0]+AD[1]) > 0.2) and (AD[0] / (AD[0]+AD[1]) < 0.8):
           tmp = outfile.write(chr+"\t"+str(pos-1)+"\t"+str(pos)+"\n")