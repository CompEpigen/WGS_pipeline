import os
import argparse
import vcfpy

parser = argparse.ArgumentParser()
parser.add_argument('-i', type = str, help='Input VCF file')
parser.add_argument('-o', type = str, help='Output tsv file')
args = parser.parse_args()

#args.i="/home/e840r/Documents/AML/out/SV/manta/C010-AML-15PB19457/C010-AML-15PB19457_SV_filtered.vcf"
#args.o = "/home/e840r/Documents/WGS/SVs/shatterseek/15PB19457_shatterseekSV.tsv"

reader = vcfpy.Reader.from_path(args.i)


def get_breakpoint_info(record):
    chr = record.CHROM
    pos = record.POS
    if record.INFO["SVTYPE"] in ["DEL","DUP","INS"]:
        chr2 = chr
        pos2 = record.INFO["END"]
        if record.INFO["SVTYPE"]=="DEL":
            orientation = "-"
            orientation2 = "+"
        elif record.INFO["SVTYPE"]=="DUP":
            orientation = "+"
            orientation2 = "-"
        else:
            print("warning!!")
            orientation = "-"
            orientation2 = "+"
    else:
        chr2 = record.ALT[0].mate_chrom
        pos2 = record.ALT[0].mate_pos
        orientation = record.ALT[0].orientation
        orientation2 = record.ALT[0].mate_orientation
    return ( (chr,pos,orientation) , (chr2,pos2,orientation2) )

def invert_strand(strand): # shatterseek seems to use a different convention than vcfpy for the strand
    if strand=="+": return "-"
    else: return "+"

def chr2int(chr):
    chr = chr.lstrip("chr")
    if chr.isdigit():
        return int(chr)
    elif chr=="X": return 23
    else: return 24

with open(args.o,"w") as outfile:
    tmp = outfile.write("\t".join(["chr1","pos1","chr2","pos2","svclass","strand1","strand2"]) + "\n")
    for record in reader:
        ( (chr1,pos1,orientation1) , (chr2,pos2,orientation2) ) = get_breakpoint_info(record)
        if chr1!=chr2:
            svclass="TRA"
        else:
            if orientation1=="-" and orientation2=="+":
                svclass = "DEL"
            elif orientation1=="+" and orientation2=="-":
                svclass = "DUP"
            elif orientation1=="+" and orientation2=="+":
                svclass= "t2tINV"
            else:
                svclass = "h2hINV"
        if chr2int(chr1)< chr2int(chr2) or (chr1==chr2 and int(pos1)<=int(pos2)):
            if chr1!="Y" and chr2!="Y":
                tmp = outfile.write("\t".join([str(chr1),str(pos1),str(chr2),str(pos2),str(svclass),invert_strand(orientation1),invert_strand(orientation2)]) + "\n")