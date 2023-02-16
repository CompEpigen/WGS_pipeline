

import os
import argparse
import vcfpy
from varcode import Variant
from pyensembl import ensembl_grch37, ensembl_grch38
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-i', type = str, help='Input VCF file')
parser.add_argument('-o', type = str, help='Output tsv file')
args = parser.parse_args()



genes_CH = ["DNMT3A","TET2","ASXL1"] # genes often mutated in clonal hematopoiesis. Some somatic mutations in those genes might be found in gnomAD.

ensembl_ref= ensembl_grch37
reader = vcfpy.Reader.from_path(args.i)
results={"chr":[],"pos":[],"ref":[],"alt":[],"gene":[],"effect":[],"VAF":[],"filters":[]}
results_filtered={"chr":[],"pos":[],"ref":[],"alt":[],"gene":[],"effect":[],"VAF":[],"filters":[]}
for record in reader:
    chr = str(record.CHROM)
    pos = str(record.POS)
    ref = str(record.REF)
    alt = str(record.ALT[0].value)
    VAF_normal = float(record.calls[0].data.get("AF")[0])
    VAF_tumor = float(record.calls[1].data.get("AF")[0])
    variant = Variant(contig=chr,start=int(pos),ref=ref,alt=alt,ensembl=ensembl_ref)
    effects = variant.effects()
    topPriorityEffect = effects.top_priority_effect()
    nonsynonymous_variant = (topPriorityEffect.gene_name is not None) and (not topPriorityEffect.short_description in ["silent","intronic","3' UTR","5' UTR","incomplete"])
    if topPriorityEffect.gene_name is not None:
        gene_name = topPriorityEffect.gene_name
    else:
        gene_name = "intergenic"
    if topPriorityEffect.short_description is not None:
        description = topPriorityEffect.short_description
    else:
        description = ""
    #TODO: filter intronic

    pop_freq = 0
    if "AF" in record.INFO: pop_freq = float(record.INFO["AF"][0])
    AC=0
    if "AC" in record.INFO: AC= int(record.INFO["AC"][0])
    filters = record.FILTER
    if "PASS" in filters: filters.remove("PASS")
    variant_id = chr+":"+pos+":"+ref+":"+alt
    if description in ["intronic","5' UTR","3' UTR","non-coding-transcript","incomplete","silent"]: filters.append("no-effect")
    if VAF_normal>0.1 or VAF_tumor < VAF_normal*2:
        filters.append("germline")
    filters_str = "_".join(filters)
    results["chr"].append(chr)
    results["pos"].append(pos)
    results["ref"].append(ref)
    results["alt"].append(alt)
    results["gene"].append(gene_name)
    results["effect"].append(description)
    results["VAF"].append(VAF_tumor)
    results["filters"].append(filters_str)
    if len(filters)==0:
        results_filtered["chr"].append(chr)
        results_filtered["pos"].append(pos)
        results_filtered["ref"].append(ref)
        results_filtered["alt"].append(alt)
        results_filtered["gene"].append(gene_name)
        results_filtered["effect"].append(description)
        results_filtered["VAF"].append(VAF_tumor)
        results_filtered["filters"].append(filters_str)

def chr_to_int(s):
    res =[]
    for x in s:
        if x.isdigit():
            res.append(int(x))
        else:
            if x=="X":res.append(23)
            else: res.append(24)
    return res
df = pd.DataFrame(results)
df.sort_values(by="chr",key=chr_to_int,inplace=True)
df.to_csv(args.o,sep="\t",index=None,header=True)

df_filtered = pd.DataFrame(results_filtered)
df_filtered.sort_values(by="chr",key=chr_to_int,inplace=True)
df_filtered.to_csv(args.o[:-4]+"_filtered.tsv",sep="\t",index=None,header=True)