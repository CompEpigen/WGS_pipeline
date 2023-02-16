import pandas as pd 
import pysam

chromosomes = [str(x) for x in range(1,23)] + ["X","Y"]
chromosomes_noY = [str(x) for x in range(1,23)] + ["X"]


df_samplesheet = pd.read_csv(config["samplesheet"],sep=",",index_col=0)

if not "sex" in df_samplesheet.columns:
	df_samplesheet["sex"] = ["XX" for x in df_samplesheet.index] # default: XX
else:
	df_samplesheet["sex"] = ["XX" if x in ["F","female","Female","XX"] else "XY" for x in df_samplesheet["sex"]]

samples_control=[]
samples_nocontrol=[]
if config["use_control"]:
	for sample in df_samplesheet.index:
		if "path_bam_control" in df_samplesheet.columns and df_samplesheet.loc[sample,"path_bam_control"]!="":
			samples_control.append(sample)
		else:
			samples_nocontrol.append(sample)
else:
	for sample in df_samplesheet.index:
		samples_nocontrol.append(sample)
		if "path_bam_control" in df_samplesheet.columns and df_samplesheet.loc[sample,"path_bam_control"]!="":
			samples_nocontrol.append(sample+"_control")


# Mutect2 needs the sample name of the control in the header of the bam file.
if "path_bam_control" in df_samplesheet.columns:
	control_samplenames=[]
	for sample in df_samplesheet.index:
		if df_samplesheet.loc[sample,"path_bam_control"]!="":
			pysam_file = pysam.AlignmentFile(df_samplesheet.loc[sample,"path_bam_control"],"rb")
			control_samplenames.append(pysam_file.header["RG"][0]["SM"])
		else:
			control_samplenames.append("")
	df_samplesheet["control_samplename"] = control_samplenames


samples_RNA=[]
if "path_BAM_RNA" in df_samplesheet.columns:
	for sample in df_samplesheet.index:
		if df_samplesheet.loc[sample,"path_BAM_RNA"]!="":
			samples_RNA.append(sample)



rule all:
	input:
		expand(config["working_dir"]+"/out_nocontrol/plots/chrplots_png/{sample}/{sample}_chr{chr}.png",sample=samples_nocontrol,chr=chromosomes),
		expand(config["working_dir"]+"/out_nocontrol/plots/circos/circos_{sample}.svg",sample=samples_nocontrol),
		expand(config["working_dir"]+"/out_nocontrol/SNV/mutect2_genePanel/{sample}/{sample}.tsv",sample=samples_nocontrol),
		expand(config["working_dir"]+"/out_control/plots/chrplots_png/{sample}/{sample}_chr{chr}.png",sample=samples_control,chr=chromosomes),
		expand(config["working_dir"]+"/out_control/plots/circos/circos_{sample}.svg",sample=samples_control),
		expand(config["working_dir"]+"/out_control/SNV/mutect2_wholegenome/{sample}/{sample}.tsv",sample=samples_control),
		expand(config["working_dir"]+"/out_control/SNV/freebayes_RNA/{sample}/{sample}.vcf.gz",sample=samples_control)


# Filter the BAM files by removing reads which align to several locations
def sample2bamfile(wcs):
	sample = wcs["sample"]
	if sample.endswith("control"):
		sample = sample[:-8]
		return df_samplesheet.loc[sample,"path_bam_control"]
	else:
		return df_samplesheet.loc[sample,"path_bam"]

rule filter_BAM:
	params:
		threads="8",
		runtime="20:00",
		memory="2G"
	input:
		sample2bamfile
	output:
		BAM = config["working_dir"]+"/BAM/{sample}.bam",
		BAM_i = config["working_dir"]+"/BAM/{sample}.bam.bai"
	shell:
		"sambamba view -t 8 -F 'mapping_quality >5' {input} -f bam -o {output.BAM}"
		#"sambamba view -t 8 -F '[XA]==null and mapping_quality >35' {input} -f bam -o {output.BAM}"
"""

rule symlink_BAM:
	params:
		threads="8",
		runtime="20:00",
		memory="2G"
	input:
		sample2bamfile
	output:
		BAM = config["working_dir"]+"/BAM/{sample}.bam",
		BAM_i = config["working_dir"]+"/BAM/{sample}.bam.bai"
	shell:
		"ln -s {input} {output.BAM}; ln -s {input}.bai {output.BAM_i}"

"""


def sample2bamfile_RNA(wcs):
	return df_samplesheet.loc[wcs["sample"],"path_bam_RNA"]

rule symmlink_BAM_RNA:
	params:
		threads="1",
		runtime="20",
		memory="1000"
	input:
		sample2bamfile_RNA
	output:
		BAM = config["working_dir"]+"/BAM_RNA/{sample}.bam",
		BAM_i = config["working_dir"]+"/BAM_RNA/{sample}.bam.bai"
	shell:
		"ln -s {input} {output.BAM}; ln -s {input}.bai {output.BAM_i}"


include: "rules/FREEC.smk"
include: "rules/manta.smk"
include: "rules/SNV.smk"
include: "rules/plots.smk"


rule run_shatterseek:
	params:
		threads="1",
		runtime="40",
		memory="16000"
	input:
		CNA=config["working_dir"]+"/out/CNA/FREEC/{sample}/{sample}.bam_ratio.txt",
		SV=config["working_dir"]+"/out/SV/manta/{sample}/{sample}_SV_filtered.vcf"
	output:
		CNA=config["working_dir"]+"/out/CNA/FREEC/{sample}/{sample}_shatterseekCNV",
		SV=config["working_dir"]+"/out/SV/manta/{sample}/{sample}_shatterseekSV",
		shatterseek=config["working_dir"]+"/out/chromothripsis/shatterseek/{sample}/{sample}_shatterseek.tsv"

	conda:
		"envs/WGS.yaml"
	shell:
		"python scripts/freec2shatterseek.py -i {input.CNA} -o {output.CNA};"\
		"python scripts/manta2shatterseek.py -i {input.SV} -o {output.SV};"\
		"Rscript scripts/shatterseek.R  {output.CNA} {output.SV} {output.shatterseek};"




"""
# Create links to the BAM files, so that the other functions do not need to take the directory structure of the original BAM files into account.
rule create_link_BAM:
	params:
		threads="1",
		runtime="20",
		memory="1000"
	output:
		BAM = config["working_dir"]+"/BAM/{sample}.bam",
		BAM_i = config["working_dir"]+"/BAM/{sample}.bam.bai"
	run:
		bam_path = config["BAM_template"].replace("_SAMPLE_",wildcards.sample)
		shell("ln -s {bam_path} {output.BAM}")
		shell("ln -s {bam_path}.bai {output.BAM_i}")
"""


"""
rule filter_BAM_template:
	params:
		threads="4",
		runtime="20:00",
		memory="1000"
	output:
		BAM = config["working_dir"]+"/BAM/{sample}.bam",
		BAM_i = config["working_dir"]+"/BAM/{sample}.bam.bai"
	run:
		bam_path = config["BAM_template"].replace("_SAMPLE_",wildcards.sample)
		shell("sambamba view -t 4 -F '[XA]==null and mapping_quality >35' {bam_path} -f bam -o {output.BAM}")
		shell("samtools index {output.BAM}")
"""