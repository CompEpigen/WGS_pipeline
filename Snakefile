chromosomes = [str(x) for x in range(1,23)] + ["X","Y"]
chromosomes_noY = [str(x) for x in range(1,23)] + ["X"]

ruleorder: run_FREEC > run_FREEC_control >  run_manta > run_manta_control > postprocess_manta > postprocess_manta_control > run_mutect2 >postprocess_mutect2 > run_mutect2_control 


if isinstance(config["samples"], str):
	samples = []
	with open(config["samples"],"r") as infile:
		for line in infile:
			samples.append(line.rstrip("\n"))
else:
	samples = config["samples"]

samples_RNA=[]
if "samples_RNA" in config:
	with open(config["samples_RNA"],"r") as infile:
		for line in infile:
			samples_RNA.append(line.rstrip("\n"))

sample_sex={}
for sample in samples:
	sample_sex[sample] = "XX"
if "metadata_file" in config.keys():
	with open(config["metadata_file"],"r") as infile:
		tmp = infile.readline()
		for line in infile:
			linesplit = line.rstrip("\n").split("\t")
			if linesplit[1] in ["F","Female","female"]:
				sample_sex[linesplit[0]] = "XX"
			else:
				sample_sex[linesplit[0]] = "XY"

if not "normal_suffix" in config:
	config["normal_suffix"] = "normal"
if not "tumor_suffix" in config:
	config["tumor_suffix"] = "tumor"
if not "normal_sample_name" in config:
	config["normal_sample_name"] = "normal"
if not "template_FREEC" in config:
	config["template_FREEC"] = "data/config_template_FREEC.txt"
if not "template_FREEC_BAF" in config:
	config["template_FREEC_BAF"] = "data/config_template_FREEC_BAF.txt"



rule all:
	input:
		expand(config["output_dir"]+"/out/plots/chrplots_png/{sample}/{sample}_chr{chr}.png",sample=samples,chr=chromosomes),
		expand(config["output_dir"]+"/out/SNV/mutect2/{sample}/{sample}.tsv",sample=samples)




# Filter the BAM files by removing reads which align to several locations
rule filter_BAM:
	params:
		threads="8",
		runtime="20:00",
		memory="2G"
	input:
		BAM = config["output_dir"]+"/BAM_unfiltered/{sample}.bam"
	output:
		BAM = config["output_dir"]+"/BAM/{sample}.bam",
		BAM_i = config["output_dir"]+"/BAM/{sample}.bam.bai"
	shell:
		"sambamba view -t 8 -F '[XA]==null and mapping_quality >35' {bam_path} -f bam -o {output.BAM}"


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
		CNA=config["output_dir"]+"/out/CNA/FREEC/{sample}/{sample}.bam_ratio.txt",
		SV=config["output_dir"]+"/out/SV/manta/{sample}/{sample}_SV_filtered.vcf"
	output:
		CNA=config["output_dir"]+"/out/CNA/FREEC/{sample}/{sample}_shatterseekCNV",
		SV=config["output_dir"]+"/out/SV/manta/{sample}/{sample}_shatterseekSV",
		shatterseek=config["output_dir"]+"/out/chromothripsis/shatterseek/{sample}/{sample}_shatterseek.tsv"

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
		BAM = config["output_dir"]+"/BAM/{sample}.bam",
		BAM_i = config["output_dir"]+"/BAM/{sample}.bam.bai"
	run:
		bam_path = config["BAM_template"].replace("_SAMPLE_",wildcards.sample)
		shell("ln -s {bam_path} {output.BAM}")
		shell("ln -s {bam_path}.bai {output.BAM_i}")
"""

rule create_link_RNA:
	params:
		threads="1",
		runtime="20",
		memory="1000"
	output:
		BAM = config["output_dir"]+"/BAM_RNA/{sample}.bam",
		BAM_i = config["output_dir"]+"/BAM_RNA/{sample}.bam.bai"
	run:
		bam_path = config["BAM_template_RNA"].replace("_SAMPLE_",wildcards.sample)
		shell("ln -s {bam_path} {output.BAM}")
		shell("ln -s {bam_path}.bai {output.BAM_i}")


"""
rule filter_BAM_template:
	params:
		threads="4",
		runtime="20:00",
		memory="1000"
	output:
		BAM = config["output_dir"]+"/BAM/{sample}.bam",
		BAM_i = config["output_dir"]+"/BAM/{sample}.bam.bai"
	run:
		bam_path = config["BAM_template"].replace("_SAMPLE_",wildcards.sample)
		shell("sambamba view -t 4 -F '[XA]==null and mapping_quality >35' {bam_path} -f bam -o {output.BAM}")
		shell("samtools index {output.BAM}")
"""