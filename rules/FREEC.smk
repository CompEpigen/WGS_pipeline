rule run_FREEC:
	params:
		threads="8",
		runtime="4:00",
		memory="16000",
		sex = lambda wcs: sample_sex[wcs.sample]
	input:
		BAM= config["output_dir"]+"/BAM/{sample}.bam"
	output:
		ratios=config["output_dir"]+"/out/CNA/FREEC/{sample}/{sample}.bam_ratio.txt",
		CNAs=config["output_dir"]+"/out/CNA/FREEC/{sample}/{sample}.bam_CNVs"
	conda:
		"../envs/WGS.yaml"
	shell:
		"cp {config[template_FREEC]} {config[output_dir]}/out/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{chrLenFile}}#{config[chrLenFile]}#g' {config[output_dir]}/out/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{outdir}}#{config[output_dir]}/out/CNA/FREEC/{wildcards.sample}#g' {config[output_dir]}/out/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{BAMFILE}}#{input.BAM}#g' {config[output_dir]}/out/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{SEX}}#{params.sex}#g' {config[output_dir]}/out/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"{config[freec_bin]} -config {config[output_dir]}/out/CNA/FREEC/{wildcards.sample}/config.txt"


rule run_FREEC_BAF:
	params:
		threads="8",
		runtime="4:00",
		memory="16000",
		sex = lambda wcs: sample_sex[wcs.sample]
	input:
		BAM= config["output_dir"]+"/BAM/{sample}.bam"
	output:
		ratios=config["output_dir"]+"/out/CNA/FREEC_BAF/{sample}/{sample}.bam_ratio.txt",
		CNAs=config["output_dir"]+"/out/CNA/FREEC_BAF/{sample}/{sample}.bam_CNVs"
	conda:
		"../envs/WGS.yaml"
	shell:
		"cp {config[template_FREEC_BAF]} {config[output_dir]}/out/CNA/FREEC_BAF/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{chrLenFile}}#{config[chrLenFile]}#g' {config[output_dir]}/out/CNA/FREEC_BAF/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{outdir}}#{config[output_dir]}/out/CNA/FREEC_BAF/{wildcards.sample}#g' {config[output_dir]}/out/CNA/FREEC_BAF/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{BAMFILE}}#{input.BAM}#g' {config[output_dir]}/out/CNA/FREEC_BAF/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{SEX}}#{params.sex}#g' {config[output_dir]}/out/CNA/FREEC_BAF/{wildcards.sample}/config.txt;"\
		"{config[freec_bin]} -config {config[output_dir]}/out/CNA/FREEC_BAF/{wildcards.sample}/config.txt"

rule run_FREEC_control:
	params:
		threads="8",
		runtime="4:00",
		memory="16000",
		sex = lambda wcs: sample_sex[wcs.sample]
	input:
		BAM_normal= config["output_dir"]+"/BAM/{sample}"+config["normal_suffix"]+".bam",
		BAM_tumor= config["output_dir"]+"/BAM/{sample}"+config["tumor_suffix"]+".bam"
	output:
		ratios=config["output_dir"]+"/out/CNA/FREEC/{sample}/{sample}.bam_ratio.txt",
		CNAs=config["output_dir"]+"/out/CNA/FREEC/{sample}/{sample}.bam_CNVs"
	conda:
		"../envs/WGS.yaml"
	shell:
		"cp {config[template_FREEC]} {config[output_dir]}/out/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{chrLenFile}}#{config[chrLenFile]}#g' {config[output_dir]}/out/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{outdir}}#{config[output_dir]}/out/CNA/FREEC/{wildcards.sample}#g' {config[output_dir]}/out/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{BAMFILE_NORMAL}}#{input.BAM_normal}#g' {config[output_dir]}/out/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{BAMFILE_TUMOR}}#{input.BAM_tumor}#g' {config[output_dir]}/out/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{SEX}}#{params.sex}#g' {config[output_dir]}/out/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"{config[freec_bin]} -config {config[output_dir]}/out/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"cp {config[output_dir]}/out/CNA/FREEC/{wildcards.sample}/{wildcards.sample}{config[tumor_suffix]}.bam_ratio.txt {output.ratios}; "\
		"cp {config[output_dir]}/out/CNA/FREEC/{wildcards.sample}/{wildcards.sample}{config[tumor_suffix]}.bam_CNVs {output.CNAs}; "\

rule create_vizualization_datafiles_FREEC:
	params:
		threads="1",
		runtime="40",
		memory="8000"
	input:
		ratios = config["output_dir"]+"/out/CNA/FREEC/{sample}/{sample}.bam_ratio.txt",
		CNAs = config["output_dir"]+"/out/CNA/FREEC/{sample}/{sample}.bam_CNVs",
	output:
		circos = config["output_dir"]+"/out/CNA/FREEC/{sample}/{sample}_CNA_circos.txt",
		chrplot = config["output_dir"]+"/out/CNA/FREEC/{sample}/{sample}_CNA_chrplots.tsv"
	conda:
		"../envs/WGS.yaml"
	shell:
		"python scripts/freec2circos.py -i {input.ratios} -o {output.circos} --cnv {input.CNAs};"\
		"python scripts/freec2chrplots.py -i {input.ratios} -o {output.chrplot} --cnv {input.CNAs};"

rule FREEC_to_shatterseek:
	params:
		threads="1",
		runtime="20",
		memory="8000"
	input:
		config["output_dir"]+"/out/CNA/FREEC/{sample}/{sample}.bam_ratio.txt"
	output:
		config["output_dir"]+"/out/CNA/FREEC/{sample}/{sample}_shatterseekCNA.tsv"

	conda:
		"../envs/WGS.yaml"
	shell:
		"python scripts/freec2shatterseek.py -i {input} -o {output}"