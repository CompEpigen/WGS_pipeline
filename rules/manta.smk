rule run_manta:
	params:
		threads="8",
		runtime="10:00",
		memory="16000"
	input:
		BAM= config["output_dir"]+"/BAM/{sample}.bam"
	output:
		config["output_dir"]+"/out/SV/manta/{sample}/{sample}_SV_unfiltered.vcf"
	conda:
		"../envs/manta.yaml"
	shell:
		"python {config[manta_bin]}/configManta.py --tumorBam {input.BAM} "\
		"--referenceFasta {config[reference_fasta]} "\
		"--runDir {config[output_dir]}/tmp/SV/manta/{wildcards.sample} "\
		"--callRegions data/hs37d5_chr.bed.gz --exome;"\
		"python {config[output_dir]}/tmp/SV/manta/{wildcards.sample}/runWorkflow.py -j {params.threads} -g 16; "\
		"cp {config[output_dir]}/tmp/SV/manta/{wildcards.sample}/results/variants/tumorSV.vcf.gz {config[output_dir]}/out/SV/manta/{wildcards.sample}/{wildcards.sample}_SV_unfiltered.vcf.gz; "\
		"gunzip -d {config[output_dir]}/out/SV/manta/{wildcards.sample}/{wildcards.sample}_SV_unfiltered.vcf.gz; "\
		"rm -r {config[output_dir]}/tmp/SV/manta/{wildcards.sample}"

rule run_manta_control:
	params:
		threads="8",
		runtime="20:00",
		memory="16000"
	input:
		BAM_normal= config["output_dir"]+"/BAM/{sample}"+config["normal_suffix"]+".bam",
		BAM_tumor= config["output_dir"]+"/BAM/{sample}"+config["tumor_suffix"]+".bam",
	output:
		config["output_dir"]+"/out/SV/manta/{sample}/{sample}_SV_somatic.vcf"
	conda:
		"../envs/manta.yaml"
	shell:
		"python {config[manta_bin]}/configManta.py --normalBam {input.BAM_normal} --tumorBam {input.BAM_tumor} "\
		"--referenceFasta {config[reference_fasta]} "\
		"--runDir ${{SCRATCHDIR}}/${{LSB_JOBID}}/manta/{wildcards.sample} "\
		"--callRegions data/hs37d5_chr.bed.gz --exome;"\
		"python ${{SCRATCHDIR}}/${{LSB_JOBID}}/manta/{wildcards.sample}/runWorkflow.py -j {params.threads} -g 16; "\
		"cp ${{SCRATCHDIR}}/${{LSB_JOBID}}/manta/{wildcards.sample}/results/variants/somaticSV.vcf.gz {output}.gz; "\
		"gunzip -d {output}.gz; "\


rule postprocess_manta:
	params:
		threads="1",
		runtime="5:45",
		memory="16000"
	input:
		VCF = config["output_dir"]+"/out/SV/manta/{sample}/{sample}_SV_unfiltered.vcf",
		BAM = config["output_dir"]+"/BAM/{sample}.bam"
	output:
		config["output_dir"]+"/out/SV/manta/{sample}/{sample}_SV_filtered.vcf"
	conda:
		"../envs/WGS.yaml"
	shell:
		"python scripts/postprocess_mantaVCF.py -i {input.VCF} -o {output} "\
		"--minPR {config[manta_minPR]} --minSR {config[manta_minSR]} --minLen {config[manta_minLen]} "\
		"--pon {config[manta_pon]} --bam {input.BAM} --mappability {config[manta_mappability]} "

rule postprocess_manta_control:
	params:
		threads="1",
		runtime="5:45",
		memory="16000"
	input:
		VCF = config["output_dir"]+"/out/SV/manta/{sample}/{sample}_SV_somatic.vcf",
		BAM = config["output_dir"]+"/BAM/{sample}"+config["tumor_suffix"]+".bam"
	output:
		config["output_dir"]+"/out/SV/manta/{sample}/{sample}_SV_filtered.vcf"
	conda:
		"../envs/WGS.yaml"
	shell:
		"python scripts/postprocess_mantaVCF.py -i {input.VCF} -o {output} "\
		"--minPR {config[manta_minPR]} --minSR {config[manta_minSR]} --minLen {config[manta_minLen]} "\
		"--bam {input.BAM} --tumorindex 1"

rule create_vizualization_datafiles_manta:
	params:
		threads="1",
		runtime="20",
		memory="8000"
	input:
		config["output_dir"]+"/out/SV/manta/{sample}/{sample}_SV_filtered.vcf"
	output:
		circos = config["output_dir"]+"/out/SV/manta/{sample}/{sample}_SV_circos.txt",
		chrplot = config["output_dir"]+"/out/SV/manta/{sample}/{sample}_SV_chrplots.tsv",
	conda:
		"../envs/WGS.yaml"
	shell:
		"python scripts/manta2circos.py -i {input} -o {output.circos};"\
		"python scripts/manta2chrplots.py -i {input} -o {output.chrplot};"
