

rule run_mutect2:
	params:
		threads="6",
		runtime="4:00",
		memory="16000"
	input:
		BAM= config["output_dir"]+"/BAM/{sample}.bam"
	output:
		unfiltered = config["output_dir"]+"/tmp/SNV/mutect2/{sample}/{sample}_unfiltered.vcf.gz",
		filtertagged = config["output_dir"]+"/tmp/SNV/mutect2/{sample}/{sample}_filtertagged.vcf.gz",
		filteredVCF = config["output_dir"]+"/tmp/SNV/mutect2/{sample}/{sample}_filtered.vcf.gz",
		annotatedVCF_pop = config["output_dir"]+"/tmp/SNV/mutect2/{sample}/{sample}_annotatedPop.vcf.gz",
		annotatedVCF_COSMIC = config["output_dir"]+"/tmp/SNV/mutect2/{sample}/{sample}_annotatedCOSMIC.vcf.gz"
	conda:
		"../envs/WGS.yaml"
	shell:
		"gatk --java-options '-Xmx4G' Mutect2 -R {config[reference_fasta]} "\
		"-I {input.BAM} "\
		"--germline-resource {config[germline_resource]} --af-of-alleles-not-in-resource 0.0000000025 "\
		"--panel-of-normals {config[SNV_pon]} -O {output.unfiltered} -L {config[target_genes_intervals]};"\
		"gatk FilterMutectCalls -V {output.unfiltered} -R {config[reference_fasta]} -O {output.filtertagged};"\
		"bcftools view -f .,PASS {output.filtertagged} -o {output.filteredVCF} -Oz;"\
		"bcftools index {output.filteredVCF}; "\
		"bcftools annotate -a {config[germline_full]} -c AF,AC {output.filteredVCF} -Oz -o {output.annotatedVCF_pop}; "\
		"bcftools index {output.annotatedVCF_pop}; "\
		"bcftools annotate -a {config[SNV_COSMIC]} -c ID {output.annotatedVCF_pop} -Oz -o {output.annotatedVCF_COSMIC}; "\
		"bcftools index {output.annotatedVCF_COSMIC}; "

rule postprocess_mutect2:
	params:
		threads="6",
		runtime="4:00",
		memory="16000"
	input:
		annotatedVCF_COSMIC = config["output_dir"]+"/tmp/SNV/mutect2/{sample}/{sample}_annotatedCOSMIC.vcf.gz"
	output:
		filteredTSV = config["output_dir"]+"/out/SNV/mutect2/{sample}/{sample}.tsv"
	conda:
		"../envs/WGS.yaml"

	shell:
		"python scripts/postprocess_VCF.py -i {input.annotatedVCF_COSMIC} -o {output.filteredTSV}"


rule run_mutect2_control:
	params:
		threads="6",
		runtime="20:00",
		memory="16000",
		normal_samplename = lambda wcs: config["normal_samplename"][wcs.sample]
	input:
		BAM_normal= config["output_dir"]+"/BAM/{sample}"+config["normal_suffix"]+".bam",
		BAM_tumor= config["output_dir"]+"/BAM/{sample}"+config["tumor_suffix"]+".bam"
	output:
		unfiltered = config["output_dir"]+"/out/SNV/mutect2/{sample}/{sample}_unfiltered.vcf.gz",
		filtertagged = config["output_dir"]+"/out/SNV/mutect2/{sample}/{sample}_filtertagged.vcf.gz",
		filteredVCF = config["output_dir"]+"/out/SNV/mutect2/{sample}/{sample}_filtered.vcf",
		filteredTSV = config["output_dir"]+"/out/SNV/mutect2/{sample}/{sample}.tsv"
	conda:
		"../envs/WGS.yaml"
	shell:
		"gatk --java-options '-Xmx4G' Mutect2 -R {config[reference_fasta]} "\
		"-I {input.BAM_tumor} -I {input.BAM_normal} -normal {params.normal_samplename} "\
		"--germline-resource {config[germline_resource]} "\
		"--af-of-alleles-not-in-resource 0.0000000025 -O {output.unfiltered} -L {config[target_genes_intervals]};"\
		"gatk FilterMutectCalls -V {output.unfiltered} -R /omics/groups/OE0219/internal/Etienne/data/reference/hs37d5_PhiX.fa -O {output.filtertagged};"\
		"bcftools view -f .,PASS {output.filtertagged} -o {output.filteredVCF} -O v;"\
		"python scripts/postprocess_VCF_control.py -i {output.filteredVCF} -o {output.filteredTSV}"


rule run_mutect2_control_full:
	params:
		threads="4",
		runtime="46:00",
		memory="8000",
		normal_samplename = lambda wcs: config["normal_samplename"][wcs.sample]
	input:
		BAM_normal= config["output_dir"]+"/BAM/{sample}"+config["normal_suffix"]+".bam",
		BAM_tumor= config["output_dir"]+"/BAM/{sample}"+config["tumor_suffix"]+".bam"
	output:
		unfiltered = config["output_dir"]+"/out/SNV/mutect2_full/{sample}/{sample}_{chr}_unfiltered.vcf.gz",
		filtertagged = config["output_dir"]+"/out/SNV/mutect2_full/{sample}/{sample}_{chr}_filtertagged.vcf.gz",
		filteredVCF = config["output_dir"]+"/out/SNV/mutect2_full/{sample}/{sample}_{chr}_filtered.vcf.gz"
	conda:
		"../envs/WGS.yaml"
	shell:
		"gatk --java-options '-Xmx4G' Mutect2 -R {config[reference_fasta]} "\
		"-I {input.BAM_tumor} -I {input.BAM_normal} -normal {params.normal_samplename} "\
		"--germline-resource {config[germline_resource]} "\
		"--af-of-alleles-not-in-resource 0.0000000025 -O {output.unfiltered} -L {wildcards.chr};"\
		"gatk FilterMutectCalls -V {output.unfiltered} -R {config[reference_fasta]} -O {output.filtertagged};"\
		"bcftools view -f .,PASS {output.filtertagged} -o {output.filteredVCF} -O z;"\
		"bcftools index {output.filteredVCF}"


rule concat_mutect2_control_full:
	params:
		threads="1",
		runtime="2:00",
		memory="8000"
	input:
		expand(config["output_dir"]+"/out/SNV/mutect2_full/{{sample}}/{{sample}}_{chr}_filtered.vcf.gz",chr=chromosomes)
	output:
		vcf = config["output_dir"]+"/out/SNV/mutect2_full/{sample}/{sample}.vcf.gz",
		tsv = config["output_dir"]+"/out/SNV/mutect2_full/{sample}/{sample}.tsv"
	conda:
		"../envs/WGS.yaml"
	shell:
		"bcftools concat -o {output.vcf} -O z {input};"\
		"python scripts/postprocess_VCF_control.py -i {output.vcf} -o {output.tsv}"



rule run_freebayes_SNPs_chr:
	params:
		threads="1",
		runtime="46:00",
		memory="8000"
	input:
		BAM= config["output_dir"]+"/BAM/{sample}.bam"
	output:
		config["output_dir"]+"/tmp/SNV/freebayes_chr/{sample}/{sample}_{chr}.vcf.gz"
	conda:
		"../envs/WGS.yaml"
	shell:
		"freebayes  -f {config[reference_fasta]} -r {wildcards.chr} -q 20 -@ {config[SNPs_vcf]} -l {input.BAM} | bcftools view -i 'INFO/AB>0.25 && INFO/AB<0.75 && INFO/DP>15' -o {output} -O z;"\
		"bcftools index {output}"
	# Note: freebayes -@ option only works when the -r option is provided.

rule concat_freebayes_SNPs_chr:
	params:
		threads="1",
		runtime="4:00",
		memory="8000"
	input:
		expand(config["output_dir"]+"/tmp/SNV/freebayes_chr/{{sample}}/{{sample}}_{chr}.vcf.gz",chr=chromosomes)
	output:
		config["output_dir"]+"/out/SNV/freebayes/{sample}/{sample}.vcf.gz"
	conda:
		"../envs/WGS.yaml"
	shell:
		"bcftools concat -o {output} -O z {input}; bcftools index {output}"

rule liftover:
	params:
		threads="1",
		runtime="8:00",
		memory="16000"
	input:
		"{sample}.vcf.gz"
	output:
		"{sample}_hg38.vcf.gz"
	conda:
		"../envs/WGS.yaml"
	shell:
		"bcftools annotate --rename-chrs data/chrom_map.txt -o {wildcards.sample}_hg19.vcf.gz -O z {input};"\
		"gatk LiftoverVcf I={wildcards.sample}_hg19.vcf.gz O={output} CHAIN=/omics/groups/OE0219/internal/Etienne/data/reference/hg19ToHg38.over.chain.gz REJECT={wildcards.sample}_failedLiftover.vcf.gz R=/omics/groups/OE0219/internal/Etienne/data/reference/hg38.fa;"\
		"bcftools index {output}"



# Genotype a list of SNPs in genes. This is used to then detect monoallelic expression. 
"""
rule run_freebayes_SNPs:
	params:
		threads="1",
		runtime="16:00",
		memory="16000"
	input:
		BAM= config["output_dir"]+"/BAM/{sample}.bam"
	output:
		config["output_dir"]+"/out/SNV/freebayes/{sample}/{sample}.vcf.gz"
	conda:
		"../envs/WGS.yaml"
	shell:
		"freebayes  -f {config[reference_fasta]} -q 20 -F 0.25 -C 6 -t {config[SNPs_gene_bed]} -v {output} {input.BAM}"
"""

# Look for variants in a list of genes, including germline ones
rule run_freebayes_genelist:
	params:
		threads="2",
		runtime="3:00",
		memory="6000"
	input:
		BAM= config["output_dir"]+"/BAM/{sample}.bam"
	output:
		vcf=config["output_dir"]+"/out/SNV/freebayes_genelist/{sample}/{sample}.vcf",
		annotatedVCF_pop = config["output_dir"]+"/tmp/SNV/freebayes_genelist/{sample}/{sample}_annotatedPop.vcf",
		tsv=config["output_dir"]+"/out/SNV/freebayes_genelist/{sample}/{sample}.tsv"
	conda:
		"../envs/WGS.yaml"
	shell:
		"freebayes -f {config[reference_fasta]} -q 20 -C 5 -t {config[target_genes_intervals_bed]} {input.BAM} -v {output.vcf}; "\
		"bcftools view {output.vcf} -o {output.vcf}.gz -Oz;"\
		"bcftools index {output.vcf}.gz; "\
		"bcftools annotate -a {config[germline_full]} -c AF,AC {output.vcf}.gz -Ov -o {output.annotatedVCF_pop}; "\
		"python scripts/postprocess_VCF_germline.py -i {output.annotatedVCF_pop} -o {output.tsv}; "

# Genotype only in genes. This is used to detect monoallelic expression
rule run_freebayes_genes_chr:
	params:
		threads="2",
		runtime="3:00",
		memory="6000"
	input:
		BAM= config["output_dir"]+"/BAM/{sample}.bam"
	output:
		vcf=config["output_dir"]+"/tmp/SNV/freebayes_genes_chr/{sample}/{sample}_chr{chr}.vcf",
		bed=config["output_dir"]+"/tmp/SNV/freebayes_genes_chr/{sample}/{sample}_chr{chr}.bed",
		vcf_gz = config["output_dir"]+"/tmp/SNV/freebayes_genes_chr/{sample}/{sample}_chr{chr}.vcf.gz"
	conda:
		"../envs/WGS.yaml"
	shell:
		"freebayes -f {config[reference_fasta]} -q 20 -C 5 -t {config[genes_intervals]}_chr{wildcards.chr}.bed {input.BAM} -v {output.vcf}; "\
		"python scripts/select_heterozygousSNPs.py -i {output.vcf} -o {output.bed}; "\
		"bcftools view -O z -o {output.vcf_gz} {output.vcf}; bcftools index {output.vcf_gz}"

# | bcftools view -i 'INFO/AB>0.25 && INFO/AB<0.75 && INFO/DP>12' -o {output} -O v

# For the prostate cancer dataset, we already have VCF files of germline SNPs
rule split_vcf_chr:
	params:
		threads="2",
		runtime="3:00",
		memory="6000"
	input:
		vcf_full= config["output_dir"]+"/out/SNV/germline_SNPs_all/{sample}.vcf"
	output:
		bed=config["output_dir"]+"/tmp/SNV/freebayes_genes_chr/{sample}/{sample}_chr{chr}.bed",
	conda:
		"../envs/WGS.yaml"
	shell:
		"python scripts/select_heterozygousSNPs.py -i {input.vcf_full} --chr {wildcards.chr} -o {output.bed}; "

rule run_freebayes_RNA:
	params:
		threads="3",
		runtime="23:00",
		memory="6000"
	input:
		BAM= config["output_dir"]+"/BAM_RNA/{sample}.bam",
		bed_DNA= config["output_dir"]+"/tmp/SNV/freebayes_genes_chr/{sample}/{sample}_chr{chr}.bed"
	output:
		config["output_dir"]+"/tmp/SNV/freebayes_RNA_chr/{sample}/{sample}_chr{chr}.vcf.gz"
	conda:
		"../envs/WGS.yaml"
	shell:
		"freebayes -f {config[reference_fasta]} -F 0.0 -C 0 -G 0 --min-coverage 6 --limit-coverage 10000 -q 20 -t {input.bed_DNA} {input.BAM} | bcftools norm -f {config[reference_fasta]} | bcftools view -i 'INFO/DP>6' -o {output} -O z; "\
		"bcftools index {output}"
		#TODO: might want to use a lower threshold than 6.


rule concat_freebayes_RNA_chr:
	params:
		threads="1",
		runtime="4:00",
		memory="8000"
	input:
		expand(config["output_dir"]+"/tmp/SNV/freebayes_RNA_chr/{{sample}}/{{sample}}_chr{chr}.vcf.gz",chr=chromosomes_noY)
	output:
		tmpout= config["output_dir"]+"/tmp/SNV/freebayes_RNA/{sample}/{sample}.vcf.gz",
		annotated_out = config["output_dir"]+"/out/SNV/freebayes_RNA/{sample}/{sample}.vcf.gz"
	conda:
		"../envs/WGS.yaml"
	shell:
		"bcftools concat -o {output.tmpout} -O z {input}; bcftools index {output.tmpout}; "
		"bcftools annotate -a {config[germline_full]} -c popAF:=AF,popAC:=AC {output.tmpout} -Oz -o {output.annotated_out}; "\
		"tabix -p vcf {output.annotated_out}; "


rule concat_freebayes_genes_chr:
	params:
		threads="1",
		runtime="4:00",
		memory="8000"
	input:
		expand(config["output_dir"]+"/tmp/SNV/freebayes_genes_chr/{{sample}}/{{sample}}_chr{chr}.vcf.gz",chr=chromosomes_noY)
	output:
		tmpout= config["output_dir"]+"/tmp/SNV/freebayes_genes/{sample}/{sample}.vcf.gz",
		annotated_out = config["output_dir"]+"/out/SNV/freebayes_genes/{sample}/{sample}.vcf.gz"
	conda:
		"../envs/WGS.yaml"
	shell:
		"bcftools concat -o {output.tmpout} -O z {input}; bcftools index {output.tmpout}; "
		"bcftools annotate -a {config[germline_full]} -c popAF:=AF,popAC:=AC {output.tmpout} -Oz -o {output.annotated_out}; "\
		"tabix -p vcf {output.annotated_out}; "


