module load samtools/1.15.1                                                                         
module load bcftools/1.10.2                                                                         
module load htslib/1.12                                                                             
module load sambamba/0.7.1
module load gatk/4.2.0.0

bsub -o lsf.o%J -W 40:00 -R "rusage[mem=16G]" "snakemake --configfile $1 --use-conda --cluster 'bsub -o lsf.o%J -W {params.runtime} -n {params.threads} -R \"rusage[mem={params.memory}]\"' --jobs 100 --latency-wait 60"