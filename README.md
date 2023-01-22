# WGS_pipeline

This is a pipeline developed by Etienne Sollier for processing WGS data, either with or without matched controls. It starts from aligned bam files and generates CNA, SV and SNV calls, as well as plots (circos and chromosome plots). For now, the pipeline can only be used from the DKFZ cluster. 

## Quick start
You need a conda environment with snakemake and pysam installed. You then have to adapt the config.yaml and samplesheet.tsv files (see input), and you can finally run `./run_pipeline.sh config.yaml`

## Input
- samplesheet: csv file containing the sample names (column "sample") and the path to the aligned bam files ("path_bam"). Optionally, one can also provide a control bam ("path_bam_control") and the sex of each sample ("sex"), which determines the germline copy number of sex chromosomes. See samplesheet.tsv and samplesheet_AML.tsv for examples.
- config.yaml: yaml files where options are specified. Most can be left to their default values, but one always has to specify a working directory and the samplesheet.

## Output
The directory structure will be:

```
working_dir/
 |-BAM/
    |-<sample1>.bam/
    |-<sample1>.bam.bai/
 |-out_control/
        |-CNA
          |-FREEC
            |-<sample1>
              |-...
        |-SV
          |-manta
            |-<sample1>
              |-...
        |-SNV
          |-mutect2_genomewide
            |-<sample>
        |-plots
          |-chrplots
          |-circosplots
 |-out_nocontrol/
      |-...
```

## Tools used in the pipeline

### CNA calling with [Control-FREEC](http://boevalab.inf.ethz.ch/FREEC/)
The pipeline calls CNA using Control-FREEC. To avoid false somatic CNA calls in case a matched normal samples is not available, some regions are masked from the CNA inference. The problematic regions were identified as regions where samples from the 1000Genomes project had CNAs. The file data/GC_profile_FREEC_PoN-1000G.cnp is used for this filtering.

### SV calling with [manta](https://github.com/Illumina/manta)
SVs are called using manta. SVs are filtered using a minimum number of split reads and split pairs, defined in the config file. When no matched control is available, we use a Panel of Normal generated using samples from the 1000Genomes project (PoN_SV_1000G.bed) to filter out putative germline SVs.

### SNV calling using [mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)
Somatic SNVs are called using mutect2. When no matched control is available, we only look for somatic mutations in a list of genes (target_genes_intervals in the config file). We also use the gnomAD and COSMIC databases, and select non-synonymous variants, in order to filter out germline variants.


