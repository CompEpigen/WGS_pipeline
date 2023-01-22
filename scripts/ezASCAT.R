library("ezASCAT")
library(ASCAT)                                                                                       
                                                                                                     
args = commandArgs(trailingOnly=TRUE)
bam_file = args[1]
sample_name = args[2]
output_dir = args[3]

dir.create(output_dir)
setwd(output_dir)

counts = ezASCAT::get_counts(t_bam = bam_file ,build = "hg19",op=sample_name,nthreads=8)
ascat.bc = ezASCAT::prep_ascat_t(t_counts = paste(sample_name,"_nucleotide_counts.tsv",sep=""), sample_name = sample_name)

#ascat.plotRawData(ascat.bc)
ascat.gg = ascat.predictGermlineGenotypes(ascat.bc, "AffySNP6")

ascat.bc = ascat.aspcf(ascat.bc,ascat.gg=ascat.gg,penalty=30)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc,gamma=1.0)
write.table(ascat.output$segments,file= paste(sample_name,"_segments.tsv",sep=""),quote=FALSE,sep="\t",row.names=FALSE)