


ls --color=none -1 ../01_RNAseq_samples_readmapping/*.SplitNCigarReads.split.bam > bam.filelist
sed -e 's@../01_RNAseq_samples_readmapping/@@' -e 's@.SplitNCigarReads.split.bam@@' bam.filelist > bam.filelist.labels



