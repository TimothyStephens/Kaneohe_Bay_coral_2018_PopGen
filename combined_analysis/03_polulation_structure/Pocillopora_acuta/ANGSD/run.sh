


ls -1 ../../../Pocillopora_acuta/01_RNAseq_samples_readmapping/*.bam ../../../samples_from_SRA/Pocillopora/01_RNAseq_samples_readmapping/*.bam > bam.filelist
sed -e 's@../../../Pocillopora_acuta/01_RNAseq_samples_readmapping/@@' -e 's@../../../samples_from_SRA/Pocillopora/01_RNAseq_samples_readmapping/@@' -e 's@.SplitNCigarReads.split.bam@@' bam.filelist > bam.filelist.labels



