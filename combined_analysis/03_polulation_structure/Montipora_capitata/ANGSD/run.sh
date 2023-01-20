

ls -1 ../../../Montipora_capitata/01_RNAseq_samples_readmapping/*.bam ../../../samples_from_SRA/Montipora/01_RNAseq_samples_readmapping/*.bam > bam.filelist
sed -e 's@../../../Montipora_capitata/01_RNAseq_samples_readmapping/@@' -e 's@../../../samples_from_SRA/Montipora/01_RNAseq_samples_readmapping/@@' -e 's@.SplitNCigarReads.split.bam@@' bam.filelist > bam.filelist.labels


