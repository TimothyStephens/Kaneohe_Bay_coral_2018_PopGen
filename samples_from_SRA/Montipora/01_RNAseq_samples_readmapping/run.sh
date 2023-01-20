


while read SRR; do wget https://sra-pub-run-odp.s3.amazonaws.com/sra/$SRR/$SRR -O $SRR.sra; done < Montipora_selected_Runs.txt

while read SRR; do ~/programs/sratoolkit.2.9.6-1-centos_linux64/bin/fasterq-dump $SRR.sra; done < Montipora_selected_Runs.txt

for F in *.fastq; do echo "gzip $F"; done | parallel -j 20

while read SRR; do echo -e "$SRR\t$SRR\t$PWD/prefetch_sra_files/${SRR}_R1_Trimmed.fastq.gz\t$PWD/prefetch_sra_files/${SRR}_R2_Trimmed.fastq.gz"; done < Montipora_selected_Runs.txt > Montipora_selected_Runs.samples.txt



