# Map RNA-seq reads to *P. acuta* reference genome

## Setup analysis directory

Link data files.

```bash
ln -s ../00_databases/Pocillopora_acuta_KBHIv2.assembly.dict
ln -s ../00_databases/Pocillopora_acuta_KBHIv2.assembly.fasta
ln -s ../00_databases/Pocillopora_acuta_KBHIv2.assembly.fasta.fai
ln -s ../../Pocillopora_selected_Runs.txt
```

Setup bash environment.

```bash
conda activate py27
```

## Align *P. acuta* RNA-seq reads downloaded from SRA against the *P. acute* reference genome

**Download reads**

We have to use `prefetch` to download the SRA files (into `/scratch/timothy/tmp/ncbi/sra`) before running the genotyping process. 

```bash
PREFETCH="/home/timothy/programs/sratoolkit.2.9.6-1-centos_linux64/bin/prefetch"
OUT="prefetch_sra_files"
while read SRR;
do
  if [ ! -f "$OUT/$SRR.sra" ];
  then
    echo "$PREFETCH --output-directory $OUT --max-size 1000GB $SRR 1>$OUT/$SRR.prefetch.log 2>&1"
  fi
done < Pocillopora_selected_Runs.txt | parallel -j 1

# Check that all log files indicate the the download was a success. 
grep -L 'was downloaded successfully' $OUT/*.log
```

The above command didn't appear to work for a lot of the samples for some unknown reason. So the below workaround was used.

```bash
# Download directly using wget
while read SRR;
do
  wget https://sra-pub-run-odp.s3.amazonaws.com/sra/$SRR/$SRR -O $SRR.sra
done < Pocillopora_selected_Runs.txt

# Extract reads in fastq format using fasterq-dump
while read SRR; 
do
  ~/programs/sratoolkit.2.9.6-1-centos_linux64/bin/fasterq-dump $SRR.sra
done < Pocillopora_selected_Runs.txt

# Compress fastq files so they dont take up so much space
for F in *.fastq; do echo "gzip $F"; done | parallel -j 20
```

Output `fasta` files will be in `prefetch_sra_files/`

**Trim reads**

Quality trim the RNA-seq reads using similar thresholds as we did for the samples generated in this study.

```bash
./run_Cutadapt.sh
```

Generate a file listing the sample ID and the location (absolute paths) to the trimmed read files.

```bash
while read SRR;
do 
  echo -e "$SRR\t$SRR\t$PWD/prefetch_sra_files/${SRR}_R1_Trimmed.fastq.gz\t$PWD/prefetch_sra_files/${SRR}_R2_Trimmed.fastq.gz"; 
done < Pocillopora_selected_Runs.txt > Pocillopora_selected_Runs.samples.txt
```

**Gather read stats**

Run `statswrapper.sh` on each file so that we get the number of reads and total bases in each of the raw and trimmed `fastq` files. 

```bash
for F in *.fastq.gz;
do
  ~/programs/bbmap/statswrapper.sh $F 1>$F.stats.txt
done
```

Collect read stats from the separate `*.stats.txt` files into a single table.

```bash
while read ID;
do
  RAW_R1_N=$(awk 'NR==2{print $1}' ${ID}.sra_1.fastq.gz.stats.txt)
  RAW_R1_B=$(awk 'NR==2{print $3}' ${ID}.sra_1.fastq.gz.stats.txt)
  RAW_R2_N=$(awk 'NR==2{print $1}' ${ID}.sra_2.fastq.gz.stats.txt)
  RAW_R2_B=$(awk 'NR==2{print $3}' ${ID}.sra_2.fastq.gz.stats.txt)
  TRIM_R1_N=$(awk 'NR==2{print $1}' ${ID}_R1_Trimmed.fastq.gz.stats.txt)
  TRIM_R1_B=$(awk 'NR==2{print $3}' ${ID}_R1_Trimmed.fastq.gz.stats.txt)
  TRIM_R2_N=$(awk 'NR==2{print $1}' ${ID}_R2_Trimmed.fastq.gz.stats.txt)
  TRIM_R2_B=$(awk 'NR==2{print $3}' ${ID}_R2_Trimmed.fastq.gz.stats.txt)
  echo -e "${ID}\t${RAW_R1_N}\t${RAW_R1_B}\t${RAW_R2_N}\t${RAW_R2_B}\t${TRIM_R1_N}\t${TRIM_R1_B}\t${TRIM_R2_N}\t${TRIM_R2_B}"
done < ../Pocillopora_selected_Runs.txt > read_stats.tsv
```

**Align trimmed reads**

Once all of the SRA files have been trimmed we can run the genotyping.

```bash
./run_RNAseq_samples_readmapping.sh
```

**Call haplotypes**

Call haplotypes in the aligned RNA-seq data using `gatk HaplotypeCaller`.

```bash
./run_HaplotypeCaller-Diploid.sh
```

### Download results from server

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0022_Coral_Genotype_Analysis/03_Analysis/2022-02-05/samples_from_SRA/Pocillopora/01_RNAseq_samples_readmapping/"
rsync -zarv --delete --prune-empty-dirs --relative \
 --include="*/" \
 --include="*.tsv" --include="*.txt" \
 --include="*.sh*" --include="*.log*" \
 --exclude="*" \
 ${WD}/./* \
 . --dry-run
```







