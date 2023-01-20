# Map RNA-seq reads to *P. acuta* reference genome

## Setup analysis directory

Link data files.

```bash
ln -s ../00_databases/Pocillopora_acuta_KBHIv2.assembly.dict
ln -s ../00_databases/Pocillopora_acuta_KBHIv2.assembly.fasta
ln -s ../00_databases/Pocillopora_acuta_KBHIv2.assembly.fasta.fai
ln -s ../samples_Pacuta.txt
```

Setup bash environment.

```bash
conda activate py27
```

## Align *P. acuta* RNA-seq reads from this study against the *P. acute* reference genome

Run the below script (which take the info in `samples_Pacuta.txt` and performs the alignment, aligned read processing, duplicate removal, and splitting of reads over intron-exon boundaries) to align *P. acuta* RNA-seq reads from this study against the *P. acute* reference genome.

```bash
./run_RNAseq_samples_readmapping.sh
```

### Download results from server

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0022_Coral_Genotype_Analysis/03_Analysis/2022-02-05/Pocillopora_acuta/01_RNAseq_samples_readmapping/"
rsync -zarv --delete --prune-empty-dirs --relative \
 --include="*/" \
 --include="*.tsv" --include="*.txt" \
 --include="*.sh*" --include="*.log*" \
 --exclude="*" \
 ${WD}/./* \
 . --dry-run
```





