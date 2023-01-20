

~/programs/vcftools/bin/vcftools \
--remove-indv SRR7059699 \
--remove-indv SRR7062295 \
--remove-indv SRR7062350 \
--remove-indv SRR7039808 \
--remove-indv SRR7040514 \
--remove-indv SRR7041301 \
--remove-indv SRR7042978 \
--remove-indv SRR7043704 \
--remove-indv SRR7055829 \
--remove-indv SRR7058378 \
--remove-indv SRR7058566 \
--remove-indv SRR7058606 \
--remove-indv SRR7058616 \
--remove-indv SRR6951423 \
--remove-indv SRR6951744 \
--remove-indv SRR6952431 \
--remove-indv SRR6963586 \
--remove-indv SRR6963878 \
--remove-indv SRR6963891 \
--remove-indv SRR6964364 \
--remove-indv SRR6986864 \
--remove-indv SRR6987146 \
--remove-indv SRR6914609 \
--remove-indv SRR6914151 \
--remove-indv SRR6914908 \
--remove-indv SRR6934388 \
--remove-indv SRR6934542 \
--remove-indv SRR6935629 \
--remove-indv SRR6942678 \
--remove-indv SRR6942729 \
--remove-indv SRR7043013 \
--remove-indv SRR7046161 \
--remove-indels \
--maf 0.01 --min-meanDP 10 --max-missing 1 \
--gzvcf GVCFall.vcf.gz --out vcftools \
--recode-bcf --recode-INFO-all


