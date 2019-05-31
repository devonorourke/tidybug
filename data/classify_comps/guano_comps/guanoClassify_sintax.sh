#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/Data/tidy_paper/classifer_comps/guano_comps/sintax
#SBATCH --job-name="guano_st"
#SBATCH --output="guano_st.log"
#SBATCH --cpus-per-task=12

module purge
module load linuxbrew/colsa

VSEARCH=/mnt/lustre/macmaneslab/devon/pkgs/vsearch-2.13.4/bin/vsearch
SINFASTA=/mnt/lustre/macmaneslab/devon/guano/BOLDdb/allBOLD/tmp_foramptk.fasta
REPSEQS=/mnt/lustre/macmaneslab/devon/guano/Data/tidy_paper/qiime_tables/dada2/dada2.arthASVs.fasta/dada2.arthASVs.fasta

$VSEARCH \
--sintax $REPSEQS \
--db $SINFASTA \
--tabbedout vsintax_guano_out.txt \
--sintax_cutoff 0.9 \
--threads 12
