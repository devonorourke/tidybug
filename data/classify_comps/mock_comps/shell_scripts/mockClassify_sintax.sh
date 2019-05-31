#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/Data/tidy_paper/classifer_comps/mock_comps/classified_output/sintax
#SBATCH --job-name="mock_st"
#SBATCH --output="mock_st.log"
#SBATCH --cpus-per-task=12

module purge
module load linuxbrew/colsa

VSEARCH=/mnt/lustre/macmaneslab/devon/pkgs/vsearch-2.13.4/bin/vsearch
SINFASTA=/mnt/lustre/macmaneslab/devon/guano/BOLDdb/allBOLD/tmp_foramptk.fasta
REPSEQS=/mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4.fa

$VSEARCH \
--sintax $REPSEQS \
--db $SINFASTA \
--tabbedout vsintax_out_sedited.txt \
--sintax_cutoff 0.9 \
--threads 12
