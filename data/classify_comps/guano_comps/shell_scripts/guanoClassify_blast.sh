#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/Data/tidy_paper/classifer_comps/guano_comps/blast
#SBATCH --job-name="guano_bl"
#SBATCH --output="guano_bl.log"

module purge
module load anaconda/colsa
source activate qiime2-2019.1

REFPATH=/mnt/lustre/macmaneslab/devon/guano/BOLDdb/allBOLD
FASTAPATH=/mnt/lustre/macmaneslab/devon/guano/Data/tidy_paper/qiime_tables/dada2

qiime feature-classifier classify-consensus-blast \
  --i-query "$FASTAPATH"/dada2.arthseqs.qza \
  --i-reference-reads "$REFPATH"/boldCOI.derep.seqs.qza \
  --i-reference-taxonomy "$REFPATH"/boldCOI.derep.tax.qza \
  --p-maxaccepts 1000 \
  --p-perc-identity 0.97 \
  --p-query-cov 0.89 \
  --p-strand both \
  --o-classification guano.blast_out.qza
