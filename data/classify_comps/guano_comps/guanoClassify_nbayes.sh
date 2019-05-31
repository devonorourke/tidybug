#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/Data/tidy_paper/classifer_comps/guano_comps/nbayes
#SBATCH --job-name="guano_nb"
#SBATCH --output="guano_nb.log"

module purge
module load anaconda/colsa
source activate qiime2-2019.1

READS=/mnt/lustre/macmaneslab/devon/guano/Data/tidy_paper/qiime_tables/dada2/dada2.arthseqs.qza
CLSSFYR=/mnt/lustre/macmaneslab/devon/guano/BOLDdb/classifier/nbClassifer_allboldDB.qza

qiime feature-classifier classify-sklearn \
  --i-reads "$READS" --i-classifier "$CLSSFYR" \
  --p-n-jobs 1 --p-reads-per-batch 2000 \
  --o-classification guano.nbayes_out.qza \
  --verbose
