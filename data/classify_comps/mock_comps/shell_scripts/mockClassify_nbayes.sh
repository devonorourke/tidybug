#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/Data/tidy_paper/classifer_comps/mock_comps/nbayes
#SBATCH --job-name="mock_sk"
#SBATCH --output="mock_sk.log"

module purge
module load anaconda/colsa
source activate qiime2-2019.1

READS=/mnt/lustre/macmaneslab/devon/guano/Data/tidy_paper/classifer_comps/mock_comps/mockIM4_seqs.qza
CLSSFYR=/mnt/lustre/macmaneslab/devon/guano/BOLDdb/classifier/nbClassifer_allboldDB.qza

qiime feature-classifier classify-sklearn \
  --i-reads "$READS" --i-classifier "$CLSSFYR" \
  --p-n-jobs 1 --p-reads-per-batch 2000 \
  --o-classification mockIM4.sklearn_out.qza \
  --verbose
