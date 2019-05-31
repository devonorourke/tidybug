#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/Data/tidy_paper/classifer_comps/mock_comps/classified_output/vsearch
#SBATCH --job-name="mock_vs"
#SBATCH --cpus-per-task=12
#SBATCH --output="mock_vs.log"

module purge
module load anaconda/colsa
source activate qiime2-2019.1

qiime feature-classifier classify-consensus-vsearch \
  --i-query /mnt/lustre/macmaneslab/devon/guano/Data/tidy_paper/classifer_comps/mock_comps/mockIM4_seqs.qza \
  --i-reference-reads /mnt/lustre/macmaneslab/devon/guano/BOLDdb/allBOLD/boldCOI.derep.seqs.qza \
  --i-reference-taxonomy /mnt/lustre/macmaneslab/devon/guano/BOLDdb/allBOLD/boldCOI.derep.tax.qza \
  --p-maxaccepts 1000 \
  --p-perc-identity 0.97 \
  --p-query-cov 0.89 \
  --p-strand both \
  --p-threads 12 \
  --o-classification mockIM4.vsearch_out.qza
