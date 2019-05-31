#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/Data/tidy_paper/classifer_comps/mock_comps/classified_output/blast
#SBATCH --job-name="mock_bl"
#SBATCH --output="mock_bl.log"

module purge
module load anaconda/colsa
source activate qiime2-2019.1

  qiime feature-classifier classify-consensus-blast \
  --i-query /mnt/lustre/macmaneslab/devon/guano/Data/tidy_paper/classifer_comps/mock_comps/mockIM4_seqs.qza \
  --i-reference-reads /mnt/lustre/macmaneslab/devon/guano/BOLDdb/allBOLD/boldCOI.derep.seqs.qza \
  --i-reference-taxonomy /mnt/lustre/macmaneslab/devon/guano/BOLDdb/allBOLD/boldCOI.derep.tax.qza \
  --p-maxaccepts 1000 \
  --p-perc-identity 0.97 \
  --p-query-cov 0.89 \
  --p-strand both \
  --o-classification mockIM4.blast_out.qza
