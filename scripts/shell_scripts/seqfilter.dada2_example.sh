## example shell script used to create the Deblur-filtered output used in the analyses of quality filtering pipelines
## note this script was modified for each library

#!/bin/bash

#SBATCH --cpus-per-task=24
#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/Data/p41/qiime
#SBATCH --partition=macmanes,shared
#SBATCH --job-name="qp41"

module purge
module load anaconda/colsa
source activate qiime2-2018.11
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$LIB".trimd.qza \
  --p-trunc-len-f 175 \
  --p-trunc-len-r 175 \
  --p-n-threads 24 \
  --o-table "$LIB".table.qza \
  --o-representative-sequences "$LIB".repSeqs.qza \
  --o-denoising-stats "$LIB".denoisingStats.qza
