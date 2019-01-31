## example shell script used to create the Cutadapt-trimmed QIIME paired-end artifact by importing raw .fastq paired reads
## note this script was modified for each library

#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p21/fastq
#SBATCH --partition=macmanes,shared
#SBATCH --job-name="p21.qimport"
#SBATCH --output=p21.qimport.log
#SBATCH --cpus-per-task=24

module purge
module load anaconda/colsa
source activate qiime2-2018.11

LIB=$(pwd | cut -d '/' -f 9)

## Importing data into QIIME - from raw data!
## We don't have a typical input file type, so we'll create a manifest file and import as paired-end data.
## See their import tutorial for details: https://docs.qiime2.org/2018.11/tutorials/importing/#sequence-data-with-sequence-quality-information-i-e-fastq

## Create the manifest file:
## find and remove any file pairs with no sequence information
## note that other directory `qiime` previously created to store output of this script
find . -type f -name "*.gz" -size 1k -delete
pwd > ../qiime/pwd.tmptxt
ls -1 | cut -f 1 -d '_' | uniq -d > ../qiime/keeplist.tmptxt
find . -type f | grep -f ../qiime/keeplist.tmptxt | sort -u | cut -d '/' -f 2 > ../qiime/filenames.tmptxt
cd ../qiime
## create sample names to paste (col1)
cut -f 1 -d "_" filenames.tmptxt > col1.tmptxt
wc -l col1.tmptxt | cut -f 1 -d ' ' > lines.tmptxt
## create directory col to paste (col2)
seq $(echo $(cat lines.tmptxt)) | xargs -Iz echo $(cat pwd.tmptxt) > dirpath.tmptxt
paste dirpath.tmptxt filenames.tmptxt -d "/" > col2.tmptxt
## create read direction col to paste (col3)
for i in $(cat filenames.tmptxt); do
if [[ $i == *_R1_* ]];
then
  echo "forward"
else
  echo "reverse"
fi;
done > col3.tmptxt
paste col1.tmptxt col2.tmptxt col3.tmptxt -d "," > manifest.tmptxt
head -1 col2.tmptxt | cut -f 9 -d '/' > libname.tmptxt
echo 'sample-id,absolute-filepath,direction' | cat - manifest.tmptxt > $(cat libname.tmptxt).manifest.file
rm *.tmptxt
```

## Now we can import the data:
mkdir reads
cd /.reads

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$LIB".manifest.file \
  --output-path "$LIB".demux.qza \
  --input-format PairedEndFastqManifestPhred33


## Primer and barcode trimming with Cutadapt
## Using Cutadapt plugin within QIIME to remove primer sequences at both 5' and 3' end.
## See their documentation: https://docs.qiime2.org/2018.11/plugins/available/cutadapt/?highlight=cutadapt
## Note we're using the trim-paired option because our data is already demultiplexed.

  qiime cutadapt trim-paired \
  --i-demultiplexed-sequences "$LIB".demux.qza \
  --p-cores 24 \
  --p-adapter-f GGTCAACAAATCATAAAGATATTGG \
  --p-adapter-r GGATTTGGAAATTGATTAGTWCC \
  --o-trimmed-sequences "$LIB".trimd.qza

  qiime demux summarize \
    --i-data "$LIB".trimd.qza \
    --p-n 500000 \
    --o-visualization "$LIB".trimd.qzv
