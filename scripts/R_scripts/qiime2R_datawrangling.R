## Data created from `seqfilter.combineNfilter.sh` output
## Modifications described here produce the 'df' data.frame object used as input for subsequent analyses/plots

library(tidyverse)
library(reshape2)
library(stringi)

## import data:
## set working directory to your path of choice...
#notrun: setwd("~/Documents/dissertation/methods_paper/datasets/")
## download data
download.file("https://github.com/devonorourke/pzero/raw/master/data/vsearch.arthtable.qza", "vsearch.table.qza")
download.file("https://github.com/devonorourke/pzero/raw/master/data/deblur.arthtable.qza", "deblur.table.qza")
download.file("https://github.com/devonorourke/pzero/raw/master/data/dada2.arthtable.qza", "dada2.table.qza")

## convert .qza to matrix, then convert wide-format matrix to long-format data.frame object
## first for dada2
featuretable <- read_qza("dada2.table.qza")
mat.tmp <- featuretable$data
rm(featuretable)
df.tmp <- as.data.frame(mat.tmp)
rm(mat.tmp)
df.tmp$OTUid <- rownames(df.tmp)
rownames(df.tmp) <- NULL
dada2.tmp <- melt(df.tmp, id = "OTUid") %>% filter(value != 0) %>% mutate(Method = "dada2")
rm(df.tmp)

## repeat for deblur
featuretable <- read_qza("deblur.table.qza")
mat.tmp <- featuretable$data
rm(featuretable)
df.tmp <- as.data.frame(mat.tmp)
rm(mat.tmp)
df.tmp$OTUid <- rownames(df.tmp)
rownames(df.tmp) <- NULL
deblur.tmp <- melt(df.tmp, id = "OTUid") %>% filter(value != 0) %>% mutate(Method = "deblur")
rm(df.tmp)

## repeat for vsearch (this takes a while longer - it's a much bigger matrix...)
featuretable <- read_qza("vsearch.table.qza")
mat.tmp <- featuretable$data
rm(featuretable)
df.tmp <- as.data.frame(mat.tmp)
rm(mat.tmp)
df.tmp$OTUid <- rownames(df.tmp)
rownames(df.tmp) <- NULL
vsearch.tmp <- melt(df.tmp, id = "OTUid") %>% filter(value != 0) %>% mutate(Method = "vsearch")
rm(df.tmp)

## combine data
df <- rbind(dada2.tmp, deblur.tmp, vsearch.tmp)
rm(dada2.tmp, deblur.tmp, vsearch.tmp)
colnames(df) <- c("HashID", "SeqID", "Reads", "Method")

## import metadata, combine with df object:
download.file("https://github.com/devonorourke/pzero/raw/master/data/metadata.tsv", "metadata.tsv")
meta <- read.table("metadata.tsv", sep="", header=TRUE)
df <- merge(df, meta, by = "SeqID")   ## note we drop a few observations from sequences present in samples known to have been contaminated
rm(meta)

## it was noticed mid-workflow that the vsearch hashing algorithm differed from the denoise/dada2 processes..
## ..so we used the following trick to generate both the sha1 and md5 hash IDs in the vsearch files:
## 1) create the new fasta with both labels: vsearch -derep_fulllength dna-sequences.fasta -relabel_md5 -relabel_keep -xsize -output p*.vsearch.md5.seqs.fasta -fasta_width 0
## 2) create the list of header names: grep "^>" p*.vsearch.md5.seqs.fasta | sed 's/>//' > p*.headers.txt
## 3) concatenate those headers into a single list: p*.headers.txt >> allheaders.txt
## 4) those header names were uplaoded to github repo and are imported here:

## import hashID lists and resolve sha1 vs MD5-hash'd strings in the headers
download.file("https://github.com/devonorourke/pzero/raw/master/data/allheaders.txt.gz", "allheaders.txt.gz")
hashtable <- read.table(gzfile("allheaders.txt.gz"), header = FALSE)
colnames(hashtable) <- c("md5", "HashID")
## merge with dataframe 'df' object
df <- merge(df, hashtable, all.x = TRUE)
## substitute any 'N/A' values with appropriate term to get just the md5 list right
df$HashID <- as.character(df$HashID)
df$md5 <- as.character(df$md5)
df$md5[is.na(df$md5)] <- df$HashID[is.na(df$md5)]
df$HashID <- NULL
colnames(df)[7] <- "HashID"

## transform the df$Reads to log2 scale:
df$log2reads <- round(log2(df$Reads)) 

## create list of distinct HashIDs among all Filtering methods
set.seed(100)
uniqHashs <- df %>% distinct(HashID)
uniqHashs$Alias <- stri_rand_strings(nrow(uniqHashs), 6)
## add back into dataset
df <- merge(df, uniqHashs)  

## write file to disk
write.csv(df, "all.filtmethods.df.csv", quote = FALSE, row.names = FALSE)
## this file is gzipped and uploaded to Github repo 'github.com/devonorourke/pzero'