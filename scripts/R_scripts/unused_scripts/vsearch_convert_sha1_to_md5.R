## Input data derived from concatenating each library from Vsearch-filtered .qza feature tables
## Output of this script is a .tsv feature table with sha1 hashIDs substituted for proper md5 labels
## The .tsv file is then imported back into a .qza artifact via a .biom intermediate
# required: library(devtools)
# runonce: install_github("jbisanz/qiime2R", force=TRUE)
library(tidyverse)
library(reshape2)
library(qiime2R)
library(phyloseq)
library(biomformat)

## convert .qza to matrix, then convert wide-format matrix to long-format data.frame object
featuretable <- read_qza("vsearch.all.raw.table.qza")
mat.tmp <- featuretable$data
df.tmp <- as.data.frame(mat.tmp)
df.tmp$OTUid <- rownames(df.tmp)
rownames(df.tmp) <- NULL
vsearch.tmp <- melt(df.tmp, id = "OTUid") %>% filter(value != 0) %>% mutate(Method = "vsearch")
colnames(vsearch.tmp) <- c("sha1", "SeqID", "Reads", "Method")

hashtable <- read.table(gzfile("vsearch.all.headers.txt"), header = TRUE)   ## import hashID lists and resolve sha1 vs MD5-hash'd strings in the headers
colnames(hashtable) <- c("md5", "sha1")
df <- merge(vsearch.tmp, hashtable, all.x = TRUE)       ## merge with dataframe 'df' object
df <- df[c(2,3,5)]
df$md5 <- as.character(df$md5)
df.mat <- dcast(df, md5 ~ SeqID, value.var = "Reads", fill = 0)   ## reformat to matrix
colnames(df.mat)[1] <- "#OTU ID"
write.table(df.mat, "vsearch.all.raw.md5.table.tsv", quote = FALSE, row.names = FALSE)    ## write file to disk
