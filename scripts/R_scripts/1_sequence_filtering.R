## Input data derived from bat-host filtered .qza files, one per pipeline (aka. "basic" filter)
## Output of this script is a single .csv file grouping all pipeline (dada2/deblur/vsearch) methods ..
## ..and filtering methods (ex. "basic", "standard", "exact") together
## This file can then be parsed from one object into whatever number of objects are needed for plots or statistical tests

# required: library(devtools)
# runonce: install_github("jbisanz/qiime2R", force=TRUE)
library(tidyverse)
library(reshape2)
library(qiime2R)
library(phyloseq)
library(stringi)

## import data:
## download data
# notrun: download.file("https://github.com/devonorourke/tidybug/raw/master/data/qiime/vsearch.arthtable.qza", "vsearch.arthtable.qza")
# notrun: download.file("https://github.com/devonorourke/tidybug/raw/master/data/qiime/deblur.arthtable.qza", "deblur.arthtable.qza")
# notrun: download.file("https://github.com/devonorourke/tidybug/raw/master/data/qiime/dada2.arthtable.qza", "dada2.arthtable.qza")

##########   /\^._.^/\   ######################################################################
## part 1: generating the "basic" data.frame object
## this object is used as input to then further filter for the "standard" and "extra" data frames
###############################################################################################

## convert .qza to matrix, then convert wide-format matrix to long-format data.frame object
## first for dada2
featuretable <- read_qza("dada2.arthtable.qza")
mat.tmp <- featuretable$data
rm(featuretable)
df.tmp <- as.data.frame(mat.tmp)
rm(mat.tmp)
df.tmp$OTUid <- rownames(df.tmp)
rownames(df.tmp) <- NULL
dada2.tmp <- melt(df.tmp, id = "OTUid") %>% filter(value != 0) %>% mutate(Method = "dada2")
rm(df.tmp)

## repeat for deblur
featuretable <- read_qza("deblur.arthtable.qza")
mat.tmp <- featuretable$data
rm(featuretable)
df.tmp <- as.data.frame(mat.tmp)
rm(mat.tmp)
df.tmp$OTUid <- rownames(df.tmp)
rownames(df.tmp) <- NULL
deblur.tmp <- melt(df.tmp, id = "OTUid") %>% filter(value != 0) %>% mutate(Method = "deblur")
rm(df.tmp)

## repeat for vsearch (this takes a while longer - it's a much bigger matrix...)
featuretable <- read_qza("vsearch.arthtable.qza")
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
df$Filt <- "basic"

## import metadata, combine with df object:
meta <- read_delim("https://github.com/devonorourke/tidybug/raw/master/data/metadata/small_meta.txt", delim = " ")
df <- merge(df, meta, by = "SeqID")   ## note we drop a few observations from sequences present in samples known to have been contaminated
StudyNames <- c("mock", "negoro", "oro16")  ## we remove any BRI-associated data because it's not involved in this study (this was from a separate NH bat project with guano captured from individual bats; not passively collected as described in methods)
df_basic <- filter(df, StudyID %in% StudyNames)
rm(meta, StudyNames, df)

##########   /\^._.^/\   ######################################################################
## part 2: generating the "standard" data.frame object
## this object reqiures that we retain only samples that contain at least 5000 reads per sample
## and then we retain only sequence variants present in at least 2 samples (with min number of reads)
###############################################################################################

df_std_ReadCounts <- df_basic %>% group_by(Method, SeqID) %>% summarise(ReadSums=sum(Reads)) %>% filter(ReadSums >= 5000)
df_std <- df_basic %>% filter(SeqID %in% df_std_ReadCounts$SeqID) %>% mutate(Filt="standard")
df_std_HashCounts <- df_std %>% group_by(Method, HashID) %>% summarise(HashCounts=n()) %>% filter(HashCounts > 1)
df_std <- df_std %>% filter(HashID %in% df_std_HashCounts$HashID)
rm(df_std_ReadCounts, df_std_HashCounts)


##########   /\^._.^/\   ######################################################################
## part 3: generating the "extra" data.frame object
## this object enacts the same filters as "standard" but also subtracts each element of the count table..
## by a value determined by the maximum count observed for a "miss" sequence variant in the mock community sample..
## associated in that paritcular sequencing run (library)
###############################################################################################

## import the "exact", "partial", and "miss" strings, per pipeline
## for dada2
da.exact <- read_csv('https://github.com/devonorourke/tidybug/raw/master/data/qiime/dada2.exactseqs.txt', col_names = FALSE)
colnames(da.exact) <- "HashID"
da.exact$MockAlign <- "exact"
da.exact$Method <- "dada2"
da.partial <- read_csv('https://github.com/devonorourke/tidybug/raw/master/data/qiime/dada2.partialseqs.txt', col_names = FALSE)
colnames(da.partial) <- "HashID"
da.partial$MockAlign <- "partial"
da.partial$Method <- "dada2"
da.miss <- read_csv('https://github.com/devonorourke/tidybug/raw/master/data/qiime/dada2.missseqs.txt', col_names = FALSE)
colnames(da.miss) <- "HashID"
da.miss$MockAlign <- "miss"
da.miss$Method <- "dada2"
## for deblur
db.exact <- read_csv('https://github.com/devonorourke/tidybug/raw/master/data/qiime/deblur.exactseqs.txt', col_names = FALSE)
colnames(db.exact) <- "HashID"
db.exact$MockAlign <- "exact"
db.exact$Method <- "deblur"
db.partial <- read_csv('https://github.com/devonorourke/tidybug/raw/master/data/qiime/deblur.partialseqs.txt', col_names = FALSE)
colnames(db.partial) <- "HashID"
db.partial$MockAlign <- "partial"
db.partial$Method <- "deblur"
db.miss <- read_csv('https://github.com/devonorourke/tidybug/raw/master/data/qiime/deblur.missseqs.txt', col_names = FALSE)
colnames(db.miss) <- "HashID"
db.miss$MockAlign <- "miss"
db.miss$Method <- "deblur"
## for vsearch
vs.exact <- read_csv('https://github.com/devonorourke/tidybug/raw/master/data/qiime/vsearch.exactseqs.txt', col_names = FALSE)
colnames(vs.exact) <- "HashID"
vs.exact$MockAlign <- "exact"
vs.exact$Method <- "vsearch"
vs.partial <- read_csv('https://github.com/devonorourke/tidybug/raw/master/data/qiime/vsearch.partialseqs.txt', col_names = FALSE)
colnames(vs.partial) <- "HashID"
vs.partial$MockAlign <- "partial"
vs.partial$Method <- "vsearch"
vs.miss <- read_csv('https://github.com/devonorourke/tidybug/raw/master/data/qiime/vsearch.missseqs.txt', col_names = FALSE)
colnames(vs.miss) <- "HashID"
vs.miss$MockAlign <- "miss"
vs.miss$Method <- "vsearch"

## combine datasets
all.mockalign <- rbind(da.exact, da.partial, da.miss, db.exact, db.partial, db.miss, vs.exact, vs.partial, vs.miss)
rm(da.exact, da.partial, db.exact, db.partial, vs.exact, vs.partial)
## create per-pipeline datasets too

## Create a mock only dataframe to determine the maximum count value associated with a particular sequence variant..
## assigned as "miss", per library, per pipeline
da.mock_basic <- df_basic %>% filter(SampleType == "mock" & Method == "dada2")
db.mock_basic <- df_basic %>% filter(SampleType == "mock" & Method == "deblur")
vs.mock_basic <- df_basic %>% filter(SampleType == "mock" & Method == "vsearch")
all.mock_basic <- rbind(da.mock_basic, db.mock_basic, vs.mock_basic)

da.mock_Contams <- da.mock_basic %>% 
  filter(HashID %in% da.miss$HashID) %>% 
  group_by(Library) %>%
  summarise(MaxRead=max(Reads)) %>%  ## the output of this shows the per-library value to filter with
  mutate(Method="dada2")

db.mock_Contams <- db.mock_basic %>% 
  filter(HashID %in% db.miss$HashID) %>% 
  group_by(Library) %>%
  summarise(MaxRead=max(Reads)) %>%
  mutate(Method="deblur")

vs.mock_Contams <- vs.mock_basic %>% 
  filter(HashID %in% vs.miss$HashID) %>% 
  group_by(Library) %>%
  summarise(MaxRead=max(Reads)) %>%
  mutate(Method="vsearch")

all.mockContam <- rbind(da.mock_Contams, db.mock_Contams, vs.mock_Contams)

## write to disk: 
write.csv(all.mock_basic, file="~/Repos/tidybug/data/text_tables/mockReadCountsPerPipelinePerLibrary.csv", row.names = FALSE, quote=FALSE)
write.csv(all.mockContam, file="~/Repos/tidybug/data/text_tables/extraSubtractValues.csv", row.names = FALSE, quote=FALSE)

rm(da.mock_Contams, db.mock_Contams, vs.mock_Contams, all.mock_basic)
rm(da.mock_basic, db.mock_basic, vs.mock_basic)
rm(da.miss, db.miss, vs.miss)

## To create the "extra" datasets, subtract N reads from the `df_std` datasets, ..
## for the matching Pipeline and Filtering method in `all.mockContam` object
cut.da.libA <- all.mockContam %>% filter(Method=="dada2" & Library=="libA") %>% pull(MaxRead)
cut.da.libB <- all.mockContam %>% filter(Method=="dada2" & Library=="libB") %>% pull(MaxRead)
cut.da.libC <- all.mockContam %>% filter(Method=="dada2" & Library=="libC") %>% pull(MaxRead)
cut.da.libD <- all.mockContam %>% filter(Method=="dada2" & Library=="libD") %>% pull(MaxRead)
cut.db.libA <- all.mockContam %>% filter(Method=="deblur" & Library=="libA") %>% pull(MaxRead)
cut.db.libB <- all.mockContam %>% filter(Method=="deblur" & Library=="libB") %>% pull(MaxRead)
cut.db.libC <- all.mockContam %>% filter(Method=="deblur" & Library=="libC") %>% pull(MaxRead)
cut.db.libD <- all.mockContam %>% filter(Method=="deblur" & Library=="libD") %>% pull(MaxRead)
cut.vs.libA <- all.mockContam %>% filter(Method=="vsearch" & Library=="libA") %>% pull(MaxRead)
cut.vs.libB <- all.mockContam %>% filter(Method=="vsearch" & Library=="libB") %>% pull(MaxRead)
cut.vs.libC <- all.mockContam %>% filter(Method=="vsearch" & Library=="libC") %>% pull(MaxRead)
cut.vs.libD <- all.mockContam %>% filter(Method=="vsearch" & Library=="libD") %>% pull(MaxRead)

da_extra_libA <- df_std %>% filter(Library=="libA" & Method == "dada2") %>% mutate(Reads=Reads-cut.da.libA)
da_extra_libB <- df_std %>% filter(Library=="libB" & Method == "dada2") %>% mutate(Reads=Reads-cut.da.libB)
da_extra_libC <- df_std %>% filter(Library=="libC" & Method == "dada2") %>% mutate(Reads=Reads-cut.da.libC)
da_extra_libD <- df_std %>% filter(Library=="libD" & Method == "dada2") %>% mutate(Reads=Reads-cut.da.libD)
db_extra_libA <- df_std %>% filter(Library=="libA" & Method == "deblur") %>% mutate(Reads=Reads-cut.db.libA)
db_extra_libB <- df_std %>% filter(Library=="libB" & Method == "deblur") %>% mutate(Reads=Reads-cut.db.libB)
db_extra_libC <- df_std %>% filter(Library=="libC" & Method == "deblur") %>% mutate(Reads=Reads-cut.db.libC)
db_extra_libD <- df_std %>% filter(Library=="libD" & Method == "deblur") %>% mutate(Reads=Reads-cut.db.libD)
vs_extra_libA <- df_std %>% filter(Library=="libA" & Method == "vsearch") %>% mutate(Reads=Reads-cut.vs.libA)
vs_extra_libB <- df_std %>% filter(Library=="libB" & Method == "vsearch") %>% mutate(Reads=Reads-cut.vs.libB)
vs_extra_libC <- df_std %>% filter(Library=="libC" & Method == "vsearch") %>% mutate(Reads=Reads-cut.vs.libC)
vs_extra_libD <- df_std %>% filter(Library=="libD" & Method == "vsearch") %>% mutate(Reads=Reads-cut.vs.libD)

df_extra <- rbind(da_extra_libA, da_extra_libB, da_extra_libC, da_extra_libD, db_extra_libA, db_extra_libB, db_extra_libC, db_extra_libD, vs_extra_libA, vs_extra_libB, vs_extra_libC, vs_extra_libD)
rm(da_extra_libA, da_extra_libB, da_extra_libC, da_extra_libD, db_extra_libA, db_extra_libB, db_extra_libC, db_extra_libD, vs_extra_libA, vs_extra_libB, vs_extra_libC, vs_extra_libD)
rm(cut.da.libA, cut.da.libB, cut.da.libC, cut.da.libD, cut.db.libA, cut.db.libB, cut.db.libC, cut.db.libD, cut.vs.libA, cut.vs.libB, cut.vs.libC, cut.vs.libD)
df_extra$Filt <- "extra"
df_extra <- df_extra %>% filter(Reads > 0)

## combine these three $Filt objects into a single master file:
all.df <- rbind(df_basic, df_std, df_extra)

## create list of distinct HashIDs among all Filtering methods; used for plotting as text_labels (instead of long HashIDs)
set.seed(100)
uniqHashs <- all.df %>% distinct(HashID)
uniqHashs$Alias <- stri_rand_strings(nrow(uniqHashs), 6)
## add back into dataset
all.df <- merge(all.df, uniqHashs)

## write to disk
write.csv(all.df, file="~/Repos/tidybug/data/text_tables/all.filtmethods.df.csv", quote = FALSE, row.names = FALSE)

## this file is gzipped and uploaded to Github repo 'github.com/devonorourke/tidybug/data/text_tables'
