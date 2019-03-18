library(tidyverse)
library(scales)
library(reshape2)
library(qiime2R)
library(viridis)

## This script is modified from the 'classifier_analyses_PalmerBayesSplit.R' script
## The goal is to further understand how an ASV assigned in amptk may differ from Vsearch
## The earlier script showed that the Bayesian methods in amptk assign many of the unassigned-Vsearch ASVs...
## ..but here we find that UTAX and SINTAX methods do not work on the same kinds of ASVs.
## In particular, UTAX appears to mostly classify sequences in cases where the Vsearch alignment is < 90%...
## ..while SINTAX assignes the data between 90-97%

# theme function for custom plot style:
theme_devon <- function () { 
  theme_bw(base_size=12, base_family="Avenir") %+replace% 
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill="transparent", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA)
    )
}


## import and wrangle classification data
baseurl <- "https://github.com/devonorourke/tidybug/raw/master/data/databases/"
vs.taxa <- read_delim(file = paste0(baseurl, "derep.vsearch.pid97.gz"), delim = "\t", col_names = TRUE)

## get HashIDs for unassigned taxa (Vsearch/Blast) or UTAX/SINTAX-assigned taxa (Palmer)
vs.taxa <- rename(vs.taxa, HashID = `Feature ID`)
vs.taxa <- vs.taxa %>%  
  select(HashID, Taxon) %>% 
  mutate(Classifier="vsearch") %>%
  mutate(Method="vsearch") %>%
  separate(data = ., col = Taxon, into = c("kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name"), sep = ";")
vs.taxa <- as.data.frame(apply(vs.taxa, 2, function(y) (gsub(".__", "", y)))) ## remove the prefixes for taxa strings
vs.taxa <- as.data.frame(apply(vs.taxa, 2, function(y) (gsub("Unassigned", NA, y)))) ## "Unassigned" listed as 'NA' now...
vs.taxa <- as.data.frame(apply(vs.taxa, 2, function(y) (gsub("^$", NA, y))))  ## convert blank records to NA
vs.taxa$kingdom_name <- NULL

## import dada2 dataset with read counts and samples ####
dada2_df <- read_qza("~/Repos/tidybug/data/qiime/dada2.arthtable.qza")
mat.tmp <- dada2_df$data
rm(dada2_df)
df.tmp <- as.data.frame(mat.tmp)
rm(mat.tmp)
df.tmp$OTUid <- rownames(df.tmp)
rownames(df.tmp) <- NULL
dada2_df <- melt(df.tmp, id = "OTUid") %>% filter(value != 0)
rm(df.tmp)
## get summary of ASV read depths and sample counts
dada2.asv.sumry <- dada2_df %>%
  group_by(OTUid) %>%
  summarise(ReadCounts=sum(value), SampleCounts=n())

## import chimera list:
baseurl <- "https://github.com/devonorourke/tidybug/raw/master/data/databases/"
chimera_df <- read_delim(file = paste0(baseurl, "chimera.ref.headers.txt"), delim = ";", col_names = FALSE)
colnames(chimera_df) <- c("OTUid", "drop", "Uchime_value")
chimera_df <- chimera_df %>% select(-drop)
chimera_df$Uchime_value[grepl("uchime_ref", chimera_df$Uchime_value)] <- "chimera_detected"

### import alignment stats ####
## collect the alignments that would typically pass a Vsearch program
align_p97c94_df <- read_delim(file = paste0(baseurl, "vsearch.p97c94out.txt"), delim = "\t", col_names = FALSE)
colnames(align_p97c94_df) <- c("query", "target", "pid", "alnlen", "mism", "opens", "qstart", "qend", "targetstart", "targetend", "drop1", "drop2")
align_p97c94_tmp <- align_p97c94_df %>% 
  select(-drop1, -drop2) %>% 
  filter(target != "*") %>% 
  mutate(AlignType="p97c94")
## collect the alignments that are below our typical 94% coverage (just 10%) - these are hits that AMPTK should pick up in the global aligner as either GS or GSL (or GDL)
## we're filtering OUT the overlapping records that are good hits in the first alignment dataframe here!
align_p97c10_df <- read_delim(file = paste0(baseurl, "vsearch.p97c10out.txt"), delim = "\t", col_names = FALSE)
colnames(align_p97c10_df) <- c("query", "target", "pid", "alnlen", "mism", "opens", "qstart", "qend", "targetstart", "targetend", "drop1", "drop2")
align_p97c10_df <- align_p97c10_df %>% 
  select(-drop1, -drop2) %>% 
  filter(target != "*") %>% 
  filter(!query %in% align_p97c94_tmp$query) %>% ## filters out duplicate records in p97c94
  mutate(AlignType="p97c10")

## combine two datasets together
align_p97_tmp <- rbind(align_p97c94_tmp, align_p97c10_df)

## collect remaining alignments that fell above 90% identity that weren't picked up by either of the 97% coverage estimates above
align_p90c10_df <- read_delim(file = paste0(baseurl, "vsearch.p90c10out.txt"), delim = "\t", col_names = FALSE)
colnames(align_p90c10_df) <- c("query", "target", "pid", "alnlen", "mism", "opens", "qstart", "qend", "targetstart", "targetend", "drop1", "drop2")
align_p90c10_df <- align_p90c10_df %>% 
  select(-drop1, -drop2) %>% 
  filter(target != "*") %>% 
  filter(!query %in% align_p97_tmp$query) %>% 
  mutate(AlignType="p90c10")

## combine three datasets together
align_true_tmp <- rbind(align_p97_tmp, align_p90c10_df)

## gather remaining records -- these should be those that don't align
nohits <- align_p97c94_df %>% 
  filter(!query %in% align_true_tmp$query) %>% 
  select(-drop1, -drop2) %>% 
  mutate(AlignType="notAligned")

## combine all datasets together
align_df <- rbind(align_true_tmp, nohits)

## cleanup:
rm(align_true_tmp, align_p97_tmp, nohits, align_p90c10_df, align_p97c94_df, align_p97c94_tmp, dada2_df, align_p97c10_df)


### full vsearch dataset with all the information for alignment, chimera, and read/sample abundances
vs_df <- merge(vs.taxa, align_df, by.x = 'HashID', by.y = 'query')
rm(vs.taxa)
vs_df <- merge(vs_df, chimera_df, all=TRUE, by.x = 'HashID', by.y = 'OTUid')
vs_df$Uchime_value[is.na(vs_df$Uchime_value)] <- "no_chimera_detected"
vs_df <- merge(vs_df, dada2.asv.sumry, by.x='HashID', by.y = 'OTUid')

## load in Jon's amptk dataset for some comparisons...
am.taxa <- read_delim(file = paste0(baseurl, "amptk.taxa.txt.gz"), delim = "\t", col_names = FALSE)
colnames(am.taxa) <- c("HashID", "splitter")
am.taxa <- separate(am.taxa, col = "splitter", into = c("splitagain", "taxa"), sep = ";")
am.taxa <- separate(am.taxa, col = "splitagain", into = c("Method", "delete"), sep = "\\|") %>% select(-delete)
am.taxa <- separate(am.taxa, col = "taxa",
                    into = c("kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name"),
                    sep=",") %>% 
  select(-kingdom_name) %>%
  mutate(Classifier="Palmer")
am.taxa <- as.data.frame(apply(am.taxa, 2, function(y) (gsub(".:", "", y)))) ## remove the prefixes for taxa strings
am.taxa <- as.data.frame(apply(am.taxa, 2, function(y) (gsub("^$", NA, y))))  ## convert blank records to NA
am_df <- merge(am.taxa, dada2.asv.sumry, all.x = TRUE, by.x='HashID', by.y = 'OTUid')
rm(am.taxa)
am_df <- merge(am_df, align_df, by.x='HashID', by.y='query')
am_df <- merge(am_df, chimera_df, by.x='HashID', by.y='OTUid', all.x = TRUE)  ## add in chimera info
am_df$Uchime_value[is.na(am_df$Uchime_value)] <- "no_chimera_detected"


## plots?
## 0) Alignments have two properties to consider: identity and coverage
## AMPTK ignores coverage, but that can be very dangeraous...
## These data represent what happens if you ignore query coverage and just pick a winner based on percent identity
#1. Plot 1 - the histogram shows that most of our "winners" have high percent identities
## plot; save as db_13a_pidHist; export at 600x500
ggplot(data = vs_df %>% filter(target != "*"), aes(pid, fill=Method)) + 
  geom_histogram(bins = 20, color="black") + scale_fill_manual(values = "gray70") +
  theme_devon() + theme(legend.position = "none") + labs(x="alignment % identity", y="number of ASVs")
#2. Plot 2 - this histogram shows that those "winners" actually have a super low coverage for many of our ASVs!
## plot; save as db_13b_alnlenHist; export at 600x500
ggplot(data = vs_df %>% filter(target != "*"), aes(alnlen, fill=Method)) + 
  geom_histogram(bins = 20, color="black") + scale_fill_manual(values = "gray70") +
  theme_devon() + theme(legend.position = "none") + labs(x="alignment coverage length", y="number of ASVs")
  ## why is this happening?
  ## The vsearch method for picking a winner requires a coverage parameter; jon's AMPTK doesn't nor does BOLD API...

#3. Plot 3 - How does alignment length, sampling depth, ASV detection frequency relate to the Palmer Method of classifier?
## plot; save as db_13c_alnlenByASVcountsAndDepth; export at 
ggplot(am_df %>% filter(target != "*") %>% filter(alnlen < 181),
       aes(y=ReadCounts, x=Method, color=alnlen, size=SampleCounts)) + 
  geom_jitter(alpha = 0.7, width = 0.3) +
  scale_y_continuous(labels = comma, trans = "log10") +
  scale_color_viridis_c(option = "viridis", end = 0.9, begin = 0.2,
                        breaks = c(0, 30, 60, 90, 120, 150, 180)) + 
  theme_devon() + theme(legend.position = "top") +
  labs(y="sequences per ASV", x="Samples ASV detected", color="Alignment\nlength (bp)", size="Samples\nper ASV") +
  guides(color = guide_colorbar(barwidth = 10)) +
  scale_size(breaks = c(1, 10, 100))

## alternatively, could just plot percent identity and alignment length directly, ignoring sampling depth/frequency of ASVs...
ggplot(am_df %>% filter(target != "*") %>% filter(alnlen < 181),
       aes(x=pid, y=Method, color=alnlen)) + 
  geom_jitter(alpha = 0.6, width = 0.3) +
  scale_color_viridis_c(option = "viridis", end = 0.9, begin = 0.2,
                        breaks = c(0, 30, 60, 90, 120, 150, 180)) + 
  theme_devon() + theme(legend.position = "top") +
  labs(x="percent identity", y="", color="Alignment\nlength (bp)") +
  guides(color = guide_colorbar(barwidth = 10)) +
  scale_size(breaks = c(1, 10, 100))


  ## major takewawy: Palmer data generates most of it's SS and US data using very short aligned sequences
  ## how many more?
am_df %>% group_by(Method, AlignType) %>% summarise(count=n()) %>% arrange(Method, -count)


## how many of the differnt Palmer Methods are missing information at each level?
am_isNAtable <- am_df %>% 
  group_by(Method) %>% 
  select(phylum_name, class_name, order_name, family_name, genus_name, species_name) %>% 
  summarise_all(funs(sum(is.na(.)))) %>%
  mutate(InfoType="missing")

tmp_isNA <- am_isNAtable %>%
select(phylum_name, class_name, order_name, family_name, genus_name, species_name) %>% 
  summarise_each(funs(sum)) %>%
  mutate(Method="Total") %>%
  mutate(InfoType="missing")

am_isNA_df <- rbind(am_isNAtable, tmp_isNA)

am_notNAtable <- am_df %>% 
  group_by(Method) %>% 
  select(phylum_name, class_name, order_name, family_name, genus_name, species_name) %>% 
  summarise_all(funs(sum(!is.na(.)))) %>%
  mutate(InfoType="present")

tmp_notNA <- am_notNAtable %>%
  select(phylum_name, class_name, order_name, family_name, genus_name, species_name) %>% 
  summarise_each(funs(sum)) %>%
  mutate(Method="Total") %>%
  mutate(InfoType="present")

am_notNA_df <- rbind(am_notNAtable, tmp_notNA)


  ## takeaway: US is missing family/genus/species in almost every sample, and SS is missing family and/or genus in most

## 1) Chimera behavior 
## For data that DID align in Vsearch, whether that be at 10, 50, or 100% identity... just  non-zero...
## how does sampling frequency and total abundances of reads related to whether 
      ## A. the ASV was marked as a chimera
      ## B. the ASV aligned with a high/low percent identity?
      ## C. the ASV aligned with a high/low query coverage?
ggplot(vs_df %>% filter(target != "*") %>% filter(alnlen < 181),
       aes(y=ReadCounts, x=SampleCounts, color=alnlen)) + 
  geom_point() +
  scale_y_continuous(labels = comma, trans = "log10") +
  scale_color_viridis_c(option = "plasma", end = 0.85,
                        breaks = c(0, 30, 60, 90, 120, 150, 180)) + 
  facet_grid(~ AlignType, scales="free") +
  theme_devon() + theme(legend.position = "top") +
  labs(y="sequences per ASV", x="nSamples ASV detected", color="Alignment\nlength (bp)") +
  guides(color = guide_colorbar(barwidth = 10))

## takeaways:
  #0. this plot illustrates why you HAVE TO FILTER for query alignment coverage... all the purples ..
  #   .. might get called in AMPTK, but won't be called in Vsearch because of the 'qcov' parameter filtering those purple dots out
  #1. most ASVs aren't chimeras, and those that are have low sample counts .. but...
  #2. the number of chimeras that align with high percent identity might exist in many samples, but ...
  #   .. don't necessarily have high alignment values

##for those that DO NOT align to any reference in VSEARCH (even with just 10% coverage here!!)..
  ## most still aren't considered chimeric... and most that are chimerica exist in very low abundances 
ggplot(vs_df %>% filter(target == "*"),
       aes(x=ReadCounts, y=SampleCounts, shape=Uchime_value)) + 
  geom_point() +
  facet_grid(AlignType ~ Uchime_value) +
  scale_x_continuous(labels = comma, trans="log10") +
  theme_devon() + theme(legend.position = "none") +
  labs(x="sequences per ASV", y="samples ASV detected")
## the 'chimera_detected' outlier (HashID == 'cf49b9c5d2a1f9cb458a9f3edcc495ff') apparently comes from a moth (Phyllophaga longispina)...
  ## 160/180 bp are aligned... might be marked as chimeric because multiple records of this species exist?
  ## the 'no_chimera_detected' outliers are things like mites, but their alignment scores are low in BLAST ...

## Does Jon's data indicate that most of these NOHIT chimeras are just chordates that we're missing in our ...
  ## .. vsearch approach that contains only arthropod records?
ggplot(am_df %>% filter(target == "*"), 
       aes(x=ReadCounts, y=SampleCounts, color=Uchime_value, shape=Uchime_value)) + 
  geom_point(alpha=0.75) +
  facet_grid(phylum_name ~ Method) +
  scale_x_continuous(labels = comma, trans = "log10") +
  scale_color_manual(values=c("red", "black")) +
  theme_devon() + theme(legend.position = "top") +
  labs(x="sequences per ASV", y="samples ASV detected", color="", shape="")
  ## not even a bit. the only chordates here are "no_chimera_detected", and there's almost never any reads there...
  ## but this does show something important: all the "NOHIT"'s getting picked up by Jon's method ...
  ## are being picked up by UTAX or SINTAX - almost every by the global alignment approach


## 2) For Jon's different classification methods, let's compare how each method relates to the ..
## .. kinds of data that is aligned > 97%, less than 97%, or not at all (by VSearch standards)
am_df$AlignType <- factor(am_df$AlignType, levels = c("p97c94", "p97c10", "p90c10", "notAligned"))
ggplot(data=am_df %>% filter(Method != "GDL"), 
       aes(x=AlignType, y=ReadCounts)) +
  geom_jitter(width=0.4, alpha=0.3) +
  scale_y_continuous(labels=comma, trans = "log10") +
  labs(x="", y="sequences per ASV") +
  theme_devon() +
  theme(legend.position = "none") +
  facet_wrap(~ Method, nrow=2)

## plot to illustrate that msot of the data that is missed from Vsearch with 97% identity but > 90%
## ..ultimately get's picked up by Jon's "Sintax" algorithm
## ..while several of the data that isn't even at 90% may continue to be classified by UTAX
## plot; save as db_9f_PalmerVsearchMethodcomparison_withAbundances; exporrt at 750x650
## those aligned with > 97% tend to have greates number of seqs

## can generate the values for these too in a table:
dt <- am_df %>% group_by(Method, AlignType) %>% summarise(counts=n()) %>% spread(Method, counts)
dt$AlignType <- gsub("\n", "", dt$AlignType)
dt
write.csv(dt, file="~/Repos/tidybug/data/databases/bayes_comparisons_palmer.vsearch_counts.csv", row.names = FALSE, quote = FALSE)

## alt plot, using colors from previous figures
## plot; save as db_9falt_PalmerVsearchMethodcomparison_withAbundances; exporrt at 750x650
am_df$AlignType <- factor(am_df$AlignType, levels=c("p97c94", "p97c10", "p90c10", "notAligned"))

ggplot(data=am_df %>% filter(Method != "GDL"), aes(x=AlignType, y=ReadCounts, color=Method)) + geom_jitter(width=0.4, alpha=0.3) +
  scale_y_continuous(labels=comma, trans = "log10") + scale_color_manual(values = c("black", "gray40", "purple", "gold")) +
  labs(x="", y="sequences per ASV") + theme_devon() + theme(legend.position = "none") + facet_wrap(~ Method, nrow=2)

## alt plot, using colors to flag for Chimeras:
## plot; save as db_9falt2_PalmerVsearchMethodcomparison_withChimeras; exporrt at 1000x600
ggplot(data=am_df %>% filter(Method != "GDL"), aes(x=AlignType, y=ReadCounts, color=Uchime_value)) +
  geom_jitter(width=0.3, alpha=0.4) + scale_y_continuous(labels=comma, trans = "log10") + scale_color_manual(values = c("red", "black")) +
  labs(x="", y="sequences per ASV") + theme_devon() + theme(legend.position = "none") + facet_grid(Uchime_value ~ Method)

dt_chimera <- am_df %>% group_by(Method, AlignType, Uchime_value) %>% summarise(counts=n()) %>% spread(Method, counts)
dt_chimera$AlignType <- gsub("\n", "", dt_chimera$AlignType)
dt_chimera
write.csv(dt_chimera, file="~/Repos/tidybug/data/databases/bayes_comparisons_palmer.vsearch_CHIMERAcounts.csv", row.names = FALSE, quote = FALSE)


## Any trends among these ASVs by taxonomic Class? What about the Chordates that are assigned in amptk?
## looks like most of the data getting classified that isn't an insect is either:
#1 - aligned with usearch and > 97%
#2 - < 90% aligned and classified with UTAX
ggplot(data=am_df %>% filter(Method != "GDL") %>% filter(phylum_name != "Arthropoda"), 
#ggplot(data=am_df %>% filter(Method != "GDL") %>% filter(class_name != "Insecta"), 
  aes(x=AlignType, y=ReadCounts, color=class_name)) +
  geom_jitter(width=0.2) +
  scale_y_continuous(labels=comma, trans = "log10") +
  labs(x="", y="sequences per ASV") +
  theme_devon() +
  #theme(legend.position = "none") +
  facet_wrap(~ Method, nrow=3)

## data table for summary like before
## how many of these are non insects?
dt_noInsects <- am_df %>% filter(class_name != "Insecta") %>% 
  group_by(Method, AlignType) %>% summarise(counts=n()) %>% spread(Method, counts)
write.csv(dt_noInsects, file="~/Repos/tidybug/data/databases/bayes_comparisons_palmer.vsearch_counts-noInsects.csv", row.names = FALSE, quote = FALSE)
## how many of these are non arthropods?
dt_noArthropods <- am_df %>% filter(phylum_name != "Arthropoda") %>% 
  group_by(Method, AlignType) %>% summarise(counts=n()) %>% spread(Method, counts)
write.csv(dt_noArthropods, file="~/Repos/tidybug/data/databases/bayes_comparisons_palmer.vsearch_counts-noArthropods.csv", row.names = FALSE, quote = FALSE)

## how many of the Unclassified data from Vsearch (standard Qiime2 parameters) have data in the Palmer method?
## which Taxa Classes are more represented?
## How many of these are assigned via global alignment vs. Bayesian classifiers?
## How many of these were flagged as chimeras?
am_nohits <- am_df %>% filter(HashID %in% nohits$HashID)
am_nohits <- merge(am_nohits, chimera_df, by.x="HashID", by.y = "OTUid", all.x = TRUE)
am_nohits$Uchime_value[is.na(am_nohits$Uchime_value)] <- "no_chimera_detected"

## how many of these ASVs are there, by Class_name, by Uchime_value?
am_nohit_grouptable <- am_nohits %>% group_by(Uchime_value, phylum_name, class_name) %>% summarise(nASVs=n())
ggplot(am_nohits %>% filter(phylum_name!="Arthropoda"), aes(x=Method, y=ReadCounts, color=class_name, size=SampleCounts)) +
  geom_jitter(width=0.2) +
  theme_devon() +
  facet_wrap(~ Uchime_value, ncol=2)


## how many of the ASVs identified by Palmer's methods have a low percent coverage?
ggplot(data = am_df %>% filter(Method != "GDL") %>% filter(target != "*"), 
       aes(y=pid, x=alnlen)) + 
  geom_point() + 
  facet_wrap(~ Method, ncol=2)


##### the final plot:
## save as db_13d_alnlen_pid_byTaxaLevel_andMethod; export at 1000x800
tester <- am_df %>%
  select(HashID, Method, phylum_name, class_name, order_name, family_name, genus_name, species_name) %>%
  gather(key, val, -HashID, -Method)

hash_phylum <- am_df %>% select(HashID, phylum_name) %>% filter(complete.cases(.))
hash_class <- am_df %>% select(HashID, phylum_name, class_name) %>% filter(complete.cases(.))
hash_order <- am_df %>% select(HashID, phylum_name, class_name, order_name) %>% filter(complete.cases(.))
hash_family <- am_df %>% select(HashID, phylum_name, class_name, order_name, family_name) %>% filter(complete.cases(.))
hash_genus <- am_df %>% select(HashID, phylum_name, class_name, order_name, family_name, genus_name) %>% filter(complete.cases(.))
hash_species <- am_df %>% select(HashID, phylum_name, class_name, order_name, family_name, genus_name, species_name) %>% filter(complete.cases(.))

am_df_species <- am_df %>% filter(HashID %in% hash_species$HashID) %>%
  mutate(Group="Species")
am_df_genus <- am_df %>% filter(HashID %in% hash_genus$HashID) %>% 
  filter(!HashID %in% am_df_species$HashID) %>%
  mutate(Group="Genus")
am_df_family <- am_df %>% filter(HashID %in% hash_family$HashID) %>% 
  filter(!HashID %in% am_df_species$HashID) %>% filter(!HashID %in% am_df_genus$HashID) %>%
  mutate(Group="Family")
am_df_order <- am_df %>% filter(HashID %in% hash_order$HashID) %>% 
  filter(!HashID %in% am_df_species$HashID) %>% filter(!HashID %in% am_df_genus$HashID) %>% filter(!HashID %in% am_df_family$HashID) %>%
  mutate(Group="Order")
am_df_class <- am_df %>% filter(HashID %in% hash_class$HashID) %>% 
  filter(!HashID %in% am_df_species$HashID) %>% filter(!HashID %in% am_df_genus$HashID) %>% filter(!HashID %in% am_df_family$HashID) %>%
  filter(!HashID %in% am_df_order$HashID) %>%
  mutate(Group="Class")
am_df_phylum <- am_df %>% filter(HashID %in% hash_phylum$HashID) %>% 
  filter(!HashID %in% am_df_species$HashID) %>% filter(!HashID %in% am_df_genus$HashID) %>% filter(!HashID %in% am_df_family$HashID) %>%
  filter(!HashID %in% am_df_order$HashID) %>% filter(!HashID %in% am_df_class$HashID) %>%
  mutate(Group="Phylum")


tester <- rbind(am_df_species, am_df_genus, am_df_family, am_df_order, am_df_class, am_df_phylum)

tester$Group <- factor(tester$Group, levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"))


ggplot(tester %>% filter(Method != "GDL") %>% filter(Group != "Phylum") %>% filter(pid > 1),
       aes(x=Group, y=pid, color=alnlen)) + 
  geom_jitter(width = 0.3) +
  facet_grid(Method ~ .) +
  scale_color_viridis_c(option = "viridis", end = 0.9, begin = 0.2,
                        breaks = c(0, 30, 60, 90, 120, 150, 180)) + 
  scale_y_continuous(breaks = c(90, 95, 100)) +
  theme_devon() + theme(legend.position = "top") +
  facet_wrap(~ Method, ncol = 2) +
  labs(y="percent identity", x="", color="Alignment\nlength (bp)") +
  guides(color = guide_colorbar(barwidth = 10))


## can also include the AlignType info to partition table further...
## plot; save as 'db_13d-alt_alnlen_pid_byTaxaLevel_andMethod'; export at
ggplot(tester %>% filter(Method != "GDL") %>% filter(Group != "Phylum") %>% filter(pid > 1),
       aes(x=Group, y=pid, color=alnlen)) + 
  geom_jitter(width = 0.3) +
  facet_grid(Method ~ .) +
  scale_color_viridis_c(option = "viridis", end = 0.9, begin = 0.2,
                        breaks = c(0, 30, 60, 90, 120, 150, 180)) + 
  scale_y_continuous(breaks = c(90, 95, 100)) +
  theme_devon() + theme(legend.position = "top") +
  facet_grid(Method ~ AlignType) +
  labs(y="percent identity", x="", color="Alignment\nlength (bp)") +
  guides(color = guide_colorbar(barwidth = 10))

