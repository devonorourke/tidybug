library(tidyverse)
library(reshape2)
library(qiime2R)
library(scales)

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

## import classification data
baseurl <- "https://github.com/devonorourke/tidybug/raw/master/data/databases/"
vs.taxa <- read_delim(file = paste0(baseurl, "derep.vsearch.pid97.gz"), delim = "\t", col_names = TRUE)
bl.taxa <- read_delim(file = paste0(baseurl, "derep.blast_pid97.tsv.gz"), delim = "\t", col_names = TRUE)
am.taxa <- read_delim(file = paste0(baseurl, "amptk.taxa.txt.gz"), delim = "\t", col_names = FALSE)

## wrangle classification data to get HashIDs for unassigned taxa (Vsearch/Blast) or UTAX/SINTAX-assigned taxa (Palmer)
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
vs.hashID <- vs.taxa %>% filter(is.na(class_name)) %>% select(HashID)

bl.taxa <- rename(bl.taxa, HashID = `Feature ID`)
bl.taxa <- bl.taxa %>%  
  select(HashID, Taxon) %>% 
  mutate(Classifier="blast") %>%
  mutate(Method="blast") %>%
  separate(data = ., col = Taxon, into = c("kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name"), sep = ";")
bl.taxa <- as.data.frame(apply(bl.taxa, 2, function(y) (gsub(".__", "", y)))) ## remove the prefixes for taxa strings
bl.taxa <- as.data.frame(apply(bl.taxa, 2, function(y) (gsub("Unassigned", NA, y)))) ## "Unassigned" listed as 'NA' now...
bl.taxa <- as.data.frame(apply(bl.taxa, 2, function(y) (gsub("^$", NA, y))))  ## convert blank records to NA
bl.taxa$kingdom_name <- NULL
bl.hashID <- bl.taxa %>% filter(is.na(class_name)) %>% select(HashID)

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
am.methodlist <- c("US", "SS")
am.hashID <- am.taxa %>% filter(Method %in% am.methodlist)

## import dada2 dataset with read counts and samples
dada2_df <- read_qza("~/Repos/tidybug/data/qiime/dada2.arthtable.qza")
mat.tmp <- dada2_df$data
rm(dada2_df)
df.tmp <- as.data.frame(mat.tmp)
rm(mat.tmp)
df.tmp$OTUid <- rownames(df.tmp)
rownames(df.tmp) <- NULL
dada2_df <- melt(df.tmp, id = "OTUid") %>% filter(value != 0)
rm(df.tmp)
## remove mock samples:
removenames <- c("mockIM4p4L1", "mockIM4p4L2", "mockIM4p7L1", "mockIM4p7L2")
dada2_df <- dada2_df %>% filter(!variable %in% removenames)

## calculate the read sums and frequency of detection for each HashID listed for each Method 
vs.unassigned.sumry <- dada2_df %>% filter(OTUid %in% vs.hashID$HashID) %>%
  group_by(OTUid) %>%
  summarise(ReadCounts=sum(value), SampleCounts=n()) %>%
  mutate(Labeler="Vsearch\nUnassigned")
bl.unassigned.sumry <- dada2_df %>% filter(OTUid %in% bl.hashID$HashID) %>%
  group_by(OTUid) %>%
  summarise(ReadCounts=sum(value), SampleCounts=n()) %>%
  mutate(Labeler="Blast\nUnassigned")

## calculate for those that ARE assigned
vs.assigned.sumry <- dada2_df %>% filter(!OTUid %in% vs.hashID$HashID) %>%
  group_by(OTUid) %>%
  summarise(ReadCounts=sum(value), SampleCounts=n()) %>%
  mutate(Labeler="Vsearch\nAssigned")
bl.assigned.sumry <- dada2_df %>% filter(!OTUid %in% bl.hashID$HashID) %>%
  group_by(OTUid) %>%
  summarise(ReadCounts=sum(value), SampleCounts=n()) %>%
  mutate(Labeler="Blast\nAssigned")

vs_bl.sumry <- rbind(vs.unassigned.sumry, vs.assigned.sumry, bl.unassigned.sumry, bl.assigned.sumry)
## plot distributions of sequence depths amd counts of observed ASVS for unclassified ASVS for Vsearch and Blast
## plot sample Counts; save as db_10a_unassignedSamples; export at 500x350
ggplot(vs_bl.sumry, aes(x=Labeler, y=SampleCounts, color=Labeler)) +
  geom_jitter(alpha=0.5) +
  geom_boxplot(outlier.shape = NA, color="black") +
  #scale_y_continuous(limits = c(0,50)) +
  scale_color_manual(values = c("mediumpurple4", "mediumpurple1", "olivedrab4", "olivedrab3")) +
  labs(x="", y="number of samples") +
  theme_devon() +
  theme(legend.position = "none")

## plot sequence Counts; save as db_10b_unassignedSequences; export at 500x500
ggplot(vs_bl.sumry, aes(x=Labeler, y=ReadCounts, color=Labeler)) +
  geom_jitter(alpha=0.5) +
  geom_boxplot(outlier.shape = NA, color="black", alpha=0.75) +
  scale_y_continuous(trans = "log2", labels = comma) +
  scale_color_manual(values = c("mediumpurple4", "mediumpurple1", "olivedrab4", "olivedrab3")) +
  labs(x="", y="number of sequences") +
  theme_devon() +
  theme(legend.position = "none")

## are any of these distributions different?
library(FSA)
Summarize(ReadCounts ~ Labeler, data = vs_bl.sumry)
  ## means are always higher in Assigned, but variance is huge...
Summarize(SampleCounts ~ Labeler, data = vs_bl.sumry)
  ## similar trend here.

## Kruskal-Wallis test
vs_bl.sumry$Labeler <- as.factor(vs_bl.sumry$Labeler)
kruskal.test(ReadCounts ~ Labeler, data=vs_bl.sumry)
  # Kruskal-Wallis chi-squared = 415.69, df = 3, p-value < 2.2e-16
kruskal.test(SampleCounts ~ Labeler, data=vs_bl.sumry)
  # Kruskal-Wallis chi-squared = 47.022, df = 3, p-value = 3.439e-10

## How many of these Assigned or Unassigned ASVs are considered chimeras?
## import chimera list:
baseurl <- "https://github.com/devonorourke/tidybug/raw/master/data/databases/"
chimera_df <- read_delim(file = paste0(baseurl, "chimera.ref.headers.txt"), delim = ";", col_names = FALSE)
colnames(chimera_df) <- c("OTUid", "drop", "Uchime_value")
chimera_df <- chimera_df %>% select(-drop)
chimera_df$Uchime_value[grepl("uchime_ref", chimera_df$Uchime_value)] <- "chimera_detected"
chimera_df <- merge(vs_bl.sumry, chimera_df, all=TRUE)
chimera_df <- chimera_df %>% drop_na(ReadCounts)
chimera_df <- chimera_df %>% replace(is.na(.), "no_chimera")
chimera_df <- chimera_df %>% separate(., into = c("Method", "Information"), col = Labeler, sep = "\n")

## summary table:
# vsearchonly %>% group_by(Information) %>% summarise(ASVs=n(), fracSum=sum(ReadCounts)) %>% mutate(fraction=fracSum/allsums)
#    Information  ASVs  fracSum fraction
#  1 Assigned     8173 47103813    0.835
#  2 Unassigned   5234  9274609    0.165

## plot; save as db_11_Chimeras; export at 600x600
ggplot(chimera_df,
       aes(x=ReadCounts, y=SampleCounts, color=Uchime_value)) + 
  geom_point(data=chimera_df %>% filter(Uchime_value=="no_chimera"), alpha=0.2) +
  geom_point(data=chimera_df %>% filter(Uchime_value=="chimera_detected"), alpha=0.85) +
  scale_x_continuous(labels = comma, limits = c(0, 210000)) +
  scale_y_continuous(limits = c(0, 65)) +
  scale_color_manual(values = c("red", "black"), labels = c("yes", "no")) +
  facet_grid(Method ~ Information) +
  labs(x="sequences per ASV", y="samples with ASV detected", color="Chimera detected") +
  theme_devon() +
  theme(legend.position = "top")

