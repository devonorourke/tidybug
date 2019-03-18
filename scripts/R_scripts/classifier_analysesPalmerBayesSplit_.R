library(tidyverse)
library(scales)

## This script is modified from the 'classifier_analyses.R' script
## The goal is to identify how many shared/distinct taxa exist among Vsearch, Blast, and Palmer classifiers
## Unlike the earlier script, thsi one splits the Palmer records by UTAX, SINTAX, and USEARCH hits, ...
## ... rather than keeping them as single records. We find that more of Jon's distinct records come from ...
## ... UTAX or SINTAX; so the global aligners do similar jobs; the difference is his Bayesian rescues

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

## import data
baseurl <- "https://github.com/devonorourke/tidybug/raw/master/data/databases/"
vs.taxa <- read_delim(file = paste0(baseurl, "derep.vsearch.pid97.gz"), delim = "\t", col_names = TRUE)
bl.taxa <- read_delim(file = paste0(baseurl, "derep.blast_pid97.tsv.gz"), delim = "\t", col_names = TRUE)
am.taxa <- read_delim(file = paste0(baseurl, "amptk.taxa.txt.gz"), delim = "\t", col_names = FALSE)
bo.taxa <- read_csv(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/boldAPI_lca97.csv")

## wrangle data
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
vs.taxa <- as.data.frame(apply(vs.taxa, 2, function(y) (gsub("Ambiguous_taxa", NA, y))))  ## For this analysis, we're modifying the dataset so that "Ambiguous taxa" is assigned NA

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
bl.taxa <- as.data.frame(apply(bl.taxa, 2, function(y) (gsub("Ambiguous_taxa", NA, y))))  ## For this analysis, we're modifying the dataset so that "Ambiguous taxa" is assigned NA


colnames(am.taxa) <- c("HashID", "splitter")
am.taxa <- separate(am.taxa, col = "splitter", into = c("splitagain", "taxa"), sep = ";")
am.taxa <- separate(am.taxa, col = "splitagain", into = c("Method", "delete"), sep = "\\|") %>% select(-delete)
am.taxa <- separate(am.taxa, col = "taxa",
                    into = c("kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name"),
                    sep=",") %>% 
  select(-kingdom_name) %>%
  mutate(Classifier="AMPTK")
am.taxa <- as.data.frame(apply(am.taxa, 2, function(y) (gsub(".:", "", y)))) ## remove the prefixes for taxa strings
am.taxa <- as.data.frame(apply(am.taxa, 2, function(y) (gsub("^$", NA, y))))  ## convert blank records to NA

colnames(bo.taxa)[2] <- "HashID"
bo.taxa <- bo.taxa %>% 
  select(HashID, location_output) %>% 
  separate(., col = location_output, sep = ";",
           into = c("phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")) %>%
  mutate(Classifier="boldAPI") %>%
  mutate(Method="boldAPI")
miss.bo <- vs.taxa %>% filter(!HashID %in% bo.taxa$HashID)   ## pulling out query names that bold API didn't match
miss.bo[,2:9] <- NA
miss.bo$Method <- "boldAPI"
miss.bo$Classifier <- "boldAPI"
bo.taxa <- rbind(miss.bo, bo.taxa)
rm(miss.bo)

## 3) For ASVs with distinct information per Level, what fraction of Palmer Classifiers are Usearch vs. Baysean?
## comparisons for Class, Order, Family, Genus, and Species levels

## Class level comparisons
class.bayescompfunction <- function(data1, data2, name1, name2, level){
  pair1_df <- data1 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>%
    mutate(name1.taxstring=paste(phylum_name, class_name, sep=";")) %>% select(HashID, name1.taxstring, Method)
  colnames(pair1_df)[3] <- paste0("Method", name1)
  pair2_df <- data2 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>%
    mutate(name2.taxstring=paste(phylum_name, class_name, sep=";")) %>% select(HashID, name2.taxstring, Method)
  colnames(pair2_df)[3] <- paste0("Method", name2)
  tmp_match <- merge(pair1_df, pair2_df, all=TRUE)
  grouper1 <- colnames(pair1_df)[3]
  grouper2 <- colnames(pair2_df)[3]
  pair <- tmp_match %>% 
    mutate(matchtest = name1.taxstring==name2.taxstring) %>%
    group_by(matchtest, MethodAMPTK) %>%
    summarise(nMatches=n()) %>%
    mutate(Pair=paste0(name1,":",name2)) %>%
    mutate(Level=level) %>%
    ungroup(matchtest) %>% 
    mutate(., matchtest = ifelse(is.na(matchtest), "AMPTK_only", matchtest)) %>%
    mutate(MethodAMPTK=as.character(MethodAMPTK)) %>%
    mutate(., MethodAMPTK = ifelse(is.na(MethodAMPTK), "no_data", MethodAMPTK)) %>%
    mutate(., matchtest = ifelse(grepl("no_data", MethodAMPTK), gsub("AMPTK_only", "no_data", matchtest), matchtest)) %>%
    mutate(., matchtest = gsub(FALSE, "Distinct", matchtest)) %>%
    mutate(., matchtest = gsub(TRUE, "Shared", matchtest))
}

p.v.class <- class.bayescompfunction(am.taxa, vs.taxa, "AMPTK", "Vsearch", "Class")
p.b.class <- class.bayescompfunction(am.taxa, bl.taxa, "AMPTK", "Blast", "Class")
p.a.class <- class.bayescompfunction(am.taxa, bo.taxa, "AMPTK", "boldAPI", "Class")
all.class <- rbind(p.v.class, p.b.class, p.a.class)
rm(p.v.class, p.b.class, p.a.class)

## Order level comparisons
order.bayescompfunction <- function(data1, data2, name1, name2, level){
  pair1_df <- data1 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>%
    mutate(name1.taxstring=paste(phylum_name, class_name, order_name, sep=";")) %>% select(HashID, name1.taxstring, Method)
  colnames(pair1_df)[3] <- paste0("Method", name1)
  pair2_df <- data2 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>%
    mutate(name2.taxstring=paste(phylum_name, class_name, order_name, sep=";")) %>% select(HashID, name2.taxstring, Method)
  colnames(pair2_df)[3] <- paste0("Method", name2)
  tmp_match <- merge(pair1_df, pair2_df, all=TRUE)
  grouper1 <- colnames(pair1_df)[3]
  grouper2 <- colnames(pair2_df)[3]
  pair <- tmp_match %>% 
    mutate(matchtest = name1.taxstring==name2.taxstring) %>%
    group_by(matchtest, MethodAMPTK) %>%
    summarise(nMatches=n()) %>%
    mutate(Pair=paste0(name1,":",name2)) %>%
    mutate(Level=level) %>%
    ungroup(matchtest) %>% 
    mutate(., matchtest = ifelse(is.na(matchtest), "AMPTK_only", matchtest)) %>%
    mutate(MethodAMPTK=as.character(MethodAMPTK)) %>%
    mutate(., MethodAMPTK = ifelse(is.na(MethodAMPTK), "no_data", MethodAMPTK)) %>%
    mutate(., matchtest = ifelse(grepl("no_data", MethodAMPTK), gsub("AMPTK_only", "no_data", matchtest), matchtest)) %>%
    mutate(., matchtest = gsub(FALSE, "Distinct", matchtest)) %>%
    mutate(., matchtest = gsub(TRUE, "Shared", matchtest))
}

p.v.order <- order.bayescompfunction(am.taxa, vs.taxa, "AMPTK", "Vsearch", "Order")
p.b.order <- order.bayescompfunction(am.taxa, bl.taxa, "AMPTK", "Blast", "Order")
p.a.order <- order.bayescompfunction(am.taxa, bo.taxa, "AMPTK", "boldAPI", "Order")
all.order <- rbind(p.v.order, p.b.order, p.a.order)
rm(p.v.order, p.b.order, p.a.order)


## Family level comparisons
family.bayescompfunction <- function(data1, data2, name1, name2, level){
  pair1_df <- data1 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(family_name)) %>%
    mutate(name1.taxstring=paste(phylum_name, class_name, order_name, family_name, sep=";")) %>% select(HashID, name1.taxstring, Method)
  colnames(pair1_df)[3] <- paste0("Method", name1)
  pair2_df <- data2 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(family_name)) %>%
    mutate(name2.taxstring=paste(phylum_name, class_name, order_name, family_name, sep=";")) %>% select(HashID, name2.taxstring, Method)
  colnames(pair2_df)[3] <- paste0("Method", name2)
  tmp_match <- merge(pair1_df, pair2_df, all=TRUE)
  grouper1 <- colnames(pair1_df)[3]
  grouper2 <- colnames(pair2_df)[3]
  pair <- tmp_match %>% 
    mutate(matchtest = name1.taxstring==name2.taxstring) %>%
    group_by(matchtest, MethodAMPTK) %>%
    summarise(nMatches=n()) %>%
    mutate(Pair=paste0(name1,":",name2)) %>%
    mutate(Level=level) %>%
    ungroup(matchtest) %>% 
    mutate(., matchtest = ifelse(is.na(matchtest), "AMPTK_only", matchtest)) %>%
    mutate(MethodAMPTK=as.character(MethodAMPTK)) %>%
    mutate(., MethodAMPTK = ifelse(is.na(MethodAMPTK), "no_data", MethodAMPTK)) %>%
    mutate(., matchtest = ifelse(grepl("no_data", MethodAMPTK), gsub("AMPTK_only", "no_data", matchtest), matchtest)) %>%
    mutate(., matchtest = gsub(FALSE, "Distinct", matchtest)) %>%
    mutate(., matchtest = gsub(TRUE, "Shared", matchtest))
}

p.v.family <- family.bayescompfunction(am.taxa, vs.taxa, "AMPTK", "Vsearch", "Family")
p.b.family <- family.bayescompfunction(am.taxa, bl.taxa, "AMPTK", "Blast", "Family")
p.a.family <- family.bayescompfunction(am.taxa, bo.taxa, "AMPTK", "boldAPI", "Family")
all.family <- rbind(p.v.family, p.b.family, p.a.family)
rm(p.v.family, p.b.family, p.a.family)


## Genus level comparisons
genus.bayescompfunction <- function(data1, data2, name1, name2, level){
  pair1_df <- data1 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(family_name)) %>% filter(!is.na(genus_name)) %>%
    mutate(name1.taxstring=paste(phylum_name, class_name, order_name, family_name, genus_name, sep=";")) %>% select(HashID, name1.taxstring, Method)
  colnames(pair1_df)[3] <- paste0("Method", name1)
  pair2_df <- data2 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(family_name)) %>% filter(!is.na(genus_name)) %>%
    mutate(name2.taxstring=paste(phylum_name, class_name, order_name, family_name, genus_name, sep=";")) %>% select(HashID, name2.taxstring, Method)
  colnames(pair2_df)[3] <- paste0("Method", name2)
  tmp_match <- merge(pair1_df, pair2_df, all=TRUE)
  grouper1 <- colnames(pair1_df)[3]
  grouper2 <- colnames(pair2_df)[3]
  pair <- tmp_match %>% 
    mutate(matchtest = name1.taxstring==name2.taxstring) %>%
    group_by(matchtest, MethodAMPTK) %>%
    summarise(nMatches=n()) %>%
    mutate(Pair=paste0(name1,":",name2)) %>%
    mutate(Level=level) %>%
    ungroup(matchtest) %>% 
    mutate(., matchtest = ifelse(is.na(matchtest), "AMPTK_only", matchtest)) %>%
    mutate(MethodAMPTK=as.character(MethodAMPTK)) %>%
    mutate(., MethodAMPTK = ifelse(is.na(MethodAMPTK), "no_data", MethodAMPTK)) %>%
    mutate(., matchtest = ifelse(grepl("no_data", MethodAMPTK), gsub("AMPTK_only", "no_data", matchtest), matchtest)) %>%
    mutate(., matchtest = gsub(FALSE, "Distinct", matchtest)) %>%
    mutate(., matchtest = gsub(TRUE, "Shared", matchtest))
}

p.v.genus <- genus.bayescompfunction(am.taxa, vs.taxa, "AMPTK", "Vsearch", "Genus")
p.b.genus <- genus.bayescompfunction(am.taxa, bl.taxa, "AMPTK", "Blast", "Genus")
p.a.genus <- genus.bayescompfunction(am.taxa, bo.taxa, "AMPTK", "boldAPI", "Genus")
all.genus <- rbind(p.v.genus, p.b.genus, p.a.genus)
rm(p.v.genus, p.b.genus)

## Species level comparisons
species.bayescompfunction <- function(data1, data2, name1, name2, level){
  pair1_df <- data1 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(family_name)) %>% filter(!is.na(genus_name)) %>% filter(!is.na(species_name)) %>%
    mutate(name1.taxstring=paste(phylum_name, class_name, order_name, family_name, genus_name, species_name, sep=";")) %>% select(HashID, name1.taxstring, Method)
  colnames(pair1_df)[3] <- paste0("Method", name1)
  pair2_df <- data2 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(family_name)) %>% filter(!is.na(genus_name)) %>% filter(!is.na(species_name)) %>%
    mutate(name2.taxstring=paste(phylum_name, class_name, order_name, family_name, genus_name, species_name, sep=";")) %>% select(HashID, name2.taxstring, Method)
  colnames(pair2_df)[3] <- paste0("Method", name2)
  tmp_match <- merge(pair1_df, pair2_df, all=TRUE)
  grouper1 <- colnames(pair1_df)[3]
  grouper2 <- colnames(pair2_df)[3]
  pair <- tmp_match %>% 
    mutate(matchtest = name1.taxstring==name2.taxstring) %>%
    group_by(matchtest, MethodAMPTK) %>%
    summarise(nMatches=n()) %>%
    mutate(Pair=paste0(name1,":",name2)) %>%
    mutate(Level=level) %>%
    ungroup(matchtest) %>% 
    mutate(., matchtest = ifelse(is.na(matchtest), "AMPTK_only", matchtest)) %>%
    mutate(MethodAMPTK=as.character(MethodAMPTK)) %>%
    mutate(., MethodAMPTK = ifelse(is.na(MethodAMPTK), "no_data", MethodAMPTK)) %>%
    mutate(., matchtest = ifelse(grepl("no_data", MethodAMPTK), gsub("AMPTK_only", "no_data", matchtest), matchtest)) %>%
    mutate(., matchtest = gsub(FALSE, "Distinct", matchtest)) %>%
    mutate(., matchtest = gsub(TRUE, "Shared", matchtest))
}

p.v.species <- species.bayescompfunction(am.taxa, vs.taxa, "AMPTK", "Vsearch", "Species")
p.b.species <- species.bayescompfunction(am.taxa, bl.taxa, "AMPTK", "Blast", "Species")
p.o.species <- species.bayescompfunction(am.taxa, bo.taxa, "AMPTK", "boldAPI", "Species")
all.species <- rbind(p.v.species, p.b.species, p.o.species)
rm(p.v.species, p.b.species, p.o.species)

## merge data sets
all_comps <- rbind(all.species, all.genus, all.family, all.order, all.class)
rm(all.species, all.genus, all.family, all.order, all.class)

## plot; save as db_9d_classifierPalmerCompByTypeClassified; export at 800x600
all_comps$Level <- factor(all_comps$Level, levels = c("Class", "Order", "Family", "Genus", "Species"))
all_comps$MethodAMPTK <- factor(all_comps$MethodAMPTK, levels = c("GS", "GSL", "GDL", "SS", "US"))
all_comps$Pair <- factor(all_comps$Pair, levels = c("AMPTK:Blast", "AMPTK:Vsearch", "AMPTK:boldAPI"))

## note that we're dropping any "no_data" values, and also dropping the one instance in which there is a GDL assignment
## The GDL hit was where both the Bayesian string and USEARCH match were > 97% but differed in full taxa string, and LCA was invoked to fix
## this happened once in 13,407 times!
ggplot(data = all_comps %>% filter(MethodAMPTK != "no_data") %>% filter(MethodAMPTK !="GDL"), 
       aes(y=nMatches, x=matchtest, fill=MethodAMPTK)) + 
  geom_bar(stat="identity", color="black") + 
  facet_grid(Level ~ Pair) +
  coord_flip() +
  theme_devon() +
  labs(x="", y="distinct sequence variants", fill="AMPTK classifier method") +
  theme(legend.position="top") +
  scale_fill_manual(values = c("gray50", "gray80", "purple", "gold")) +
  scale_y_continuous(labels = comma)
  #scale_fill_manual(values = c("gray50", "gray75", "black", "white")) 