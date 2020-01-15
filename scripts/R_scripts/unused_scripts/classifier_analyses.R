library(tidyverse)
library(viridis)
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

## 1) How many unique taxa are retained for each classifier  at each Level?
## function:
uniqfunction <- function(data,dataset) {
  ClassTable <- data %>% filter(!is.na(class_name) )%>% summarise(nUniq = n_distinct(class_name)) %>% mutate(Level = "Class")
  OrderTable <- data %>% filter(!is.na(order_name) )%>% summarise(nUniq = n_distinct(order_name)) %>% mutate(Level = "Order")
  FamilyTable <- data %>% filter(!is.na(family_name) )%>% summarise(nUniq = n_distinct(family_name)) %>% mutate(Level = "Family")
  GenusTable <- data %>% filter(!is.na(genus_name) )%>% summarise(nUniq = n_distinct(genus_name)) %>% mutate(Level = "Genus")
  SpeciesTable <- data %>% filter(!is.na(species_name) )%>% summarise(nUniq = n_distinct(species_name)) %>% mutate(Level = "Species")
  fullTable <- rbind(ClassTable, OrderTable, FamilyTable, GenusTable, SpeciesTable)
  fullTable %>% mutate(Method=dataset)
}

palmer_uniq <- uniqfunction(am.taxa, "AMPTK")
vsearch_uniq <- uniqfunction(vs.taxa, "Vsearch")
blast_uniq <- uniqfunction(bl.taxa, "Blast")
api_uniq <- uniqfunction(bo.taxa, "boldAPI")

all_uniq <- rbind(palmer_uniq, vsearch_uniq, blast_uniq, api_uniq)
rm(palmer_uniq, vsearch_uniq, blast_uniq, api_uniq)

## plot; save as db_8_classifierUniq; export at 500x350
pal4 <- c("green4", "tan4", "#DA5B6AFF", "dodgerblue3")

all_uniq$Method <- factor(all_uniq$Method, levels = c("Blast", "Vsearch", "AMPTK", "boldAPI"))
all_uniq$Level <- factor(all_uniq$Level, levels = c("Class", "Order", "Family", "Genus", "Species"))
ggplot(all_uniq, aes(x=Level, y=nUniq, group=Method, color=Method, shape=Method)) +
  geom_point() +
  scale_shape_manual(values = c(0,1,2,4)) +
  scale_y_continuous(labels = comma) +
  scale_color_manual(values = pal4) +
  geom_line() +
  labs(x="", y="Distinct taxa", title="") +
  theme_devon()

## 2) What fraction of ASVs are assigned taxa information for each classifier at each Level?
## For this analysis, we're modifying the Vsearch and Blast datasets so that "Ambiguous taxa" is assigned NA
## These were converted in the first plot above

## function: counts of present/absent information per taxa Level
missfunction <- function(data,dataset) {
  Class.no <- data %>% select(class_name) %>% filter(!is.na(class_name)) %>% summarise(counts = n()) %>% mutate(Level = "Class", Value = "Present")
  Class.na <- data %>% filter(is.na(class_name)) %>% summarise(counts = n()) %>% mutate(Level = "Class", Value = "Missing")
  Class <- rbind(Class.no, Class.na) %>% spread(key=Value, value=counts)
  Order.no <- data %>% filter(!is.na(order_name)) %>% summarise(counts = n()) %>% mutate(Level = "Order", Value = "Present")
  Order.na <- data %>% filter(is.na(order_name)) %>% summarise(counts = n()) %>% mutate(Level = "Order", Value = "Missing")
  Order <- rbind(Order.no, Order.na) %>% spread(key=Value, value=counts)
  Family.no <- data %>% filter(!is.na(family_name)) %>% summarise(counts = n()) %>% mutate(Level = "Family", Value = "Present")
  Family.na <- data %>% filter(is.na(family_name)) %>% summarise(counts = n()) %>% mutate(Level = "Family", Value = "Missing")
  Family <- rbind(Family.no, Family.na) %>% spread(key=Value, value=counts)
  Genus.no <- data %>% filter(!is.na(genus_name)) %>% summarise(counts = n()) %>% mutate(Level = "Genus", Value = "Present")
  Genus.na <- data %>% filter(is.na(genus_name)) %>% summarise(counts = n()) %>% mutate(Level = "Genus", Value = "Missing")
  Genus <- rbind(Genus.no, Genus.na) %>% spread(key=Value, value=counts)
  Species.no <- data %>% filter(!is.na(species_name)) %>% summarise(counts = n()) %>% mutate(Level = "Species", Value = "Present")
  Species.na <- data %>% filter(is.na(species_name)) %>% summarise(counts = n()) %>% mutate(Level = "Species", Value = "Missing")
  Species <- rbind(Species.no, Species.na) %>% spread(key=Value, value=counts)
  tmp <- rbind(Class, Order, Family, Genus, Species)
  tmp <- tmp %>% mutate(Database=dataset)
  tmp$pPresent <- tmp$Present / (tmp$Present + tmp$Missing)
  tmp
}

palmer_miss <- missfunction(am.taxa, "AMPTK")
vsearch_miss <- missfunction(vs.taxa, "Vsearch")
blast_miss <- missfunction(bl.taxa, "Blast")
api_miss <- missfunction(bo.taxa, "boldAPI")

all_miss <- rbind(palmer_miss, vsearch_miss, blast_miss, api_miss)
rm(palmer_miss, vsearch_miss, blast_miss, api_miss)

alt_miss_df <- all_miss %>% gather(key = Type, value = nTaxa, Present, Missing)
alt_miss_df$Level <- factor(alt_miss_df$Level, levels = c("Class", "Order", "Family", "Genus", "Species"))
alt_miss_df$Database <- factor(alt_miss_df$Database, levels = c("Blast", "Vsearch", "AMPTK", "boldAPI"))
alt_miss_df$Type <- factor(alt_miss_df$Type, levels = c("Missing", "Present"))

## plot; save as db_9a_classifierMissingness; export at 600x600
ggplot(alt_miss_df, aes(y=nTaxa, fill=Type, x=Level)) +
  geom_bar(stat="identity", color="gray30") +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = c("gray85", "gray45")) +
  facet_wrap(Database ~ ., nrow = 4) +
  labs(x="", y="distinct sequences", title="", fill="Taxa\ninformation") +
  theme_devon() +
  coord_flip() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle=22.5, hjust=1))


## parse the AMPTK dataset at each Level, assigning a secondary Fill to account for the $Method used to classifiy
## similar function but uses Labeler to create distinct Fill terms for plot

altmissfunction <- function(data,dataset) {
  Class.no <- data%>% filter(!is.na(class_name)) %>% mutate(Grouper=paste(Classifier, Method, sep=".")) %>% group_by(Grouper) %>% summarise(counts = n()) %>% mutate(Level = "Class", Value = "information\nPresent")
  Class.na <- data%>% filter(is.na(class_name)) %>% mutate(Grouper=paste(Classifier, Method, sep=".")) %>% group_by(Grouper) %>% summarise(counts = n()) %>% mutate(Level = "Class", Value = "information\nMissing")
  Class <- rbind(Class.no, Class.na)
  Order.no <- data %>% filter(!is.na(order_name)) %>% mutate(Grouper=paste(Classifier, Method, sep=".")) %>% group_by(Grouper) %>% summarise(counts = n()) %>% mutate(Level = "Order", Value = "information\nPresent")
  Order.na <- data %>% filter(is.na(order_name)) %>% mutate(Grouper=paste(Classifier, Method, sep=".")) %>% group_by(Grouper) %>% summarise(counts = n()) %>% mutate(Level = "Order", Value = "information\nMissing")
  Order <- rbind(Order.no, Order.na)
  Family.no <- data %>% filter(!is.na(family_name)) %>% mutate(Grouper=paste(Classifier, Method, sep=".")) %>% group_by(Grouper) %>% summarise(counts = n()) %>% mutate(Level = "Family", Value = "information\nPresent")
  Family.na <- data %>% filter(is.na(family_name)) %>% mutate(Grouper=paste(Classifier, Method, sep=".")) %>% group_by(Grouper) %>% summarise(counts = n()) %>% mutate(Level = "Family", Value = "information\nMissing")
  Family <- rbind(Family.no, Family.na)
  Genus.no <- data %>% filter(!is.na(genus_name)) %>% mutate(Grouper=paste(Classifier, Method, sep=".")) %>% group_by(Grouper) %>% summarise(counts = n()) %>% mutate(Level = "Genus", Value = "information\nPresent")
  Genus.na <- data %>% filter(is.na(genus_name)) %>% mutate(Grouper=paste(Classifier, Method, sep=".")) %>% group_by(Grouper) %>% summarise(counts = n()) %>% mutate(Level = "Genus", Value = "information\nMissing")
  Genus <- rbind(Genus.no, Genus.na)
  Species.no <- data %>% filter(!is.na(species_name)) %>% mutate(Grouper=paste(Classifier, Method, sep=".")) %>% group_by(Grouper) %>% summarise(counts = n()) %>% mutate(Level = "Species", Value = "information\nPresent")
  Species.na <- data %>% filter(is.na(species_name)) %>% mutate(Grouper=paste(Classifier, Method, sep=".")) %>% group_by(Grouper) %>% summarise(counts = n()) %>% mutate(Level = "Species", Value = "information\nMissing")
  Species <- rbind(Species.no, Species.na)
  tmp <- rbind(Class, Order, Family, Genus, Species)
  tmp <- tmp %>% mutate(Database=dataset)
  tmp
}


palmer_altmiss <- altmissfunction(am.taxa, "AMPTK")
vsearch_altmiss <- altmissfunction(vs.taxa, "Vsearch")
blast_altmiss <- altmissfunction(bl.taxa, "Blast")
api_altmiss <- altmissfunction(bo.taxa, "boldAPI")
all_altmiss <- rbind(palmer_altmiss, vsearch_altmiss, blast_altmiss, api_altmiss)
rm(palmer_altmiss, vsearch_altmiss, blast_altmiss, api_altmiss)
all_altmiss <- all_altmiss %>% separate(., col = Grouper, into = c("Classifier", "Method"), sep="\\.")

# set levels:
all_altmiss$Level <- factor(all_altmiss$Level, levels = c("Class", "Order", "Family", "Genus", "Species"))
all_altmiss$Database <- factor(all_altmiss$Database, levels = c("Blast", "Vsearch", "AMPTK", "boldAPI"))
#all_altmiss$Method <- factor(all_altmiss$Method, levels = c("blast", "vsearch", "GS", "GSL", "SS", "US", "boldAPI", "GDL"))
all_altmiss$Method <- factor(all_altmiss$Method, levels = c("GS", "GSL", "SS", "US", "GDL", "blast", "vsearch", 'boldAPI'))

## plot; save as db_9b_classifierMissingness; export at 800x600
ggplot(all_altmiss %>% filter(Method != "GDL"), 
       aes(x=Level, y=counts, fill=Method)) +
  geom_bar(stat="identity", color="black", size=0.75) +
  facet_grid(Database ~ Value) +
  coord_flip() +
  scale_fill_manual(values = c("gray35", "gray75", "purple", "gold", "green4", "tan4", "dodgerblue3"),
                    labels = c("amptk\nGS", "amptk\nGSL", "amptk\nSS", "amptk\nUS", "blast", "vsearch", "boldAPI")) +
  theme_devon() +
  labs(x="", y="distinct sequences", color="Taxa information", fill="Classifier\nMethod") +
  theme(legend.position = "top") +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))


## 3) What fraction of ASVs are identical among Classifiers, per Level?
## we'll want to remove the instances in which there is something compared to an NA...
## comparisons for Class, Order, Family, Genus, and Species levels

## Class level comparisons
class.compfunction <- function(data1, data2, name1, name2, level){
  pair1_df <- data1 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>%
    mutate(name1.taxstring=paste(phylum_name, class_name, sep=";")) %>% select(HashID, name1.taxstring)
  pair2_df <- data2 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>%
    mutate(name2.taxstring=paste(phylum_name, class_name, sep=";")) %>% select(HashID, name2.taxstring)
  tmp_match <- merge(pair1_df, pair2_df, all=TRUE) %>% filter(complete.cases(.))
  pair <- tmp_match %>% 
    mutate(matchtest = name1.taxstring==name2.taxstring) %>%
    group_by(matchtest) %>% 
    summarise(nMatches=n()) %>%
    mutate(Pair=paste0(name1,":",name2)) %>%
    mutate(Level=level)
}

p.v.class <- class.compfunction(am.taxa, vs.taxa, "AMPTK", "Vsearch", "Class")
p.b.class <- class.compfunction(am.taxa, bl.taxa, "AMPTK", "Blast", "Class")
p.a.class <- class.compfunction(am.taxa, bo.taxa, "AMPTK", "boldAPI", "Class")
v.b.class <- class.compfunction(vs.taxa, bl.taxa, "Vsearch", "Blast", "Class")
v.a.class <- class.compfunction(vs.taxa, bl.taxa, "Vsearch", "boldAPI", "Class")
b.a.class <- class.compfunction(vs.taxa, bl.taxa, "Blast", "boldAPI", "Class")

all.class <- rbind(p.v.class, p.b.class, p.a.class, v.b.class, v.a.class, b.a.class)
rm(p.v.class, p.b.class, p.a.class, v.b.class, v.a.class, b.a.class)

## Order level comparisons
order.compfunction <- function(data1, data2, name1, name2, level){
  pair1_df <- data1 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>%
    mutate(name1.taxstring=paste(phylum_name, class_name, order_name, sep=";")) %>% select(HashID, name1.taxstring)
  pair2_df <- data2 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>%
    mutate(name2.taxstring=paste(phylum_name, class_name, order_name, sep=";")) %>% select(HashID, name2.taxstring)
  tmp_match <- merge(pair1_df, pair2_df, all=TRUE) %>% filter(complete.cases(.))
  pair <- tmp_match %>% 
    mutate(matchtest = name1.taxstring==name2.taxstring) %>%
    group_by(matchtest) %>% 
    summarise(nMatches=n()) %>%
    mutate(Pair=paste0(name1,":",name2)) %>%
    mutate(Level=level)
}


p.v.order <- order.compfunction(am.taxa, vs.taxa, "AMPTK", "Vsearch", "Order")
p.b.order <- order.compfunction(am.taxa, bl.taxa, "AMPTK", "Blast", "Order")
p.a.order <- order.compfunction(am.taxa, bo.taxa, "AMPTK", "boldAPI", "Order")
v.b.order <- order.compfunction(vs.taxa, bl.taxa, "Vsearch", "Blast", "Order")
v.a.order <- order.compfunction(vs.taxa, bl.taxa, "Vsearch", "boldAPI", "Order")
b.a.order <- order.compfunction(vs.taxa, bl.taxa, "Blast", "boldAPI", "Order")

all.order <- rbind(p.v.order, p.b.order, p.a.order, v.b.order, v.a.order, b.a.order)
rm(p.v.order, p.b.order, p.a.order, v.b.order, v.a.order, b.a.order)

## Family level comparisons
family.compfunction <- function(data1, data2, name1, name2, level){
  pair1_df <- data1 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(family_name)) %>%
    mutate(name1.taxstring=paste(phylum_name, class_name, order_name, family_name, sep=";")) %>% select(HashID, name1.taxstring)
  pair2_df <- data2 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(order_name)) %>%
    mutate(name2.taxstring=paste(phylum_name, class_name, order_name, family_name, sep=";")) %>% select(HashID, name2.taxstring)
  tmp_match <- merge(pair1_df, pair2_df, all=TRUE) %>% filter(complete.cases(.))
  pair <- tmp_match %>% 
    mutate(matchtest = name1.taxstring==name2.taxstring) %>%
    group_by(matchtest) %>% 
    summarise(nMatches=n()) %>%
    mutate(Pair=paste0(name1,":",name2)) %>%
    mutate(Level=level)
}

p.v.family <- family.compfunction(am.taxa, vs.taxa, "AMPTK", "Vsearch", "Family")
p.b.family <- family.compfunction(am.taxa, bl.taxa, "AMPTK", "Blast", "Family")
p.a.family <- family.compfunction(am.taxa, bo.taxa, "AMPTK", "boldAPI", "Family")
v.b.family <- family.compfunction(vs.taxa, bl.taxa, "Vsearch", "Blast", "Family")
v.a.family <- family.compfunction(vs.taxa, bl.taxa, "Vsearch", "boldAPI", "Family")
b.a.family <- family.compfunction(vs.taxa, bl.taxa, "Blast", "boldAPI", "Family")

all.family <- rbind(p.v.family, p.b.family, p.a.family, v.b.family, v.a.family, b.a.family)
rm(p.v.family, p.b.family, p.a.family, v.b.family, v.a.family, b.a.family)

## Genus level comparisons
genus.compfunction <- function(data1, data2, name1, name2, level){
  pair1_df <- data1 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(family_name)) %>% filter(!is.na(genus_name)) %>%
    mutate(name1.taxstring=paste(phylum_name, class_name, order_name, family_name, genus_name, sep=";")) %>% select(HashID, name1.taxstring)
  pair2_df <- data2 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(genus_name)) %>%
    mutate(name2.taxstring=paste(phylum_name, class_name, order_name, family_name, genus_name, sep=";")) %>% select(HashID, name2.taxstring)
  tmp_match <- merge(pair1_df, pair2_df, all=TRUE) %>% filter(complete.cases(.))
  pair <- tmp_match %>% 
    mutate(matchtest = name1.taxstring==name2.taxstring) %>%
    group_by(matchtest) %>% 
    summarise(nMatches=n()) %>%
    mutate(Pair=paste0(name1,":",name2)) %>%
    mutate(Level=level)
}

p.v.genus <- genus.compfunction(am.taxa, vs.taxa, "AMPTK", "Vsearch", "Genus")
p.b.genus <- genus.compfunction(am.taxa, bl.taxa, "AMPTK", "Blast", "Genus")
p.a.genus <- genus.compfunction(am.taxa, bo.taxa, "AMPTK", "boldAPI", "Genus")
v.b.genus <- genus.compfunction(vs.taxa, bl.taxa, "Vsearch", "Blast", "Genus")
v.a.genus <- genus.compfunction(vs.taxa, bl.taxa, "Vsearch", "boldAPI", "Genus")
b.a.genus <- genus.compfunction(vs.taxa, bl.taxa, "Blast", "boldAPI", "Genus")

all.genus <- rbind(p.v.genus, p.b.genus, p.a.genus, v.b.genus, v.a.genus, b.a.genus)
rm(p.v.genus, p.b.genus, p.a.genus, v.b.genus, v.a.genus, b.a.genus)


## Species level comparisons
species.compfunction <- function(data1, data2, name1, name2, level){
  pair1_df <- data1 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(family_name)) %>% filter(!is.na(genus_name)) %>% filter(!is.na(species_name)) %>%
    mutate(name1.taxstring=paste(phylum_name, class_name, order_name, family_name, genus_name, species_name, sep=";")) %>% select(HashID, name1.taxstring)
  pair2_df <- data2 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(genus_name)) %>% filter(!is.na(species_name)) %>%
    mutate(name2.taxstring=paste(phylum_name, class_name, order_name, family_name, genus_name, species_name, sep=";")) %>% select(HashID, name2.taxstring)
  tmp_match <- merge(pair1_df, pair2_df, all=TRUE) %>% filter(complete.cases(.))
  pair <- tmp_match %>% 
    mutate(matchtest = name1.taxstring==name2.taxstring) %>%
    group_by(matchtest) %>% 
    summarise(nMatches=n()) %>%
    mutate(Pair=paste0(name1,":",name2)) %>%
    mutate(Level=level)
}

p.v.species <- species.compfunction(am.taxa, vs.taxa, "AMPTK", "Vsearch", "Species")
p.b.species <- species.compfunction(am.taxa, bl.taxa, "AMPTK", "Blast", "Species")
p.a.species <- species.compfunction(am.taxa, bo.taxa, "AMPTK", "boldAPI", "Species")
v.b.species <- species.compfunction(vs.taxa, bl.taxa, "Vsearch", "Blast", "Species")
v.a.species <- species.compfunction(vs.taxa, bl.taxa, "Vsearch", "boldAPI", "Species")
b.a.species <- species.compfunction(vs.taxa, bl.taxa, "Blast", "boldAPI", "Species")

all.species <- rbind(p.v.species, p.b.species, p.a.species, v.b.species, v.a.species, b.a.species)
rm(p.v.species, p.b.species, p.a.species, v.b.species, v.a.species, b.a.species)

## combine all comparisons at all Levels
all.comps <- rbind(all.class, all.order, all.family, all.genus, all.species)
rm(all.class, all.order, all.family, all.genus, all.species)
## write to disk:
# write.csv(all.comps, file="~/Repos/tidybug/data/databases/classifier_taxacomparisons.csv", quote = FALSE, row.names = FALSE)

## plot; save as db_9c_classificationDistinctTaxa; export at 600x600
all.comps$matchtest <- gsub("TRUE", "Shared", all.comps$matchtest)
all.comps$matchtest <- gsub("FALSE", "Distinct", all.comps$matchtest)
all.comps$Level <- factor(all.comps$Level, levels = c("Class", "Order", "Family", "Genus", "Species"))
all.comps$Pair <- as.factor(all.comps$Pair)
all.comps$Pair <- factor(all.comps$Pair, 
                         levels = c("AMPTK:Blast", "AMPTK:Vsearch", "AMPTK:boldAPI", 
                                    "Blast:boldAPI", "Vsearch:Blast", "Vsearch:boldAPI"))

ggplot(all.comps, aes(y=nMatches, x=Pair, fill=matchtest)) +
  geom_bar(stat="identity") +
  facet_grid(Level ~ .) +
  coord_flip() +
  labs(x="", y="distinct sequence variants", fill="Taxa Information") +
  scale_fill_manual(values=c("gray30", "gray70")) +
  theme_devon() +
  theme(legend.position = "top")

rm(all.comps)

## create data.farme of all taxa that are NOT SHARED at each level:

## Class level
class.distinctTaxafunction <- function(data1, data2, name1, name2, level){
  pair1_df <- data1 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>%
    mutate(name1.taxstring=paste(phylum_name, class_name, sep=";")) %>% select(HashID, name1.taxstring)
  pair2_df <- data2 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>%
    mutate(name2.taxstring=paste(phylum_name, class_name, sep=";")) %>% select(HashID, name2.taxstring)
  tmp_match <- merge(pair1_df, pair2_df, all=TRUE) %>% filter(complete.cases(.))
  pair <- tmp_match %>% 
    mutate(matchtest = name1.taxstring==name2.taxstring) %>%
    filter(matchtest==FALSE) %>% 
    mutate(Pair=paste0(name1,":",name2)) %>%
    mutate(Level=level)
  colnames(pair)[2] <- "pair1Taxa"
  colnames(pair)[3] <- "pair2Taxa"
  pair
}

p.v.class.distinct <- class.distinctTaxafunction(am.taxa, vs.taxa, "AMPTK", "Vsearch", "Class")
p.b.class.distinct <- class.distinctTaxafunction(am.taxa, bl.taxa, "AMPTK", "Blast", "Class")
p.a.class.distinct <- class.distinctTaxafunction(am.taxa, bl.taxa, "AMPTK", "boldAPI", "Class")
v.b.class.distinct <- class.distinctTaxafunction(bl.taxa, vs.taxa, "Vsearch", "Blast", "Class")
v.a.class.distinct <- class.distinctTaxafunction(bl.taxa, vs.taxa, "Vsearch", "boldAPI", "Class")
b.a.class.distinct <- class.distinctTaxafunction(bl.taxa, vs.taxa, "Blast", "boldAPI", "Class")
all.class.distinct <- rbind(p.v.class.distinct, p.b.class.distinct, p.a.class.distinct, v.b.class.distinct, v.a.class.distinct, b.a.class.distinct)
rm(p.v.class.distinct, p.b.class.distinct, p.a.class.distinct, v.b.class.distinct, v.a.class.distinct, b.a.class.distinct)


## Order level
order.distinctTaxafunction <- function(data1, data2, name1, name2, level){
  pair1_df <- data1 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>%
    mutate(name1.taxstring=paste(phylum_name, class_name, order_name, sep=";")) %>% select(HashID, name1.taxstring)
  pair2_df <- data2 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>%
    mutate(name2.taxstring=paste(phylum_name, class_name, order_name, sep=";")) %>% select(HashID, name2.taxstring)
  tmp_match <- merge(pair1_df, pair2_df, all=TRUE) %>% filter(complete.cases(.))
  pair <- tmp_match %>% 
    mutate(matchtest = name1.taxstring==name2.taxstring) %>%
    filter(matchtest==FALSE) %>% 
    mutate(Pair=paste0(name1,":",name2)) %>%
    mutate(Level=level)
  colnames(pair)[2] <- "pair1Taxa"
  colnames(pair)[3] <- "pair2Taxa"
  pair
}

p.v.order.distinct <- order.distinctTaxafunction(am.taxa, vs.taxa, "AMPTK", "Vsearch", "Order")
p.b.order.distinct <- order.distinctTaxafunction(am.taxa, bl.taxa, "AMPTK", "Blast", "Order")
p.a.order.distinct <- order.distinctTaxafunction(am.taxa, bl.taxa, "AMPTK", "boldAPI", "Order")
v.b.order.distinct <- order.distinctTaxafunction(bl.taxa, vs.taxa, "Vsearch", "Blast", "Order")
v.a.order.distinct <- order.distinctTaxafunction(bl.taxa, vs.taxa, "Vsearch", "boldAPI", "Order")
b.a.order.distinct <- order.distinctTaxafunction(bl.taxa, vs.taxa, "Blast", "boldAPI", "Order")
all.order.distinct <- rbind(p.v.order.distinct, p.b.order.distinct, p.a.order.distinct, v.b.order.distinct, v.a.order.distinct, b.a.order.distinct)
rm(p.v.order.distinct, p.b.order.distinct, p.a.order.distinct, v.b.order.distinct, v.a.order.distinct, b.a.order.distinct)

## Family level
family.distinctTaxafunction <- function(data1, data2, name1, name2, level){
  pair1_df <- data1 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(family_name)) %>%
    mutate(name1.taxstring=paste(phylum_name, class_name, order_name, family_name, sep=";")) %>% select(HashID, name1.taxstring)
  pair2_df <- data2 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(family_name)) %>%
    mutate(name2.taxstring=paste(phylum_name, class_name, order_name, family_name, sep=";")) %>% select(HashID, name2.taxstring)
  tmp_match <- merge(pair1_df, pair2_df, all=TRUE) %>% filter(complete.cases(.))
  pair <- tmp_match %>% 
    mutate(matchtest = name1.taxstring==name2.taxstring) %>%
    filter(matchtest==FALSE) %>% 
    mutate(Pair=paste0(name1,":",name2)) %>%
    mutate(Level=level)
  colnames(pair)[2] <- "pair1Taxa"
  colnames(pair)[3] <- "pair2Taxa"
  pair
}

p.v.family.distinct <- family.distinctTaxafunction(am.taxa, vs.taxa, "AMPTK", "Vsearch", "Family")
p.b.family.distinct <- family.distinctTaxafunction(am.taxa, bl.taxa, "AMPTK", "Blast", "Family")
p.a.family.distinct <- family.distinctTaxafunction(am.taxa, bl.taxa, "AMPTK", "boldAPI", "Family")
v.b.family.distinct <- family.distinctTaxafunction(bl.taxa, vs.taxa, "Vsearch", "Blast", "Family")
v.a.family.distinct <- family.distinctTaxafunction(bl.taxa, vs.taxa, "Vsearch", "boldAPI", "Family")
b.a.family.distinct <- family.distinctTaxafunction(bl.taxa, vs.taxa, "Blast", "boldAPI", "Family")
all.family.distinct <- rbind(p.v.family.distinct, p.b.family.distinct, p.a.family.distinct, v.b.family.distinct, v.a.family.distinct, b.a.family.distinct)
rm(p.v.family.distinct, p.b.family.distinct, p.a.family.distinct, v.b.family.distinct, v.a.family.distinct, b.a.family.distinct)


## Genus level
genus.distinctTaxafunction <- function(data1, data2, name1, name2, level){
  pair1_df <- data1 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(family_name)) %>% filter(!is.na(genus_name)) %>%
    mutate(name1.taxstring=paste(phylum_name, class_name, order_name, family_name, genus_name, sep=";")) %>% select(HashID, name1.taxstring)
  pair2_df <- data2 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(family_name)) %>% filter(!is.na(genus_name)) %>%
    mutate(name2.taxstring=paste(phylum_name, class_name, order_name, family_name, genus_name, sep=";")) %>% select(HashID, name2.taxstring)
  tmp_match <- merge(pair1_df, pair2_df, all=TRUE) %>% filter(complete.cases(.))
  pair <- tmp_match %>% 
    mutate(matchtest = name1.taxstring==name2.taxstring) %>%
    filter(matchtest==FALSE) %>% 
    mutate(Pair=paste0(name1,":",name2)) %>%
    mutate(Level=level)
  colnames(pair)[2] <- "pair1Taxa"
  colnames(pair)[3] <- "pair2Taxa"
  pair
}

p.v.genus.distinct <- genus.distinctTaxafunction(am.taxa, vs.taxa, "AMPTK", "Vsearch", "Genus")
p.b.genus.distinct <- genus.distinctTaxafunction(am.taxa, bl.taxa, "AMPTK", "Blast", "Genus")
p.a.genus.distinct <- genus.distinctTaxafunction(am.taxa, bl.taxa, "AMPTK", "boldAPI", "Genus")
v.b.genus.distinct <- genus.distinctTaxafunction(bl.taxa, vs.taxa, "Vsearch", "Blast", "Genus")
v.a.genus.distinct <- genus.distinctTaxafunction(bl.taxa, vs.taxa, "Vsearch", "boldAPI", "Genus")
b.a.genus.distinct <- genus.distinctTaxafunction(bl.taxa, vs.taxa, "Blast", "boldAPI", "Genus")
all.genus.distinct <- rbind(p.v.genus.distinct, p.b.genus.distinct, p.a.genus.distinct, v.b.genus.distinct, v.a.genus.distinct, b.a.genus.distinct)
rm(p.v.genus.distinct, p.b.genus.distinct, p.a.genus.distinct, v.b.genus.distinct, v.a.genus.distinct, b.a.genus.distinct)

## Species level
species.distinctTaxafunction <- function(data1, data2, name1, name2, level){
  pair1_df <- data1 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(family_name)) %>% filter(!is.na(genus_name)) %>% filter(!is.na(species_name)) %>%
    mutate(name1.taxstring=paste(phylum_name, class_name, order_name, family_name, genus_name, species_name, sep=";")) %>% select(HashID, name1.taxstring)
  pair2_df <- data2 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(family_name)) %>% filter(!is.na(genus_name)) %>% filter(!is.na(species_name)) %>%
    mutate(name2.taxstring=paste(phylum_name, class_name, order_name, family_name, genus_name, species_name, sep=";")) %>% select(HashID, name2.taxstring)
  tmp_match <- merge(pair1_df, pair2_df, all=TRUE) %>% filter(complete.cases(.))
  pair <- tmp_match %>% 
    mutate(matchtest = name1.taxstring==name2.taxstring) %>%
    filter(matchtest==FALSE) %>% 
    mutate(Pair=paste0(name1,":",name2)) %>%
    mutate(Level=level)
  colnames(pair)[2] <- "pair1Taxa"
  colnames(pair)[3] <- "pair2Taxa"
  pair
}


p.v.species.distinct <- species.distinctTaxafunction(am.taxa, vs.taxa, "AMPTK", "Vsearch", "Species")
p.b.species.distinct <- species.distinctTaxafunction(am.taxa, bl.taxa, "AMPTK", "Blast", "Species")
p.a.species.distinct <- species.distinctTaxafunction(am.taxa, bl.taxa, "AMPTK", "boldAPI", "Species")
v.b.species.distinct <- species.distinctTaxafunction(bl.taxa, vs.taxa, "Vsearch", "Blast", "Species")
v.a.species.distinct <- species.distinctTaxafunction(bl.taxa, vs.taxa, "Vsearch", "boldAPI", "Species")
b.a.species.distinct <- species.distinctTaxafunction(bl.taxa, vs.taxa, "Blast", "boldAPI", "Species")
all.species.distinct <- rbind(p.v.species.distinct, p.b.species.distinct, p.a.species.distinct, v.b.species.distinct, v.a.species.distinct, b.a.species.distinct)
rm(p.v.species.distinct, p.b.species.distinct, p.a.species.distinct, v.b.species.distinct, v.a.species.distinct, b.a.species.distinct)

## combine all datasets and write to disk:
all.distinct <- rbind(all.class.distinct, all.order.distinct, all.family.distinct, all.genus.distinct, all.species.distinct)
rm(all.class.distinct, all.order.distinct, all.family.distinct, all.genus.distinct, all.species.distinct)
# write.csv(all.distinct, file="~/Repos/tidybug/data/databases/classifier_allDistinctTaxa.csv", row.names = FALSE, quote = TRUE)
