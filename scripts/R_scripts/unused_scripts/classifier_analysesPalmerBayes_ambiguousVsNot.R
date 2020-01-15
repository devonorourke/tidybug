library(tidyverse)
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


## For ASVs that are distinct (but information is retained for both Classifiers), how many of these are 'Ambiguous_taxa' in our data?
## set up the five functions
class.distinctTaxafunction <- function(data1, data2, name1, name2, level){
  pair1_df <- data1 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>%
    mutate(name1.taxstring=paste(phylum_name, class_name, sep=";")) %>% select(HashID, name1.taxstring, Method)
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

order.distinctTaxafunction <- function(data1, data2, name1, name2, level){
  pair1_df <- data1 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>%
    mutate(name1.taxstring=paste(phylum_name, class_name, order_name, sep=";")) %>% select(HashID, name1.taxstring, Method)
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


family.distinctTaxafunction <- function(data1, data2, name1, name2, level){
  pair1_df <- data1 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(family_name)) %>%
    mutate(name1.taxstring=paste(phylum_name, class_name, order_name, family_name, sep=";")) %>% select(HashID, name1.taxstring, Method)
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

genus.distinctTaxafunction <- function(data1, data2, name1, name2, level){
  pair1_df <- data1 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(family_name)) %>% filter(!is.na(genus_name)) %>%
    mutate(name1.taxstring=paste(phylum_name, class_name, order_name, family_name, genus_name, sep=";")) %>% select(HashID, name1.taxstring, Method)
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


species.distinctTaxafunction <- function(data1, data2, name1, name2, level){
  pair1_df <- data1 %>%
    filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(family_name)) %>% filter(!is.na(genus_name)) %>% filter(!is.na(species_name)) %>%
    mutate(name1.taxstring=paste(phylum_name, class_name, order_name, family_name, genus_name, species_name, sep=";")) %>% select(HashID, name1.taxstring, Method)
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


## gather data and combine
p.v.class.distinct <- class.distinctTaxafunction(am.taxa, vs.taxa, "Palmer", "Vsearch", "Class")
p.b.class.distinct <- class.distinctTaxafunction(am.taxa, vs.taxa, "Palmer", "Blast", "Class")
p.v.order.distinct <- order.distinctTaxafunction(am.taxa, vs.taxa, "Palmer", "Vsearch", "Order")
p.b.order.distinct <- order.distinctTaxafunction(am.taxa, vs.taxa, "Palmer", "Blast", "Order")
p.v.family.distinct <- family.distinctTaxafunction(am.taxa, vs.taxa, "Palmer", "Vsearch", "Family")
p.b.family.distinct <- family.distinctTaxafunction(am.taxa, vs.taxa, "Palmer", "Blast", "Family")
p.v.genus.distinct <- genus.distinctTaxafunction(am.taxa, vs.taxa, "Palmer", "Vsearch", "Genus")
p.b.genus.distinct <- genus.distinctTaxafunction(am.taxa, vs.taxa, "Palmer", "Blast", "Genus")
p.v.species.distinct <- species.distinctTaxafunction(am.taxa, vs.taxa, "Palmer", "Vsearch", "Species")
p.b.species.distinct <- species.distinctTaxafunction(am.taxa, vs.taxa, "Palmer", "Blast", "Species")

all.distinct <- rbind(p.v.class.distinct, p.b.class.distinct, p.v.order.distinct, p.b.order.distinct, p.v.family.distinct, p.b.family.distinct, p.v.genus.distinct, p.b.genus.distinct, p.v.species.distinct, p.b.species.distinct)
rm(p.v.class.distinct, p.b.class.distinct, p.v.order.distinct, p.b.order.distinct, p.v.family.distinct, p.b.family.distinct, p.v.genus.distinct, p.b.genus.distinct, p.v.species.distinct, p.b.species.distinct)

## Which Palmer classifier method was invoked when the "match" was to our Ambiguous taxa?
## we find that most of the time his LCA method wasn't invoked - it was invoked after the match...
## So our database things that it's an ambiguous level, but Palmer's database thinks there is just one representative
  ## we know this because there are only "GS" or "GSL" matches, and not "US", or "SS" matches...
backup <- all.distinct

all.distinct$AmbigTest <- grepl("Ambiguous_taxa", all.distinct$name2.taxstring)

distinct_table <- all.distinct %>% group_by(pair2Taxa, Level, AmbigTest, Pair) %>% summarise(counts=n())
distinct_table$AmbigTest <- gsub(TRUE, "Ambiguous Taxa", distinct_table$AmbigTest)
distinct_table$AmbigTest <- gsub(FALSE, "No information", distinct_table$AmbigTest)
distinct_table$Level <- factor(distinct_table$Level, levels = c("Class", "Order", "Family", "Genus", "Species"))

## plot; save as db_9e_ambiguousDataMatch; export at 550x550
ggplot(distinct_table, aes(x=Level, y=counts, fill=pair2Taxa)) +
  geom_bar(stat="identity", color="black") +
  facet_grid(AmbigTest ~ Pair) +
  coord_flip() +
  labs(x="", y="distinct sequence variants", fill="Palmer classifier method") +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = c("gray50", "gray80", "gold")) +
  theme_devon() +
  theme(legend.position = "top")

## alternative plot; save as db_9ealt_ambiguousDataMatch; export at 650x500
ggplot(distinct_table, aes(x=AmbigTest, y=counts, fill=pair2Taxa)) +
  geom_bar(stat="identity", color="black") +
  facet_grid(Level ~ Pair) +
  coord_flip() +
  labs(x="", y="distinct sequence variants", fill="Palmer classifier method") +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = c("gray50", "gray80", "gold")) +
  theme_devon() +
  theme(legend.position = "top")

