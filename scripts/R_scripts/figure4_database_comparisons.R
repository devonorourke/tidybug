library(tidyverse)
library(scales)
library(ggpubr)
library(viridis)
library(UpSetR)
library(grid)
library(formattable)

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

## import and process datasets one at a time to conserve disk space
## 1) tidy taxa data as needed
## 2) calculate number of unique taxa, per level
## 3) calculate number of most frequent taxa per Level (Class through Species)
## 4) calculate missingness per taxa Level; 

## functions applied:
## 2) calculating unique taxa:
uniqfunction <- function(data,dataset) {
  ClassTable <- data %>% filter(!is.na(class_name) )%>% summarise(nUniq = n_distinct(class_name)) %>% mutate(Level = "Class")
  OrderTable <- data %>% filter(!is.na(order_name) )%>% summarise(nUniq = n_distinct(order_name)) %>% mutate(Level = "Order")
  FamilyTable <- data %>% filter(!is.na(family_name) )%>% summarise(nUniq = n_distinct(family_name)) %>% mutate(Level = "Family")
  GenusTable <- data %>% filter(!is.na(genus_name) )%>% summarise(nUniq = n_distinct(genus_name)) %>% mutate(Level = "Genus")
  SpeciesTable <- data %>% filter(!is.na(species_name) )%>% summarise(nUniq = n_distinct(species_name)) %>% mutate(Level = "Species")
  fullTable <- rbind(ClassTable, OrderTable, FamilyTable, GenusTable, SpeciesTable)
  fullTable %>% mutate(Database=dataset)
}


## 3) calculating most abundant taxa:
uniqCountfunction <- function(data, database) {
  class_counts <- data %>% select(class_name) %>% filter(!is.na(class_name)) %>% group_by(class_name) %>% summarise(counts=n()) %>% arrange(-counts)
  class_top20 <- class_counts %>% top_n(20) %>% slice(., 1:20)    ## ADDED 'slice' to avoid instances where a tie existed in 20th row (with >20th row)
  class_filter <- as.character(class_top20$class_name)
  class_bottom <- class_counts %>% filter(!class_name %in% class_filter) %>% summarise(counts=sum(counts)) %>% mutate(class_name="others")
  class_uniq <- rbind(class_top20, class_bottom) %>% filter(counts != 0) %>% mutate(Level="Class") %>% mutate(Database=database)
  colnames(class_uniq)[1] <- "TaxaName"
  order_counts <- data %>% select(order_name) %>% filter(!is.na(order_name)) %>% group_by(order_name) %>% summarise(counts=n()) %>% arrange(-counts)
  order_top20 <- order_counts %>% top_n(20) %>% slice(., 1:20)
  order_filter <- as.character(order_top20$order_name)
  order_bottom <- order_counts %>% filter(!order_name %in% order_filter) %>% summarise(counts=sum(counts)) %>% mutate(order_name="others")
  order_uniq <- rbind(order_top20, order_bottom) %>% mutate(Level="Order") %>% mutate(Database=database)
  colnames(order_uniq)[1] <- "TaxaName"
  family_counts <- data %>% select(family_name) %>% filter(!is.na(family_name)) %>% group_by(family_name) %>% summarise(counts=n()) %>% arrange(-counts)
  family_top20 <- family_counts %>% top_n(20) %>% slice(., 1:20)
  family_filter <- as.character(family_top20$family_name)
  family_bottom <- family_counts %>% filter(!family_name %in% family_filter) %>% summarise(counts=sum(counts)) %>% mutate(family_name="others")
  family_uniq <- rbind(family_top20, family_bottom) %>% mutate(Level="Family") %>% mutate(Database=database)
  colnames(family_uniq)[1] <- "TaxaName"
  genus_counts <- data %>% select(genus_name) %>% filter(!is.na(genus_name)) %>% group_by(genus_name) %>% summarise(counts=n()) %>% arrange(-counts)
  genus_top20 <- genus_counts %>% top_n(20) %>% slice(., 1:20)
  genus_filter <- as.character(genus_top20$genus_name)
  genus_bottom <- genus_counts %>% filter(!genus_name %in% genus_filter) %>% summarise(counts=sum(counts)) %>% mutate(genus_name="others")
  genus_uniq <- rbind(genus_top20, genus_bottom) %>% mutate(Level="Genus") %>% mutate(Database=database)
  colnames(genus_uniq)[1] <- "TaxaName"
  species_counts <- data %>% select(species_name) %>% filter(!is.na(species_name)) %>% group_by(species_name) %>% summarise(counts=n()) %>% arrange(-counts)
  species_top20 <- species_counts %>% top_n(20)
  species_filter <- as.character(species_top20$species_name)
  species_bottom <- species_counts %>% filter(!species_name %in% species_filter) %>% summarise(counts=sum(counts)) %>% mutate(species_name="others")
  species_uniq <- rbind(species_top20, species_bottom) %>% mutate(Level="Species") %>% mutate(Database=database)
  colnames(species_uniq)[1] <- "TaxaName"
  rbind(class_uniq, order_uniq, family_uniq, genus_uniq, species_uniq)
}

## 4) counts of present/absent information per taxa Level
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


## ----  palmer data:  ----  ##
palmer_taxa <- read_csv(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/palmer.taxa.txt.gz", col_names = FALSE)
colnames(palmer_taxa) <- c("class_name", "order_name", "family_name", "genus_name", "species_name")
palmer_taxa <- as.data.frame(apply(palmer_taxa, 2, function(y) (gsub(".:", "", y))))  ## remove prefixes for each taxa (ex. 'p:', or 'c:', etc.)
## apply three function
palmer_uniq <- uniqfunction(palmer_taxa, "Palmer")
palmer_abund <- uniqCountfunction(palmer_taxa, "Palmer")
palmer_miss <- missfunction(palmer_taxa, "Palmer")
#rm(palmer_taxa)

## ----  porter data:  ----  ##
porter_taxa <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/porter.taxa.txt.gz", col_names = FALSE, delim = ';')
colnames(porter_taxa) <- c("sequenceID", "class_name", "order_name", "family_name", "genus_name", "species_name")
porter_taxa <- as.data.frame(apply(porter_taxa, 2, function(y) (gsub("_", " ", y))))    ## replace underscore in species names with spaces to match other database naming conventions
porter_taxa <- as.data.frame(apply(porter_taxa, 2, function(y) (gsub("Ambiguous taxa", NA, y))))    ## substitute the Ambiguous Taxa designation with NA (it's unknown!)
porter_uniq <- uniqfunction(porter_taxa, "Porter")
porter_abund <- uniqCountfunction(porter_taxa, "Porter")
porter_miss <- missfunction(porter_taxa, "Porter")
#rm(porter_taxa)

## ----  (our) tidybug data:  ----  ##
tidybug_taxa <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/derep.taxa.txt.gz", delim = ";", col_names = FALSE)
colnames(tidybug_taxa) <- c("class_name", "order_name", "family_name", "genus_name", "species_name")
tidybug_taxa <- as.data.frame(apply(tidybug_taxa, 2, function(y) (gsub("Ambiguous_taxa", NA, y))))    ## substitute the Ambiguous Taxa designation with NA (it's unknown!)
tidybug_uniq <- uniqfunction(tidybug_taxa, "tidybug")
tidybug_abund <- uniqCountfunction(tidybug_taxa, "tidybug")
tidybug_miss <- missfunction(tidybug_taxa, "tidybug")
#rm(tidybug_taxa)

## combine the three sets of data
uniq_df <- rbind(palmer_uniq, porter_uniq, tidybug_uniq)
abund_df <- rbind(palmer_abund, porter_abund, tidybug_abund)
miss_df <- rbind(palmer_miss, porter_miss, tidybug_miss)

## color palette for plot
vpal3 <- viridis(3, option = "plasma", end=0.85)
dpal3 <- c("#7ab7bf", "#ffcc33", "#e57375")

# set levels for plot
uniq_df$Database <- factor(uniq_df$Database, levels = c("tidybug", "Palmer", "Porter"))
uniq_df$Level <- factor(uniq_df$Level, levels = c("Species", "Genus", "Family", "Order", "Class"))

## set up comma-separated labels for bar val text
uniq_df <- uniq_df %>% 
  mutate(LabelVals = comma(uniq_df$nUniq, digits = 0))

## generate nUniq table for supplementary material:
tmpuniq <- uniq_df %>% spread(data = ., key = Level, value = nUniq)
write.csv(tmpuniq, file="~/Repos/tidybug/data/text_tables/database_evals/uniqTaxa.csv", 
          row.names = FALSE, quote = FALSE)

## generate nSeqs table for supplementary material:
write.csv(miss_df, file="~/Repos/tidybug/data/text_tables/database_evals/uniqSeqs.csv",
          row.names = FALSE, quote = FALSE)

## reshape data for plot
alt_miss_df <- miss_df %>% 
  gather(key = Type, value = nTaxa, Present) %>% 
  mutate(LabelVals = comma(nTaxa, digits = 0))

## set levels for plot
alt_miss_df$Level <- factor(alt_miss_df$Level, levels = c("Species", "Genus", "Family", "Order", "Class"))
alt_miss_df$Database <- factor(alt_miss_df$Database, levels = c("tidybug", "Palmer", "Porter"))

## second plot assesses number of distinct sequences (note even same species can have more than one non-identical seq in DB)
## plot; save as figure4_distinctseqs; export at 800x400
## this figure facet labels were manually edited in Illustrator for publication
ggplot(alt_miss_df, aes(x=Database, y=nTaxa, fill=Database)) +
  geom_bar(stat="identity") +
  facet_wrap(Level ~ ., nrow = 5) +
  scale_y_continuous(breaks = scales::pretty_breaks(4), 
                     limits = c(0, NA), 
                     labels = scales::number_format(accuracy = 1, big.mark = ",")) +
  scale_fill_manual(values = dpal3) +
  geom_line() +
  labs(x="", y="unique Arthropod sequences\n", title="") +
  theme_devon() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "top", legend.text = element_text(size=14),
        strip.text.x = element_text(size=13)) +
  geom_text(aes(x=Database, y=nTaxa, label=LabelVals),
            data = alt_miss_df, 
            vjust=1)


################################################################################
## upset plots to compare overlaps in shared taxa
################################################################################
## functions to gather distinct Family, Genus, or Species names
uSpecies.function <- function(data, DBname) {
  data %>% 
    filter(!is.na(species_name)) %>% 
    mutate(taxa = paste(class_name, order_name, family_name, genus_name, species_name, sep = '_')) %>% 
    distinct(taxa) %>% 
    mutate(DBname = DBname)
}

uGenus.function <- function(data, DBname) {
  data %>% 
    filter(!is.na(genus_name)) %>% 
    mutate(taxa = paste(class_name, order_name, family_name, genus_name, sep = '_')) %>% 
    distinct(taxa) %>% 
    mutate(DBname = DBname)
}

uFamily.function <- function(data, DBname) {
  data %>% 
    filter(!is.na(family_name)) %>% 
    mutate(taxa = paste(class_name, order_name, family_name, sep = '_')) %>% 
    distinct(taxa) %>% 
    mutate(DBname = DBname)
}

uOrder.function <- function(data, DBname) {
  data %>% 
    filter(!is.na(order_name)) %>% 
    mutate(taxa = paste(class_name, order_name, sep = '_')) %>% 
    distinct(taxa) %>% 
    mutate(DBname = DBname)
}

uClass.function <- function(data, DBname) {
  data %>% 
    filter(!is.na(class_name)) %>% 
    mutate(taxa = paste(class_name)) %>% 
    distinct(taxa) %>% 
    mutate(DBname = DBname)
}

# generate unique taxa strings for Palmer data
palmer_uSpecies <- uSpecies.function(palmer_taxa, 'Palmer')
palmer_uGenus <- uGenus.function(palmer_taxa, 'Palmer')
palmer_uFamily <- uFamily.function(palmer_taxa, 'Palmer')
palmer_uOrder <- uOrder.function(palmer_taxa, 'Palmer')
palmer_uClass <- uClass.function(palmer_taxa, 'Palmer')
rm(palmer_taxa)

# generate unique taxa strings for Porter data
porter_uSpecies <- uSpecies.function(porter_taxa, 'Porter')
porter_uGenus <- uGenus.function(porter_taxa, 'Porter')
porter_uFamily <- uFamily.function(porter_taxa, 'Porter')
porter_uOrder <- uOrder.function(porter_taxa, 'Porter')
porter_uClass <- uClass.function(porter_taxa, 'Porter')
rm(porter_taxa)

## generate unique taxa strings for Tidybug data
tidybug_uSpecies <- uSpecies.function(tidybug_taxa, 'tidybug')
tidybug_uGenus <- uGenus.function(tidybug_taxa, 'tidybug')
tidybug_uFamily <- uFamily.function(tidybug_taxa, 'tidybug')
tidybug_uOrder <- uOrder.function(tidybug_taxa, 'tidybug')
tidybug_uClass <- uClass.function(tidybug_taxa, 'tidybug')
rm(tidybug_taxa)


## combine datasets per taxa Rank and convert to Upset format with binary presence/absence
## the resulting individual files were imported into Illustrator only to group into a grid; 
## no modifications to the images themselves were made otherwise

#########
#species-data combined here:
#########

tmp1 <- merge(porter_uSpecies, palmer_uSpecies, by = 'taxa', all = TRUE)
upset_species <- merge(tmp1, tidybug_uSpecies, by = 'taxa', all=TRUE)
rm(tmp1)
colnames(upset_species)[2:4] <- c("Porter", "Palmer", "tidybug")
## data tidying
upset_species$tidybug <- ifelse(grepl("tidybug", upset_species$tidybug), gsub("tidybug", 1, upset_species$tidybug), upset_species$tidybug)
upset_species$Porter <- ifelse(grepl("Porter", upset_species$Porter), gsub("Porter", 1, upset_species$Porter), upset_species$Porter)
upset_species$Palmer <- ifelse(grepl("Palmer", upset_species$Palmer), gsub("Palmer", 1, upset_species$Palmer), upset_species$Palmer)
upset_species[is.na(upset_species)] <- 0
upset_species$tidybug <- as.numeric(upset_species$tidybug)
upset_species$Porter <- as.numeric(upset_species$Porter)
upset_species$Palmer <- as.numeric(upset_species$Palmer)

## upset plot for taxa with Species-rank information
## export figure as 'figure4_upset_species'; export at 500x490
upset(upset_species, sets = c("tidybug", "Palmer", "Porter"), 
      sets.bar.color = dpal3,
      keep.order = TRUE,
      order.by = "freq", empty.intersections = "on",
      mainbar.y.label = "Unique Species names in group\n\n",
      sets.x.label = "Unique species names", scale.sets = "identity",
      show.numbers = "yes",
      text.scale = c(1.4, 1.4, 1, 1, 1, 1.4))  # follows this order: (intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)


#########
#Genus-data combined here:
#########

tmp1 <- merge(porter_uGenus, palmer_uGenus, by = 'taxa', all = TRUE)
upset_Genus <- merge(tmp1, tidybug_uGenus, by = 'taxa', all=TRUE)
rm(tmp1)
colnames(upset_Genus)[2:4] <- c("Porter", "Palmer", "tidybug")
## data tidying
upset_Genus$tidybug <- ifelse(grepl("tidybug", upset_Genus$tidybug), gsub("tidybug", 1, upset_Genus$tidybug), upset_Genus$tidybug)
upset_Genus$Porter <- ifelse(grepl("Porter", upset_Genus$Porter), gsub("Porter", 1, upset_Genus$Porter), upset_Genus$Porter)
upset_Genus$Palmer <- ifelse(grepl("Palmer", upset_Genus$Palmer), gsub("Palmer", 1, upset_Genus$Palmer), upset_Genus$Palmer)
upset_Genus[is.na(upset_Genus)] <- 0
upset_Genus$tidybug <- as.numeric(upset_Genus$tidybug)
upset_Genus$Porter <- as.numeric(upset_Genus$Porter)
upset_Genus$Palmer <- as.numeric(upset_Genus$Palmer)

## upset plot for taxa with at least Genus-rank information
## export figure as 'figure4_upset_genus'; export at 500x490
upset(upset_Genus, sets = c("tidybug", "Palmer", "Porter"), 
      sets.bar.color = dpal3,
      keep.order = TRUE,
      order.by = "freq", empty.intersections = "on",
      mainbar.y.label = "Unique Genus names in group\n\n",
      sets.x.label = "Unique genus names", scale.sets = "identity",
      show.numbers = "yes",
      text.scale = c(1.4, 1.4, 1, 1, 1, 1.4))  # follows this order: (intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)


#########
#Family-data combined here:
#########

tmp1 <- merge(porter_uFamily, palmer_uFamily, by = 'taxa', all = TRUE)
upset_Family <- merge(tmp1, tidybug_uFamily, by = 'taxa', all=TRUE)
rm(tmp1)
colnames(upset_Family)[2:4] <- c("Porter", "Palmer", "tidybug")
## data tidying
upset_Family$tidybug <- ifelse(grepl("tidybug", upset_Family$tidybug), gsub("tidybug", 1, upset_Family$tidybug), upset_Family$tidybug)
upset_Family$Porter <- ifelse(grepl("Porter", upset_Family$Porter), gsub("Porter", 1, upset_Family$Porter), upset_Family$Porter)
upset_Family$Palmer <- ifelse(grepl("Palmer", upset_Family$Palmer), gsub("Palmer", 1, upset_Family$Palmer), upset_Family$Palmer)
upset_Family[is.na(upset_Family)] <- 0
upset_Family$tidybug <- as.numeric(upset_Family$tidybug)
upset_Family$Porter <- as.numeric(upset_Family$Porter)
upset_Family$Palmer <- as.numeric(upset_Family$Palmer)

## upset plot for taxa with at least Family-rank information
## export figure as 'figure4_upset_family'; export at 500x490
upset(upset_Family, sets = c("tidybug", "Palmer", "Porter"), 
      sets.bar.color = dpal3,
      keep.order = TRUE,
      order.by = "freq", empty.intersections = "on",
      mainbar.y.label = "Unique Family names in group\n\n",
      sets.x.label = "Unique Family names", scale.sets = "identity",
      show.numbers = "yes",
      text.scale = c(1.4, 1.4, 1, 1, 1, 1.4))  # follows this order: (intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)


#########
#Order-data combined here:
#########

tmp1 <- merge(porter_uOrder, palmer_uOrder, by = 'taxa', all = TRUE)
upset_Order <- merge(tmp1, tidybug_uOrder, by = 'taxa', all=TRUE)
rm(tmp1)
colnames(upset_Order)[2:4] <- c("Porter", "Palmer", "tidybug")
## data tidying
upset_Order$tidybug <- ifelse(grepl("tidybug", upset_Order$tidybug), gsub("tidybug", 1, upset_Order$tidybug), upset_Order$tidybug)
upset_Order$Porter <- ifelse(grepl("Porter", upset_Order$Porter), gsub("Porter", 1, upset_Order$Porter), upset_Order$Porter)
upset_Order$Palmer <- ifelse(grepl("Palmer", upset_Order$Palmer), gsub("Palmer", 1, upset_Order$Palmer), upset_Order$Palmer)
upset_Order[is.na(upset_Order)] <- 0
upset_Order$tidybug <- as.numeric(upset_Order$tidybug)
upset_Order$Porter <- as.numeric(upset_Order$Porter)
upset_Order$Palmer <- as.numeric(upset_Order$Palmer)

## upset plot for taxa with at least Order-rank information
## export figure as 'figure4_upset_order'; export at 500x490
upset(upset_Order, sets = c("tidybug", "Palmer", "Porter"), 
      sets.bar.color = dpal3,
      keep.order = TRUE,
      order.by = "freq", empty.intersections = "on",
      mainbar.y.label = "Unique Order names in group\n\n",
      sets.x.label = "Unique Order names", scale.sets = "identity",
      show.numbers = "yes",
      nintersects = 6, ## invoke this value to remove the empty intersection otherwise displayed
      text.scale = c(1.4, 1.4, 1, 1, 1, 1.4))



#########
#Class-data combined here:
#########

tmp1 <- merge(porter_uClass, palmer_uClass, by = 'taxa', all = TRUE)
upset_Class <- merge(tmp1, tidybug_uClass, by = 'taxa', all=TRUE)
rm(tmp1)
colnames(upset_Class)[2:4] <- c("Porter", "Palmer", "tidybug")
## data tidying
upset_Class$tidybug <- ifelse(grepl("tidybug", upset_Class$tidybug), gsub("tidybug", 1, upset_Class$tidybug), upset_Class$tidybug)
upset_Class$Porter <- ifelse(grepl("Porter", upset_Class$Porter), gsub("Porter", 1, upset_Class$Porter), upset_Class$Porter)
upset_Class$Palmer <- ifelse(grepl("Palmer", upset_Class$Palmer), gsub("Palmer", 1, upset_Class$Palmer), upset_Class$Palmer)
upset_Class[is.na(upset_Class)] <- 0
upset_Class$tidybug <- as.numeric(upset_Class$tidybug)
upset_Class$Porter <- as.numeric(upset_Class$Porter)
upset_Class$Palmer <- as.numeric(upset_Class$Palmer)

## upset plot for taxa with at least Class-rank information
## export figure as 'figure4_upset_class'; export at 500x490
upset(upset_Class, sets = c("tidybug", "Palmer", "Porter"), 
      sets.bar.color = dpal3,
      keep.order = TRUE,
      order.by = "freq", empty.intersections = "on",
      mainbar.y.label = "Unique Class names in group\n\n",
      sets.x.label = "Unique Class names", scale.sets = "identity",
      show.numbers = "yes",
      nintersects = 5,
      text.scale = c(1.4, 1.4, 1, 1, 1, 1.4))


################################################################################
## unused code
################################################################################

## grab screenshots to convert into format needed for ggarrange master graphic
upset(upset_species, sets = c("tidybug", "Palmer", "Porter"), sets.bar.color = dpal3, keep.order = TRUE, order.by = "freq", 
      empty.intersections = "on", mainbar.y.label = "Unique Species names in group\n\n", sets.x.label = "Unique species names", scale.sets = "identity", 
      show.numbers = "yes", text.scale = c(1.4, 1.4, 1, 1, 1, 1.4))
grid.edit('arrange',name='arrange2')
u1 = grid.grab()

upset(upset_Genus, sets = c("tidybug", "Palmer", "Porter"), sets.bar.color = dpal3, keep.order = TRUE, order.by = "freq", empty.intersections = "on",
      mainbar.y.label = "Unique Genus names in group\n\n", sets.x.label = "Unique genus names", scale.sets = "identity",
      show.numbers = "yes", text.scale = c(1.4, 1.4, 1, 1, 1, 1.4))  
grid.edit('arrange',name='arrange2')
u2 = grid.grab()

upset(upset_Family, sets = c("tidybug", "Palmer", "Porter"), sets.bar.color = dpal3, keep.order = TRUE, order.by = "freq", 
      empty.intersections = "on", mainbar.y.label = "Unique Family names in group\n\n", sets.x.label = "Unique Family names", scale.sets = "identity",
      show.numbers = "yes", text.scale = c(1.4, 1.4, 1, 1, 1, 1.4))
grid.edit('arrange',name='arrange2')
u3 = grid.grab()

upset(upset_Order, sets = c("tidybug", "Palmer", "Porter"), sets.bar.color = dpal3, keep.order = TRUE, order.by = "freq", empty.intersections = "on",
      mainbar.y.label = "Unique Order names in group\n\n", sets.x.label = "Unique Order names", scale.sets = "identity",
      show.numbers = "yes", nintersects = 6, text.scale = c(1.4, 1.4, 1, 1, 1, 1.4))
grid.edit('arrange',name='arrange2')
u4 = grid.grab()

upset(upset_Class, sets = c("tidybug", "Palmer", "Porter"),  sets.bar.color = dpal3, keep.order = TRUE, order.by = "freq", empty.intersections = "on", 
      mainbar.y.label = "Unique Class names in group\n\n", sets.x.label = "Unique Class names", scale.sets = "identity",
      show.numbers = "yes", nintersects = 5, text.scale = c(1.4, 1.4, 1, 1, 1, 1.4))
grid.edit('arrange',name='arrange2')
u5 = grid.grab()

ggarrange(u1,u2,u3,u4,u5, nrow = 5)
p1
png("~/Desktop/fig4_plot.png", width = 900, height = 1200)
ggarrange(p1, p2, p3, ncol=3, labels="A", "B", "C")
dev.off()
