library(tidyverse)
library(scales)
library(ggpubr)

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
#notrun: palmer_taxa <- fread(file = "~/Repos/tidybug/data/databases/palmer.taxa.txt.gz", fill = TRUE)
palmer_taxa <- read_csv(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/palmer.taxa.txt.gz", col_names = FALSE)
colnames(palmer_taxa) <- c("class_name", "order_name", "family_name", "genus_name", "species_name")
palmer_taxa <- as.data.frame(apply(palmer_taxa, 2, function(y) (gsub(".:", "", y))))  ## remove prefixes for each taxa (ex. 'p:', or 'c:', etc.)

## apply three function
palmer_uniq <- uniqfunction(palmer_taxa, "Palmer")
palmer_abund <- uniqCountfunction(palmer_taxa, "Palmer")
palmer_miss <- missfunction(palmer_taxa, "Palmer")
rm(palmer_taxa)

## ----  porter data:  ----  ##
## porter data:
porter_taxa <- read_csv(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/porter.derep.taxa.txt.gz", col_names = FALSE)
colnames(porter_taxa) <- c("class_name", "order_name", "family_name", "genus_name", "species_name")
porter_taxa <- as.data.frame(apply(porter_taxa, 2, function(y) (gsub("_", " ", y))))    ## replace underscore in species names with spaces to match other database naming conventions
porter_uniq <- uniqfunction(porter_taxa, "Porter")
porter_abund <- uniqCountfunction(porter_taxa, "Porter")
porter_miss <- missfunction(porter_taxa, "Porter")
rm(porter_taxa)

## ----  derep data:  ----  ##
derep_taxa <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/derep.taxa.txt.gz", delim = ";", col_names = FALSE)
colnames(derep_taxa) <- c("class_name", "order_name", "family_name", "genus_name", "species_name")
derep_taxa <- as.data.frame(apply(derep_taxa, 2, function(y) (gsub("Ambiguous_taxa", NA, y))))    ## substitute the Ambiguous Taxa designation with NA (it's unknown!)
derep_uniq <- uniqfunction(derep_taxa, "derep")
derep_abund <- uniqCountfunction(derep_taxa, "derep")
derep_miss <- missfunction(derep_taxa, "derep")
rm(derep_taxa)

## ----  raw data:  ----  ##
raw_taxa <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/raw.taxa.txt.gz", delim = ";", col_names = TRUE)
raw_uniq <- uniqfunction(raw_taxa, "raw")
#notrun: raw_abund <- uniqCountfunction(raw_taxa, "raw")
raw_miss <- missfunction(raw_taxa, "raw")
rm(raw_taxa)

## combine the three sets of data
uniq_df <- rbind(palmer_uniq, porter_uniq, derep_uniq, raw_uniq)
rm(palmer_uniq, porter_uniq, derep_uniq, raw_uniq)
abund_df <- rbind(palmer_abund, porter_abund, derep_abund)
rm(palmer_abund, porter_abund, derep_abund)
miss_df <- rbind(palmer_miss, porter_miss, derep_miss, raw_miss)
rm(palmer_miss, porter_miss, derep_miss, raw_miss)


## color palette for plot
vpal4 <- viridis(4, option = "plasma", end=0.85)
pal4 <- c("#8A09A5FF", "#0D0887FF", "#DA5B6AFF", "#FEBA2CFF")

   
## plot; save as db_2_uniqueness; export at 700x525
uniq_df$Database <- factor(uniq_df$Database, levels = c("raw", "derep", "Palmer", "Porter"))
uniq_df$Level <- factor(uniq_df$Level, levels = c("Class", "Order", "Family", "Genus", "Species"))
ggplot(uniq_df, aes(x=Level, y=nUniq, group=Database, color=Database, shape=Database)) +
  geom_point() +
  scale_shape_manual(values = c(0,1,2,3)) +
  scale_y_continuous(labels = comma) +
  scale_color_manual(values = pal4) +
  geom_line() +
  labs(x="", y="Distinct taxa", title="") +
  theme_devon()

## plot; save as db_3_missingness; export at 1000x600
alt_miss_df <- miss_df %>% gather(key = Type, value = nTaxa, Present, Missing)
alt_miss_df$Level <- factor(alt_miss_df$Level, levels = c("Class", "Order", "Family", "Genus", "Species"))
alt_miss_df$Database <- factor(alt_miss_df$Database, levels = c("raw", "derep", "Palmer", "Porter"))
alt_miss_df$Type <- factor(alt_miss_df$Type, levels = c("Missing", "Present"))
ggplot(alt_miss_df, aes(y=nTaxa, fill=Type, x=Level)) +
  geom_bar(stat="identity", color="gray30") +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = c("gray85", "gray45")) +
  facet_wrap(Database ~ ., nrow = 4) +
  labs(x="", y="distinct taxa", title="", fill="Taxa information") +
  theme_devon() +
  coord_flip() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle=22.5, hjust=1))

## plot; save as db_3b_completeness; export at 600x475
miss_df$Database <- factor(miss_df$Database, levels = c("raw", "derep", "Palmer", "Porter"))
miss_df$Level <- factor(miss_df$Level, levels = c("Class", "Order", "Family", "Genus", "Species"))
ggplot(miss_df, aes(x=Level, y=pPresent, group=Database, color=Database, shape=Database)) +
  geom_point(size = 3, stroke=1.2) +
  scale_shape_manual(values = c(0,1,2,3)) +
  scale_color_manual(values = pal4) +
  geom_line() +
  labs(x="", y="Fraction taxa with information", title="") +
  theme_devon()



## add in abundance freuquencies
alt_abund_df <- abund_df %>% group_by(Level, Database) %>% mutate(nObs=sum(counts)) %>% mutate(pCounts=counts/nObs)
## export entire data.frame as single .csv file:
write.csv(alt_abund_df, file="~/Repos/tidybug/data/text_tables/db_abundances.csv", row.names = FALSE)

## plots separated by Levels, then grouped after creating unique x-axes; 
class_abund <- alt_abund_df %>% filter(Level == "Class") %>% group_by(Database) %>% mutate(cumFrac=cumsum(pCounts))
order_abund <- alt_abund_df %>% filter(Level == "Order") %>% group_by(Database) %>% mutate(cumFrac=cumsum(pCounts))
family_abund <- alt_abund_df %>% filter(Level == "Family") %>% group_by(Database) %>% mutate(cumFrac=cumsum(pCounts))
genus_abund <- alt_abund_df %>% filter(Level == "Genus") %>% group_by(Database) %>% mutate(cumFrac=cumsum(pCounts))

## select just the TaxaNames with common Levels from top hits (not the "others")
hitsfunction <- function(data, filter1, filter2){
  filter1_enq <- enquo(filter1)
  filter2_enq <- enquo(filter2)
  tmp.df <- data %>% filter(!!filter1_enq) %>% filter(!!filter2_enq) %>% arrange(-counts) %>% group_by(Level) %>% select(TaxaName)
  as.character(tmp.df$TaxaName)
}

derep.ordhits <- hitsfunction(order_abund, Database=="derep", TaxaName!="others")
palmer.ordhits <- hitsfunction(order_abund, Database=="Palmer", TaxaName!="others")
porter.ordhits <- hitsfunction(order_abund, Database=="Porter", TaxaName!="others")
ordhits <- intersect(intersect(derep.ordhits, palmer.ordhits),porter.ordhits)

derep.famhits <- hitsfunction(family_abund, Database=="derep", TaxaName!="others")
palmer.famhits <- hitsfunction(family_abund, Database=="Palmer", TaxaName!="others")
porter.famhits <- hitsfunction(family_abund, Database=="Porter", TaxaName!="others")
famhits <- intersect(intersect(derep.famhits, palmer.famhits),porter.famhits)

derep.genhits <- hitsfunction(genus_abund, Database=="derep", TaxaName!="others")
palmer.genhits <- hitsfunction(genus_abund, Database=="Palmer", TaxaName!="others")
porter.genhits <- hitsfunction(genus_abund, Database=="Porter", TaxaName!="others")
genhits <- intersect(intersect(derep.genhits, palmer.genhits),porter.genhits)

## 3 color palette from 4pal
pal3 <- c("#0D0887FF", "#DA5B6AFF", "#FEBA2CFF")

## plots - creating a 2x3 plot with total counts (left column) and percentages (right column) for Order through Genus Levels (rows)
## sort the plot order by descending $count for the "derep" Database
order.order <- order_abund %>% filter(Database=="derep") %>% filter(TaxaName %in% ordhits) %>% arrange(-counts) %>% group_by(Level) %>% select(TaxaName)
order.order <- as.character(order.order$TaxaName)
order_abund$TaxaName <- factor(order_abund$TaxaName, levels = order.order)

family.order <- family_abund %>% filter(Database=="derep") %>% filter(TaxaName %in% famhits) %>% arrange(-counts) %>% group_by(Level) %>% select(TaxaName)
family.order <- as.character(family.order$TaxaName)
family_abund$TaxaName <- factor(family_abund$TaxaName, levels = family.order)

genus.order <- genus_abund %>% filter(Database=="derep") %>% filter(TaxaName %in% genhits) %>% arrange(-counts) %>% group_by(Level) %>% select(TaxaName)
genus.order <- as.character(genus.order$TaxaName)
genus_abund$TaxaName <- factor(genus_abund$TaxaName, levels = genus.order)


pop <- ggplot(data=order_abund %>% filter(TaxaName %in% ordhits), aes(fill=Database, y=pCounts, x=TaxaName)) +
  geom_bar(stat = "identity", color="gray30", position="dodge") +
  scale_fill_manual(values=pal3) + 
  scale_y_continuous(limits = c(0,0.45), breaks = c(0, 0.2, 0.4)) +
  labs(x="", y="", subtitle = "common abundant Orders") +
  theme_devon() +
  theme(axis.text.x = element_text(angle=22.5, hjust = 1),
        axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
        plot.subtitle = element_text(size = 10),
        legend.position = "top")

poc <- ggplot(data=order_abund %>% filter(TaxaName %in% ordhits), aes(fill=Database, y=counts, x=TaxaName)) +
  geom_bar(stat = "identity", color="gray30", position="dodge") +
  scale_fill_manual(values=pal3) + 
  scale_y_continuous(labels = comma, breaks = c(0, 300000, 600000)) +
  labs(x="", y="", subtitle = "common abundant Orders") +
  theme_devon() +
  theme(axis.text.x = element_text(angle=22.5, hjust = 1),
        axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
        plot.subtitle = element_text(size = 10),
        legend.position = "top")

pfp <- ggplot(data=family_abund %>% filter(TaxaName %in% famhits), aes(fill=Database, y=pCounts, x=TaxaName)) +
  geom_bar(stat = "identity", color="gray30", position="dodge") +
  scale_fill_manual(values=pal3) + 
  labs(x="", y="fraction distinct sequences", subtitle = "common abundant Families") +
  scale_y_continuous(limits = c(0,0.1), breaks = c(0, 0.05, 0.1)) +
  theme_devon() +
  theme(axis.text.x = element_text(angle=22.5, hjust = 1),
        axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
        plot.subtitle = element_text(size = 10),
        legend.position = "top")

pfc <- ggplot(data=family_abund %>% filter(TaxaName %in% famhits), aes(fill=Database, y=counts, x=TaxaName)) +
  geom_bar(stat = "identity", color="gray30", position="dodge") +
  scale_fill_manual(values=pal3) + 
  scale_y_continuous(labels = comma, breaks = c(0, 75000, 150000)) +
  labs(x="", y="distinct sequences", subtitle = "common abundant Families") +
  theme_devon() +
  theme(axis.text.x = element_text(angle=22.5, hjust = 1),
        axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
        plot.subtitle = element_text(size = 10),
        legend.position = "top")

pgp <- ggplot(data=genus_abund %>% filter(TaxaName %in% genhits), aes(fill=Database, y=pCounts, x=TaxaName)) +
  geom_bar(stat = "identity", color="gray30", position="dodge") +
  scale_fill_manual(values=pal3) +
  labs(y="", x="", subtitle = "common abundant Genera") +
  scale_y_continuous(limits = c(0,0.04), breaks = c(0, 0.02, 0.04)) +
  theme_devon() +
  theme(axis.text.x = element_text(angle=22.5, hjust = 1),
        axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
        plot.subtitle = element_text(size = 10),
        legend.position = "top")

pgc <- ggplot(data=genus_abund %>% filter(TaxaName %in% genhits), aes(fill=Database, y=counts, x=TaxaName)) +
  geom_bar(stat = "identity", color="gray30", position="dodge") +
  scale_fill_manual(values=pal3) +
  scale_y_continuous(labels = comma, breaks = c(0, 20000, 40000)) +
  labs(y="", x="", subtitle = "common abundant Genera") +
  theme_devon() +
  theme(axis.text.x = element_text(angle=22.5, hjust = 1),
        axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
        plot.subtitle = element_text(size = 10),
        legend.position = "top")


## plot all 6 graphics together; save as db_4_abundantTaxa; export at 1000 x 800
ggarrange(pop, poc, pfp, pfc, pgp, pgc,
          common.legend = TRUE, legend = "top",
          ncol=2, nrow = 3, 
          heights = c(1.2, 1, 1))
