library(tidyverse)
library(scales)
library(viridis)
library(ggpubr)

## theme function:
theme_devon <- function () { 
  theme_bw(base_size=12, base_family="Avenir") %+replace% 
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill="transparent", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA)
    )
}


## import data:
raw_meta <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/reduced.allArth.meta.txt.gz", delim = ";", col_names = TRUE)
derep_taxa <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/boldCOI.derep.txt.gz", delim = ";", col_names = FALSE)
colnames(derep_taxa) <- c("SeqID", "class_name", "order_name", "family_name", "genus_name", "species_name")
#notrun: we aren't going to substitute ambiguous taxa as NA because we're quantifying how often they ambiguous taxa are present intitially!
#notrun: derep_taxa <- as.data.frame(apply(derep_taxa, 2, function(y) (gsub("Ambiguous_taxa", NA, y))))
derep_IDs <- derep_taxa %>% select(SeqID)
rm(derep_taxa)


## calculating number of records by institution
countrylist <- c("United States", "Canada")
## raw data, filtered by country

raw_unfilt_df <- raw_meta %>% group_by(institution_storing) %>% 
  summarise(nRecords = n()) %>% mutate(pRecords = nRecords / sum(nRecords)) %>% arrange(., -nRecords) %>% 
  mutate(RecordType="raw") %>%
  mutate(FilterType="Any location")

raw_filtd_df <- raw_meta %>% filter(country %in% countrylist) %>% group_by(institution_storing) %>% 
  summarise(nRecords = n()) %>% mutate(pRecords = nRecords / sum(nRecords)) %>% arrange(., -nRecords) %>% 
  mutate(RecordType="raw") %>%
  mutate(FilterType="US|Canada")

derep_unfilt_df <- raw_meta %>% filter(sequenceID %in% derep_IDs$SeqID) %>% group_by(institution_storing) %>% 
  summarise(nRecords = n()) %>%  mutate(pRecords = nRecords / sum(nRecords)) %>% arrange(., -nRecords) %>% 
  mutate(RecordType="derep") %>%
  mutate(FilterType="Any location")

derep_filtd_df <- raw_meta %>% filter(sequenceID %in% derep_IDs$SeqID) %>% filter(country %in% countrylist) %>% group_by(institution_storing) %>% 
  summarise(nRecords = n()) %>% mutate(pRecords = nRecords / sum(nRecords)) %>% arrange(., -nRecords) %>% 
  mutate(RecordType="derep") %>%
  mutate(FilterType="US|Canada")

rm(derep_IDs, raw_meta)

# combine datasets
df <- rbind(raw_unfilt_df, raw_filtd_df, derep_unfilt_df, derep_filtd_df)
rm(raw_unfilt_df, raw_filtd_df, derep_unfilt_df, derep_filtd_df)
## export dataset:
write.csv(df, file="~/Repos/tidybug/data/text_tables/database_geography.csv", row.names = FALSE, quote = TRUE)

# pull out just the top 2 datasets of interest, then group all others into a single "other" bin
selectInstitution <- c("Centre for Biodiversity Genomics", "Mined from GenBank, NCBI")
tmp1 <- df %>% group_by(FilterType, RecordType) %>% 
  filter(institution_storing %in% selectInstitution) %>% 
  select(FilterType, RecordType, institution_storing, pRecords)
tmp2 <- df %>% group_by(FilterType, RecordType) %>% 
  filter(!institution_storing %in% selectInstitution)%>% 
  select(FilterType, RecordType, institution_storing, pRecords) %>% 
  group_by(FilterType, RecordType) %>% 
  summarise(pRecords=sum(pRecords)) %>% 
  mutate(institution_storing="other")
df_plot <- rbind(tmp1, tmp2)

# set levels
df_plot$RecordType <- factor(df_plot$RecordType, levels = c("raw", "derep"))

viridis(3, option="cividis", begin=0.2)
pal3 <- c("#31446BFF","#FFEA46FF", "#958F78FF")

# plot; save as db_6_geography; export at 650x557
ggplot(data=df_plot, aes(x=FilterType, y=pRecords, label=institution_storing, fill=institution_storing)) +
  geom_bar(stat="identity") +
  facet_grid(RecordType ~ .) +
  scale_fill_manual(values=pal3) + 
  labs(x="", y="fraction of sequences", fill="Institution") +
  theme_devon() +
  theme(legend.position = "top", 
        strip.text = element_text(size=12))


## to calculate the proportion of taxa represented we'll download the data again but filter it differently:
## we're going to explore just the dereplicated data...
## import data:
countrylist <- c("United States", "Canada")
raw_meta <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/reduced.allArth.meta.txt.gz", delim = ";", col_names = TRUE)
raw_NORTH <- raw_meta %>% select(sequenceID, country) %>% filter(country %in% countrylist)
derep_taxa <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/boldCOI.derep.txt.gz", delim = ";", col_names = FALSE)
colnames(derep_taxa) <- c("SeqID", "class_name", "order_name", "family_name", "genus_name", "species_name")

## create filtered dataset by country requirement 
derep_taxa_wLocation <- merge(derep_taxa, raw_NORTH, by.x='SeqID', by.y='sequenceID')

## substitute ambiguous taxa as NA this time!
derep_taxa_wLocation <- as.data.frame(apply(derep_taxa_wLocation, 2, function(y) (gsub("Ambiguous_taxa", NA, y))))
derep_taxa <- as.data.frame(apply(derep_taxa, 2, function(y) (gsub("Ambiguous_taxa", NA, y))))
derep_taxa$country <- "Any location"

## composition function
compfunction <- function(data,dataset) {
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

## plot; save as db_7a_MissingnessByGeography; export at 450x400
## we can see that geography-based filtering does NOT change the overall fraction of information...
derep_allSites <- compfunction(derep_taxa, "Any location")
derep_NAsites <- compfunction(derep_taxa_wLocation, "US | Canada")
df1 <- rbind(derep_allSites, derep_NAsites)
rm(derep_allSites, derep_NAsites, raw_NORTH, raw_meta)
df1$Level <- factor(df1$Level, levels = c("Class", "Order", "Family", "Genus", "Species"))
ggplot(df1, aes(x=Level, y=pPresent, fill=Database)) + 
  geom_bar(stat="identity", position = "dodge", color="black", size=0.25) +
  scale_fill_manual(values = c("gray40", "gray80")) +
  labs(x="", y="fraction records with information", fill="Country of record") +
  theme_devon() +
  theme(legend.position = "top")

## can also present this by showing missing vs. present info:
df1alt <- df1 %>% gather(key=InfoType, value=Value, Present, Missing)
## plot; save as db_7a-alt_MissingnessByGeography; export at 500x500
ggplot(df1alt, aes(x=Level, y=Value, fill=InfoType)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("gray40", "gray80")) +
  scale_y_continuous(labels=comma) + 
  facet_grid(Database ~ .) +
  labs(x="", y="number records", fill="Information content") +
  theme_devon() +
  theme(legend.position = "top") +
  coord_flip()


## calculating most abundant taxa:
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
  tmp <- rbind(class_uniq, order_uniq, family_uniq, genus_uniq, species_uniq)
  tmp %>% group_by(Level) %>% mutate(pCounts=(counts/sum(counts)))
}

## generate the two datasets and merge
uniq_allSites <- uniqCountfunction(derep_taxa, "Any location")
uniq_allSites <- as.data.frame(uniq_allSites)
uniq_NAsites <- uniqCountfunction(derep_taxa_wLocation, "US | Canada")
uniq_NAsites <- as.data.frame(uniq_NAsites)
uniq_all <- rbind(uniq_allSites, uniq_NAsites)
rm(uniq_allSites, uniq_NAsites, derep_taxa_wLocation, derep_taxa)

## plots separated by Levels, then grouped after creating unique x-axes; analyse same 3 Levels (Order through Genus) as earlier plots from `database_composition.R`
order_abund <- uniq_all %>% filter(Level == "Order") %>% group_by(Database) %>% mutate(cumFrac=cumsum(pCounts))
family_abund <- uniq_all %>% filter(Level == "Family") %>% group_by(Database) %>% mutate(cumFrac=cumsum(pCounts))
genus_abund <- uniq_all %>% filter(Level == "Genus") %>% group_by(Database) %>% mutate(cumFrac=cumsum(pCounts))

## select just the TaxaNames with common Levels from top hits (not the "others")
hitsfunction <- function(data, filter1){
  filter1_enq <- enquo(filter1)
  tmp.df <- data %>% filter(!!filter1_enq) %>% arrange(-counts) %>% group_by(Level) %>% select(TaxaName)
  as.character(tmp.df$TaxaName)
}

allSites.ordhits <- hitsfunction(order_abund, Database=="Any location")
NAsites.ordhits <- hitsfunction(order_abund, Database=="US | Canada")
ordhits <- intersect(allSites.ordhits, NAsites.ordhits)

allSites.famhits <- hitsfunction(family_abund, Database=="Any location")
NAsites.famhits <- hitsfunction(family_abund, Database=="US | Canada")
famhits <- intersect(allSites.famhits, NAsites.famhits)

allSites.genhits <- hitsfunction(genus_abund, Database=="Any location")
NAsites.genhits <- hitsfunction(genus_abund, Database=="US | Canada")
genhits <- intersect(allSites.genhits, NAsites.genhits)


## sort the data by most abundant taxa at 3 Levels (Order, Family, Genus); descending for "Any location" $Database
order.order <- order_abund %>% filter(Database=="Any location") %>% filter(TaxaName %in% ordhits) %>% arrange(-counts) %>% group_by(Level) %>% select(TaxaName)
order.order <- as.character(order.order$TaxaName)
order_abund$TaxaName <- factor(order_abund$TaxaName, levels = order.order)

family.order <- family_abund %>% filter(Database=="Any location") %>% filter(TaxaName %in% famhits) %>% arrange(-counts) %>% group_by(Level) %>% select(TaxaName)
family.order <- as.character(family.order$TaxaName)
family_abund$TaxaName <- factor(family_abund$TaxaName, levels = family.order)

genus.order <- genus_abund %>% filter(Database=="Any location") %>% filter(TaxaName %in% genhits) %>% arrange(-counts) %>% group_by(Level) %>% select(TaxaName)
genus.order <- as.character(genus.order$TaxaName)
genus_abund$TaxaName <- factor(genus_abund$TaxaName, levels = genus.order)


pal2 <- c("#0D0887FF", "light blue")    ## 2 color palette

## plots stitched together like before in a 3x2 grid for % composition and total number of reads
pop <- ggplot(data=order_abund %>% filter(TaxaName %in% ordhits) %>% filter(TaxaName != "others"), aes(fill=Database, y=pCounts, x=TaxaName)) +
  geom_bar(stat = "identity", color="gray30", position="dodge") +
  scale_fill_manual(values=pal2) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4)) +
  labs(x="", y="", subtitle = "common abundant Orders", fill="") +
  theme_devon() +
  theme(axis.text.x = element_text(angle=22.5, hjust = 1),
        axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
        plot.subtitle = element_text(size = 10),
        legend.position = "top")

poc <- ggplot(data=order_abund %>% filter(TaxaName %in% ordhits) %>% filter(TaxaName != "others"), aes(fill=Database, y=counts, x=TaxaName)) +
  geom_bar(stat = "identity", color="gray30", position="dodge") +
  scale_fill_manual(values=pal2) + 
  scale_y_continuous(labels = comma, breaks = c(0, 300000, 600000)) +
  labs(x="", y="", subtitle = "common abundant Orders", fill="") +
  theme_devon() +
  theme(axis.text.x = element_text(angle=22.5, hjust = 1),
        axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
        plot.subtitle = element_text(size = 10),
        legend.position = "top")

pfp <- ggplot(data=family_abund %>% filter(TaxaName %in% famhits) %>% filter(TaxaName != "others"), aes(fill=Database, y=pCounts, x=TaxaName)) +
  geom_bar(stat = "identity", color="gray30", position="dodge") +
  scale_fill_manual(values=pal2) + 
  labs(x="", y="fraction distinct sequences", subtitle = "common abundant Families", fill = "") +
  scale_y_continuous(limits = c(0,0.15), breaks = c(0, 0.05, 0.1)) +
  theme_devon() +
  theme(axis.text.x = element_text(angle=22.5, hjust = 1),
        axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
        plot.subtitle = element_text(size = 10),
        legend.position = "top")

pfc <- ggplot(data=family_abund %>% filter(TaxaName %in% famhits) %>% filter(TaxaName != "others"), aes(fill=Database, y=counts, x=TaxaName)) +
  geom_bar(stat = "identity", color="gray30", position="dodge") +
  scale_fill_manual(values=pal2) + 
  scale_y_continuous(labels = comma, breaks = c(0, 75000, 150000)) +
  labs(x="", y="distinct sequences", subtitle = "common abundant Families", fill="") +
  theme_devon() +
  theme(axis.text.x = element_text(angle=22.5, hjust = 1),
        axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
        plot.subtitle = element_text(size = 10),
        legend.position = "top")

pgp <- ggplot(data=genus_abund %>% filter(TaxaName %in% genhits) %>% filter(TaxaName != "others"), aes(fill=Database, y=pCounts, x=TaxaName)) +
  geom_bar(stat = "identity", color="gray30", position="dodge") +
  scale_fill_manual(values=pal2) +
  labs(y="", x="", subtitle = "common abundant Genera", fill="") +
  #scale_y_continuous(breaks = c(0, 0.015, 0.03)) +
  theme_devon() +
  theme(axis.text.x = element_text(angle=22.5, hjust = 1),
        axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
        plot.subtitle = element_text(size = 10),
        legend.position = "top")

pgc <- ggplot(data=genus_abund %>% filter(TaxaName %in% genhits) %>% filter(TaxaName != "others"), aes(fill=Database, y=counts, x=TaxaName)) +
  geom_bar(stat = "identity", color="gray30", position="dodge") +
  scale_fill_manual(values=pal2) +
  #scale_y_continuous(labels = comma, breaks = c(0, 20000, 40000)) +
  labs(y="", x="", subtitle = "common abundant Genera", fill="") +
  theme_devon() +
  theme(axis.text.x = element_text(angle=22.5, hjust = 1),
        axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
        plot.subtitle = element_text(size = 10),
        legend.position = "top")


## plot all 6 graphics together; save as db_7b_CompositionByGeography; export at 1200 x 900
ggarrange(pop, poc, pfp, pfc, pgp, pgc,
          common.legend = TRUE, legend = "top",
          ncol=2, nrow = 3, 
          heights = c(1.2, 1, 1))
