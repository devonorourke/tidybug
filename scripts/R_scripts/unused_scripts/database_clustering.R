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


## function to calculate missingness and incorporate ambiguous data information
missfunction <- function(data,dataset) {
  Class.no <- data %>% filter(!is.na(class_name)) %>% filter(class_name != "Ambiguous_taxa") %>% summarise(counts = n()) %>% mutate(Level = "Class", Value = "Present")
  Class.na <- data %>% filter(is.na(class_name)) %>% summarise(counts = n()) %>% mutate(Level = "Class", Value = "Missing")
  Class.am <- data %>% filter(class_name=="Ambiguous_taxa") %>% summarise(counts = n()) %>% mutate(Level = "Class", Value = "ambiguous")
  Class <- rbind(Class.no, Class.na, Class.am) %>% spread(key=Value, value=counts)
  Order.no <- data %>% filter(!is.na(order_name)) %>% filter(order_name != "Ambiguous_taxa") %>% summarise(counts = n()) %>% mutate(Level = "Order", Value = "Present")
  Order.na <- data %>% filter(is.na(order_name)) %>% summarise(counts = n()) %>% mutate(Level = "Order", Value = "Missing")
  Order.am <- data %>% filter(order_name=="Ambiguous_taxa") %>% summarise(counts = n()) %>% mutate(Level = "Order", Value = "ambiguous")
  Order <- rbind(Order.no, Order.na, Order.am) %>% spread(key=Value, value=counts)
  Family.no <- data %>% filter(!is.na(family_name)) %>% filter(family_name != "Ambiguous_taxa") %>% summarise(counts = n()) %>% mutate(Level = "Family", Value = "Present")
  Family.na <- data %>% filter(is.na(family_name)) %>% summarise(counts = n()) %>% mutate(Level = "Family", Value = "Missing")
  Family.am <- data %>% filter(family_name=="Ambiguous_taxa") %>% summarise(counts = n()) %>% mutate(Level = "Family", Value = "ambiguous")
  Family <- rbind(Family.no, Family.na, Family.am) %>% spread(key=Value, value=counts)
  Genus.no <- data %>% filter(!is.na(genus_name)) %>% filter(genus_name != "Ambiguous_taxa") %>% summarise(counts = n()) %>% mutate(Level = "Genus", Value = "Present")
  Genus.na <- data %>% filter(is.na(genus_name)) %>% summarise(counts = n()) %>% mutate(Level = "Genus", Value = "Missing")
  Genus.am <- data %>% filter(genus_name=="Ambiguous_taxa") %>% summarise(counts = n()) %>% mutate(Level = "Genus", Value = "ambiguous")
  Genus <- rbind(Genus.no, Genus.na, Genus.am) %>% spread(key=Value, value=counts)
  Species.no <- data %>% filter(!is.na(species_name)) %>% filter(species_name != "Ambiguous_taxa") %>% summarise(counts = n()) %>% mutate(Level = "Species", Value = "Present")
  Species.na <- data %>% filter(is.na(species_name)) %>% summarise(counts = n()) %>% mutate(Level = "Species", Value = "Missing")
  Species.am <- data %>% filter(species_name=="Ambiguous_taxa") %>% summarise(counts = n()) %>% mutate(Level = "Species", Value = "ambiguous")
  Species <- rbind(Species.no, Species.na, Species.am) %>% spread(key=Value, value=counts)
  tmp <- rbind(Class, Order, Family, Genus, Species)
  tmp <- tmp %>% mutate(Database=dataset)
  tmp$pPresent <- tmp$Present / (tmp$Present + tmp$Missing + tmp$ambiguous)
  tmp$pMissing <- tmp$Missing / (tmp$Present + tmp$Missing + tmp$ambiguous)
  tmp$pAmbiguous <- tmp$ambiguous / (tmp$Present + tmp$Missing + tmp$ambiguous)
  tmp
}


## import data:
derep_taxa <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/boldCOI.derep.txt.gz", delim = ";", col_names = FALSE)
colnames(derep_taxa) <- c("SeqID", "class_name", "order_name", "family_name", "genus_name", "species_name")

clust99_taxa <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/clust99.names.txt.gz", delim = ";", col_names = FALSE)
colnames(clust99_taxa) <- c("SeqID", "kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")
clust99_taxa <- clust99_taxa %>% select(-kingdom_name, -phylum_name)

clust97_taxa <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/clust97.names.txt.gz", delim = ";", col_names = FALSE)
colnames(clust97_taxa) <- c("SeqID", "kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")
clust97_taxa <- clust97_taxa %>% select(-kingdom_name, -phylum_name)

clust95_taxa <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/clust95.names.txt.gz", delim = ";", col_names = FALSE)
colnames(clust95_taxa) <- c("SeqID", "kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")
clust95_taxa <- clust95_taxa %>% select(-kingdom_name, -phylum_name)

## generate missingness per level
derep_miss <- missfunction(derep_taxa, "clust 100%")
clust99_miss <- missfunction(clust99_taxa, "clust 99%")
clust97_miss <- missfunction(clust97_taxa, "clust 97%")
clust95_miss <- missfunction(clust95_taxa, "clust 95%")

rm(clust99_taxa, clust97_taxa, clust95_taxa)
all_miss <- rbind(derep_miss, clust99_miss, clust97_miss, clust95_miss)
#rm(derep_miss, clust99_miss, clust97_miss, clust95_miss)

## levels for plot
all_miss$Level <- factor(all_miss$Level, levels = c("Class", "Order", "Family", "Genus", "Species"))
all_miss$Database <- factor(all_miss$Database, levels = c("clust 100%", "clust 99%", "clust 97%", "clust 95%"))

## plot - total number of records; save as db5_clustRecords; export at 800x600
alt_miss <- all_miss %>% select(Level, Database, Present, Missing, ambiguous) %>% gather(key = InfoType, value = Value, Present, Missing, ambiguous)
alt_miss$InfoType <- gsub("ambiguous", "Ambiguous", alt_miss$InfoType)
alt_miss$InfoType <- gsub("ambiguous", "Ambiguous", alt_miss$InfoType)
alt_miss$Labeler <- paste(alt_miss$Database, alt_miss$InfoType, sep=" - ")
alt_miss$Labeler <- factor(alt_miss$Labeler, levels = c(
  "clust 100% - Present", "clust 100% - Missing", "clust 100% - Ambiguous","clust 99% - Present", "clust 99% - Missing", "clust 99% - Ambiguous",
  "clust 97% - Present", "clust 97% - Missing", "clust 97% - Ambiguous","clust 95% - Present", "clust 95% - Missing", "clust 95% - Ambiguous"))

## plot; save as 
p1 <- ggplot(alt_miss, aes(x=Level, y=Value, fill=InfoType)) +  
  geom_bar(stat="identity") +
  scale_y_continuous(labels = comma, breaks = c(0, 750000, 1500000)) +
  coord_flip() +
  facet_wrap(~ Database, nrow = 4) +
  #scale_fill_manual(values=c("ambiguous", "missing", "present")) +
  scale_fill_manual(values=c("red", "gray20", "gray70"),
                    labels=c("Present", "Missing", "Ambiguous"),
                    breaks=c("Present", "Missing", "Ambiguous")) +
  labs(x="", y="number of sequences", fill="Taxa Information") +
  theme_devon() +
  theme(legend.position = "top")


## notrun: pal12 <- c(rep('#ca0020',3), rep('#f4a582',3), rep('#92c5de',3), rep('#0571b0', 3))
## notrun: shape12 <- rep(c(16, 17, 15),4)

## notrun: ggplot(alt_miss, aes(x=Level, y=Value, color=Labeler, group=Labeler, shape=Labeler)) +
##   geom_point(size=4, alpha=0.7) +
##  scale_y_continuous(labels = comma) +
##  scale_color_manual(values = pal12) +
##  scale_shape_manual(values = shape12) +
##  geom_line(data=alt_miss %>% filter (InfoType=="Present")) +
##  geom_line(data=alt_miss %>% filter(InfoType=="Missing"), linetype="dotted") +
##  geom_line(data=alt_miss %>% filter(InfoType=="Ambiguous"), linetype="dashed") +
##  labs(x="", y="distinct sequences", color="", shape="") +
##  theme_devon() +
##  theme(legend.position = "top") +
##  guides(fill=guide_legend(ncol=3))

## plot - fraction of information by clustering; save as db5b_clustFractions; export at 800x600
alt_miss2 <- all_miss %>% select(Level, Database, pPresent, pMissing, pAmbiguous) %>% 
  gather(key = InfoType, value = Value, pPresent, pMissing, pAmbiguous)

alt_miss2$InfoType <- factor(alt_miss2$InfoType, levels = c("pPresent", "pMissing", "pAmbiguous"))

p2 <- ggplot(alt_miss2, aes(x=Level, y=Value, group=Database, fill=InfoType)) +  
  geom_bar(stat="identity") +
  facet_wrap(~ Database, nrow = 4) +
  coord_flip() +
  scale_fill_manual(values=c("gray70", "gray20", "red"),
                    labels=c("Present", "Missing", "Ambiguous"),
                    breaks=c("pPresent", "pMissing", "pAmbiguous")) +
  labs(x="", y="fraction of sequences", fill="Taxa Information") +
  theme_devon() +
  theme(legend.position = "top")

## plot together; save as 7_figure_dbClusteringAmbiguity; export at 
ggarrange(p1, p2, ncol = 2, common.legend = TRUE, labels = c("A", "B"))

### Examining how clustering impacts the proportion of information for specific arthropod Orders:
# 1) What are the most frequently observed arthropod Orders in our dereplicated dataset?
derep.topArthOrders <- derep_taxa %>% group_by(order_name) %>% summarise(OrderCounts=n()) %>% arrange(-OrderCounts)
## we're going to select the top 6 Orders, ad well as 3 others in that list we see in our bat samples a lot:
SelectOrders <- c("Diptera", "Lepidoptera", "Hymenoptera", "Coleoptera", "Hemiptera", "Araneae", "Trichoptera", "Psocodea", "Ephemeroptera")

## import data again, but this time convert any 'Ambiguous_name' value to NA
## import data:
derep_taxa2 <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/boldCOI.derep.txt.gz", delim = ";", col_names = FALSE)
colnames(derep_taxa2) <- c("SeqID", "class_name", "order_name", "family_name", "genus_name", "species_name")
derep_taxa2 <- as.data.frame(apply(derep_taxa2, 2, function(y) (gsub("Ambiguous_taxa", NA, y))))

clust99_taxa2 <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/clust99.names.txt.gz", delim = ";", col_names = FALSE)
colnames(clust99_taxa2) <- c("SeqID", "kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")
clust99_taxa2 <- clust99_taxa2 %>% select(-kingdom_name, -phylum_name)
clust99_taxa2 <- as.data.frame(apply(clust99_taxa2, 2, function(y) (gsub("Ambiguous_taxa", NA, y))))

clust97_taxa2 <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/clust97.names.txt.gz", delim = ";", col_names = FALSE)
colnames(clust97_taxa2) <- c("SeqID", "kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")
clust97_taxa2 <- clust97_taxa2 %>% select(-kingdom_name, -phylum_name)
clust97_taxa2 <- as.data.frame(apply(clust97_taxa2, 2, function(y) (gsub("Ambiguous_taxa", NA, y))))

clust95_taxa2 <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/clust95.names.txt.gz", delim = ";", col_names = FALSE)
colnames(clust95_taxa2) <- c("SeqID", "kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")
clust95_taxa2 <- clust95_taxa2 %>% select(-kingdom_name, -phylum_name)
clust95_taxa2 <- as.data.frame(apply(clust95_taxa2, 2, function(y) (gsub("Ambiguous_taxa", NA, y))))

## function to calculate proportion of information retained at Order, Family, Genus, and Species level depending on clustering %
taxOrderComp.function <- function(data, pid) {
  tmp1<- data %>% filter(order_name %in% SelectOrders) %>% 
    filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(family_name)) %>% filter(!is.na(genus_name)) %>% filter(!is.na(species_name)) %>%
    group_by(order_name) %>% summarise(nRecords = n()) %>% mutate(Rank = "Includes Species")
  tmp2 <- data %>% filter(order_name %in% SelectOrders) %>% 
    filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(family_name)) %>% filter(!is.na(genus_name)) %>% filter(is.na(species_name)) %>%
    group_by(order_name) %>% summarise(nRecords = n()) %>% mutate(Rank = "Missing Species")
  tmp3 <- data %>% filter(order_name %in% SelectOrders) %>% 
    filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(family_name)) %>% filter(is.na(genus_name)) %>% filter(is.na(species_name)) %>%
    group_by(order_name) %>% summarise(nRecords = n()) %>% mutate(Rank = "Missing Genus")
  tmp4 <- data %>% filter(order_name %in% SelectOrders) %>% 
    filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(is.na(family_name)) %>% filter(is.na(genus_name)) %>% filter(is.na(species_name)) %>%
    group_by(order_name) %>% summarise(nRecords = n()) %>% mutate(Rank = "Missing Family")
  tmp <- rbind(tmp1, tmp2, tmp3, tmp4) %>% mutate(Cluster=pid)
  Totals <- data %>% filter(order_name %in% SelectOrders) %>% group_by(order_name) %>% summarise(nTotal=n())
  out <- merge(tmp, Totals) %>% mutate(pRecords=nRecords/nTotal)
  out
}

derep.selectOrder_df <- taxOrderComp.function(derep_taxa2, "clust 100%")
clust99.selectOrder_df <- taxOrderComp.function(clust99_taxa2, "clust 99%")
clust97.selectOrder_df <- taxOrderComp.function(clust97_taxa2, "clust 97%")
clust95.selectOrder_df <- taxOrderComp.function(clust95_taxa2, "clust 95%")

all.selectOrder_df <- rbind(derep.selectOrder_df, clust99.selectOrder_df, clust97.selectOrder_df, clust95.selectOrder_df)
#rm(derep.selectOrder_df, clust99.selectOrder_df, clust97.selectOrder_df, clust95.selectOrder_df)

## plot
pal4 <- c("#d7191c", '#fdae61', '#abd9e9', '#2c7bb6')
all.selectOrder_df$Rank <- factor(all.selectOrder_df$Rank, levels=c("Includes Species", "Missing Species", "Missing Genus", "Missing Family"))
all.selectOrder_df$Cluster <- factor(all.selectOrder_df$Cluster, levels = c("clust 100%", "clust 99%", "clust 97%", "clust 95%"))

## plotting counts of records first
## save as db_13a_clusteringTopOrders_counts; export at 1000x800
p3 <- ggplot(all.selectOrder_df, 
       aes(y=nRecords, x=reorder(order_name, -nRecords), fill=Rank)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=pal4) +
  facet_grid(. ~ Cluster, scales = "free") +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  theme_devon() + theme(legend.position="top") +
  labs(x="", y="number of sequence records", fill="Record\ninformation")

## plotting percentage of records second
## save as db_13b_clusteringTopOrders_percent; export at 
p4 <- ggplot(all.selectOrder_df, 
       aes(y=pRecords, x=reorder(order_name, -nRecords), fill=Rank)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=pal4) +
  facet_grid(. ~ Cluster, scales = "free") +
  coord_flip() +
  theme_devon() + theme(legend.position="top") +
  labs(x="", y="fraction of sequence records", fill="Record\ninformation")

## plot; save as s10_figure_selectArths_clusteringReducesTaxaInfo; export at 1000x500
ggarrange(p3, NULL, p4, nrow=3, common.legend = TRUE, labels = c("A", "", "B"), heights = c(1.1, .1, 1.1))
