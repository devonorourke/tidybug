library(tidyverse)
library(scales)
library(viridis)
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
#notrun: derep_taxa <- as.data.frame(apply(derep_taxa, 2, function(y) (gsub("Ambiguous_taxa", NA, y))))
clust99_names <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/clust99.names.txt.gz", delim = ";", col_names = FALSE)
colnames(clust99_names) <- "SeqID"
clust97_names <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/clust97.names.txt.gz", delim = ";", col_names = FALSE)
colnames(clust97_names) <- "SeqID"
clust95_names <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/clust95.names.txt.gz", delim = ";", col_names = FALSE)
colnames(clust95_names) <- "SeqID"

## generate missingness per level
derep_miss <- missfunction(derep_taxa, "derep")
clust99_miss <- missfunction(derep_taxa %>% filter(SeqID %in% clust99_names$SeqID), "clust99")
clust97_miss <- missfunction(derep_taxa %>% filter(SeqID %in% clust97_names$SeqID), "clust97")
clust95_miss <- missfunction(derep_taxa %>% filter(SeqID %in% clust95_names$SeqID), "clust95")

rm(clust95_names, clust97_names, clust99_names, derep_taxa)
all_miss <- rbind(derep_miss, clust99_miss, clust97_miss, clust95_miss)

## levels for plot
all_miss$Level <- factor(all_miss$Level, levels = c("Class", "Order", "Family", "Genus", "Species"))
all_miss$Database <- factor(all_miss$Database, levels = c("derep", "clust99", "clust97", "clust95"))

## plot - total number of records; save as db5_clustRecords; export at 800x600
alt_miss <- all_miss %>% select(Level, Database, Present, Missing, ambiguous) %>% gather(key = InfoType, value = Value, Present, Missing, ambiguous)
alt_miss$InfoType <- gsub("ambiguous", "Ambiguous", alt_miss$InfoType)
alt_miss$InfoType <- gsub("ambiguous", "Ambiguous", alt_miss$InfoType)
alt_miss$Labeler <- paste(alt_miss$Database, alt_miss$InfoType, sep=" - ")
alt_miss$Labeler <- factor(alt_miss$Labeler, levels = c(
  "derep - Present", "derep - Missing", "derep - Ambiguous","clust99 - Present", "clust99 - Missing", "clust99 - Ambiguous",
  "clust97 - Present", "clust97 - Missing", "clust97 - Ambiguous","clust95 - Present", "clust95 - Missing", "clust95 - Ambiguous"))

pal12 <- c(rep('#ca0020',3), rep('#f4a582',3), rep('#92c5de',3), rep('#0571b0', 3))
shape12 <- rep(c(16, 17, 15),4)

ggplot(alt_miss, aes(x=Level, y=Value, color=Labeler, group=Labeler, shape=Labeler)) +
  geom_point(size=4, alpha=0.7) +
  scale_y_continuous(labels = comma) +
  scale_color_manual(values = pal12) +
  scale_shape_manual(values = shape12) +
  geom_line(data=alt_miss %>% filter (InfoType=="Present")) +
  geom_line(data=alt_miss %>% filter(InfoType=="Missing"), linetype="dotted") +
  geom_line(data=alt_miss %>% filter(InfoType=="Ambiguous"), linetype="dashed") +
  labs(x="", y="distinct sequences", color="", shape="") +
  theme_devon() +
  theme(legend.position = "top") +
  guides(fill=guide_legend(ncol=3))

## plot - fraction of information by clustering; save as db5b_clustFractions; export at 800x600
alt_miss2 <- all_miss %>% select(Level, Database, pPresent, pMissing, pAmbiguous) %>% gather(key = InfoType, value = Value, pPresent, pMissing, pAmbiguous)

ggplot(alt_miss2, aes(x=Level, y=Value, group=Database, fill=InfoType)) +  
  geom_bar(stat="identity", position="dodge") +
  facet_wrap(~ Database, nrow = 4) +
  coord_flip() +
  scale_fill_manual(values=c("gray20", "gray50", "gray80"),
                    labels=c("Ambigous", "Missing", "Present"),
                    breaks=c("pAmbiguous", "pMissing", "pPresent")) +
  labs(x="", y="fraction of sequences", fill="Taxa Information") +
  theme_devon() +
  theme(legend.position = "top", legend.text = )
