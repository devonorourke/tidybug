library(tidyverse)

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


## import data:
## gather metadata from raw data, subset from list of IDs present in dereplicated dataset
raw_meta <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/reduced.allArth.meta.txt.gz", delim = ";", col_names = TRUE)
derep_taxa <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/boldCOI.derep.txt.gz", delim = ";", col_names = FALSE)
colnames(derep_taxa) <- c("sequenceID", "class_name", "order_name", "family_name", "genus_name", "species_name")
derep_meta <- raw_meta %>% filter(sequenceID %in% derep_taxa$sequenceID) %>% select(sequenceID, country, institution_storing)
rm(raw_meta)
derep_df <- merge(derep_meta, derep_taxa)
rm(derep_meta, derep_taxa)

#notrun: we aren't going to substitute ambiguous taxa as NA because we're quantifying how often they ambiguous taxa are present intitially!
#notrun: derep_taxa <- as.data.frame(apply(derep_taxa, 2, function(y) (gsub("Ambiguous_taxa", NA, y))))


## list of countries?
CountryCounts <- derep_df %>% group_by(country) %>% summarise(CountryCounts=n())
CountryCounts_withTaxa <- derep_df %>% group_by(country, order_name) %>% summarise(CountryCounts=n())

## Selecting pairs of countries that share neighboring country where one country is 10x as sampled as other
CountryList <- c("Canada", "United States", "Costa Rica", "Panama", "Germany", "Austria", "Finland", "Sweden", NA)
select_df <- derep_df %>% filter(country %in% CountryList)
rm(derep_df)

# which of these orders are in common for each country pair?
intersect.function <- function(data, location, pairname){
  data %>%
    filter(country %in% location) %>%
    group_by(country, order_name) %>%
    summarise(nTaxa=n()) %>%
    mutate(Pair=pairname) %>%
    group_by(country) %>%
    top_n(., 10, nTaxa) %>%
    slice(., 1:10) %>%    ## picking first winner in whatever ties exist in 10th value...
    mutate(pTaxa=nTaxa/sum(nTaxa)) 
}

## four pairs of countries to compare:
northamerica <- c("Canada", "United States")
scandanavia <- c("Finland", "Sweden")
bavaria <- c("Germany", "Austria")
centralamerica <- c("Costa Rica", "Panama")

## gather data
topOrders.northamerica <- intersect.function(select_df, northamerica, "United States | Canada")
topOrders.scandanavia <- intersect.function(select_df, scandanavia, "Finland | Sweden")
topOrders.bavaria <- intersect.function(select_df, bavaria, "Germany | Austria")
topOrders.centralamerica <- intersect.function(select_df, centralamerica, "Costa Rica | Panama")
topOrders.unlisted <- intersect.function(select_df, NA, "Unlisted")
topOrders.unlisted <- topOrders.unlisted %>% ungroup(.) %>% mutate(country = replace_na(country, "unlisted"))

## merge:
topOrders.all <- rbind(topOrders.northamerica, topOrders.scandanavia, topOrders.bavaria, topOrders.centralamerica,topOrders.unlisted)

## identify private Orders from each country pair
uniqOrder.function <- function(data, filter_exp){
  filter_exp_enq <- enquo(filter_exp)
  tmp <- data %>% 
    filter(!!filter_exp_enq) %>% 
    ungroup(.) %>% select(order_name)
  as.character(tmp$order_name)
}

us.orders <- uniqOrder.function(topOrders.all, country=="United States")
cdn.orders <- uniqOrder.function(topOrders.all, country=="Canada")
fin.orders <- uniqOrder.function(topOrders.all, country=="Finland")
swe.orders <- uniqOrder.function(topOrders.all, country=="Sweden")
ger.orders <- uniqOrder.function(topOrders.all, country=="Germany")
aus.orders <- uniqOrder.function(topOrders.all, country=="Austria")
cos.orders <- uniqOrder.function(topOrders.all, country=="Costa Rica")
pan.orders <- uniqOrder.function(topOrders.all, country=="Panama")

intersect.northamerica <- intersect(us.orders, cdn.orders)
intersect.scandanavia <- intersect(fin.orders, swe.orders)
intersect.bavaria <- intersect(ger.orders, aus.orders)
intersect.centralamerica <- intersect(cos.orders, pan.orders)

## filter dataset to retain just those intersecting taxa:
select.northamerica <- topOrders.northamerica %>% filter(order_name %in% intersect.northamerica)
select.centralamerica <- topOrders.centralamerica %>% filter(order_name %in% intersect.centralamerica)
select.bavaria <- topOrders.bavaria %>% filter(order_name %in% intersect.bavaria)
select.scandanavia <- topOrders.scandanavia %>% filter(order_name %in% intersect.scandanavia)

select.all <- rbind(select.northamerica, select.centralamerica, select.bavaria, select.scandanavia)

## plot
select.all$Pair <- factor(select.all$Pair, levels = c("Costa Rica | Panama", "Finland | Sweden", "Germany | Austria", "United States | Canada"))
select.all$country <- factor(select.all$country, levels = c("Costa Rica", "Panama", "Finland", "Sweden", "Germany", "Austria", "United States", "Canada"))

pal10 <- c("gray65", "#f58231", "#808000", "#469990", "#000075",
           "#e6194B", "#9A6324", "white", "#42d4f4", "#f032e6")

## plot; save as db_12_databaseCompositionByCountry; export at 800x600
ggplot(data = select.all, aes(x=country, y=pTaxa, fill=order_name)) + 
  geom_bar(stat="identity", color="black", size=0.5) + 
  facet_grid(. ~ Pair, scales = "free_x") +
  scale_fill_manual(values=pal10) +
  theme_devon() +
  labs(x="", y="fraction of records", fill="Taxonomic\nOrder") +
  theme(legend.position = "top")


## can also plot the "Unlisted" ones, but note that we're going ot have a few different colors
altpal10 <- c("#fabebe", "gray65", "#f58231","#ffd8b1", "#808000", "#000075",
              "#e6194B", "#9A6324", "#fffac8", "#aaffc3")

## plot; save as db_12b_databaseCompositionByCountry-Unlisted; export at 800x600
ggplot(data = topOrders.unlisted, aes(x=country, y=pTaxa, fill=order_name)) + 
  geom_bar(stat="identity", color="black", size=0.5) + 
  facet_grid(. ~ Pair, scales = "free_x") +
  scale_fill_manual(values=altpal10) +
  theme_devon() +
  labs(x="", y="fraction of records", fill="Taxonomic\nOrder") +
  theme(plot.margin = unit(c(0,11,0,11), "cm"))
