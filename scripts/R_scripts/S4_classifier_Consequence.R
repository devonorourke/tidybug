library(tidyverse)
library(reshape2)
library(ggpubr)
library(ggrepel)
library(scales)

################################################################################
## 1. read in data
################################################################################

readerfunction <- function(urlpath, classifier) {
  tmp <- read_delim(file = urlpath, delim = "\t", col_names = FALSE) %>%
    rename(., HashID = X1, Taxon = X2) %>%
    separate(.,
             col = Taxon,
             into = c("kingdom", "phylum", "Class", "Order", "Family", "Genus", "Species"),
             sep = ";") %>%
    select(-kingdom, -phylum)
  tmp <- as.data.frame(apply(tmp, 2, function(y) (gsub(".__", "", y))))
  tmp <- as.data.frame(apply(tmp, 2, function(y) (gsub("^$|^ $", NA, y))))
  tmp <- as.data.frame(apply(tmp, 2, function(y) (gsub("Ambiguous taxa", NA, y))))
  tmp <- as.data.frame(apply(tmp, 2, function(y) (gsub("^NA$", NA, y))))
  tmp[] <- lapply(tmp, as.character)
  tmp %>%
    gather(key = level, value = Taxon, -HashID) %>%
    mutate(Classifier=classifier)
}

bl_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/guano_comps/blast/guano_blast_taxonomy.txt"
vs_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/guano_comps/vsearch/guano_vsearch_taxonomy.txt"
st_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/guano_comps/sintax/guano_sintax_taxonomy.txt"
nb_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/guano_comps/nbayes/guano_nbayes_taxonomy.txt"
bo_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/guano_comps/bold/guano_bold_taxonomy.txt"

bl_df <- readerfunction(bl_url, "blast")
vs_df <- readerfunction(vs_url, "vsearch")
st_df <- readerfunction(st_url, "sintax")
nb_df <- readerfunction(nb_url, "nbayes")
bo_df <- readerfunction(bo_url, "bold")

## combine df's
all_df <- rbind(bl_df, vs_df, st_df, nb_df, bo_df)
rm(bl_df, vs_df, st_df, nb_df, bo_df)
rm(bl_url, nb_url, st_url, vs_url, bo_url, readerfunction)

################################################################################
## 2. get instances where BOLD species missing, but present in all nonBOLD
################################################################################

boldMiss_df <- all_df %>%
  filter(level == 'Species') %>%
  spread(key = Classifier, value = Taxon) %>%
  mutate(missBOLDonly = case_when(
    !is.na(blast) & !is.na(vsearch) & !is.na(sintax) & !is.na(nbayes) & is.na(bold) ~ TRUE
  )) %>% 
  filter(missBOLDonly == TRUE)
  ## 1,641 instances where this is true
    ## sanity check: are all these HashID's distinct? length(unique(boldMiss_df$HashID)) ... yes

## get list of HashIDs for instances where just BOLD failed to provide a named species:
bold_miss_species_hash <- boldMiss_df %>% select(HashID) %>% pull()

## use this list to pull out the Order-level information of those missing taxa
## ignore the differences between BOLD and nonBOLD values - the difference occurs because of instances ...
## ...where BOLD has the Order but not the Species (+1 to count) or lacks the Species and the Order (0 to count)
boldMiss_taxa <- all_df %>%
  filter(HashID %in% bold_miss_species_hash) %>%
  filter(level == "Order") %>%
  group_by(Classifier, Taxon) %>%
  summarise(nTaxa = n()) %>%
  spread(key = Classifier, value = nTaxa, fill = 0) %>% 
  filter(!is.na(Taxon))
  ## write as "TableS7_missingBOLDspeciesNames_ComparedToNonBOLD"
#  write_csv(boldMiss_taxa, quote_escape = FALSE,
#            path = "~/Documents/nau_projects/guano/EcoEvo_methods_paper/TableS7_missingBOLDspeciesNames_ComparedToNonBOLD.csv")

## Most missing fall into the Orders we generally see repeatedly: Diptera, Coleoptera, Lepidoptera...
## but we're missing a ton of Megaloptera HashIDs... what are those ones??

## get species names of missing Megalopteras - how many distinct species do those missing ASVs represent?
## first, get HashIDs from BOLD-missing Megalopterans:
bold_miss_Megalopter_hashID <- all_df %>%
  filter(HashID %in% bold_miss_species_hash) %>%
  filter(level == "Order" & Taxon == "Megaloptera") %>%
  select(HashID) %>%
  distinct(HashID) %>%
  pull()
  ## there are 283 Megaloptera missing species names from BOLD that are named in all nonBOLD

## next, pull the species names from each Classifier and summarize distinct species counts
bold_miss_Megaloptera_df <- all_df %>%
  filter(HashID %in% bold_miss_Megalopter_hashID) %>%
  filter(level == "Species") %>%
  group_by(Classifier, Taxon) %>%
  summarise(nTaxa = n()) %>%
  spread(key = Classifier, value = nTaxa)
  ## Amazingly, those 283 missing ASVs represent just 4 named species!
    ## In our nonBOLD datasets, these compress into:
      # Chauliodes pectinicornis (249)
      # Chauliodes rastricornis (20)
      # Sialis vagans (1), and 
      # Sialis velata (13)
    ## those values were consistent among all 4 nonBOLD classifiers, by the way...

## but, are these missing species names from the 283 BOLD ASVs absent in all the other BOLD ASVs?
## in other words, are those 4 species named in other ASVs by BOLD?
bold_miss_Megalopter_speciesNames <- all_df %>%
  filter(HashID %in% bold_miss_Megalopter_hashID) %>%
  filter(level == "Species") %>%
  distinct(Taxon) %>% 
  filter(!is.na(Taxon)) %>% 
  pull()

bold_Megalopter_namedSpecies <- all_df %>%
  filter(Classifier == "bold" & Taxon %in% bold_miss_Megalopter_speciesNames) %>% # query just 4 Megalopteran species in bold
  group_by(Taxon) %>% 
  tally()
  ## here we see that just 1 of the 4 species is absent from our dataset: Silas veleta
  ## We find that BOLD has named other ASVs with these species names in these amounts:
    ## Chauliodes pectinicornis (16); 
    ## Chauliodes rastricornis (28); 
    ## Sialis vagans (1) 
  ## In other words, it clearly is missing a lot of our ASVs, but in terms of all fishfly's, we're really just missing a single species

#####################################################################3

## What about the really common bugs we see in guano?
## Repeat the same analysis for the missing Lepidoptera, Coleoptera, and Diptera:
ColDipLepTarget <- c("Coleoptera", "Diptera", "Lepidoptera")
bold_miss_ColDipLep_hashID <- all_df %>%
  filter(HashID %in% bold_miss_species_hash) %>%
  filter(level == "Order" & Taxon %in% ColDipLepTarget) %>%
  select(HashID) %>%
  distinct(HashID) %>%
  pull()
  ## there are 1,158 missing ASVs with that were assigned Coleoptera, Diptera, and Lepidoptera species names...
  ## in all nonBOLD classifiers, but missing from BOLD classifier

bold_miss_ColDipLep_df <- all_df %>%
  filter(HashID %in% bold_miss_ColDipLep_hashID) %>%
  filter(level == "Species") %>%
  group_by(Classifier, Taxon) %>%
  summarise(nTaxa = n()) %>%
  spread(key = Classifier, value = nTaxa) %>%
  select(-bold) %>%
  filter(complete.cases(.))
  ## in this case, we're missing hundreds of distinct species, so unlike Megaloptera, these ASVs don't necessarily all collapse to a few species
  ## many species named are rare, with only 36 of the species being assigned to more than 5 ASVs
  ## among the most frequently named species among the unassigned (by BOLD) ASVs:
    ## Pseudolimnophila luteipennis (182)
    ## Psilocorsis reflexella (57)
    ## Dendroides canadensis (52)
    ## Phyllophaga anxia (48)
    ## Malacosoma disstria (18)

## Do these Col/Dip/Lep ASVs that are frequently named a species name by nonBOLD exist in other BOLD ASVs?
bold_miss_ColDipLep_speciesNames <- all_df %>%
  filter(HashID %in% bold_miss_ColDipLep_hashID) %>%
  filter(level == "Species") %>%
  distinct(Taxon) %>% 
  filter(!is.na(Taxon)) %>% 
  pull()
  ## 341 distinct species names
  ## note this value differs from the number of rows in the 'bold_miss_ColDipLep_df' object because...
  ##...some Classifiers may assign different Species names to the same ASV!

bold_ColDipLep_namedSpecies <- all_df %>%
  filter(Classifier == "bold" & Taxon %in% bold_miss_ColDipLep_speciesNames) %>%
  group_by(Taxon) %>% 
  tally()
  ## interesting that our top hits aren't the ones assigned the most ASVs by BOLD necessarily:
    ## Coquillettidia perturbans (27) is top hit, but this was named in just (8) of the nonBOLD classifier ASVs
    ## Pseudolimnophila luteipennis (was 182, now just 5 additional)
    ## Psilocorsis reflexella (was 57, 0 additional) -- this is a species that is completely absent from BOLD
    ## Dendroides canadensis (was 52, now just 19 additional)
    ## Phyllophaga anxia (was 48, now just 18 additional)
    ## Malacosoma disstria (was 18, now just 6 additional)

## The main question:
# For those ASVs that were missing a BOLD species name, but had a nonBOLD species in all other classifiers,
  # (1) How often did the same species get assigned by BOLD to a different ASV? 
  # (2) How often were they just absent?

bold_miss_ALL_speciesNames <- all_df %>%
  filter(HashID %in% bold_miss_species_hash) %>%
  filter(level == "Species") %>%
  distinct(Taxon) %>% 
  filter(!is.na(Taxon)) %>% 
  pull()
  ## 432 unique species names among ASVs not assigned a species name by BOLD, but assigned by all nonBOLD

bold_ALL_namedSpecies <- all_df %>%
  filter(Classifier == "bold" & Taxon %in% bold_miss_ALL_speciesNames) %>%
  group_by(Taxon) %>% 
  tally()
  ## there are 131 unique species that BOLD does in fact assign to alternative ASVs
  ## this is the answer to: (1) How often did the same species get assigned by BOLD to a different ASV? 
  ## Given that there were 432 unique species names in the 'bold_miss_ALL_speciesNames' object...
  ## and 131 of these were found in alternative BOLD ASVs, there must be 301 instances where the ASVs...
  ## are unnamed in BOLD and lack any alternative ASV. 
    ## This is the answer to: (2) How often were they just absent?

## Can we plot the frequency of occurrence and sum of the reads among these two groups of named/notnamed species?

## get list of species names that fall into each group:
MissingSpeciesNameHasAlternativeASVNamedInBOLD <- bold_ALL_namedSpecies %>% select(Taxon) %>% pull()
NoMissingSpeciesNameAlternativeASVNamedInBOLD <- setdiff(bold_miss_ALL_speciesNames, MissingSpeciesNameHasAlternativeASVNamedInBOLD)

## get the list of all HashIDs associated with these samples in the BOLD dataset
MissingSpeciesNameHasAlternativeASVNamedInBOLD_hash <- all_df %>% 
  filter(Taxon %in% MissingSpeciesNameHasAlternativeASVNamedInBOLD) %>% 
  distinct(HashID) %>% 
  filter(HashID %in% boldMiss_df$HashID)
  ## 1,065 of 1,641 ASVs are those with named species names in BOLD, but aren't naming other ASVs that were named the same species by nonBOLD classifiers

NoMissingSpeciesNameAlternativeASVNamedInBOLD_hash <- all_df %>% 
  filter(Taxon %in% NoMissingSpeciesNameAlternativeASVNamedInBOLD) %>% 
  distinct(HashID) %>% 
  filter(HashID %in% boldMiss_df$HashID)
  ## 598 of 1,641 ASVs are those with no named species name in BOLD, and no alternative ASVs have the missing species name that was assigned by nonBOLD classifiers

## pull in the read data and select from the HashIDs present in the missing BOLD datasets but found in all nonBOLD
read_data <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz") %>%
  filter(Filt == "basic", Method == "dada2", SampleType == "sample", HashID %in% boldMiss_df$HashID) %>%
  group_by(HashID) %>%
  summarise(SumReads = sum(Reads), nObs=n())

## merge species names into plot as necessary:
TaxMap <- all_df %>%
  filter(level == "Species") %>%
  filter(complete.cases(.)) %>%
  select(HashID, Taxon) %>%
  distinct()

plotdat <- merge(read_data, TaxMap, all.x = TRUE) %>% 
  filter(Taxon != "Ambiguous_taxa")
  ## this value is greater than the original because we can have multiple species names...
  ## assigned to the same HashID by different classifiers (happens in just 54 instances of 1,544)

## to simplify, we'll just drop these instances
tmp <- plotdat %>% 
  group_by(HashID) %>% 
  tally() %>% 
  filter(n > 1) %>% 
  select(HashID) %>% pull()

## final data to plot includes a grouping variable for color/shape
plotdat <- plotdat %>% 
  filter(!HashID %in% disagreeTaxa_hash) %>% 
  mutate(Grouper = case_when(HashID %in% MissingSpeciesNameHasAlternativeASVNamedInBOLD_hash$HashID ~ "YesAlt",
                             HashID %in% NoMissingSpeciesNameAlternativeASVNamedInBOLD_hash$HashID ~ "NoAlt"))

  
## now plot!
ggplot(plotdat, aes(x=SumReads,y=nObs,shape=Grouper,colour=Grouper,label=Taxon)) +
  geom_point(size=3) +
  scale_x_continuous(trans = log10_trans(), labels = comma, limits=c(1,100000000), breaks = c(100,10000,1000000)) +
  geom_text_repel(data=plotdat %>%
                    filter(Grouper == 'NoAlt') %>%
                    filter(nObs > 19 & SumReads >= 10000),
                  size=3.5, segment.size  = 0.2, segment.color = "grey50",
                  nudge_y = 100, nudge_x = +100, force = 5) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), legend.position = "none") +
  scale_color_manual(values=c("dodgerblue", "chocolate3")) +
  labs(x="\nRead abundance per ASV", y="Number of ASVs detected in dataset\n", colour="",shape="")

## save as "S4_classifier_Consequence"; export at w950 x h600
