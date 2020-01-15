## specialthanks to Scott Chamberlain for the R 'bold' library
## runonce: install.packages("remotes")
#  library(remotes)
## runonce: remotes::install_github("ropensci/bold")
# can also run the async branch
## runonce: remotes::install_github("ropensci/bold@async")

library(readr)
library(dplyr)
library(bold)
library(tidyr)

df <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/mock_community/CFMR_insect_mock4.tsv", delim = "\t", col_names = FALSE)
colnames(df) <- c("id", "seq")

in_list <- as.list(setNames(df$seq, df$id))
out_list <- bold_identify(in_list, db = "COX1")
out_df <- do.call("rbind", lapply(out_list, data.frame))
out_df$query <- row.names(out_df)
out_df$query <- gsub("\\..*","", out_df$query)
colnames(out_df)[1] <- "processid"

## create three data.frames to generate winners from: 95%, 97%, and 99% identity
out_95_df <- out_df %>% select(processid, similarity, query) %>% filter(similarity >= 0.95)
out_97_df <- out_df %>% select(processid, similarity, query) %>% filter(similarity >= 0.97)
out_99_df <- out_df %>% select(processid, similarity, query) %>% filter(similarity >= 0.99) 

rm(out_list, in_list)

## load meta data to merge BOLD taxa records
##runfromlocal: meta <- read_delim(file = "~/Desktop/bold.meta.tmp.gz", delim = ";", col_names = FALSE)

## download the file from OSF account here: https://osf.io/k3eh6/files/
colnames(meta) <- c("sequenceID", "processid", "bin_uri", "genbank_accession", "country", "institution_storing", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")

## merge data frames
data_95 <- merge(out_95_df, meta, by='processid', all.x=TRUE) %>%
  select(query, sequenceID, phylum_name, class_name, order_name, family_name, genus_name, species_name) %>%
  filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(genus_name))

data_97 <- merge(out_97_df, meta, by='processid', all.x=TRUE) %>%
  select(query, sequenceID, phylum_name, class_name, order_name, family_name, genus_name, species_name) %>%
  filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(genus_name))

data_99 <- merge(out_99_df, meta, by='processid', all.x=TRUE) %>%
  select(query, sequenceID, phylum_name, class_name, order_name, family_name, genus_name, species_name) %>%
  filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(genus_name))

rm(meta)

## collect LCA results:

lca_95_tmp <- data_95 %>%
  gather(key, val, -query, -sequenceID) %>%
  select(-sequenceID, -key) %>%
  group_by(query, val) %>%
  add_count() %>%
  group_by(query) %>%
  filter(n == max(n)) %>%
  summarise(location_output = paste0(unique(val[!is.na(val)]), collapse = ";"))

lca_97_tmp <- data_97 %>%
  gather(key, val, -query, -sequenceID) %>%
  select(-sequenceID, -key) %>%
  group_by(query, val) %>%
  add_count() %>%
  group_by(query) %>%
  filter(n == max(n)) %>%
  summarise(location_output = paste0(unique(val[!is.na(val)]), collapse = ";"))

lca_99_tmp <- data_99 %>%
  gather(key, val, -query, -sequenceID) %>%
  select(-sequenceID, -key) %>%
  group_by(query, val) %>%
  add_count() %>%
  group_by(query) %>%
  filter(n == max(n)) %>%
  summarise(location_output = paste0(unique(val[!is.na(val)]), collapse = ";"))

standardformatfunction <- function(data) {
  data %>% 
  separate(location_output, into=c("phylum", "class", "order", "family", "genus", "species"),sep = ";") %>%
  mutate(phylum = gsub("^", "p__", phylum)) %>% 
  mutate(class = gsub("^", "c__", class)) %>% 
  mutate(order = gsub("^", "o__", order)) %>% 
  mutate(family = gsub("^", "f__", family)) %>% 
  mutate(genus = gsub("^", "g__", genus)) %>% 
  mutate(species = gsub("^", "s__", species)) %>% 
  mutate(Taxon = paste("k__Animalia", phylum, class, order, family, genus, species, sep = ";")) %>%
  select(query, Taxon)
}

lca_99_out <- standardformatfunction(lca_99_tmp)
lca_97_out <- standardformatfunction(lca_97_tmp)
lca_95_out <- standardformatfunction(lca_95_tmp)

## write files to disk; this is what was shared in 'tidybug' GitHub repo
write.table(lca_95_out, file="~/Repos/tidybug/data/classify_comps/mock_comps/observed_mockData/bold/boldAPI_mock_lca95.txt", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(lca_97_out, file="~/Repos/tidybug/data/classify_comps/mock_comps/observed_mockData/bold/boldAPI_mock_lca97.txt", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(lca_99_out, file="~/Repos/tidybug/data/classify_comps/mock_comps/observed_mockData/bold/boldAPI_mock_lca99.txt", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

## write QIIME-formatted files for evaluating recall and precision
qiime.function <- function(data) {
  data %>% 
    separate(., 
             col=location_output, 
             into=c('phylum', 'class', 'order', 'family', 'genus', 'species'),
             sep = ';') %>% 
    mutate(kingdom='Animalia') %>% 
    mutate(keeper=paste('k__',kingdom,';', 'p__',phylum, ';', 'c__',class, ';', 'o__',order, ';',
                        'f__',family, ';', 'g__',genus, ';', 's__',species, ';', sep = "")) %>% 
    select(query, keeper) %>% 
    mutate(keeper=gsub('*..__NA.*', '', .$keeper))
}

lca95_qiime <- qiime.function(lca_95_out)
lca97_qiime <- qiime.function(lca_97_out)
lca99_qiime <- qiime.function(lca_99_out)

write.csv(lca95_qiime, file="~/Repos/tidybug/data/classify_comps/mock_comps/observed_mockData/bold/boldAPI_mock_lca95_qiimeFormat.csv", row.names = FALSE, quote = FALSE)
write.csv(lca97_qiime, file="~/Repos/tidybug/data/classify_comps/mock_comps/observed_mockData/bold/boldAPI_mock_lca97_qiimeFormat.csv", row.names = FALSE, quote = FALSE)
write.csv(lca99_qiime, file="~/Repos/tidybug/data/classify_comps/mock_comps/observed_mockData/bold/boldAPI_mock_lca99_qiimeFormat.csv", row.names = FALSE, quote = FALSE)
