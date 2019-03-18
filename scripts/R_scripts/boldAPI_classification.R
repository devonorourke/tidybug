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

df <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/qiime/dada2.arthASVs.txt.gz", delim = "\t", col_names = FALSE)
colnames(df) <- c("id", "seq")

#### needed to break up data.frame with:
#   df1 <- df[1:3000,]
#   df2 <- df[3001:6500,]
#   df3 <- df[6501:10000,]
#   df4 <- df[10001:13407,]

## then iterated by substituting these
# in_list <- as.list(setNames(df1$seq, df1$id))
# in_list <- as.list(setNames(df2$seq, df2$id))
# in_list <- as.list(setNames(df3$seq, df3$id))
# in_list <- as.list(setNames(df4$seq, df4$id))


out_list <- bold_identify(in_list, db = "COX1")
out_df <- do.call("rbind", lapply(out_list, data.frame))
out_df$query <- row.names(out_df)
out_df$query <- gsub("\\..*","", out_df$query)
colnames(out_df)[1] <- "processid"

## create two data.frames to generate winners from: 97% identity and 99% identity
out_97_df <- out_df %>% select(processid, similarity, query) %>% filter(similarity >= 0.97)
out_99_df <- out_df %>% select(processid, similarity, query) %>% filter(similarity >= 0.99) 

rm(out_list, in_list)

## load meta data to merge BOLD taxa records
meta <- read_delim(file = "~/Desktop/bold.meta.tmp.gz", delim = ";", col_names = FALSE)
colnames(meta) <- c("sequenceID", "processid", "bin_uri", "genbank_accession", "country", "institution_storing", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")

## merge data frames
data_97 <- merge(out_97_df, meta, by='processid', all.x=TRUE) %>%
  select(query, sequenceID, phylum_name, class_name, order_name, family_name, genus_name, species_name) %>%
  filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(genus_name))

data_99 <- merge(out_99_df, meta, by='processid', all.x=TRUE) %>%
  select(query, sequenceID, phylum_name, class_name, order_name, family_name, genus_name, species_name) %>%
  filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>% filter(!is.na(order_name)) %>% filter(!is.na(genus_name))

## collect LCA results:
lca_97_out <- data_97 %>%
  gather(key, val, -query, -sequenceID) %>%
  select(-sequenceID, -key) %>%
  group_by(query, val) %>%
  add_count() %>%
  group_by(query) %>%
  filter(n == max(n)) %>%
  summarise(location_output = paste0(unique(val[!is.na(val)]), collapse = ";"))

lca_99_out <- data_99 %>%
  gather(key, val, -query, -sequenceID) %>%
  select(-sequenceID, -key) %>%
  group_by(query, val) %>%
  add_count() %>%
  group_by(query) %>%
  filter(n == max(n)) %>%
  summarise(location_output = paste0(unique(val[!is.na(val)]), collapse = ";"))

## wrote individual parts (df1 through df4) to .csv files, then reimported
## example for p97: 
     #    write.csv(lca_97_out, file="~/Desktop/boldAPI_p97-part4_taxa.csv", row.names = FALSE, quote = TRUE)
     #    write.csv(lca_99_out, file="~/Desktop/boldAPI_p99_taxa-part4.csv", row.names = FALSE, quote = TRUE)

## import individual files:
lca_97_1 <- read.csv("~/Desktop/boldAPI_p97-part1_taxa.csv")
lca_97_2 <- read.csv("~/Desktop/boldAPI_p97-part2_taxa.csv")
lca_97_3 <- read.csv("~/Desktop/boldAPI_p97-part3_taxa.csv")
lca_97_4 <- read.csv("~/Desktop/boldAPI_p97-part4_taxa.csv")
lca_97 <- rbind(lca_97_1, lca_97_2, lca_97_3, lca_97_4)
rm(lca_97_1, lca_97_2, lca_97_3, lca_97_4)

lca_99_1 <- read.csv("~/Desktop/boldAPI_p99_taxa-part1.csv")
lca_99_2 <- read.csv("~/Desktop/boldAPI_p99_taxa-part2.csv")
lca_99_3 <- read.csv("~/Desktop/boldAPI_p99_taxa-part3.csv")
lca_99_4 <- read.csv("~/Desktop/boldAPI_p99_taxa-part4.csv")
lca_99 <- rbind(lca_99_1, lca_99_2, lca_99_3, lca_99_4)
rm(lca_99_1, lca_99_2, lca_99_3, lca_99_4)

## write combined files to disk; this is what was shared in 'tidybug' GitHub repo
write.csv(lca_97, file="~/Repos/tidybug/data/databases/boldAPI_lca97.csv")
write.csv(lca_99, file="~/Repos/tidybug/data/databases/boldAPI_lca99.csv")
