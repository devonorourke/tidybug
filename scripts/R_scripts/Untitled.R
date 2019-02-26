#runonce: install.packages('bold')
#runonce: install.packages('dplyr')
#runonce: install.packages('taxize')
library(bold)
library(dplyr)
library(taxize)

## getting list of all Arthropod names in BOLD
x <- downstream("Arthropoda", db = "ncbi", downto = "class")
x.nms <- x$Arthropoda$childtaxa_name
x.checks <- bold_tax_name(x.nms)

## generating two data pulls: one with just Insects, one with all non-Insects
nonInsect_names <- checks %>% filter(!is.na(taxonrep) & taxonrep != "Insecta") %>% select(taxon)

#nonInsect_names2 <- checks %>% filter(!is.na(taxonrep) & taxonrep != "Insecta" & taxonrep != "Arachnida") %>% select(taxon)

nonInsects_list <- lapply(nonInsect_names, bold_seqspec)
nonInsects_df <- do.call(rbind.data.frame, nonInsects_list) %>% 
  filter(markercode == "COI-5P") %>% 
  select(sequenceID, processid, bin_uri, genbank_accession, nucleotides, 
         country, institution_storing, 
         phylum_name, class_name, order_name, family_name, genus_name, species_name)
write.csv(nonInsects_df, file='allBOLD.nonInsects.csv', row.names = FALSE, quote = FALSE)


y <- downstream("Insecta", db = "ncbi", downto = "order")
y.nms <- y$Insecta$childtaxa_name
y.checks <- bold_tax_name(y.nms)
## these are missing a few names, oddly...
## generate the insect names with a bit of web scraping:
library('rvest')
library('stringr')
boldurl <- html("http://v4.boldsystems.org/index.php/Taxbrowser_Taxonpage?taxid=82")
boldtext <- boldurl %>% 
  html_nodes("div.col-md-6") %>%
  #html_nodes("div.ibox") %>%
  html_text()
tmptext <- substr(boldtext[7], start=71, stop=nchar(boldtext[7]))
tmptext
tmptext2 <- gsub('[[:digit:]]+', '', tmptext)
str(tmptext2)
tmptext3 <- unlist(strsplit(tmptext2, '\\['))
tmptext4 <- gsub('\\]', '', tmptext3)
tmptext5 <- str_trim(tmptext4)
insectNames <- tmptext5[tmptext5 != ""]
  rm(list=ls(pattern = "tmptext"))
## use the `insectNames` vector as the input for our list
z.checks <- bold_tax_name(insectNames)  
  ## now they're all identified!

Insects_list <- lapply(nonInsect_names, bold_seqspec)
Insects_df <- do.call(rbind.data.frame, Insects_df) %>% 
  filter(markercode == "COI-5P") %>% 
  select(sequenceID, processid, bin_uri, genbank_accession, nucleotides, 
         country, institution_storing, 
         phylum_name, class_name, order_name, family_name, genus_name, species_name)
write.csv(Insects_df, file='allBOLD.Insects.csv', row.names = FALSE, quote = FALSE)

## combine both data.frame objects into a full-arthropod list, and generate the taxonomy file and fasta file

