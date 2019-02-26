#runonce: install.packages('bold')
#runonce: install.packages('dplyr')
#runonce: install.packages('taxize')
library(bold)
library(dplyr)
library(taxize)
library(rvest)
library(stringr)

## getting list of all Arthropod names in BOLD
x <- downstream("Arthropoda", db = "ncbi", downto = "class")
x.nms <- x$Arthropoda$childtaxa_name
x.checks <- bold_tax_name(x.nms)

## generating two data pulls: one with just Insects, one with all non-Insects
nonInsect_names <- x.checks %>% filter(!is.na(taxon) & taxon != "Insecta") %>% select(taxon)

nonInsects_list <- lapply(nonInsect_names, bold_seqspec)
nonInsects_df <- do.call(rbind.data.frame, nonInsects_list) %>% 
  filter(markercode == "COI-5P") %>% 
  select(sequenceID, processid, bin_uri, genbank_accession, nucleotides, 
         country, institution_storing, 
         phylum_name, class_name, order_name, family_name, genus_name, species_name)

## generate the insect names with a bit of web scraping:
boldurl <- read_html("http://v4.boldsystems.org/index.php/Taxbrowser_Taxonpage?taxid=82")
boldtext <- boldurl %>% 
  html_nodes("div.col-md-6") %>%
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

## diptera only
tmpurl <- 'http://www.boldsystems.org/index.php/Taxbrowser_Taxonpage?taxon=diptera&searchTax=Search+Taxonomy'
boldurl <- read_html(tmpurl)
boldtext <- boldurl %>% html_nodes("div.col-md-6") %>% html_text()
tmptext <- substr(boldtext[7], start=71, stop=nchar(boldtext[7]))
tmptext
tmptext2 <- gsub('[[:digit:]]+', '', tmptext)
str(tmptext2)
tmptext3 <- unlist(strsplit(tmptext2, '\\['))
tmptext4 <- gsub('\\]', '', tmptext3)
tmptext5 <- gsub(')', '', tmptext4)
tmptext6 <- str_trim(tmptext5)
DipNames <- tmptext6[tmptext6 != ""]
rm(list=ls(pattern = "tmptext"))


Dip_list <- lapply(DipNames, bold_seqspec)
Dip_df <- do.call(rbind.data.frame, Dip_list) %>% 
  filter(markercode == "COI-5P") %>% 
  select(sequenceID, processid, bin_uri, genbank_accession, nucleotides, 
         country, institution_storing, 
         phylum_name, class_name, order_name, family_name, genus_name, species_name)

## Coloeopterans only
tmpurl <- 'http://www.boldsystems.org/index.php/Taxbrowser_Taxonpage?taxon=coleoptera&searchTax=Search+Taxonomy'
boldurl <- read_html(tmpurl)
boldtext <- boldurl %>% html_nodes("ol:nth-child(2)") %>% html_text()
tmptext <- gsub('[[:digit:]]+', '', boldtext)
tmptext2 <- unlist(strsplit(tmptext, '\\['))
tmptext3 <- gsub('\\]', '', tmptext2)
tmptext4 <- gsub(')', '', tmptext3)
tmptext5 <- str_trim(tmptext4)
## there is a 'Genera' listing we have to remove here too:
tmptext6 <- gsub("Genera \\( coleopteraJanzen", "", tmptext5)
ColeoNames <- tmptext6[tmptext6 != ""]
rm(list=ls(pattern = "tmptext"))

Coleo_list <- lapply(ColeoNames, bold_seqspec)
Coleo_df <- do.call(rbind.data.frame, Coleo_list) %>% 
  filter(markercode == "COI-5P") %>% 
  select(sequenceID, processid, bin_uri, genbank_accession, nucleotides, 
         country, institution_storing, 
         phylum_name, class_name, order_name, family_name, genus_name, species_name)

## Lepidoptera
tmpurl <- 'http://www.boldsystems.org/index.php/Taxbrowser_Taxonpage?taxon=lepidoptera&searchTax=Search+Taxonomy'
boldurl <- read_html(tmpurl)
boldtext <- boldurl %>% html_nodes("ol:nth-child(2)") %>% html_text()
tmptext <- gsub('[[:digit:]]+', '', boldtext)
tmptext2 <- unlist(strsplit(tmptext, '\\['))
tmptext3 <- gsub('\\]', '', tmptext2)
tmptext4 <- gsub(')', '', tmptext3)
tmptext5 <- str_trim(tmptext4)
tmptext5
LepNames <- tmptext5[tmptext5 != ""]
rm(list=ls(pattern = "tmptext"))

Lep_list <- lapply(LepNames, bold_seqspec)
Lep_df <- do.call(rbind.data.frame, Lep_list) %>% 
  filter(markercode == "COI-5P") %>% 
  select(sequenceID, processid, bin_uri, genbank_accession, nucleotides, 
         country, institution_storing, 
         phylum_name, class_name, order_name, family_name, genus_name, species_name)

## Hymenoptera
tmpurl <- 'http://www.boldsystems.org/index.php/Taxbrowser_Taxonpage?taxon=hymenoptera&searchTax=Search+Taxonomy'
boldurl <- read_html(tmpurl)
boldtext <- boldurl %>% html_nodes("ol:nth-child(2)") %>% html_text()
tmptext <- gsub('[[:digit:]]+', '', boldtext)
tmptext2 <- unlist(strsplit(tmptext, '\\['))
tmptext3 <- gsub('\\]', '', tmptext2)
tmptext4 <- gsub(')', '', tmptext3)
tmptext5 <- str_trim(tmptext4)
tmptext5
HymnNames <- tmptext5[tmptext5 != ""]
rm(list=ls(pattern = "tmptext"))
rm(boldurl, boldtext, tmpurl)

Hymn_list <- lapply(HymnNames, bold_seqspec)
Hymn_df <- do.call(rbind.data.frame, Hymn_list) %>% 
  filter(markercode == "COI-5P") %>% 
  select(sequenceID, processid, bin_uri, genbank_accession, nucleotides, 
         country, institution_storing, 
         phylum_name, class_name, order_name, family_name, genus_name, species_name)

## for all the rest of the Insect Orders
excludeNames <- c("Coleoptera|Diptera|Hymenoptera|Lepidoptera")
otherInsectOrders <- str_remove_all(insectNames, excludeNames)
otherInsectOrders <- otherInsectOrders[otherInsectOrders != ""]

otherInsects_list <- lapply(otherInsectOrders, bold_seqspec)
otherInsects_df <- do.call(rbind.data.frame, otherInsects_list) %>% 
  filter(markercode == "COI-5P") %>% 
  select(sequenceID, processid, bin_uri, genbank_accession, nucleotides, 
         country, institution_storing, 
         phylum_name, class_name, order_name, family_name, genus_name, species_name)


## combine all data.frame objects into a full-arthropod list, and generate the taxonomy file and fasta file
x.raw <- rbind(nonInsects_df, otherInsects_df, Hymn_df, Lep_df, Coleo_df, Dip_df)
rm(nonInsects_df, otherInsects_df, Hymn_df, Lep_df, Coleo_df, Dip_df)

x.taxon <- x.raw %>% select(sequenceID, phylum_name, class_name, order_name, family_name, genus_name, species_name)
x.taxon$kingdom_name <- "k__Animalia"
x.taxon$phylum_name <- x.taxon$phylum_name %>% replace(is.na(.), "") %>% sub("^", "p__", .)
x.taxon$class_name <- x.taxon$class_name %>% replace(is.na(.), "") %>% sub("^", "c__", .)
x.taxon$order_name <- x.taxon$order_name %>% replace(is.na(.), "") %>% sub("^", "o__", .)
x.taxon$family_name <- x.taxon$family_name %>% replace(is.na(.), "") %>% sub("^", "f__", .)
x.taxon$genus_name <- x.taxon$genus_name %>% replace(is.na(.), "") %>% sub("^", "g__", .)
x.taxon$species_name <- x.taxon$species_name %>% replace(is.na(.), "") %>% sub("^", "s__", .)
x.taxon$taxon <- paste(x.taxon$kingdom_name, x.taxon$phylum_name, x.taxon$class_name, 
                       x.taxon$order_name, x.taxon$family_name, x.taxon$genus_name, 
                       x.taxon$species_name, sep = ";")
x.taxon <- x.taxon %>% select(sequenceID, taxon)

x.fasta <- x.raw %>% select(sequenceID, nucleotides)
x.fasta <- merge(x.taxon, x.fasta)
rm(x.taxon, x.seqs)
write.csv(x.fasta, file = "boldCustom.allArth.seqNtaxa.csv", quote = FALSE, row.names = FALSE)   ## file used to create taxonomy and fasta files

x.meta <- x.raw %>% select(sequenceID, processid, bin_uri, genbank_accession, country, institution_storing,
                           phylum_name, class_name, order_name, family_name, genus_name, species_name)
write.table(x.meta, file = "boldCustom.allArth.meta.txt", quote = TRUE, row.names = FALSE, sep = ";")   ## file used in further R analyses