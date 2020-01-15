#runonce: install.packages('bold')
#runonce: install.packages('dplyr')
#runonce: install.packages('taxize')
library(bold)
library(dplyr)
library(taxize)
library(rvest)
library(stringr)

## getting list of all Arthropod names in NCBI
x <- downstream("Chordata", db = "ncbi", downto = "class")
x.nms <- x$Chordata$childtaxa_name
x.checks <- bold_tax_name(x.nms)
  ## NCBI only reports 9 names, some missing from BOLD list... 
  ## this indicates we need to try another method to pull the complete set of names from BOLD directly

## generate the Chordata names with a bit of web scraping:
boldurl <- read_html("http://v4.boldsystems.org/index.php/Taxbrowser_Taxonpage?taxid=18")
boldtext <- boldurl %>% 
  html_nodes("div.col-md-6") %>%
  html_text()
tmptext <- substr(boldtext[7], start=72, stop=nchar(boldtext[7]))
tmptext2 <- gsub('[[:digit:]]+', '', tmptext)
tmptext3 <- unlist(strsplit(tmptext2, '\\['))
tmptext4 <- gsub('\\]', '', tmptext3)
tmptext5 <- str_trim(tmptext4)
ChordateNames <- tmptext5[tmptext5 != ""]
rm(list=ls(pattern = "tmptext"))


## download and filter BOLD records
Chordate_list <- lapply(ChordateNames, bold_seqspec)

Chordate_df <- do.call(rbind.data.frame, Chordate_list) %>% 
  filter(markercode == "COI-5P") %>% 
  select(sequenceID, processid, bin_uri, genbank_accession, nucleotides, country, institution_storing, 
         phylum_name, class_name, order_name, family_name, genus_name, species_name)

rm(Chordate_df)

## combine all data.frame objects into a full-chordate data.frame and generate the taxonomy file and fasta files

x.taxon <- Chordate_df %>% select(sequenceID, phylum_name, class_name, order_name, family_name, genus_name, species_name)
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

x.fasta <- Chordate_df %>% select(sequenceID, nucleotides)
x.fasta <- merge(x.taxon, x.fasta)

write.csv(x.fasta, file = "boldCustom.allChordate.seqNtaxa.csv", quote = FALSE, row.names = FALSE)   ## file used to create taxonomy and fasta files

x.meta <- chordate_df %>% select(sequenceID, processid, bin_uri, genbank_accession, country, institution_storing,
                           phylum_name, class_name, order_name, family_name, genus_name, species_name)
write.table(x.meta, file = "boldCustom.allChordate.meta.txt", quote = TRUE, row.names = FALSE, sep = ";")   ## file used in further R analyses