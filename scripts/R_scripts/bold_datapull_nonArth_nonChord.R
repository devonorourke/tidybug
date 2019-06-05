## gathering remaining non-Arthropod, non-Chordate records from BOLD
## also pulling the Fungi and Protist records too

library(bold)
library(dplyr)

## copied the non-Arthropod, non-Chordate Animal names directly from BOLD website: http://v4.boldsystems.org/index.php/TaxBrowser_Home
AnmlNames <- c('Acanthocephala', 'Acoelomorpha', 'Annelida', 'Brachiopoda', 
              'Bryozoa', 'Chaetognatha', 'Cnidaria', 'Ctenophora', 'Cycliophora', 
              'Echinodermata', 'Entoprocta', 'Gastrotricha', 'Gnathostomulida', 'Hemichordata', 
              'Kinorhyncha', 'Mollusca', 'Nematoda', 'Nematomorpha', 'Nemertea', 'Onychophora',
              'Phoronida', 'Placozoa', 'Platyhelminthes', 'Porifera', 'Priapulida', 'Rhombozoa', 
              'Rotifera','Sipuncula','Tardigrada', 'Xenacoelomorpha')


## download and filter BOLD records
Anml_list <- lapply(AnmlNames, bold_seqspec)

Anml_df <- do.call(rbind.data.frame, Anml_list) %>% 
  filter(markercode == "COI-5P") %>% 
  select(sequenceID, processid, bin_uri, genbank_accession, nucleotides, country, institution_storing, 
         phylum_name, class_name, order_name, family_name, genus_name, species_name)

x.taxon <- Anml_df %>% select(sequenceID, phylum_name, class_name, order_name, family_name, genus_name, species_name)
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

x.fasta <- Anml_df %>% select(sequenceID, nucleotides)
x.fasta <- merge(x.taxon, x.fasta)

write.csv(x.fasta, file = "boldCustom.allNonArthChordAnml.seqNtaxa.csv", quote = FALSE, row.names = FALSE)   ## file used to create taxonomy and fasta files

x.meta <- Anml_df %>% select(sequenceID, processid, bin_uri, genbank_accession, country, institution_storing,
                           phylum_name, class_name, order_name, family_name, genus_name, species_name)
write.table(x.meta, file = "boldCustom.allNonArthChordAnml.meta.txt", quote = TRUE, row.names = FALSE, sep = ";")   ## file used in further R analyses

### --------- repeat for Protist and Fungi:

nonAnmlNames <- c('Ascomycota', 'Basidiomycota', 'Chytridiomycota', 'Glomeromycota', 'Myxomycota', 'Zygomycota', 
                  'Chlorarachniophyta', 'Ciliophora', 'Heterokontophyta', 'Pyrrophycophyta')

nonAnml_list <- lapply(nonAnmlNames, bold_seqspec)

nonAnml_df <- do.call(rbind.data.frame, nonAnml_list) %>% 
  filter(markercode == "COI-5P") %>% 
  select(sequenceID, processid, bin_uri, genbank_accession, nucleotides, country, institution_storing, 
         phylum_name, class_name, order_name, family_name, genus_name, species_name)

x.taxon <- nonAnml_df %>% select(sequenceID, phylum_name, class_name, order_name, family_name, genus_name, species_name)
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

x.fasta <- nonAnml_df %>% select(sequenceID, nucleotides)
x.fasta <- merge(x.taxon, x.fasta)

write.csv(x.fasta, file = "boldCustom.allnonAnml.seqNtaxa.csv", quote = FALSE, row.names = FALSE)   ## file used to create taxonomy and fasta files

x.meta <- nonAnml_df %>% select(sequenceID, processid, bin_uri, genbank_accession, country, institution_storing,
                                 phylum_name, class_name, order_name, family_name, genus_name, species_name)
write.table(x.meta, file = "boldCustom.allnonAnml.meta.txt", quote = TRUE, row.names = FALSE, sep = ";")   ## file used in further R analyses