setwd("~/Repos/tidybug/data/databases/")
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

## import taxa data
palmer_taxa <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/palmer.taxa.txt.gz", delim = "\t")

porter_taxa <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/porter.taxa.txt.gz", delim = ";", col_names = FALSE)
colnames(porter_taxa) <- c("Phylum_name", "Class_name", "Order_name", "Family_name", "Genus_name", "Species_name")
