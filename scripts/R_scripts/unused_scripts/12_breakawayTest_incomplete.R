#runonce: install.packages("breakaway")
library(breakaway)
library(tidyverse)
library(reshape2)

## import data:
df <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz")
df <- df %>% filter(SampleType != "mock")
meta <- read_delim('https://github.com/devonorourke/pzero/raw/master/data/clean_metadata.txt', delim="\t")

## select just the data we want:
## select just the metadata we want:
meta_names <- dplyr::intersect(df$SeqID, meta$SeqID)
meta <- meta %>% filter(SeqID %in% meta_names)
meta <- meta %>% filter(Date != "unknown")
meta <- meta[c(5,8,14:15)]
meta$Site <- gsub("control", "ncontrol", meta$Site)
meta$Date <- gsub("control", "ncontrol", meta$Date)

## overwrite new $Date and $WOY columns with lubridate package to ensure we're selecting a consistent WOY
meta$Date <- as.character(lubridate::mdy(meta$Date))
meta$WOY <- as.character(lubridate::isoweek(meta$Date))
## adding in "startingMonth" to WOY for later groupings
WOYstring <- c("14",  "15",  "16",  "17",  "18",  "19",  "20",  "21",  "22",  "23",  "24",  "25",  "26",  "27",  "28",  "29",  "30",  "31",  "32",  "33",  "34",  "35",  "36",  "37",  "38",  "39",  "40",  "41",  "42",  "43", "ncontrol", "mock")
StartMonthstring <- c("April","April","April","April","April","May","May","May","May","June","June","June","June","July","July","July","July","July","August","August","August","August","September","September","September","September","October","October","October","October", "ncontrol", "mock")
tmp <- data.frame(WOYstring, StartMonthstring)
colnames(tmp) <- c("WOY", "MonthStart")
meta <- merge(meta, tmp, all.x=TRUE)
meta$MonthStart <- as.character(meta$MonthStart)
meta[is.na(meta)] <- "ncontrol"

## merge data:
df <- merge(df, meta, by='SeqID', all.x=TRUE)
rm(meta, meta_names, studytargets, tmp, StartMonthstring, WOYstring)

## Rename the Site and WOY functions for the 'ncontrol' values so that they reflect the..
## ..per-DNAplate where they came from... it's the best way I can think of to group negative control samples
tmp <- df %>% filter(SampleType=='ncontrol')
df$WOY[df$SampleType=='ncontrol'] <- paste("ncontrol", tmp$DNAplate, sep = ".")
df$Site[df$SampleType=='ncontrol'] <- paste("ncontrol", tmp$DNAplate, sep = ".")
rm(tmp)

## Select FOX site for testing:
fox.dada2.basic <- df %>% filter(Site=="FOX" & Method=="dada2" & Filt=="basic")

## create otu table
otu_data <- dcast(fox.dada2.basic, HashID ~ SeqID, value.var = "Reads", fill = 0)
row.names(otu_data) <- otu_data$HashID
otu_data$HashID <- NULL
## create metadata 
meta_data <- fox.dada2.basic %>% distinct(SeqID, SampleType, Library, WOY, Site, Date)
row.names(meta_data) <- meta_data$SeqID
meta_data$SeqID <- NULL
rm(fox.dada2.basic)

## generate frequency table
freq_df <- fox.dada2.basic %>% 
  count(SeqID, Reads) %>%
  group_by(SeqID)

## convert to list of frequency tables (one table per sample)
frequencytablelist <- split( freq_df , f = freq_df$SeqID )
## move sample names to row names and remove from column
list <- lapply(frequencytablelist,function(DF) {DF$SeqID <- NULL; DF})

## check it worked:
head(frequencytablelist[[3]]) ## the 3'rd sample
head(list[[3]]) ## the 3'rd sample
head(frequencytablelist[[32]])  ## the 32'nd sample

breakaway_nof1(frequencytablelist[[3]])
breakaway_nof1(list[[63]])

