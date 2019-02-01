# Table used in paper

library(tidyverse)
library(reshape2)

## import data and select mock samples:
df <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz")
mock <- df %>% filter(SampleType == "mock")
## generate filtered mock data
Sample_filtered <- df %>% group_by(SeqID) %>% summarise(SumReads=sum(Reads)) %>% filter(SumReads >= 5000) ## require at least 5000 reads per sample
df.filt <- df %>% filter(SeqID %in% Sample_filtered$SeqID)
rm(df)
Hash_filtered <- df.filt %>% group_by(HashID) %>% summarise(HashCounts=n()) %>% filter(HashCounts > 1)
df.filt <- df.filt %>% filter(HashID %in% Hash_filtered$HashID)
mock.filt <- df.filt %>% filter(SampleType == "mock")

## import mock matches (exact, partial, and misses)
mockExact <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/mock_community/mock.exactHits.txt", col_names = FALSE)
mockExact$Type <- "exact"
colnames(mockExact)[1] <- "HashID"
mockPartial <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/mock_community/mock.partialHits.txt", col_names = FALSE)
mockPartial$Type <- "partial"
colnames(mockPartial)[1] <- "HashID"
mockMiss <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/mock_community/mock.partialMisses.txt", col_names = FALSE)
mockMiss$Type <- "miss"
colnames(mockMiss)[1] <- "HashID"
mockHitTypes <- rbind(mockExact, mockPartial, mockMiss)
rm(mockExact, mockPartial, mockMiss)

## merge mockHitTypes $Type to mock object (by $HashID)
mock <- merge(mock, mockHitTypes)
mock.filt <- merge(mock.filt, mockHitTypes)
rm(mockHitTypes)
## add Labeler for faceting plots
mock$Labeler <- "unfiltered"
mock.filt$Labeler <- "filtered"
## merge together
mock.all <- rbind(mock, mock.filt)
rm(mock, mock.filt)

## tiny summary of mock data, per Library, per filtering/not
mock_Sumry <- mock.all %>% group_by(Library, Method, Labeler, Type) %>% summarise(CountTypes=n())
## convert into matrix of values
filt.exact.mock_Sumry <- mock_Sumry %>% filter(Labeler=="filtered" & Type=="exact")
unfilt.exact.mock_Sumry <- mock_Sumry %>% filter(Labeler=="unfiltered" & Type=="exact")
filt.exact.mock_SumryTable <- dcast(data=filt.exact.mock_Sumry, Library ~ Method, value.var = "CountTypes")
unfilt.exact.mock_SumryTable <- dcast(data=unfilt.exact.mock_Sumry, Library ~ Method, value.var = "CountTypes")
## we'll look to find how many shared ASVs occur between Libraries (within Filtering algorithm) ..
## and how many shared ASVs occur across Libraries (between Filtering algorithms)
## these values are manually added to the two tables presented above

## to find potential intersections (shared ASVs) by various parameters
## generating vectors of individual exact-match HashIDs, filtered by Method, Library and filtering Labeler
## Only focusing on unfiltered exact matches
exact.unfilt.vs.libA <- mock.all %>% filter(Type=="exact" & Labeler=="unfiltered" & Method=="vsearch" & Library=="libA") %>% select(HashID)
exact.unfilt.vs.libB <- mock.all %>% filter(Type=="exact" & Labeler=="unfiltered" & Method=="vsearch" & Library=="libB") %>% select(HashID)
exact.unfilt.vs.libC <- mock.all %>% filter(Type=="exact" & Labeler=="unfiltered" & Method=="vsearch" & Library=="libC") %>% select(HashID)
exact.unfilt.vs.libD <- mock.all %>% filter(Type=="exact" & Labeler=="unfiltered" & Method=="vsearch" & Library=="libD") %>% select(HashID)
exact.unfilt.da.libA <- mock.all %>% filter(Type=="exact" & Labeler=="unfiltered" & Method=="dada2" & Library=="libA") %>% select(HashID)
exact.unfilt.da.libB <- mock.all %>% filter(Type=="exact" & Labeler=="unfiltered" & Method=="dada2" & Library=="libB") %>% select(HashID)
exact.unfilt.da.libC <- mock.all %>% filter(Type=="exact" & Labeler=="unfiltered" & Method=="dada2" & Library=="libC") %>% select(HashID)
exact.unfilt.da.libD <- mock.all %>% filter(Type=="exact" & Labeler=="unfiltered" & Method=="dada2" & Library=="libD") %>% select(HashID)
exact.unfilt.db.libA <- mock.all %>% filter(Type=="exact" & Labeler=="unfiltered" & Method=="deblur" & Library=="libA") %>% select(HashID)
exact.unfilt.db.libB <- mock.all %>% filter(Type=="exact" & Labeler=="unfiltered" & Method=="deblur" & Library=="libB") %>% select(HashID)
exact.unfilt.db.libC <- mock.all %>% filter(Type=="exact" & Labeler=="unfiltered" & Method=="deblur" & Library=="libC") %>% select(HashID)
exact.unfilt.db.libD <- mock.all %>% filter(Type=="exact" & Labeler=="unfiltered" & Method=="deblur" & Library=="libD") %>% select(HashID)
## Intersections across Filtering methods, by Library
int.exact.unfilt.libA <- intersect(exact.unfilt.da.libA$HashID, intersect(exact.unfilt.db.libA$HashID, exact.unfilt.vs.libA$HashID))
int.exact.unfilt.libB <- intersect(exact.unfilt.da.libB$HashID, intersect(exact.unfilt.db.libB$HashID, exact.unfilt.vs.libB$HashID))
int.exact.unfilt.libC <- intersect(exact.unfilt.da.libC$HashID, intersect(exact.unfilt.db.libC$HashID, exact.unfilt.vs.libC$HashID))
int.exact.unfilt.libD <- intersect(exact.unfilt.da.libD$HashID, intersect(exact.unfilt.db.libD$HashID, exact.unfilt.vs.libD$HashID))
## Intersections across Libraries, by Filtering Method
int.exact.unfilt.da <- intersect(exact.unfilt.da.libA$HashID, intersect(exact.unfilt.da.libB$HashID, intersect(exact.unfilt.da.libC$HashID, exact.unfilt.da.libD$HashID)))
int.exact.unfilt.db <- intersect(exact.unfilt.db.libA$HashID, intersect(exact.unfilt.db.libB$HashID, intersect(exact.unfilt.db.libC$HashID, exact.unfilt.db.libD$HashID)))
int.exact.unfilt.vs <- intersect(exact.unfilt.vs.libA$HashID, intersect(exact.unfilt.vs.libB$HashID, intersect(exact.unfilt.vs.libC$HashID, exact.unfilt.vs.libD$HashID)))

## again for Filtered exact matches
exact.filt.vs.libA <- mock.all %>% filter(Type=="exact" & Labeler=="filtered" & Method=="vsearch" & Library=="libA") %>% select(HashID)
exact.filt.vs.libB <- mock.all %>% filter(Type=="exact" & Labeler=="filtered" & Method=="vsearch" & Library=="libB") %>% select(HashID)
exact.filt.vs.libC <- mock.all %>% filter(Type=="exact" & Labeler=="filtered" & Method=="vsearch" & Library=="libC") %>% select(HashID)
exact.filt.vs.libD <- mock.all %>% filter(Type=="exact" & Labeler=="filtered" & Method=="vsearch" & Library=="libD") %>% select(HashID)
exact.filt.da.libA <- mock.all %>% filter(Type=="exact" & Labeler=="filtered" & Method=="dada2" & Library=="libA") %>% select(HashID)
exact.filt.da.libB <- mock.all %>% filter(Type=="exact" & Labeler=="filtered" & Method=="dada2" & Library=="libB") %>% select(HashID)
exact.filt.da.libC <- mock.all %>% filter(Type=="exact" & Labeler=="filtered" & Method=="dada2" & Library=="libC") %>% select(HashID)
exact.filt.da.libD <- mock.all %>% filter(Type=="exact" & Labeler=="filtered" & Method=="dada2" & Library=="libD") %>% select(HashID)
exact.filt.db.libA <- mock.all %>% filter(Type=="exact" & Labeler=="filtered" & Method=="deblur" & Library=="libA") %>% select(HashID)
exact.filt.db.libB <- mock.all %>% filter(Type=="exact" & Labeler=="filtered" & Method=="deblur" & Library=="libB") %>% select(HashID)
exact.filt.db.libC <- mock.all %>% filter(Type=="exact" & Labeler=="filtered" & Method=="deblur" & Library=="libC") %>% select(HashID)
exact.filt.db.libD <- mock.all %>% filter(Type=="exact" & Labeler=="filtered" & Method=="deblur" & Library=="libD") %>% select(HashID)
## Intersections across Filtering methods, by Library
int.exact.filt.libA <- intersect(exact.filt.da.libA$HashID, intersect(exact.filt.db.libA$HashID, exact.filt.vs.libA$HashID))
int.exact.filt.libB <- intersect(exact.filt.da.libB$HashID, intersect(exact.filt.db.libB$HashID, exact.filt.vs.libB$HashID))
int.exact.filt.libC <- intersect(exact.filt.da.libC$HashID, intersect(exact.filt.db.libC$HashID, exact.filt.vs.libC$HashID))
int.exact.filt.libD <- intersect(exact.filt.da.libD$HashID, intersect(exact.filt.db.libD$HashID, exact.filt.vs.libD$HashID))
## Intersections across Libraries, by Filtering Method
int.exact.filt.da <- intersect(exact.filt.da.libA$HashID, intersect(exact.filt.da.libB$HashID, intersect(exact.filt.da.libC$HashID, exact.filt.da.libD$HashID)))
int.exact.filt.db <- intersect(exact.filt.db.libA$HashID, intersect(exact.filt.db.libB$HashID, intersect(exact.filt.db.libC$HashID, exact.filt.db.libD$HashID)))
int.exact.filt.vs <- intersect(exact.filt.vs.libA$HashID, intersect(exact.filt.vs.libB$HashID, intersect(exact.filt.vs.libC$HashID, exact.filt.vs.libD$HashID)))

### Repeat again for the "partial" values
## First, generate the table:
## convert into matrix of values
filt.partial.mock_Sumry <- mock_Sumry %>% filter(Labeler=="filtered" & Type=="partial")
unfilt.partial.mock_Sumry <- mock_Sumry %>% filter(Labeler=="unfiltered" & Type=="partial")
filt.partial.mock_SumryTable <- dcast(data=filt.partial.mock_Sumry, Library ~ Method, value.var = "CountTypes")
unfilt.partial.mock_SumryTable <- dcast(data=unfilt.partial.mock_Sumry, Library ~ Method, value.var = "CountTypes")

## and generate all the intersections
## Unfiltered
partial.unfilt.vs.libA <- mock.all %>% filter(Type=="partial" & Labeler=="unfiltered" & Method=="vsearch" & Library=="libA") %>% select(HashID)
partial.unfilt.vs.libB <- mock.all %>% filter(Type=="partial" & Labeler=="unfiltered" & Method=="vsearch" & Library=="libB") %>% select(HashID)
partial.unfilt.vs.libC <- mock.all %>% filter(Type=="partial" & Labeler=="unfiltered" & Method=="vsearch" & Library=="libC") %>% select(HashID)
partial.unfilt.vs.libD <- mock.all %>% filter(Type=="partial" & Labeler=="unfiltered" & Method=="vsearch" & Library=="libD") %>% select(HashID)
partial.unfilt.da.libA <- mock.all %>% filter(Type=="partial" & Labeler=="unfiltered" & Method=="dada2" & Library=="libA") %>% select(HashID)
partial.unfilt.da.libB <- mock.all %>% filter(Type=="partial" & Labeler=="unfiltered" & Method=="dada2" & Library=="libB") %>% select(HashID)
partial.unfilt.da.libC <- mock.all %>% filter(Type=="partial" & Labeler=="unfiltered" & Method=="dada2" & Library=="libC") %>% select(HashID)
partial.unfilt.da.libD <- mock.all %>% filter(Type=="partial" & Labeler=="unfiltered" & Method=="dada2" & Library=="libD") %>% select(HashID)
partial.unfilt.db.libA <- mock.all %>% filter(Type=="partial" & Labeler=="unfiltered" & Method=="deblur" & Library=="libA") %>% select(HashID)
partial.unfilt.db.libB <- mock.all %>% filter(Type=="partial" & Labeler=="unfiltered" & Method=="deblur" & Library=="libB") %>% select(HashID)
partial.unfilt.db.libC <- mock.all %>% filter(Type=="partial" & Labeler=="unfiltered" & Method=="deblur" & Library=="libC") %>% select(HashID)
partial.unfilt.db.libD <- mock.all %>% filter(Type=="partial" & Labeler=="unfiltered" & Method=="deblur" & Library=="libD") %>% select(HashID)
## Intersections across Filtering methods, by Library
int.partial.unfilt.libA <- intersect(partial.unfilt.da.libA$HashID, intersect(partial.unfilt.db.libA$HashID, partial.unfilt.vs.libA$HashID))
int.partial.unfilt.libB <- intersect(partial.unfilt.da.libB$HashID, intersect(partial.unfilt.db.libB$HashID, partial.unfilt.vs.libB$HashID))
int.partial.unfilt.libC <- intersect(partial.unfilt.da.libC$HashID, intersect(partial.unfilt.db.libC$HashID, partial.unfilt.vs.libC$HashID))
int.partial.unfilt.libD <- intersect(partial.unfilt.da.libD$HashID, intersect(partial.unfilt.db.libD$HashID, partial.unfilt.vs.libD$HashID))
## Intersections across Libraries, by Filtering Method
int.partial.unfilt.da <- intersect(partial.unfilt.da.libA$HashID, intersect(partial.unfilt.da.libB$HashID, intersect(partial.unfilt.da.libC$HashID, partial.unfilt.da.libD$HashID)))
int.partial.unfilt.db <- intersect(partial.unfilt.db.libA$HashID, intersect(partial.unfilt.db.libB$HashID, intersect(partial.unfilt.db.libC$HashID, partial.unfilt.db.libD$HashID)))
int.partial.unfilt.vs <- intersect(partial.unfilt.vs.libA$HashID, intersect(partial.unfilt.vs.libB$HashID, intersect(partial.unfilt.vs.libC$HashID, partial.unfilt.vs.libD$HashID)))

## again for Filtered partial matches
partial.filt.vs.libA <- mock.all %>% filter(Type=="partial" & Labeler=="filtered" & Method=="vsearch" & Library=="libA") %>% select(HashID)
partial.filt.vs.libB <- mock.all %>% filter(Type=="partial" & Labeler=="filtered" & Method=="vsearch" & Library=="libB") %>% select(HashID)
partial.filt.vs.libC <- mock.all %>% filter(Type=="partial" & Labeler=="filtered" & Method=="vsearch" & Library=="libC") %>% select(HashID)
partial.filt.vs.libD <- mock.all %>% filter(Type=="partial" & Labeler=="filtered" & Method=="vsearch" & Library=="libD") %>% select(HashID)
partial.filt.da.libA <- mock.all %>% filter(Type=="partial" & Labeler=="filtered" & Method=="dada2" & Library=="libA") %>% select(HashID)
partial.filt.da.libB <- mock.all %>% filter(Type=="partial" & Labeler=="filtered" & Method=="dada2" & Library=="libB") %>% select(HashID)
partial.filt.da.libC <- mock.all %>% filter(Type=="partial" & Labeler=="filtered" & Method=="dada2" & Library=="libC") %>% select(HashID)
partial.filt.da.libD <- mock.all %>% filter(Type=="partial" & Labeler=="filtered" & Method=="dada2" & Library=="libD") %>% select(HashID)
partial.filt.db.libA <- mock.all %>% filter(Type=="partial" & Labeler=="filtered" & Method=="deblur" & Library=="libA") %>% select(HashID)
partial.filt.db.libB <- mock.all %>% filter(Type=="partial" & Labeler=="filtered" & Method=="deblur" & Library=="libB") %>% select(HashID)
partial.filt.db.libC <- mock.all %>% filter(Type=="partial" & Labeler=="filtered" & Method=="deblur" & Library=="libC") %>% select(HashID)
partial.filt.db.libD <- mock.all %>% filter(Type=="partial" & Labeler=="filtered" & Method=="deblur" & Library=="libD") %>% select(HashID)
## Intersections across Filtering methods, by Library
int.partial.filt.libA <- intersect(partial.filt.da.libA$HashID, intersect(partial.filt.db.libA$HashID, partial.filt.vs.libA$HashID))
int.partial.filt.libB <- intersect(partial.filt.da.libB$HashID, intersect(partial.filt.db.libB$HashID, partial.filt.vs.libB$HashID))
int.partial.filt.libC <- intersect(partial.filt.da.libC$HashID, intersect(partial.filt.db.libC$HashID, partial.filt.vs.libC$HashID))
int.partial.filt.libD <- intersect(partial.filt.da.libD$HashID, intersect(partial.filt.db.libD$HashID, partial.filt.vs.libD$HashID))
## Intersections across Libraries, by Filtering Method
int.partial.filt.da <- intersect(partial.filt.da.libA$HashID, intersect(partial.filt.da.libB$HashID, intersect(partial.filt.da.libC$HashID, partial.filt.da.libD$HashID)))
int.partial.filt.db <- intersect(partial.filt.db.libA$HashID, intersect(partial.filt.db.libB$HashID, intersect(partial.filt.db.libC$HashID, partial.filt.db.libD$HashID)))
int.partial.filt.vs <- intersect(partial.filt.vs.libA$HashID, intersect(partial.filt.vs.libB$HashID, intersect(partial.filt.vs.libC$HashID, partial.filt.vs.libD$HashID)))


### Repeat again for the "miss" values
## First, generate the table:
## convert into matrix of values
filt.miss.mock_Sumry <- mock_Sumry %>% filter(Labeler=="filtered" & Type=="miss")
unfilt.miss.mock_Sumry <- mock_Sumry %>% filter(Labeler=="unfiltered" & Type=="miss")
filt.miss.mock_SumryTable <- dcast(data=filt.miss.mock_Sumry, Library ~ Method, value.var = "CountTypes")
unfilt.miss.mock_SumryTable <- dcast(data=unfilt.miss.mock_Sumry, Library ~ Method, value.var = "CountTypes")

## and generate all the intersections
## Unfiltered
miss.unfilt.vs.libA <- mock.all %>% filter(Type=="miss" & Labeler=="unfiltered" & Method=="vsearch" & Library=="libA") %>% select(HashID)
miss.unfilt.vs.libB <- mock.all %>% filter(Type=="miss" & Labeler=="unfiltered" & Method=="vsearch" & Library=="libB") %>% select(HashID)
miss.unfilt.vs.libC <- mock.all %>% filter(Type=="miss" & Labeler=="unfiltered" & Method=="vsearch" & Library=="libC") %>% select(HashID)
miss.unfilt.vs.libD <- mock.all %>% filter(Type=="miss" & Labeler=="unfiltered" & Method=="vsearch" & Library=="libD") %>% select(HashID)
miss.unfilt.da.libA <- mock.all %>% filter(Type=="miss" & Labeler=="unfiltered" & Method=="dada2" & Library=="libA") %>% select(HashID)
miss.unfilt.da.libB <- mock.all %>% filter(Type=="miss" & Labeler=="unfiltered" & Method=="dada2" & Library=="libB") %>% select(HashID)
miss.unfilt.da.libC <- mock.all %>% filter(Type=="miss" & Labeler=="unfiltered" & Method=="dada2" & Library=="libC") %>% select(HashID)
miss.unfilt.da.libD <- mock.all %>% filter(Type=="miss" & Labeler=="unfiltered" & Method=="dada2" & Library=="libD") %>% select(HashID)
miss.unfilt.db.libA <- mock.all %>% filter(Type=="miss" & Labeler=="unfiltered" & Method=="deblur" & Library=="libA") %>% select(HashID)
miss.unfilt.db.libB <- mock.all %>% filter(Type=="miss" & Labeler=="unfiltered" & Method=="deblur" & Library=="libB") %>% select(HashID)
miss.unfilt.db.libC <- mock.all %>% filter(Type=="miss" & Labeler=="unfiltered" & Method=="deblur" & Library=="libC") %>% select(HashID)
miss.unfilt.db.libD <- mock.all %>% filter(Type=="miss" & Labeler=="unfiltered" & Method=="deblur" & Library=="libD") %>% select(HashID)
## Intersections across Filtering methods, by Library
int.miss.unfilt.libA <- intersect(miss.unfilt.da.libA$HashID, intersect(miss.unfilt.db.libA$HashID, miss.unfilt.vs.libA$HashID))
int.miss.unfilt.libB <- intersect(miss.unfilt.da.libB$HashID, intersect(miss.unfilt.db.libB$HashID, miss.unfilt.vs.libB$HashID))
int.miss.unfilt.libC <- intersect(miss.unfilt.da.libC$HashID, intersect(miss.unfilt.db.libC$HashID, miss.unfilt.vs.libC$HashID))
int.miss.unfilt.libD <- intersect(miss.unfilt.da.libD$HashID, intersect(miss.unfilt.db.libD$HashID, miss.unfilt.vs.libD$HashID))
## Intersections across Libraries, by Filtering Method
int.miss.unfilt.da <- intersect(miss.unfilt.da.libA$HashID, intersect(miss.unfilt.da.libB$HashID, intersect(miss.unfilt.da.libC$HashID, miss.unfilt.da.libD$HashID)))
int.miss.unfilt.db <- intersect(miss.unfilt.db.libA$HashID, intersect(miss.unfilt.db.libB$HashID, intersect(miss.unfilt.db.libC$HashID, miss.unfilt.db.libD$HashID)))
int.miss.unfilt.vs <- intersect(miss.unfilt.vs.libA$HashID, intersect(miss.unfilt.vs.libB$HashID, intersect(miss.unfilt.vs.libC$HashID, miss.unfilt.vs.libD$HashID)))

## again for Filtered miss matches
miss.filt.vs.libA <- mock.all %>% filter(Type=="miss" & Labeler=="filtered" & Method=="vsearch" & Library=="libA") %>% select(HashID)
miss.filt.vs.libB <- mock.all %>% filter(Type=="miss" & Labeler=="filtered" & Method=="vsearch" & Library=="libB") %>% select(HashID)
miss.filt.vs.libC <- mock.all %>% filter(Type=="miss" & Labeler=="filtered" & Method=="vsearch" & Library=="libC") %>% select(HashID)
miss.filt.vs.libD <- mock.all %>% filter(Type=="miss" & Labeler=="filtered" & Method=="vsearch" & Library=="libD") %>% select(HashID)
miss.filt.da.libA <- mock.all %>% filter(Type=="miss" & Labeler=="filtered" & Method=="dada2" & Library=="libA") %>% select(HashID)
miss.filt.da.libB <- mock.all %>% filter(Type=="miss" & Labeler=="filtered" & Method=="dada2" & Library=="libB") %>% select(HashID)
miss.filt.da.libC <- mock.all %>% filter(Type=="miss" & Labeler=="filtered" & Method=="dada2" & Library=="libC") %>% select(HashID)
miss.filt.da.libD <- mock.all %>% filter(Type=="miss" & Labeler=="filtered" & Method=="dada2" & Library=="libD") %>% select(HashID)
miss.filt.db.libA <- mock.all %>% filter(Type=="miss" & Labeler=="filtered" & Method=="deblur" & Library=="libA") %>% select(HashID)
miss.filt.db.libB <- mock.all %>% filter(Type=="miss" & Labeler=="filtered" & Method=="deblur" & Library=="libB") %>% select(HashID)
miss.filt.db.libC <- mock.all %>% filter(Type=="miss" & Labeler=="filtered" & Method=="deblur" & Library=="libC") %>% select(HashID)
miss.filt.db.libD <- mock.all %>% filter(Type=="miss" & Labeler=="filtered" & Method=="deblur" & Library=="libD") %>% select(HashID)
## Intersections across Filtering methods, by Library
int.miss.filt.libA <- intersect(miss.filt.da.libA$HashID, intersect(miss.filt.db.libA$HashID, miss.filt.vs.libA$HashID))
int.miss.filt.libB <- intersect(miss.filt.da.libB$HashID, intersect(miss.filt.db.libB$HashID, miss.filt.vs.libB$HashID))
int.miss.filt.libC <- intersect(miss.filt.da.libC$HashID, intersect(miss.filt.db.libC$HashID, miss.filt.vs.libC$HashID))
int.miss.filt.libD <- intersect(miss.filt.da.libD$HashID, intersect(miss.filt.db.libD$HashID, miss.filt.vs.libD$HashID))
## Intersections across Libraries, by Filtering Method
int.miss.filt.da <- intersect(miss.filt.da.libA$HashID, intersect(miss.filt.da.libB$HashID, intersect(miss.filt.da.libC$HashID, miss.filt.da.libD$HashID)))
int.miss.filt.db <- intersect(miss.filt.db.libA$HashID, intersect(miss.filt.db.libB$HashID, intersect(miss.filt.db.libC$HashID, miss.filt.db.libD$HashID)))
int.miss.filt.vs <- intersect(miss.filt.vs.libA$HashID, intersect(miss.filt.vs.libB$HashID, intersect(miss.filt.vs.libC$HashID, miss.filt.vs.libD$HashID)))


