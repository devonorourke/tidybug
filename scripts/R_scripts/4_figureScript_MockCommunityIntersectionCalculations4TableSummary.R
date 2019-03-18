# Table used in paper

library(tidyverse)
library(reshape2)

theme_devon <- function () { 
  theme_bw(base_size=12, base_family="Avenir") %+replace% 
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill="transparent", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA)
    )
}

## import data and select mock samples:
df <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz")
mock <- df %>% filter(SampleType == "mock")
mock$Labeler <- paste(mock$HashID, mock$Method, sep="-")
rm(df)
HashFiltLabels <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/HashIDs_withFiltLabels.csv")   ## from '1_sequence_filter.R` script`
mock <- merge(mock, HashFiltLabels)
mock$Labeler <- NULL
rm(HashFiltLabels)

## tiny summary of mock data, per Library, per filtering/not
mock_Sumry <- mock %>% group_by(Library, Method, Filt, MockAlign) %>% summarise(CountTypes=n())

## convert into matrix of values
basic.exact.mock_Sumry <- mock_Sumry %>% filter(Filt=="basic" & MockAlign=="exact")
basic.partial.mock_Sumry <- mock_Sumry %>% filter(Filt=="basic" & MockAlign=="partial")
basic.miss.mock_Sumry <- mock_Sumry %>% filter(Filt=="basic" & MockAlign=="miss")
standard.exact.mock_Sumry <- mock_Sumry %>% filter(Filt=="standard" & MockAlign=="exact")
standard.partial.mock_Sumry <- mock_Sumry %>% filter(Filt=="standard" & MockAlign=="partial")
standard.miss.mock_Sumry <- mock_Sumry %>% filter(Filt=="standard" & MockAlign=="miss")
extra.exact.mock_Sumry <- mock_Sumry %>% filter(Filt=="extra" & MockAlign=="exact")
extra.partial.mock_Sumry <- mock_Sumry %>% filter(Filt=="extra" & MockAlign=="partial")
extra.miss.mock_Sumry <- mock_Sumry %>% filter(Filt=="extra" & MockAlign=="miss") ## does not produce any data (as expected!)

basic.exact.mock_Table <- dcast(basic.exact.mock_Sumry, Library ~ Method, value.var = "CountTypes")
basic.partial.mock_Table <- dcast(basic.partial.mock_Sumry, Library ~ Method, value.var = "CountTypes")
basic.miss.mock_Table <- dcast(basic.miss.mock_Sumry, Library ~ Method, value.var = "CountTypes")
standard.exact.mock_Table <- dcast(standard.exact.mock_Sumry, Library ~ Method, value.var = "CountTypes")
standard.partial.mock_Table <- dcast(standard.partial.mock_Sumry, Library ~ Method, value.var = "CountTypes")
standard.miss.mock_Table <- dcast(standard.miss.mock_Sumry, Library ~ Method, value.var = "CountTypes")
extra.exact.mock_Table <- dcast(extra.exact.mock_Sumry, Library ~ Method, value.var = "CountTypes")
extra.partial.mock_Table <- dcast(extra.partial.mock_Sumry, Library ~ Method, value.var = "CountTypes")
##not run: no data exists...   extra.miss.mock_Table <- dcast(extra.miss.mock_Sumry, Library ~ Method, value.var = "CountTypes")


## we'll look to find how many shared ASVs occur between Libraries (within Filtering algorithm) ..
## and how many shared ASVs occur across Libraries (between Filtering algorithms)
## these values are manually added to the two tables presented above

## to find potential intersections (shared ASVs) by various parameters
## generating vectors of individual exact-match HashIDs, filtered by Method, Library and filtering Labeler
## Only focusing on "basic" exact matches
exact.basic.vs.libA <- mock %>% filter(MockAlign=="exact" & Filt=="basic" & Method=="vsearch" & Library=="libA") %>% select(HashID)
exact.basic.vs.libB <- mock %>% filter(MockAlign=="exact" & Filt=="basic" & Method=="vsearch" & Library=="libB") %>% select(HashID)
exact.basic.vs.libC <- mock %>% filter(MockAlign=="exact" & Filt=="basic" & Method=="vsearch" & Library=="libC") %>% select(HashID)
exact.basic.vs.libD <- mock %>% filter(MockAlign=="exact" & Filt=="basic" & Method=="vsearch" & Library=="libD") %>% select(HashID)
exact.basic.da.libA <- mock %>% filter(MockAlign=="exact" & Filt=="basic" & Method=="dada2" & Library=="libA") %>% select(HashID)
exact.basic.da.libB <- mock %>% filter(MockAlign=="exact" & Filt=="basic" & Method=="dada2" & Library=="libB") %>% select(HashID)
exact.basic.da.libC <- mock %>% filter(MockAlign=="exact" & Filt=="basic" & Method=="dada2" & Library=="libC") %>% select(HashID)
exact.basic.da.libD <- mock %>% filter(MockAlign=="exact" & Filt=="basic" & Method=="dada2" & Library=="libD") %>% select(HashID)
exact.basic.db.libA <- mock %>% filter(MockAlign=="exact" & Filt=="basic" & Method=="deblur" & Library=="libA") %>% select(HashID)
exact.basic.db.libB <- mock %>% filter(MockAlign=="exact" & Filt=="basic" & Method=="deblur" & Library=="libB") %>% select(HashID)
exact.basic.db.libC <- mock %>% filter(MockAlign=="exact" & Filt=="basic" & Method=="deblur" & Library=="libC") %>% select(HashID)
exact.basic.db.libD <- mock %>% filter(MockAlign=="exact" & Filt=="basic" & Method=="deblur" & Library=="libD") %>% select(HashID)
## Intersections across Filtering methods, by Library
int.exact.basic.libA <- length(intersect(exact.basic.da.libA$HashID, intersect(exact.basic.db.libA$HashID, exact.basic.vs.libA$HashID)))
int.exact.basic.libB <- length(intersect(exact.basic.da.libB$HashID, intersect(exact.basic.db.libB$HashID, exact.basic.vs.libB$HashID)))
int.exact.basic.libC <- length(intersect(exact.basic.da.libC$HashID, intersect(exact.basic.db.libC$HashID, exact.basic.vs.libC$HashID)))
int.exact.basic.libD <- length(intersect(exact.basic.da.libD$HashID, intersect(exact.basic.db.libD$HashID, exact.basic.vs.libD$HashID)))
int.exact.basic.meth.vals <- c(int.exact.basic.libA, int.exact.basic.libB, int.exact.basic.libC, int.exact.basic.libD)
int.exact.basic.meth.names <- c('int.exact.basic.libA', 'int.exact.basic.libB', 'int.exact.basic.libC', 'int.exact.basic.libD')
int.exact.basic.meth.df <- data.frame(int.exact.basic.meth.vals, int.exact.basic.meth.names)
colnames(int.exact.basic.meth.df) <- c("intersections", "name")

## Intersections across Libraries, by Filtering Method
int.exact.basic.da <- length(intersect(exact.basic.da.libA$HashID, intersect(exact.basic.da.libB$HashID, intersect(exact.basic.da.libC$HashID, exact.basic.da.libD$HashID))))
int.exact.basic.db <- length(intersect(exact.basic.db.libA$HashID, intersect(exact.basic.db.libB$HashID, intersect(exact.basic.db.libC$HashID, exact.basic.db.libD$HashID))))
int.exact.basic.vs <- length(intersect(exact.basic.vs.libA$HashID, intersect(exact.basic.vs.libB$HashID, intersect(exact.basic.vs.libC$HashID, exact.basic.vs.libD$HashID))))
int.exact.basic.lib.vals <- c(int.exact.basic.da, int.exact.basic.db, int.exact.basic.vs)
int.exact.basic.lib.names <- c('int.exact.basic.da', 'int.exact.basic.db', 'int.exact.basic.vs')
int.exact.basic.lib.df <- data.frame(int.exact.basic.lib.vals, int.exact.basic.lib.names)
colnames(int.exact.basic.lib.df) <- c("intersections", "name")

rm(exact.basic.vs.libA, exact.basic.vs.libB, exact.basic.vs.libC, exact.basic.vs.libD, exact.basic.da.libA, exact.basic.da.libB, exact.basic.da.libC, exact.basic.da.libD,
   exact.basic.db.libA, exact.basic.db.libB, exact.basic.db.libC, exact.basic.db.libD, int.exact.basic.libA, int.exact.basic.libB, int.exact.basic.libC, int.exact.basic.libD,
   int.exact.basic.meth.vals, int.exact.basic.meth.names, int.exact.basic.lib.vals, int.exact.basic.lib.names)

## Again, for on "standard" exact matches
exact.standard.vs.libA <- mock %>% filter(MockAlign=="exact" & Filt=="standard" & Method=="vsearch" & Library=="libA") %>% select(HashID)
exact.standard.vs.libB <- mock %>% filter(MockAlign=="exact" & Filt=="standard" & Method=="vsearch" & Library=="libB") %>% select(HashID)
exact.standard.vs.libC <- mock %>% filter(MockAlign=="exact" & Filt=="standard" & Method=="vsearch" & Library=="libC") %>% select(HashID)
exact.standard.vs.libD <- mock %>% filter(MockAlign=="exact" & Filt=="standard" & Method=="vsearch" & Library=="libD") %>% select(HashID)
exact.standard.da.libA <- mock %>% filter(MockAlign=="exact" & Filt=="standard" & Method=="dada2" & Library=="libA") %>% select(HashID)
exact.standard.da.libB <- mock %>% filter(MockAlign=="exact" & Filt=="standard" & Method=="dada2" & Library=="libB") %>% select(HashID)
exact.standard.da.libC <- mock %>% filter(MockAlign=="exact" & Filt=="standard" & Method=="dada2" & Library=="libC") %>% select(HashID)
exact.standard.da.libD <- mock %>% filter(MockAlign=="exact" & Filt=="standard" & Method=="dada2" & Library=="libD") %>% select(HashID)
exact.standard.db.libA <- mock %>% filter(MockAlign=="exact" & Filt=="standard" & Method=="deblur" & Library=="libA") %>% select(HashID)
exact.standard.db.libB <- mock %>% filter(MockAlign=="exact" & Filt=="standard" & Method=="deblur" & Library=="libB") %>% select(HashID)
exact.standard.db.libC <- mock %>% filter(MockAlign=="exact" & Filt=="standard" & Method=="deblur" & Library=="libC") %>% select(HashID)
exact.standard.db.libD <- mock %>% filter(MockAlign=="exact" & Filt=="standard" & Method=="deblur" & Library=="libD") %>% select(HashID)
## Intersections across Filtering methods, by Library
int.exact.standard.libA <- length(intersect(exact.standard.da.libA$HashID, intersect(exact.standard.db.libA$HashID, exact.standard.vs.libA$HashID)))
int.exact.standard.libB <- length(intersect(exact.standard.da.libB$HashID, intersect(exact.standard.db.libB$HashID, exact.standard.vs.libB$HashID)))
int.exact.standard.libC <- length(intersect(exact.standard.da.libC$HashID, intersect(exact.standard.db.libC$HashID, exact.standard.vs.libC$HashID)))
int.exact.standard.libD <- length(intersect(exact.standard.da.libD$HashID, intersect(exact.standard.db.libD$HashID, exact.standard.vs.libD$HashID)))
int.exact.standard.meth.vals <- c(int.exact.standard.libA, int.exact.standard.libB, int.exact.standard.libC, int.exact.standard.libD)
int.exact.standard.meth.names <- c('int.exact.standard.libA', 'int.exact.standard.libB', 'int.exact.standard.libC', 'int.exact.standard.libD')
int.exact.standard.meth.df <- data.frame(int.exact.standard.meth.vals, int.exact.standard.meth.names)
colnames(int.exact.standard.meth.df) <- c("intersections", "name")

## Intersections across Libraries, by Filtering Method
int.exact.standard.da <- length(intersect(exact.standard.da.libA$HashID, intersect(exact.standard.da.libB$HashID, intersect(exact.standard.da.libC$HashID, exact.standard.da.libD$HashID))))
int.exact.standard.db <- length(intersect(exact.standard.db.libA$HashID, intersect(exact.standard.db.libB$HashID, intersect(exact.standard.db.libC$HashID, exact.standard.db.libD$HashID))))
int.exact.standard.vs <- length(intersect(exact.standard.vs.libA$HashID, intersect(exact.standard.vs.libB$HashID, intersect(exact.standard.vs.libC$HashID, exact.standard.vs.libD$HashID))))
int.exact.standard.lib.vals <- c(int.exact.standard.da, int.exact.standard.db, int.exact.standard.vs)
int.exact.standard.lib.names <- c('int.exact.standard.da', 'int.exact.standard.db', 'int.exact.standard.vs')
int.exact.standard.lib.df <- data.frame(int.exact.standard.lib.vals, int.exact.standard.lib.names)
colnames(int.exact.standard.lib.df) <- c("intersections", "name")

rm(exact.standard.vs.libA, exact.standard.vs.libB, exact.standard.vs.libC, exact.standard.vs.libD, exact.standard.da.libA, exact.standard.da.libB, exact.standard.da.libC, exact.standard.da.libD,
   exact.standard.db.libA, exact.standard.db.libB, exact.standard.db.libC, exact.standard.db.libD, int.exact.standard.libA, int.exact.standard.libB, int.exact.standard.libC, int.exact.standard.libD,
   int.exact.standard.meth.vals, int.exact.standard.meth.names, int.exact.standard.lib.vals, int.exact.standard.lib.names)


## Again, for on "extra" exact matches
exact.extra.vs.libA <- mock %>% filter(MockAlign=="exact" & Filt=="extra" & Method=="vsearch" & Library=="libA") %>% select(HashID)
exact.extra.vs.libB <- mock %>% filter(MockAlign=="exact" & Filt=="extra" & Method=="vsearch" & Library=="libB") %>% select(HashID)
exact.extra.vs.libC <- mock %>% filter(MockAlign=="exact" & Filt=="extra" & Method=="vsearch" & Library=="libC") %>% select(HashID)
exact.extra.vs.libD <- mock %>% filter(MockAlign=="exact" & Filt=="extra" & Method=="vsearch" & Library=="libD") %>% select(HashID)
exact.extra.da.libA <- mock %>% filter(MockAlign=="exact" & Filt=="extra" & Method=="dada2" & Library=="libA") %>% select(HashID)
exact.extra.da.libB <- mock %>% filter(MockAlign=="exact" & Filt=="extra" & Method=="dada2" & Library=="libB") %>% select(HashID)
exact.extra.da.libC <- mock %>% filter(MockAlign=="exact" & Filt=="extra" & Method=="dada2" & Library=="libC") %>% select(HashID)
exact.extra.da.libD <- mock %>% filter(MockAlign=="exact" & Filt=="extra" & Method=="dada2" & Library=="libD") %>% select(HashID)
exact.extra.db.libA <- mock %>% filter(MockAlign=="exact" & Filt=="extra" & Method=="deblur" & Library=="libA") %>% select(HashID)
exact.extra.db.libB <- mock %>% filter(MockAlign=="exact" & Filt=="extra" & Method=="deblur" & Library=="libB") %>% select(HashID)
exact.extra.db.libC <- mock %>% filter(MockAlign=="exact" & Filt=="extra" & Method=="deblur" & Library=="libC") %>% select(HashID)
exact.extra.db.libD <- mock %>% filter(MockAlign=="exact" & Filt=="extra" & Method=="deblur" & Library=="libD") %>% select(HashID)
## Intersections across Filtering methods, by Library
int.exact.extra.libA <- length(intersect(exact.extra.da.libA$HashID, intersect(exact.extra.db.libA$HashID, exact.extra.vs.libA$HashID)))
int.exact.extra.libB <- length(intersect(exact.extra.da.libB$HashID, intersect(exact.extra.db.libB$HashID, exact.extra.vs.libB$HashID)))
int.exact.extra.libC <- length(intersect(exact.extra.da.libC$HashID, intersect(exact.extra.db.libC$HashID, exact.extra.vs.libC$HashID)))
int.exact.extra.libD <- length(intersect(exact.extra.da.libD$HashID, intersect(exact.extra.db.libD$HashID, exact.extra.vs.libD$HashID)))
int.exact.extra.meth.vals <- c(int.exact.extra.libA, int.exact.extra.libB, int.exact.extra.libC, int.exact.extra.libD)
int.exact.extra.meth.names <- c('int.exact.extra.libA', 'int.exact.extra.libB', 'int.exact.extra.libC', 'int.exact.extra.libD')
int.exact.extra.meth.df <- data.frame(int.exact.extra.meth.vals, int.exact.extra.meth.names)
colnames(int.exact.extra.meth.df) <- c("intersections", "name")

## Intersections across Libraries, by Filtering Method
int.exact.extra.da <- length(intersect(exact.extra.da.libA$HashID, intersect(exact.extra.da.libB$HashID, intersect(exact.extra.da.libC$HashID, exact.extra.da.libD$HashID))))
int.exact.extra.db <- length(intersect(exact.extra.db.libA$HashID, intersect(exact.extra.db.libB$HashID, intersect(exact.extra.db.libC$HashID, exact.extra.db.libD$HashID))))
int.exact.extra.vs <- length(intersect(exact.extra.vs.libA$HashID, intersect(exact.extra.vs.libB$HashID, intersect(exact.extra.vs.libC$HashID, exact.extra.vs.libD$HashID))))
int.exact.extra.lib.vals <- c(int.exact.extra.da, int.exact.extra.db, int.exact.extra.vs)
int.exact.extra.lib.names <- c('int.exact.extra.da', 'int.exact.extra.db', 'int.exact.extra.vs')
int.exact.extra.lib.df <- data.frame(int.exact.extra.lib.vals, int.exact.extra.lib.names)
colnames(int.exact.extra.lib.df) <- c("intersections", "name")

rm(exact.extra.vs.libA, exact.extra.vs.libB, exact.extra.vs.libC, exact.extra.vs.libD, exact.extra.da.libA, exact.extra.da.libB, exact.extra.da.libC, exact.extra.da.libD,
   exact.extra.db.libA, exact.extra.db.libB, exact.extra.db.libC, exact.extra.db.libD, int.exact.extra.libA, int.exact.extra.libB, int.exact.extra.libC, int.exact.extra.libD,
   int.exact.extra.meth.vals, int.exact.extra.meth.names, int.exact.extra.lib.vals, int.exact.extra.lib.names)


## Now focusing on "basic" & partial matches
partial.basic.vs.libA <- mock %>% filter(MockAlign=="partial" & Filt=="basic" & Method=="vsearch" & Library=="libA") %>% select(HashID)
partial.basic.vs.libB <- mock %>% filter(MockAlign=="partial" & Filt=="basic" & Method=="vsearch" & Library=="libB") %>% select(HashID)
partial.basic.vs.libC <- mock %>% filter(MockAlign=="partial" & Filt=="basic" & Method=="vsearch" & Library=="libC") %>% select(HashID)
partial.basic.vs.libD <- mock %>% filter(MockAlign=="partial" & Filt=="basic" & Method=="vsearch" & Library=="libD") %>% select(HashID)
partial.basic.da.libA <- mock %>% filter(MockAlign=="partial" & Filt=="basic" & Method=="dada2" & Library=="libA") %>% select(HashID)
partial.basic.da.libB <- mock %>% filter(MockAlign=="partial" & Filt=="basic" & Method=="dada2" & Library=="libB") %>% select(HashID)
partial.basic.da.libC <- mock %>% filter(MockAlign=="partial" & Filt=="basic" & Method=="dada2" & Library=="libC") %>% select(HashID)
partial.basic.da.libD <- mock %>% filter(MockAlign=="partial" & Filt=="basic" & Method=="dada2" & Library=="libD") %>% select(HashID)
partial.basic.db.libA <- mock %>% filter(MockAlign=="partial" & Filt=="basic" & Method=="deblur" & Library=="libA") %>% select(HashID)
partial.basic.db.libB <- mock %>% filter(MockAlign=="partial" & Filt=="basic" & Method=="deblur" & Library=="libB") %>% select(HashID)
partial.basic.db.libC <- mock %>% filter(MockAlign=="partial" & Filt=="basic" & Method=="deblur" & Library=="libC") %>% select(HashID)
partial.basic.db.libD <- mock %>% filter(MockAlign=="partial" & Filt=="basic" & Method=="deblur" & Library=="libD") %>% select(HashID)
## Intersections across Filtering methods, by Library
int.partial.basic.libA <- length(intersect(partial.basic.da.libA$HashID, intersect(partial.basic.db.libA$HashID, partial.basic.vs.libA$HashID)))
int.partial.basic.libB <- length(intersect(partial.basic.da.libB$HashID, intersect(partial.basic.db.libB$HashID, partial.basic.vs.libB$HashID)))
int.partial.basic.libC <- length(intersect(partial.basic.da.libC$HashID, intersect(partial.basic.db.libC$HashID, partial.basic.vs.libC$HashID)))
int.partial.basic.libD <- length(intersect(partial.basic.da.libD$HashID, intersect(partial.basic.db.libD$HashID, partial.basic.vs.libD$HashID)))
int.partial.basic.meth.vals <- c(int.partial.basic.libA, int.partial.basic.libB, int.partial.basic.libC, int.partial.basic.libD)
int.partial.basic.meth.names <- c('int.partial.basic.libA', 'int.partial.basic.libB', 'int.partial.basic.libC', 'int.partial.basic.libD')
int.partial.basic.meth.df <- data.frame(int.partial.basic.meth.vals, int.partial.basic.meth.names)
colnames(int.partial.basic.meth.df) <- c("intersections", "name")

## Intersections across Libraries, by Filtering Method
int.partial.basic.da <- length(intersect(partial.basic.da.libA$HashID, intersect(partial.basic.da.libB$HashID, intersect(partial.basic.da.libC$HashID, partial.basic.da.libD$HashID))))
int.partial.basic.db <- length(intersect(partial.basic.db.libA$HashID, intersect(partial.basic.db.libB$HashID, intersect(partial.basic.db.libC$HashID, partial.basic.db.libD$HashID))))
int.partial.basic.vs <- length(intersect(partial.basic.vs.libA$HashID, intersect(partial.basic.vs.libB$HashID, intersect(partial.basic.vs.libC$HashID, partial.basic.vs.libD$HashID))))
int.partial.basic.lib.vals <- c(int.partial.basic.da, int.partial.basic.db, int.partial.basic.vs)
int.partial.basic.lib.names <- c('int.partial.basic.da', 'int.partial.basic.db', 'int.partial.basic.vs')
int.partial.basic.lib.df <- data.frame(int.partial.basic.lib.vals, int.partial.basic.lib.names)
colnames(int.partial.basic.lib.df) <- c("intersections", "name")

rm(partial.basic.vs.libA, partial.basic.vs.libB, partial.basic.vs.libC, partial.basic.vs.libD, partial.basic.da.libA, partial.basic.da.libB, partial.basic.da.libC, partial.basic.da.libD,
   partial.basic.db.libA, partial.basic.db.libB, partial.basic.db.libC, partial.basic.db.libD, int.partial.basic.libA, int.partial.basic.libB, int.partial.basic.libC, int.partial.basic.libD,
   int.partial.basic.meth.vals, int.partial.basic.meth.names, int.partial.basic.lib.vals, int.partial.basic.lib.names)

## Again, for on "standard" partial matches
partial.standard.vs.libA <- mock %>% filter(MockAlign=="partial" & Filt=="standard" & Method=="vsearch" & Library=="libA") %>% select(HashID)
partial.standard.vs.libB <- mock %>% filter(MockAlign=="partial" & Filt=="standard" & Method=="vsearch" & Library=="libB") %>% select(HashID)
partial.standard.vs.libC <- mock %>% filter(MockAlign=="partial" & Filt=="standard" & Method=="vsearch" & Library=="libC") %>% select(HashID)
partial.standard.vs.libD <- mock %>% filter(MockAlign=="partial" & Filt=="standard" & Method=="vsearch" & Library=="libD") %>% select(HashID)
partial.standard.da.libA <- mock %>% filter(MockAlign=="partial" & Filt=="standard" & Method=="dada2" & Library=="libA") %>% select(HashID)
partial.standard.da.libB <- mock %>% filter(MockAlign=="partial" & Filt=="standard" & Method=="dada2" & Library=="libB") %>% select(HashID)
partial.standard.da.libC <- mock %>% filter(MockAlign=="partial" & Filt=="standard" & Method=="dada2" & Library=="libC") %>% select(HashID)
partial.standard.da.libD <- mock %>% filter(MockAlign=="partial" & Filt=="standard" & Method=="dada2" & Library=="libD") %>% select(HashID)
partial.standard.db.libA <- mock %>% filter(MockAlign=="partial" & Filt=="standard" & Method=="deblur" & Library=="libA") %>% select(HashID)
partial.standard.db.libB <- mock %>% filter(MockAlign=="partial" & Filt=="standard" & Method=="deblur" & Library=="libB") %>% select(HashID)
partial.standard.db.libC <- mock %>% filter(MockAlign=="partial" & Filt=="standard" & Method=="deblur" & Library=="libC") %>% select(HashID)
partial.standard.db.libD <- mock %>% filter(MockAlign=="partial" & Filt=="standard" & Method=="deblur" & Library=="libD") %>% select(HashID)
## Intersections across Filtering methods, by Library
int.partial.standard.libA <- length(intersect(partial.standard.da.libA$HashID, intersect(partial.standard.db.libA$HashID, partial.standard.vs.libA$HashID)))
int.partial.standard.libB <- length(intersect(partial.standard.da.libB$HashID, intersect(partial.standard.db.libB$HashID, partial.standard.vs.libB$HashID)))
int.partial.standard.libC <- length(intersect(partial.standard.da.libC$HashID, intersect(partial.standard.db.libC$HashID, partial.standard.vs.libC$HashID)))
int.partial.standard.libD <- length(intersect(partial.standard.da.libD$HashID, intersect(partial.standard.db.libD$HashID, partial.standard.vs.libD$HashID)))
int.partial.standard.meth.vals <- c(int.partial.standard.libA, int.partial.standard.libB, int.partial.standard.libC, int.partial.standard.libD)
int.partial.standard.meth.names <- c('int.partial.standard.libA', 'int.partial.standard.libB', 'int.partial.standard.libC', 'int.partial.standard.libD')
int.partial.standard.meth.df <- data.frame(int.partial.standard.meth.vals, int.partial.standard.meth.names)
colnames(int.partial.standard.meth.df) <- c("intersections", "name")

## Intersections across Libraries, by Filtering Method
int.partial.standard.da <- length(intersect(partial.standard.da.libA$HashID, intersect(partial.standard.da.libB$HashID, intersect(partial.standard.da.libC$HashID, partial.standard.da.libD$HashID))))
int.partial.standard.db <- length(intersect(partial.standard.db.libA$HashID, intersect(partial.standard.db.libB$HashID, intersect(partial.standard.db.libC$HashID, partial.standard.db.libD$HashID))))
int.partial.standard.vs <- length(intersect(partial.standard.vs.libA$HashID, intersect(partial.standard.vs.libB$HashID, intersect(partial.standard.vs.libC$HashID, partial.standard.vs.libD$HashID))))
int.partial.standard.lib.vals <- c(int.partial.standard.da, int.partial.standard.db, int.partial.standard.vs)
int.partial.standard.lib.names <- c('int.partial.standard.da', 'int.partial.standard.db', 'int.partial.standard.vs')
int.partial.standard.lib.df <- data.frame(int.partial.standard.lib.vals, int.partial.standard.lib.names)
colnames(int.partial.standard.lib.df) <- c("intersections", "name")

rm(partial.standard.vs.libA, partial.standard.vs.libB, partial.standard.vs.libC, partial.standard.vs.libD, partial.standard.da.libA, partial.standard.da.libB, partial.standard.da.libC, partial.standard.da.libD,
   partial.standard.db.libA, partial.standard.db.libB, partial.standard.db.libC, partial.standard.db.libD, int.partial.standard.libA, int.partial.standard.libB, int.partial.standard.libC, int.partial.standard.libD,
   int.partial.standard.meth.vals, int.partial.standard.meth.names, int.partial.standard.lib.vals, int.partial.standard.lib.names)


## Again, for on "extra" partial matches
partial.extra.vs.libA <- mock %>% filter(MockAlign=="partial" & Filt=="extra" & Method=="vsearch" & Library=="libA") %>% select(HashID)
partial.extra.vs.libB <- mock %>% filter(MockAlign=="partial" & Filt=="extra" & Method=="vsearch" & Library=="libB") %>% select(HashID)
partial.extra.vs.libC <- mock %>% filter(MockAlign=="partial" & Filt=="extra" & Method=="vsearch" & Library=="libC") %>% select(HashID)
partial.extra.vs.libD <- mock %>% filter(MockAlign=="partial" & Filt=="extra" & Method=="vsearch" & Library=="libD") %>% select(HashID)
partial.extra.da.libA <- mock %>% filter(MockAlign=="partial" & Filt=="extra" & Method=="dada2" & Library=="libA") %>% select(HashID)
partial.extra.da.libB <- mock %>% filter(MockAlign=="partial" & Filt=="extra" & Method=="dada2" & Library=="libB") %>% select(HashID)
partial.extra.da.libC <- mock %>% filter(MockAlign=="partial" & Filt=="extra" & Method=="dada2" & Library=="libC") %>% select(HashID)
partial.extra.da.libD <- mock %>% filter(MockAlign=="partial" & Filt=="extra" & Method=="dada2" & Library=="libD") %>% select(HashID)
partial.extra.db.libA <- mock %>% filter(MockAlign=="partial" & Filt=="extra" & Method=="deblur" & Library=="libA") %>% select(HashID)
partial.extra.db.libB <- mock %>% filter(MockAlign=="partial" & Filt=="extra" & Method=="deblur" & Library=="libB") %>% select(HashID)
partial.extra.db.libC <- mock %>% filter(MockAlign=="partial" & Filt=="extra" & Method=="deblur" & Library=="libC") %>% select(HashID)
partial.extra.db.libD <- mock %>% filter(MockAlign=="partial" & Filt=="extra" & Method=="deblur" & Library=="libD") %>% select(HashID)
## Intersections across Filtering methods, by Library
int.partial.extra.libA <- length(intersect(partial.extra.da.libA$HashID, intersect(partial.extra.db.libA$HashID, partial.extra.vs.libA$HashID)))
int.partial.extra.libB <- length(intersect(partial.extra.da.libB$HashID, intersect(partial.extra.db.libB$HashID, partial.extra.vs.libB$HashID)))
int.partial.extra.libC <- length(intersect(partial.extra.da.libC$HashID, intersect(partial.extra.db.libC$HashID, partial.extra.vs.libC$HashID)))
int.partial.extra.libD <- length(intersect(partial.extra.da.libD$HashID, intersect(partial.extra.db.libD$HashID, partial.extra.vs.libD$HashID)))
int.partial.extra.meth.vals <- c(int.partial.extra.libA, int.partial.extra.libB, int.partial.extra.libC, int.partial.extra.libD)
int.partial.extra.meth.names <- c('int.partial.extra.libA', 'int.partial.extra.libB', 'int.partial.extra.libC', 'int.partial.extra.libD')
int.partial.extra.meth.df <- data.frame(int.partial.extra.meth.vals, int.partial.extra.meth.names)
colnames(int.partial.extra.meth.df) <- c("intersections", "name")

## Intersections across Libraries, by Filtering Method
int.partial.extra.da <- length(intersect(partial.extra.da.libA$HashID, intersect(partial.extra.da.libB$HashID, intersect(partial.extra.da.libC$HashID, partial.extra.da.libD$HashID))))
int.partial.extra.db <- length(intersect(partial.extra.db.libA$HashID, intersect(partial.extra.db.libB$HashID, intersect(partial.extra.db.libC$HashID, partial.extra.db.libD$HashID))))
int.partial.extra.vs <- length(intersect(partial.extra.vs.libA$HashID, intersect(partial.extra.vs.libB$HashID, intersect(partial.extra.vs.libC$HashID, partial.extra.vs.libD$HashID))))
int.partial.extra.lib.vals <- c(int.partial.extra.da, int.partial.extra.db, int.partial.extra.vs)
int.partial.extra.lib.names <- c('int.partial.extra.da', 'int.partial.extra.db', 'int.partial.extra.vs')
int.partial.extra.lib.df <- data.frame(int.partial.extra.lib.vals, int.partial.extra.lib.names)
colnames(int.partial.extra.lib.df) <- c("intersections", "name")

rm(partial.extra.vs.libA, partial.extra.vs.libB, partial.extra.vs.libC, partial.extra.vs.libD, partial.extra.da.libA, partial.extra.da.libB, partial.extra.da.libC, partial.extra.da.libD,
   partial.extra.db.libA, partial.extra.db.libB, partial.extra.db.libC, partial.extra.db.libD, int.partial.extra.libA, int.partial.extra.libB, int.partial.extra.libC, int.partial.extra.libD,
   int.partial.extra.meth.vals, int.partial.extra.meth.names, int.partial.extra.lib.vals, int.partial.extra.lib.names)


## Now focusing on "basic" & miss matches
miss.basic.vs.libA <- mock %>% filter(MockAlign=="miss" & Filt=="basic" & Method=="vsearch" & Library=="libA") %>% select(HashID)
miss.basic.vs.libB <- mock %>% filter(MockAlign=="miss" & Filt=="basic" & Method=="vsearch" & Library=="libB") %>% select(HashID)
miss.basic.vs.libC <- mock %>% filter(MockAlign=="miss" & Filt=="basic" & Method=="vsearch" & Library=="libC") %>% select(HashID)
miss.basic.vs.libD <- mock %>% filter(MockAlign=="miss" & Filt=="basic" & Method=="vsearch" & Library=="libD") %>% select(HashID)
miss.basic.da.libA <- mock %>% filter(MockAlign=="miss" & Filt=="basic" & Method=="dada2" & Library=="libA") %>% select(HashID)
miss.basic.da.libB <- mock %>% filter(MockAlign=="miss" & Filt=="basic" & Method=="dada2" & Library=="libB") %>% select(HashID)
miss.basic.da.libC <- mock %>% filter(MockAlign=="miss" & Filt=="basic" & Method=="dada2" & Library=="libC") %>% select(HashID)
miss.basic.da.libD <- mock %>% filter(MockAlign=="miss" & Filt=="basic" & Method=="dada2" & Library=="libD") %>% select(HashID)
miss.basic.db.libA <- mock %>% filter(MockAlign=="miss" & Filt=="basic" & Method=="deblur" & Library=="libA") %>% select(HashID)
miss.basic.db.libB <- mock %>% filter(MockAlign=="miss" & Filt=="basic" & Method=="deblur" & Library=="libB") %>% select(HashID)
miss.basic.db.libC <- mock %>% filter(MockAlign=="miss" & Filt=="basic" & Method=="deblur" & Library=="libC") %>% select(HashID)
miss.basic.db.libD <- mock %>% filter(MockAlign=="miss" & Filt=="basic" & Method=="deblur" & Library=="libD") %>% select(HashID)
## Intersections across Filtering methods, by Library
int.miss.basic.libA <- length(intersect(miss.basic.da.libA$HashID, intersect(miss.basic.db.libA$HashID, miss.basic.vs.libA$HashID)))
int.miss.basic.libB <- length(intersect(miss.basic.da.libB$HashID, intersect(miss.basic.db.libB$HashID, miss.basic.vs.libB$HashID)))
int.miss.basic.libC <- length(intersect(miss.basic.da.libC$HashID, intersect(miss.basic.db.libC$HashID, miss.basic.vs.libC$HashID)))
int.miss.basic.libD <- length(intersect(miss.basic.da.libD$HashID, intersect(miss.basic.db.libD$HashID, miss.basic.vs.libD$HashID)))
int.miss.basic.meth.vals <- c(int.miss.basic.libA, int.miss.basic.libB, int.miss.basic.libC, int.miss.basic.libD)
int.miss.basic.meth.names <- c('int.miss.basic.libA', 'int.miss.basic.libB', 'int.miss.basic.libC', 'int.miss.basic.libD')
int.miss.basic.meth.df <- data.frame(int.miss.basic.meth.vals, int.miss.basic.meth.names)
colnames(int.miss.basic.meth.df) <- c("intersections", "name")

## Intersections across Libraries, by Filtering Method
int.miss.basic.da <- length(intersect(miss.basic.da.libA$HashID, intersect(miss.basic.da.libB$HashID, intersect(miss.basic.da.libC$HashID, miss.basic.da.libD$HashID))))
int.miss.basic.db <- length(intersect(miss.basic.db.libA$HashID, intersect(miss.basic.db.libB$HashID, intersect(miss.basic.db.libC$HashID, miss.basic.db.libD$HashID))))
int.miss.basic.vs <- length(intersect(miss.basic.vs.libA$HashID, intersect(miss.basic.vs.libB$HashID, intersect(miss.basic.vs.libC$HashID, miss.basic.vs.libD$HashID))))
int.miss.basic.lib.vals <- c(int.miss.basic.da, int.miss.basic.db, int.miss.basic.vs)
int.miss.basic.lib.names <- c('int.miss.basic.da', 'int.miss.basic.db', 'int.miss.basic.vs')
int.miss.basic.lib.df <- data.frame(int.miss.basic.lib.vals, int.miss.basic.lib.names)
colnames(int.miss.basic.lib.df) <- c("intersections", "name")

rm(miss.basic.vs.libA, miss.basic.vs.libB, miss.basic.vs.libC, miss.basic.vs.libD, miss.basic.da.libA, miss.basic.da.libB, miss.basic.da.libC, miss.basic.da.libD,
   miss.basic.db.libA, miss.basic.db.libB, miss.basic.db.libC, miss.basic.db.libD, int.miss.basic.libA, int.miss.basic.libB, int.miss.basic.libC, int.miss.basic.libD,
   int.miss.basic.meth.vals, int.miss.basic.meth.names, int.miss.basic.lib.vals, int.miss.basic.lib.names)

## Again, for on "standard" miss matches
miss.standard.vs.libA <- mock %>% filter(MockAlign=="miss" & Filt=="standard" & Method=="vsearch" & Library=="libA") %>% select(HashID)
miss.standard.vs.libB <- mock %>% filter(MockAlign=="miss" & Filt=="standard" & Method=="vsearch" & Library=="libB") %>% select(HashID)
miss.standard.vs.libC <- mock %>% filter(MockAlign=="miss" & Filt=="standard" & Method=="vsearch" & Library=="libC") %>% select(HashID)
miss.standard.vs.libD <- mock %>% filter(MockAlign=="miss" & Filt=="standard" & Method=="vsearch" & Library=="libD") %>% select(HashID)
miss.standard.da.libA <- mock %>% filter(MockAlign=="miss" & Filt=="standard" & Method=="dada2" & Library=="libA") %>% select(HashID)
miss.standard.da.libB <- mock %>% filter(MockAlign=="miss" & Filt=="standard" & Method=="dada2" & Library=="libB") %>% select(HashID)
miss.standard.da.libC <- mock %>% filter(MockAlign=="miss" & Filt=="standard" & Method=="dada2" & Library=="libC") %>% select(HashID)
miss.standard.da.libD <- mock %>% filter(MockAlign=="miss" & Filt=="standard" & Method=="dada2" & Library=="libD") %>% select(HashID)
miss.standard.db.libA <- mock %>% filter(MockAlign=="miss" & Filt=="standard" & Method=="deblur" & Library=="libA") %>% select(HashID)
miss.standard.db.libB <- mock %>% filter(MockAlign=="miss" & Filt=="standard" & Method=="deblur" & Library=="libB") %>% select(HashID)
miss.standard.db.libC <- mock %>% filter(MockAlign=="miss" & Filt=="standard" & Method=="deblur" & Library=="libC") %>% select(HashID)
miss.standard.db.libD <- mock %>% filter(MockAlign=="miss" & Filt=="standard" & Method=="deblur" & Library=="libD") %>% select(HashID)
## Intersections across Filtering methods, by Library
int.miss.standard.libA <- length(intersect(miss.standard.da.libA$HashID, intersect(miss.standard.db.libA$HashID, miss.standard.vs.libA$HashID)))
int.miss.standard.libB <- length(intersect(miss.standard.da.libB$HashID, intersect(miss.standard.db.libB$HashID, miss.standard.vs.libB$HashID)))
int.miss.standard.libC <- length(intersect(miss.standard.da.libC$HashID, intersect(miss.standard.db.libC$HashID, miss.standard.vs.libC$HashID)))
int.miss.standard.libD <- length(intersect(miss.standard.da.libD$HashID, intersect(miss.standard.db.libD$HashID, miss.standard.vs.libD$HashID)))
int.miss.standard.meth.vals <- c(int.miss.standard.libA, int.miss.standard.libB, int.miss.standard.libC, int.miss.standard.libD)
int.miss.standard.meth.names <- c('int.miss.standard.libA', 'int.miss.standard.libB', 'int.miss.standard.libC', 'int.miss.standard.libD')
int.miss.standard.meth.df <- data.frame(int.miss.standard.meth.vals, int.miss.standard.meth.names)
colnames(int.miss.standard.meth.df) <- c("intersections", "name")

## Intersections across Libraries, by Filtering Method
int.miss.standard.da <- length(intersect(miss.standard.da.libA$HashID, intersect(miss.standard.da.libB$HashID, intersect(miss.standard.da.libC$HashID, miss.standard.da.libD$HashID))))
int.miss.standard.db <- length(intersect(miss.standard.db.libA$HashID, intersect(miss.standard.db.libB$HashID, intersect(miss.standard.db.libC$HashID, miss.standard.db.libD$HashID))))
int.miss.standard.vs <- length(intersect(miss.standard.vs.libA$HashID, intersect(miss.standard.vs.libB$HashID, intersect(miss.standard.vs.libC$HashID, miss.standard.vs.libD$HashID))))
int.miss.standard.lib.vals <- c(int.miss.standard.da, int.miss.standard.db, int.miss.standard.vs)
int.miss.standard.lib.names <- c('int.miss.standard.da', 'int.miss.standard.db', 'int.miss.standard.vs')
int.miss.standard.lib.df <- data.frame(int.miss.standard.lib.vals, int.miss.standard.lib.names)
colnames(int.miss.standard.lib.df) <- c("intersections", "name")

rm(miss.standard.vs.libA, miss.standard.vs.libB, miss.standard.vs.libC, miss.standard.vs.libD, miss.standard.da.libA, miss.standard.da.libB, miss.standard.da.libC, miss.standard.da.libD,
   miss.standard.db.libA, miss.standard.db.libB, miss.standard.db.libC, miss.standard.db.libD, int.miss.standard.libA, int.miss.standard.libB, int.miss.standard.libC, int.miss.standard.libD,
   int.miss.standard.meth.vals, int.miss.standard.meth.names, int.miss.standard.lib.vals, int.miss.standard.lib.names)


## Again, for on "extra" miss matches
miss.extra.vs.libA <- mock %>% filter(MockAlign=="miss" & Filt=="extra" & Method=="vsearch" & Library=="libA") %>% select(HashID)
miss.extra.vs.libB <- mock %>% filter(MockAlign=="miss" & Filt=="extra" & Method=="vsearch" & Library=="libB") %>% select(HashID)
miss.extra.vs.libC <- mock %>% filter(MockAlign=="miss" & Filt=="extra" & Method=="vsearch" & Library=="libC") %>% select(HashID)
miss.extra.vs.libD <- mock %>% filter(MockAlign=="miss" & Filt=="extra" & Method=="vsearch" & Library=="libD") %>% select(HashID)
miss.extra.da.libA <- mock %>% filter(MockAlign=="miss" & Filt=="extra" & Method=="dada2" & Library=="libA") %>% select(HashID)
miss.extra.da.libB <- mock %>% filter(MockAlign=="miss" & Filt=="extra" & Method=="dada2" & Library=="libB") %>% select(HashID)
miss.extra.da.libC <- mock %>% filter(MockAlign=="miss" & Filt=="extra" & Method=="dada2" & Library=="libC") %>% select(HashID)
miss.extra.da.libD <- mock %>% filter(MockAlign=="miss" & Filt=="extra" & Method=="dada2" & Library=="libD") %>% select(HashID)
miss.extra.db.libA <- mock %>% filter(MockAlign=="miss" & Filt=="extra" & Method=="deblur" & Library=="libA") %>% select(HashID)
miss.extra.db.libB <- mock %>% filter(MockAlign=="miss" & Filt=="extra" & Method=="deblur" & Library=="libB") %>% select(HashID)
miss.extra.db.libC <- mock %>% filter(MockAlign=="miss" & Filt=="extra" & Method=="deblur" & Library=="libC") %>% select(HashID)
miss.extra.db.libD <- mock %>% filter(MockAlign=="miss" & Filt=="extra" & Method=="deblur" & Library=="libD") %>% select(HashID)
## Intersections across Filtering methods, by Library
int.miss.extra.libA <- length(intersect(miss.extra.da.libA$HashID, intersect(miss.extra.db.libA$HashID, miss.extra.vs.libA$HashID)))
int.miss.extra.libB <- length(intersect(miss.extra.da.libB$HashID, intersect(miss.extra.db.libB$HashID, miss.extra.vs.libB$HashID)))
int.miss.extra.libC <- length(intersect(miss.extra.da.libC$HashID, intersect(miss.extra.db.libC$HashID, miss.extra.vs.libC$HashID)))
int.miss.extra.libD <- length(intersect(miss.extra.da.libD$HashID, intersect(miss.extra.db.libD$HashID, miss.extra.vs.libD$HashID)))
int.miss.extra.meth.vals <- c(int.miss.extra.libA, int.miss.extra.libB, int.miss.extra.libC, int.miss.extra.libD)
int.miss.extra.meth.names <- c('int.miss.extra.libA', 'int.miss.extra.libB', 'int.miss.extra.libC', 'int.miss.extra.libD')
int.miss.extra.meth.df <- data.frame(int.miss.extra.meth.vals, int.miss.extra.meth.names)
colnames(int.miss.extra.meth.df) <- c("intersections", "name")

## Intersections across Libraries, by Filtering Method
int.miss.extra.da <- length(intersect(miss.extra.da.libA$HashID, intersect(miss.extra.da.libB$HashID, intersect(miss.extra.da.libC$HashID, miss.extra.da.libD$HashID))))
int.miss.extra.db <- length(intersect(miss.extra.db.libA$HashID, intersect(miss.extra.db.libB$HashID, intersect(miss.extra.db.libC$HashID, miss.extra.db.libD$HashID))))
int.miss.extra.vs <- length(intersect(miss.extra.vs.libA$HashID, intersect(miss.extra.vs.libB$HashID, intersect(miss.extra.vs.libC$HashID, miss.extra.vs.libD$HashID))))
int.miss.extra.lib.vals <- c(int.miss.extra.da, int.miss.extra.db, int.miss.extra.vs)
int.miss.extra.lib.names <- c('int.miss.extra.da', 'int.miss.extra.db', 'int.miss.extra.vs')
int.miss.extra.lib.df <- data.frame(int.miss.extra.lib.vals, int.miss.extra.lib.names)
colnames(int.miss.extra.lib.df) <- c("intersections", "name")

rm(miss.extra.vs.libA, miss.extra.vs.libB, miss.extra.vs.libC, miss.extra.vs.libD, miss.extra.da.libA, miss.extra.da.libB, miss.extra.da.libC, miss.extra.da.libD,
   miss.extra.db.libA, miss.extra.db.libB, miss.extra.db.libC, miss.extra.db.libD, int.miss.extra.libA, int.miss.extra.libB, int.miss.extra.libC, int.miss.extra.libD,
   int.miss.extra.meth.vals, int.miss.extra.meth.names, int.miss.extra.lib.vals, int.miss.extra.lib.names)


## Now merge the intersection dataframes into a single object
all.int.df <- do.call(rbind, lapply(ls(pattern = "^int.*.df"), get))
all.int.df$IntType <- all.int.df$name
all.int.df$IntType <- gsub(".*lib.*", "Library", all.int.df$IntType)
all.int.df$IntType <- gsub(".*int.*", "Method", all.int.df$IntType)
all.int.df <- separate(data=all.int.df, col = name, into = c("tmp", "AlignType", "Filt", "Method"), sep="\\.")
all.int.df$tmp <- NULL
## this data.frame served as input to create the data
write.csv(all.int.df, file = "~/Repos/tidybug/data/text_tables/intersection_data.txt", quote = FALSE, row.names = FALSE)

## which "miss" sequences are present in either Dada2/Deblur data?
missers <- mock %>% filter(Method!="vsearch" & MockAlign=="miss" & Filt=="basic") %>% distinct(HashID)
write.table(missers, file="~/Repos/tidybug/data/mock_community/miss.nonVsearch.HashIDs.txt", row.names = FALSE, quote=FALSE, col.names = FALSE)
## this list was then used to query the `mock.partialMisses.fasta` file and then that output was used as input to NCBI blast
## 1) grep -f miss.nonVsearch.HashIDs.txt mock.partialMisses.fasta -A 1 | sed '/^--$/d' > miss.nonVsearch.fasta
## 2) use `miss.nonVsearch.fasta` for NCBI blast 

## how many times are these "miss" ASVs observed across the entire dataset?
df <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz")
nonV.miss <- df %>% filter(HashID %in% missers$HashID)
missSums <- nonV.miss %>% group_by(Method, Filt, HashID, Library) %>% summarise(Counts=n())
missSums$Labeler <- paste(missSums$Method, missSums$Library, missSums$Filt, sep=".")

basic.MissTable <- dcast(missSums %>% filter(Filt=="basic"), HashID ~ Labeler, value.var='Counts')
standard.MissTable <- dcast(missSums %>% filter(Filt=="standard"), HashID ~ Labeler, value.var='Counts')
extra.MissTable <- dcast(missSums %>% filter(Filt=="extra"), HashID ~ Labeler, value.var='Counts')

## what's the mean and median times an ASV is observed in the data (per Method per Library)?
dfSums <- df %>% 
  group_by(Method, HashID, Filt, Library) %>% 
  summarise(Counts=n()) %>% 
  group_by(Method, Filt, Library) %>%
  summarise(MeanCounts=mean(Counts), MedianCounts=median(Counts)) %>%
  mutate(Labeler=paste(Method, Library, Filt, sep="."))
## somewhere between 1-9! nothing like the 100 we're seeing frequently (in plot below)
rm(df)

## easier to visualize?
## order $Filt levels
missSums$Filt <- factor(missSums$Filt, levels=c("basic", "standard", "extra"))
## and plot; save as '4_figure_missASVoutliers-likelyIndexBleedOrContaminant.png'; exported at 800x800
ggplot(missSums, aes(x=Library, y=Counts, color=Library)) + 
  geom_jitter(width=0.3, alpha=0.7) +
  facet_grid(Filt ~ Method) +
  scale_color_manual(values=c('#1f78b4', '#b2df8a', '#a6cee3', '#33a02c')) +
  labs(title="", x="", y="Frequency ASV detected") +
  #geom_hline(mapping = TRUE, data = dfSums, yintercept = dfSums$MeanCounts, linetype="dotted", color="blue") +
  #geom_hline(mapping = TRUE, data = dfSums, yintercept = dfSums$MedianCounts, linetype="dashed", color="red") +
  theme_devon() + theme(legend.position="none")


## did not run: adding in blast pid and pcov info
## did not run: added by adding the `miss.nonVsearch.fasta` file to NCBI blast to collect percent identity (pid) and query coverage (qcov)..
## did not run:  ..values were selected from top hit 
## did not run: blastdat <- read_delim("~/Repos/tidybug/data/mock_community/blastdata.tsv", delim = "\t")
## did not run: missSums <- merge(missSums, blastdat)

## did not run: ggplot(missSums, aes(x=qcov, y=pid, color=Counts)) + 
## did not run:   geom_point() +
## did not run: facet_grid(Library ~ Method) +
## did not run:   scale_color_viridis_c(option="plasma", end = 0.95) +
## did not run:   labs(x="query coverage", y="percent identity", title="", color="frequency\ndetected") +
## did not run:   theme_devon()
