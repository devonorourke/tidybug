library(tidyverse)

meta <- read_delim('https://github.com/devonorourke/pzero/raw/master/data/clean_metadata.txt', delim = "\t")
p41 <- read.csv("~/Desktop/p41.manifest.file")
p42 <- read.csv("~/Desktop/p42.manifest.file")
p71 <- read.csv("~/Desktop/p71.manifest.file")
p72 <- read.csv("~/Desktop/p72.manifest.file")
pall <- rbind(p41, p42, p71, p72)

pSamples <- pall$sample.id
samplemeta <- intersect(meta$SeqID, pall$sample.id)
meta <- meta %>% filter(SeqID %in% samplemeta)
#rm(p41, p42, p71, p72, pSamples, samplemeta)
########

meta$Site <- gsub("control", "ncontrol", meta$Site)
meta$Date <- gsub("control", "ncontrol", meta$Date)

## overwrite new $Date and $WOY columns with lubridate package to ensure we're selecting a consistent WOY
meta$NewDate <- as.character(lubridate::mdy(meta$Date))
## any samples we'd be discarding or are NA are converted to generic date
meta$NewDate[which(meta$SampleType == "ncontrol")] = "2000-01-01"
meta$NewDate[which(meta$Date == "unknown")] = "2000-01-01"
meta$NewDate[which(meta$SampleID == "35A01")] = "2000-01-01"
meta$NewDate[which(meta$SampleType == "mock")] = "2000-01-01"

meta$Library <- meta$SeqBatch
meta$Library <- gsub("4.1", "libA", meta$Library)
meta$Library <- gsub("4.2", "libB", meta$Library)
meta$Library <- gsub("7.1", "libC", meta$Library)
meta$Library <- gsub("7.2", "libD", meta$Library)


newmeta <- meta[,c(5,11,14,23,24)]
newmeta$newSite <- newmeta$Site
newmeta$newSite <- gsub("HOP", "Hopkinton NH", newmeta$newSite)
newmeta$newSite <- gsub("ACA", "Acadia National Park ME", newmeta$newSite)
newmeta$newSite <- gsub("YRK", "Yorktown Naval Weapons Station VA", newmeta$newSite)
newmeta$newSite <- gsub("ELY", "Ely Mine VT", newmeta$newSite)
newmeta$newSite <- gsub("BRN", "Brown Lane Hollis NH", newmeta$newSite)
newmeta$newSite <- gsub("MAP", "Maple Hill Hollis NH", newmeta$newSite)
newmeta$newSite <- gsub("ROL", "Rollinsford NH", newmeta$newSite)
newmeta$newSite <- gsub("PEN", "Penacook NH", newmeta$newSite)
newmeta$newSite <- gsub("EPS", "Epsom NH", newmeta$newSite)
newmeta$newSite <- gsub("CNB", "Canterbury NH", newmeta$newSite)
newmeta$newSite <- gsub("GIL", "Gilsum NH", newmeta$newSite)
newmeta$newSite <- gsub("HOL", "Holderness NH", newmeta$newSite)
newmeta$newSite <- gsub("CHI", "Chichester NH", newmeta$newSite)
newmeta$newSite <- gsub("MAS", "Massabesic NH", newmeta$newSite)
newmeta$newSite <- gsub("FOX", "Fox State Forest NH", newmeta$newSite)
newmeta$newSite <- gsub("FAR", "Fairfield ME", newmeta$newSite)
newmeta$newSite <- gsub("MTV", "Mount Vernon NH", newmeta$newSite)
newmeta$newSite <- gsub("ALS", "Alstead NH", newmeta$newSite)
newmeta$newSite[which(newmeta$SampleType == "mock")] = "mock"

newmeta$Site <- NULL

newmeta$bioproject_accession <- "SUB5109424"
newmeta$organism <- "not collected"
newmeta$env_broad_scale <- "not collected"
newmeta$env_local_scale <- "not applicable"
newmeta$env_medium <- "not applicable"
newmeta$host <- "not collected"
newmeta$lat_lon <- "not applicable"
newmeta$AltDate <- newmeta$NewDate

newmeta$AltDate[which(newmeta$SampleType == "ncontrol")] = "not applicable"
newmeta$AltDate[which(newmeta$NewDate == "unknown")] = "not collected"
newmeta$AltDate[which(newmeta$SampleType == "mock")] = "not applicable"
newmeta$AltDate[which(newmeta$newSite == "Canterbury NH" & newmeta$Library == "libD" & newmeta$NewDate == "2000-01-01")] = "not collected"

write.csv(newmeta, file="~/Desktop/newmeta.csv", quote=FALSE, row.names = FALSE)
