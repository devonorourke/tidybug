# Plots used in paper

library(tidyverse)
library(reshape2)
library(vegan)

theme_devon <- function () { 
  theme_bw(base_size=12, base_family="Avenir") %+replace% 
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill="transparent", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA)
    )
}

## import data:
df <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz")
df <- df %>% filter(SampleType != "mock")
meta <- read_delim("https://github.com/devonorourke/tidybug/raw/master/data/metadata/large_meta.txt", delim = "\t")
meta <- meta %>% select(SeqID, StudyID, SampleType, Site, Date, DNAplate)
meta_names <- dplyr::intersect(df$SeqID, meta$SeqID)
meta <- meta %>% filter(SeqID %in% meta_names)
rm(meta_names)
meta <- meta %>% filter(Date != "unknown")
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
rm(WOYstring, StartMonthstring, tmp)
meta$MonthStart <- as.character(meta$MonthStart)
meta[is.na(meta)] <- "ncontrol"
meta <- meta %>% select(-StudyID, -SampleType)
## merge data:
df <- merge(df, meta, by='SeqID', all.x=TRUE)

######### Filtering for FOX only data 
targetLibs <- c("libA", "libD")
tmp.df <- df %>% filter(Library %in% targetLibs)
select1 <- "FOX"
dat.tmp1 <- tmp.df %>% filter(Site %in% select1)  ## data from FOX
select2 <- dat.tmp1 %>% distinct(DNAplate)
select2 <- as.character(select2$DNAplate)
dat.tmp2 <- tmp.df %>% filter(SampleType=="ncontrol" & DNAplate %in% select2)
## no negative control samples in any plates associated with FOX samples
## selecting just FOX samples from df object moving forward for diversity analyses
rm(targetLibs, tmp.df, select1, dat.tmp1, select2, dat.tmp2)

fox.df <- df %>% filter(Site=="FOX")

## calculate alpha diversity (OTUs and Simpsons) for each dataframe, ..
# ..split by $Method (dada2, deblur, vsearch) and $Filt (basic, standard, extra)

## create generic alpha function
alpha.function <- function(data, filter_exp, filter_exp2, z) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  tmp.mat$SeqID <- NULL
  tmp.mat$simpson <- diversity(tmp.mat, index=z)
  SeqID <- row.names(tmp.mat)
  AlphaVal <- tmp.mat[,ncol(tmp.mat)]
  data.frame(SeqID, AlphaVal)
}

## calculate alphas
dada2.basic.simpson <- alpha.function(fox.df, Method=="dada2", Filt=="basic", "simpson")
dada2.basic.simpson$Label <- "dada2-basic-simpson"
dada2.standard.simpson <- alpha.function(fox.df, Method=="dada2", Filt=="standard", "simpson")
dada2.standard.simpson$Label <- "dada2-standard-simpson"
dada2.extra.simpson <- alpha.function(fox.df, Method=="dada2", Filt=="extra", "simpson")
dada2.extra.simpson$Label <- "dada2-extra-simpson"
deblur.basic.simpson <- alpha.function(fox.df, Method=="deblur", Filt=="basic", "simpson")
deblur.basic.simpson$Label <- "deblur-basic-simpson"
deblur.standard.simpson <- alpha.function(fox.df, Method=="deblur", Filt=="standard", "simpson")
deblur.standard.simpson$Label <- "deblur-standard-simpson"
deblur.extra.simpson <- alpha.function(fox.df, Method=="deblur", Filt=="extra", "simpson")
deblur.extra.simpson$Label <- "deblur-extra-simpson"
vsearch.basic.simpson <- alpha.function(fox.df, Method=="vsearch", Filt=="basic", "simpson")
vsearch.basic.simpson$Label <- "vsearch-basic-simpson"
vsearch.standard.simpson <- alpha.function(fox.df, Method=="vsearch", Filt=="standard", "simpson")
vsearch.standard.simpson$Label <- "vsearch-standard-simpson"
vsearch.extra.simpson <- alpha.function(fox.df, Method=="vsearch", Filt=="extra", "simpson")
vsearch.extra.simpson$Label <- "vsearch-extra-simpson"

## merge into single dataframe
all.fox.simpson <- rbind(dada2.basic.simpson, dada2.standard.simpson, dada2.extra.simpson, deblur.basic.simpson, deblur.standard.simpson, deblur.extra.simpson, vsearch.basic.simpson, vsearch.standard.simpson, vsearch.extra.simpson)
rm(dada2.basic.simpson, dada2.standard.simpson, dada2.extra.simpson, deblur.basic.simpson, deblur.standard.simpson, deblur.extra.simpson, vsearch.basic.simpson, vsearch.standard.simpson, vsearch.extra.simpson)
rm(df)

## split the $Label
all.fox.simpson <- separate(all.fox.simpson, col=Label, into=c("Method", "Filt", "AlphaTest"), sep = "-")
## add in a the WOY and MonthStart metadata
tmp <- meta %>% select(SeqID, WOY, MonthStart)
all.fox.simpson <- merge(all.fox.simpson, tmp)

## set the levels
all.fox.simpson$Filt <- factor(all.fox.simpson$Filt,levels = c("basic", "standard", "extra"))
all.fox.simpson$WOY <- factor(all.fox.simpson$WOY,
                              levels=c("14","15","16","17","19","20","21","22","23","35","36","37","38","39","40","41","42","43"))
## set the palette to match figure 6
pipepal3 <- c('#9f9244', '#6c42b8', '#628a47', '#a9d190')

## the full plot is really noisy...
ggplot(data = all.fox.simpson, aes(x = WOY, y = AlphaVal, color=WOY)) +
  geom_boxplot() +
  #geom_jitter(alpha=0.55, width = 0.25) +
  facet_grid(Filt ~ Method) +
  #scale_color_manual(values=pipepal3) +
  labs(title="", x="", y="Simpson 1-D", color="") +
  theme_devon() +
  theme(legend.position="top", legend.text = element_text(size = 12), axis.text.x = element_text(angle=22.5, hjust = 1))


## running Kruskal-Wallis rank sum test
kruskal.test(AlphaVal ~ WOY, data=all.fox.simpson %>% filter(Method=="deblur" & Filt=="extra"))


## function per batch
kruskal.function <- function(data, filter_exp, filter_exp2) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  print(kruskal.test(AlphaVal ~ WOY, data=tmp.df))
}

dada2.basic.kruskal <- kruskal.function(all.fox.simpson, Method=="dada2", Filt=="basic")
