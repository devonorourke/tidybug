# Plots used in paper

library(tidyverse)
library(scales)
library(viridis)

# create theme function for all plots
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

## import data and select mock samples:
df <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz")
mock <- df %>% filter(SampleType == "mock")
mock$Labeler <- paste(mock$HashID, mock$Method, sep="-")
rm(df)
HashFiltLabels <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/HashIDs_withFiltLabels.csv")
mock <- merge(mock, HashFiltLabels)
mock$Labeler <- NULL
rm(HashFiltLabels)
#mock <- x
mock <- mock %>% select(Method, Reads, Filt, Library, MockAlign, Alias)
## summarize number of unique ASVs per group
uniqHashTable <- mock %>% group_by(Method,Library,Filt) %>% distinct(Alias) %>% summarise(counts=n())
uniqHashTable$Merger <- paste(uniqHashTable$Method, uniqHashTable$Library, uniqHashTable$Filt, sep="-")
uniqHashTable <- uniqHashTable[,c(4:5)]
## add these to label in plot
mock$Merger <- paste(mock$Method, mock$Library, mock$Filt, sep="-")
mock <- merge(mock, uniqHashTable)
## build facet label
mock$FacetLabel <- paste(mock$Library, paste0("ASVs =  ",mock$counts), sep=" ")

## generate 3 color palette:
pal3 <- viridis::inferno(3, begin = 0.15, end = 0.85)

## reset the levels for plot
mock$MockAlign <- factor(mock$MockAlign,levels = c("exact", "partial", "miss"))
mock$Filt <- factor(mock$Filt,levels = c("basic", "standard", "extra"))

## change facet labels?
...

## plot; saved as '3_figure_mockSeqCounts_byFilterMethod'; export at 750x750
## note the numbero of distinct samples was manually added following the additional code below producing the `uniqHashTable` object
ggplot(data = mock, aes(x = Method, y = Reads, color=MockAlign)) +
  geom_jitter(alpha=0.55, width = 0.25) +
  scale_color_manual(values = pal3, name = "alignment type") +
  scale_y_continuous(labels = comma, trans = "log2", limits = c(0.01,250000)) +
  facet_grid(Filt ~ Library) +
  labs(title="", x="", y="sequence counts", caption="") +
  theme_devon() +
  theme(legend.position="top", legend.text = element_text(size = 12),
        panel.spacing = unit(1.5, "lines"))
  

        