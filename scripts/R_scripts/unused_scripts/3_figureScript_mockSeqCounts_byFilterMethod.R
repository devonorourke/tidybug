# Plots used in paper

library(tidyverse)
library(scales)
library(viridis)
library(FSA)

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
HashFiltLabels <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/HashIDs_withFiltLabels.csv")  ## from '1_sequence_filter.R` script`
mock <- merge(mock, HashFiltLabels)
mock$Labeler <- NULL
rm(HashFiltLabels)
mock <- mock %>% select(Method, Reads, Filt, Library, MockAlign, Alias)

## generate 3 color palette:
pal3 <- viridis::inferno(3, begin = 0.15, end = 0.85)

## reset the levels for plot
mock$MockAlign <- factor(mock$MockAlign,levels = c("exact", "partial", "miss"))
mock$Filt <- factor(mock$Filt,levels = c("basic", "standard", "extra"))

## plot; saved as '3_figure_mockSeqCounts_byFilterMethod'; export at 750x750
## note the numbero of distinct samples was manually added following the additional code below producing the `uniqHashTable` object
ggplot(data = mock, aes(x = MockAlign, y = Reads, color=MockAlign)) +
  geom_jitter(alpha=0.55, width = 0.25) +
  scale_color_manual(values = pal3, name = "alignment type") +
  scale_y_continuous(labels = comma, trans = "log2", 
                     breaks = c(16, 1024, 65536),
                     limits = c(0.1, max(mock$Reads))) +
  facet_grid(Filt ~ Method) +
  labs(title="", x="", y="sequence counts", caption="") +
  theme_devon() +
  theme(legend.position="top", legend.text = element_text(size = 12),
        panel.spacing = unit(1.5, "lines"), 
        panel.grid.major.x = element_blank() )  


## summarize number of unique ASVs per group; these are manually added to image
## summarize number of unique ASVs per group
uniqHashTable <- mock %>% group_by(Method,Filt,MockAlign) %>% summarise(counts=n())
uniqHashTable_perLib <- mock %>% group_by(Method,Filt,MockAlign,Library) %>% summarise(counts=n())

## stats comparing if number of ASVs in each 'MockAlign' group differ by Method (Filt groups seperate tests)
## Kruskal-Wallis version -- group the Method + MockAlign into single factor:
## perform the KW test and post hoc Dunn's Test for 3 filtering methods (basic, standard, extra):

mock_stats_func <- function(filtmethod){
  tmp_df <- uniqHashTable_perLib %>% filter(Filt == filtmethod)
  tmp_df$labeler <- as.factor(paste(tmp_df$Method, tmp_df$MockAlign, sep="-"))
  print(kruskal.test(counts ~ labeler, data=tmp_df))
  dunn_tmp <- dunnTest(counts ~ labeler, data = tmp_df, method="bh")
  dunn_tmp_df <- dunn_tmp$res
  dunn_tmp_df
}

stats_mock_basic <- mock_stats_func("basic")
stats_mock_standard <- mock_stats_func("standard")
stats_mock_extra <- mock_stats_func("extra")



################################################################################
## second plot... summing up the individual data points by $MockAlign group
################################################################################

## summarizing the total number of reads per `$MockAlign` group PER LIBRARY:
sumReads_MethFiltLibAligngrp <- mock %>% 
  group_by(Method, Filt, Library, MockAlign) %>% 
  summarise(groupReads=sum(Reads)) %>% 
  group_by(Method, Filt, Library) %>% 
  mutate(perLibReads = sum(groupReads)) %>% 
  mutate(fracGroupReads = groupReads / perLibReads)

## plot
sumReads_MethFiltLibAligngrp$Filt <- factor(sumReads_MethFiltLibAligngrp$Filt, levels = c("basic", "standard", "extra"))
sumReads_MethFiltLibAligngrp$MockAlign <- factor(sumReads_MethFiltLibAligngrp$MockAlign, levels = c("exact", "partial", "miss"))

## plot; saved as 's1_figure_mockSeqCounts_byFilterMethodLib_stackplot'; export at 750x750
p1 <- ggplot(sumReads_MethFiltLibAligngrp, aes(x=Library, y=groupReads, fill=MockAlign)) + 
  geom_bar(stat="identity") + 
  facet_grid(Filt ~ Method) +
  scale_fill_manual(values = pal3, name = "alignment type") +
  scale_y_continuous(labels=comma) +
  labs(x="", y="sequence counts") +
  #scale_y_continuous(labels = comma, trans = 'log10', breaks = c(10, 1000, 100000)) +
  theme_devon() + theme(legend.position = "top")

## plot; saved as 's2_figure_mockSeqCounts_byFilterMethodLib_stackplotSCALED'; export at 750x750
p2 <- ggplot(sumReads_MethFiltLibAligngrp, aes(x=Library, y=fracGroupReads, fill=MockAlign)) + 
  geom_bar(stat="identity") + 
  facet_grid(Filt ~ Method) +
  scale_fill_manual(values = pal3, name = "alignment type") +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x="", y="fraction of sequence counts") +
  theme_devon() + theme(legend.position = "top")

## stitch both plots together; save as 's1_figure_mockSeqCounts_byFilterMethodLib_stackplotBOTH'; export at 1000x500
require(ggpubr)
ggarrange(p1, NULL, p2, common.legend = TRUE, widths = c(1, 0.1, 1), ncol=3)
