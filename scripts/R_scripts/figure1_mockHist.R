# author: devon o'rourke
# contact: devon@outermostlab.com
# Script used for figure 1 in paper

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
HashFiltLabels <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/HashIDs_withFiltLabels.csv")  ## from '1_sequence_filter.R` script`
mock_tmp <- df %>% filter(SampleType == "mock") %>% 
  mutate(Labeler = paste(HashID, Method, sep="-"),
         log2Reads = log2(Reads))
mock <- merge(mock_tmp, HashFiltLabels) %>% 
  select(-Labeler)
rm(HashFiltLabels, mock_tmp, df)

## summarize number of unique ASVs per group; these are manually added to image and used for stats below
uniqHashTable <- mock %>% 
  group_by(Method,Filt,MockAlign) %>% 
  summarise(counts=n()) %>% 
  mutate(log2Reads=-2)


## generate 3 color palette:
pal3 <- viridis::inferno(3, begin = 0.15, end = 0.85)

## reset the levels for plot
mock$MockAlign <- factor(mock$MockAlign,levels = c("exact", "partial", "miss"))
mock$Filt <- factor(mock$Filt,levels = c("basic", "standard", "extra"))

## plot; export at 750x750
## note the numbero of distinct samples was manually added following the additional code below producing the `uniqHashTable` object
ggplot(data = mock, aes(x = MockAlign, y = log2Reads)) +
  geom_jitter(alpha=0.55, width = 0.25) +
  facet_grid(Filt ~ Method) +
  scale_y_continuous(labels = c("1", "16", "256", "4,096", "65,536"),
                     breaks = c(0, 4, 8, 12, 16),
                     limits = c(-3, 17)) +
  labs(title="", x="\n alignment type", y="sequence counts\n") +
  theme_devon() +
  theme(legend.position="none",
        strip.text.x = element_text(size=13), strip.text.y = element_text(size=13),
        panel.spacing = unit(1.5, "lines"), panel.grid.major.x = element_blank() ) +
  geom_text(aes(x=MockAlign,y=log2Reads,label=counts),
            data=uniqHashTable)



################################################################################
## stats evaluating if mean ranks of 
################################################################################

## read data are not normally distributed...
shapiro_function1 <- function(filtmethod){
  set.seed(10)  
  da_tmp <- mock %>% filter(Filt==filtmethod, Method=="dada2")
  db_tmp <- mock %>% filter(Filt==filtmethod, Method=="deblur")
  vs_tmp <- mock %>% filter(Filt==filtmethod, Method=="vsearch")
  print(shapiro.test(da_tmp$Reads))
  print(shapiro.test(db_tmp$Reads))
  print(shapiro.test(vs_tmp$Reads))
}

shapiro_function1("basic")    ## NOT normal for any dataset
shapiro_function1("standard")    ## NOT normal for any dataset
shapiro_function1("extra")    ## NOT normal for any dataset

## run Kruskal-Wallis to test for rank ordered mean differences in # Reads between each Method (one per filter parameter)
## apply post hoc Dunn's Test for pairwise differences
mock_kw_function <- function(filtmethod){
  print(kruskal.test(Reads ~ as.factor(Method), data = mock %>% filter(Filt == filtmethod)))
}

mock_kw_function("basic")    ## kw sig diff
mock_kw_function("standard") ## kw sig diff
mock_kw_function("extra")    ## kw sig diff

mock_dunn_function <- function(filtmethod){
  dunn_tmp <- dunnTest(Reads ~ as.factor(Method), data = mock %>% filter(Filt == filtmethod), method = "bh")
  comps = dunn_tmp$res$Comparison
  stat = dunn_tmp$res$Z
  pval = round(dunn_tmp$res$P.adj, 4)
  print(data.frame(comps, stat, pval))
}

mock_dunn_function("basic")    ## 
mock_dunn_function("standard") ## 
mock_dunn_function("extra")    ## 

################################################################################
## unused code
################################################################################

## second plot... summing up the individual data points by $MockAlign group
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

