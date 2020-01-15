# Plots used in paper

library(tidyverse)
library(scales)
library(viridis)
library(ggrepel)
library(reshape2)

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

## generate 3 color palette:
pal3 <- viridis::inferno(3, begin = 0.15, end = 0.85)

## reset the levels for plot
mock.all$Type <- factor(mock.all$Type,levels = c("exact", "partial", "miss"))
mock.all$Library <- factor(mock.all$Library,levels = c("libA", "libB", "libC", "libD"))
mock.all$Labeler <- factor(mock.all$Labeler,levels = c("unfiltered", "filtered"))

## save as fig2A; export at EXACTLY 1000x504... any resizing screws with label positions
ggplot(data = mock.all, aes(x = Reads, y = Method, color=Type, label = Alias)) +
  geom_jitter(alpha=0.6, height = 0.2) +
  scale_color_manual(values = pal3, name = "") +
  scale_x_continuous(labels = comma, trans = "log2") +
  facet_grid(Labeler ~ Library) +
  #geom_text_repel(data = mock %>% filter(Type == "partial" & Reads > 2000), direction="x", size=3, vjust = 0, force = 3, nudge_y = 0.2, angle=45,segment.size = 0.2) +
  labs(title="", x="number of sequences", y="") +
  theme_devon() +
  theme(legend.position="top")