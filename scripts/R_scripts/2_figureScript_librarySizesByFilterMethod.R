library(tidyverse)
library(scales)
library(ggrepel)

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

## read in data and calculate sum of reads per Library, per Method
df <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz")
guanoSums_tmp <- df %>% 
  filter(Filt=="basic" & SampleType=="sample") %>% 
  group_by(Method, Library) %>% 
  summarize(TrimSum=sum(Reads)) %>% 
  mutate(ReadType="guano")
mockSums_tmp <- df %>% 
  filter(Filt=="basic"& SampleType=="mock") %>% 
  group_by(Method, Library) %>% 
  summarize(TrimSum=sum(Reads)) %>% 
  mutate(ReadType="mock")
rm(df)

## add in raw read sums per library; raw "Joined" values assigned above derived by exporting "$LIB".joind.seqs.qza and ..
## ..counting number of N/4 lines of all data. see "seqfilter.qfilt_example.sh" for details
## Need to subtract the values that count towards 'mock' samples per library though too!
## operations were: 1) qiime tools export --import-path "$LIB".joind.seqs.qza --output-path joind
## and  ...       : 2) cat ./joind/*.gz | zcat | awk '{s++}END{print s/4}'
## Repeated same idea for ONLY the four mock samples so that we could calculate all guano contributions separately from mock data
## 3) for example: zgrep -c '1:N:0:ACTGTGTA+CTAGTATG' mockIM4p4L1_L001_R1_001.fastq.gz

LibSums <- data.frame(Library=as.character(c("libA", "libB", "libC", "libD")),
                      LibReads=as.numeric(c(14724385, 17071355, 14279338, 14171652)),
                      MockReads=as.numeric(c(181012, 456441, 228383, 1176419)))
LibSums$GuanoReads <- LibSums$LibReads - LibSums$MockReads

## merge guano and mock data separately with the `LibSums` object
guanoLibSums <- merge(guanoSums_tmp, LibSums) %>% 
  select(-LibReads, -MockReads) %>% 
  mutate(ReadsRetained = round((TrimSum / GuanoReads), 2)) %>% 
  select(-GuanoReads)

mockLibSums <- merge(mockSums_tmp, LibSums) %>% 
  select(-LibReads, -GuanoReads) %>% 
  mutate(ReadsRetained = round((TrimSum / MockReads), 2)) %>% 
  select(-MockReads)

rm(guanoSums_tmp, mockSums_tmp)

## combine together into single object for plot:
plotdat <- rbind(guanoLibSums, mockLibSums)
rm(guanoLibSums, mockLibSums)

## set levels for plot
plotdat$Method <- factor(plotdat$Method,levels = c("vsearch", "deblur", "dada2"))

## plot; exported at 881x468; save as '2_figure_librarySizesByFilterMethod'
ggplot(data = plotdat, aes(x = ReadsRetained, y = Method, label=Library)) + 
  geom_text_repel(vjust=1, nudge_y = 0.2, segment.alpha = 0.2) +
  geom_point() +
  facet_wrap(~ReadType) +
  scale_shape_manual(values=c(65,66,67,68)) +
  scale_x_continuous(limits = c(0.40, 1.0)) +
  geom_line(aes(group = Method), color="gray50") +
  labs(title="", y="", x="fraction of retained sequences") +
  theme_devon()

## calculate percentage of reads retained per Method (aggregate all Libraries, combining Mock and Bat Guano together)
## not used for plots; just used for calculation...
MethodSums <- df %>% filter(Filt=="basic") %>% group_by(Method) %>% summarize(TrimSum=sum(Reads))
JoinedReadSums <- sum(14724385, 17071355, 14279338, 14171652)
MethodSums$JoinedReadSums <- JoinedReadSums
MethodSums <- MethodSums %>% mutate(ReadsRetained = TrimSum/JoinedReadSums) %>% select(Method, ReadsRetained)

