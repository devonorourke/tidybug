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
LibSums <- df %>% filter(Filt=="basic") %>% group_by(Method, Library) %>% summarize(TrimSum=sum(Reads))

## add in raw read sums per library; raw "Joined" values assigned above derived by exporting "$LIB".joind.seqs.qza and ..
## ..counting number of N/4 lines of all data. see "seqfilter.qfilt_example.sh" for details
## operations were: 1) qiime tools export --import-path "$LIB".joind.seqs.qza --output-path joind
## and  ...       : 2) cat ./joind/*.gz | zcat | awk '{s++}END{print s/4}'
LibSums$JoindReads <- ""
LibSums$JoindReads[LibSums$Library=="libA"] <- 14724385
LibSums$JoindReads[LibSums$Library=="libB"] <- 17071355
LibSums$JoindReads[LibSums$Library=="libC"] <- 14279338
LibSums$JoindReads[LibSums$Library=="libD"] <- 14171652
LibSums$JoindReads <- as.numeric(LibSums$JoindReads) 

## calculate percent of reads retained per Library per Method
LibSums <- LibSums %>% mutate(ReadsRetained = TrimSum/JoindReads)
LibSums$Method <- factor(LibSums$Method,levels = c("vsearch", "dada2", "deblur"))

## plot; exported at 750x300; save as '2_figure_librarySizesByFilterMethod.png'
ggplot(data = LibSums, aes(x = ReadsRetained, y = Method, label=Library)) + 
  geom_text_repel(vjust=1, nudge_y = 0.2, segment.alpha = 0.2) +
  geom_point() +
  #geom_point(aes(shape=Library), size=4, alpha=1) +
  scale_shape_manual(values=c(65,66,67,68)) +
  #scale_shape_manual(values=c(0,1,2,5)) +
  geom_line(aes(group = Method), color="gray50") +
  labs(title="", y="", x="fraction of retained sequences") +
  theme_devon()

## calculate percentage of reads retained per Method (aggregate all Libraries)
## not used for plots; just used for calculation...
MethodSums <- df %>% filter(Filt=="basic") %>% group_by(Method) %>% summarize(TrimSum=sum(Reads))
JoinedReadSums <- sum(14724385, 17071355, 14279338, 14171652)
MethodSums$JoinedReadSums <- JoinedReadSums
MethodSums <- MethodSums %>% mutate(ReadsRetained = TrimSum/JoinedReadSums) %>% select(Method, ReadsRetained)
