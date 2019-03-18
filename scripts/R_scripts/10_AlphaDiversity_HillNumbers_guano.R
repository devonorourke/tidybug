library(tidyverse)
library(reshape2)
library(vegan)
library(cowplot)

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
## filter out mock samples:
df <- df %>% filter(SampleType != "mock") %>% select(-StudyID, -Alias)
df$Labeler <- paste(df$Method, df$Filt, df$Library, sep="-")

## create similar dataframe with rarefied data:
rrarewithdrop <- function(x, sample) {
  rrarefy(x[rowSums(x) >= sample, , drop = FALSE], sample)
}

rare.function <- function(data, filter_exp, filter_exp2) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  meta.df <- tmp.df %>% distinct(SeqID, Library, Method, Filt, SampleType, Labeler)
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  SeqID <- row.names(tmp.mat)
  tmp.mat$SeqID <- NULL
  tmp.mat2 <- data.frame(rrarewithdrop(tmp.mat, 5000))
  tmp.mat2$SeqID <- row.names(tmp.mat2)
  tmp.df2 <- melt(data = tmp.mat2, id.vars = "SeqID",
                  variable.name = "HashID", 
                  value.name = "Reads") %>% 
    filter(Reads > 0) %>%
    mutate(HashID=sub('.', '', HashID))
  tmp.df2 <- merge(tmp.df2, meta.df)
}

## generate rarefied dataset:
tmp.dada2.basic <- rare.function(df, Method=="dada2", Filt=="basic")
tmp.dada2.standard <- rare.function(df, Method=="dada2", Filt=="standard")
tmp.dada2.extra <- rare.function(df, Method=="dada2", Filt=="extra")
tmp.deblur.basic <- rare.function(df, Method=="deblur", Filt=="basic")
tmp.deblur.standard <- rare.function(df, Method=="deblur", Filt=="standard")
tmp.deblur.extra <- rare.function(df, Method=="deblur", Filt=="extra")
tmp.vsearch.basic <- rare.function(df, Method=="vsearch", Filt=="basic")
tmp.vsearch.standard <- rare.function(df, Method=="vsearch", Filt=="standard")
tmp.vsearch.extra <- rare.function(df, Method=="vsearch", Filt=="extra")
rare.df <- rbind(tmp.dada2.basic, tmp.dada2.standard, tmp.dada2.extra,
                 tmp.deblur.basic, tmp.deblur.standard, tmp.deblur.extra,
                 tmp.vsearch.basic, tmp.vsearch.standard, tmp.vsearch.extra)
rm(list=ls(pattern="tmp.*"))


## combine both datasets and add Rarefy label
df$RarefyType <- "unrarefied"
rare.df$RarefyType <- "rarefied"
df <- rbind(df, rare.df)
rm(rare.df)

## Hull Numbers function applied to calculate diversity measures: observed OTUs, simpson, and shannon
hillnumber.function <- function(data, filter_exp, filter_exp2, filter_exp3) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  filter_exp_enq3 <- enquo(filter_exp3)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2) %>% filter(!!filter_exp_enq3)
  Method <- tmp.df %>% distinct(Method)
  Filt <- tmp.df %>% distinct(Filt)
  RarefyType <- tmp.df %>% distinct(RarefyType)
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  SeqID <- row.names(tmp.mat)
  tmp.mat$SeqID <- NULL
  tmp.hill <- renyi(tmp.mat, scales = c(0,1,2), hill=TRUE)
  tmp <- data.frame(tmp.hill, SeqID, Method, Filt, RarefyType, row.names = NULL)
  colnames(tmp)[1:3] <- c("q=0", "q=1", "q=2")
  gather(tmp, key="Hill_qType", value = "Hill_value", c('q=0', 'q=1', 'q=2'))
}

## calculate Hill Numbers for unrarefied and rarefied data
unrarefied.dada2.basic <- hillnumber.function(df, Method=="dada2", Filt=="basic", RarefyType=="unrarefied")
unrarefied.dada2.standard <- hillnumber.function(df, Method=="dada2", Filt=="standard", RarefyType=="unrarefied")
unrarefied.dada2.extra <- hillnumber.function(df, Method=="dada2", Filt=="extra", RarefyType=="unrarefied")
unrarefied.deblur.basic <- hillnumber.function(df, Method=="deblur", Filt=="basic", RarefyType=="unrarefied")
unrarefied.deblur.standard <- hillnumber.function(df, Method=="deblur", Filt=="standard", RarefyType=="unrarefied")
unrarefied.deblur.extra <- hillnumber.function(df, Method=="deblur", Filt=="extra", RarefyType=="unrarefied")
unrarefied.vsearch.basic <- hillnumber.function(df, Method=="vsearch", Filt=="basic", RarefyType=="unrarefied")
unrarefied.vsearch.standard <- hillnumber.function(df, Method=="vsearch", Filt=="standard", RarefyType=="unrarefied")
unrarefied.vsearch.extra <- hillnumber.function(df, Method=="vsearch", Filt=="extra", RarefyType=="unrarefied")

rarefied.dada2.basic <- hillnumber.function(df, Method=="dada2", Filt=="basic", RarefyType=="rarefied")
rarefied.dada2.standard <- hillnumber.function(df, Method=="dada2", Filt=="standard", RarefyType=="rarefied")
rarefied.dada2.extra <- hillnumber.function(df, Method=="dada2", Filt=="extra", RarefyType=="rarefied")
rarefied.deblur.basic <- hillnumber.function(df, Method=="deblur", Filt=="basic", RarefyType=="rarefied")
rarefied.deblur.standard <- hillnumber.function(df, Method=="deblur", Filt=="standard", RarefyType=="rarefied")
rarefied.deblur.extra <- hillnumber.function(df, Method=="deblur", Filt=="extra", RarefyType=="rarefied")
rarefied.vsearch.basic <- hillnumber.function(df, Method=="vsearch", Filt=="basic", RarefyType=="rarefied")
rarefied.vsearch.standard <- hillnumber.function(df, Method=="vsearch", Filt=="standard", RarefyType=="rarefied")
rarefied.vsearch.extra <- hillnumber.function(df, Method=="vsearch", Filt=="extra", RarefyType=="rarefied")

# rm(df)

## merge into single dataframe
all.guano.hill <- rbind(unrarefied.dada2.basic, unrarefied.dada2.standard, unrarefied.dada2.extra, 
                        unrarefied.deblur.basic, unrarefied.deblur.standard, unrarefied.deblur.extra, 
                        unrarefied.vsearch.basic, unrarefied.vsearch.standard, unrarefied.vsearch.extra,
                        rarefied.dada2.basic, rarefied.dada2.standard, rarefied.dada2.extra, 
                        rarefied.deblur.basic, rarefied.deblur.standard, rarefied.deblur.extra, 
                        rarefied.vsearch.basic, rarefied.vsearch.standard, rarefied.vsearch.extra)

rm(list=ls(pattern = "rarefied*"))

## set the levels
all.guano.hill$Filt <- factor(all.guano.hill$Filt, levels = c("basic", "standard", "extra"))
## add labeler to make facet labels easier to understand
all.guano.hill$Labeler <- paste(all.guano.hill$Method, all.guano.hill$Hill_qType, sep="  ")
all.guano.hill$Grouper <- paste(all.guano.hill$Filt, all.guano.hill$RarefyType, sep="-")

## split df by Hill q-value:
u.q0.hill.df <- all.guano.hill %>% filter(Hill_qType == "q=0")
u.q1.hill.df <- all.guano.hill %>% filter(Hill_qType == "q=1")
u.q2.hill.df <- all.guano.hill %>% filter(Hill_qType == "q=2")

## generate palette for 6 colors following plot5 color scheme but altering hue:
#pal3 <- c('#9f9244', '#6c42b8', '#628a47')
pal6 <- c('#9f9244', '#6c42b8', '#628a47',
          '#ebdb8e', '#c8b2e8', '#a9d190')

## order levels:
all.guano.hill$Grouper <- factor(all.guano.hill$Grouper, levels = c("basic-unrarefied", "standard-unrarefied","extra-unrarefied",
                                                                    "basic-rarefied", "standard-rarefied", "extra-rarefied"))

## make separate plots per qtype, then merge into single figure
q0 <- ggplot(data=all.guano.hill %>% filter(Hill_qType=="q=0"),
             aes(x=Grouper, y=Hill_value, color=Grouper, group=Grouper)) +
  scale_color_manual(values = pal6) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.2, position=position_jitterdodge(jitter.width = 0.3)) +
  scale_y_continuous(limits = c(0,250)) +
  facet_grid( ~ Labeler) +
  labs(title = "", x="", y="sequence variants", color="", shape="") +
  theme_devon() + theme(legend.position = "top", 
                        axis.text.x = element_blank(),
                        #axis.text.x = element_text(angle = 22.5, hjust=1), 
                        axis.ticks.x = element_blank(),
                        axis.title.y = element_text(size = 10),
                        plot.margin = margin(c(0,0.5,0,1), unit = "cm")) +
  guides(col = guide_legend(nrow = 1))


q1 <- ggplot(data=all.guano.hill %>% filter(Hill_qType=="q=1"),
             aes(x=Grouper, y=Hill_value, color=Grouper, group=Grouper)) +
  scale_color_manual(values = pal6) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.2, position=position_jitterdodge(jitter.width = 0.3)) +
  facet_grid( ~ Labeler) +
  labs(title = "", x="", y="sequence variant equivalents", color="", shape="") +
  theme_devon() + 
  theme(legend.position = "none", 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 10),
        plot.margin = margin(c(0,0.5,0,1), unit = "cm"))


q2 <- ggplot(data=all.guano.hill %>% filter(Hill_qType=="q=2"),
             aes(x=Grouper, y=Hill_value, color=Grouper, group=Grouper)) +
  scale_color_manual(values = pal6) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.2, position=position_jitterdodge(jitter.width = 0.3)) +
  facet_grid( ~ Labeler) +
  labs(title = "", x="", y="sequence variant equivalents", color="", shape="",
       caption = "21 outliers with > 250 sequence variants (all Vsearch) not shown for q=0") +
  theme_devon() + theme(legend.position = "none", 
                        axis.text.x = element_blank(),
                        #axis.text.x = element_text(angle = 22.5, hjust=1, size = 8), 
                        axis.ticks.x = element_blank(),
                        axis.title.y = element_text(size = 10),
                        plot.margin = margin(c(0,0.5,0,1), unit = "cm"))


## plot; save as 10_guanoAlpha_HillNumbers; export at 900x900
#plot_grid(q0,q1,q2, ncol=1, rel_heights = c(1.3, 0.9, 1.2))
plot_grid(q0,q1,q2, ncol=1, rel_heights = c(1.2, 1, 1))
