##runonce: install.packages("iNEXT")
library(iNEXT)
library(tidyverse)
library(reshape2)
library(scales)
library(cowplot)
## import mock data:
df <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz")
mock <- df %>% filter(SampleType == "mock")
rm(df)

## inext function
inext.function <- function(data, filter_exp, filter_exp2) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  Method <- tmp.df %>% distinct(Method)
  Filt <- tmp.df %>% distinct(Filt)
  tmp.mat <- dcast(tmp.df, HashID ~ SeqID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$HashID
  tmp.mat$HashID <- NULL
  tmp.inext <- iNEXT(tmp.mat, q=c(0,2), datatype = "abundance")
  df <- fortify(tmp.inext, type=1)
  data.frame(df, Method, Filt) %>% select(-datatype, -plottype)
}

## generate data for plot
i.dada2.basic <- inext.function(mock, Method=="dada2", Filt=="basic")
i.dada2.standard <- inext.function(mock, Method=="dada2", Filt=="standard")
i.dada2.extra <- inext.function(mock, Method=="dada2", Filt=="extra")
i.deblur.basic <- inext.function(mock, Method=="deblur", Filt=="basic")
i.deblur.standard <- inext.function(mock, Method=="deblur", Filt=="standard")
i.deblur.extra <- inext.function(mock, Method=="deblur", Filt=="extra")
i.vsearch.basic <- inext.function(mock, Method=="vsearch", Filt=="basic")
i.vsearch.standard <- inext.function(mock, Method=="vsearch", Filt=="standard")
i.vsearch.extra <- inext.function(mock, Method=="vsearch", Filt=="extra")

## combine:
mock.inext.all <- rbind(i.dada2.basic, i.dada2.standard, i.dada2.extra, i.deblur.basic, i.deblur.standard, i.deblur.extra, i.vsearch.basic, i.vsearch.standard, i.vsearch.extra)
rm(i.dada2.basic, i.dada2.standard, i.dada2.extra, i.deblur.basic, i.deblur.standard, i.deblur.extra, i.vsearch.basic, i.vsearch.standard, i.vsearch.extra)

## relabel mock names to match LibA-D like all other plots:
mock.inext.all$site <- as.character(mock.inext.all$site)
mock.inext.all$site[which(mock.inext.all$site=="mockIM4p4L1")] = "libA"
mock.inext.all$site[which(mock.inext.all$site=="mockIM4p4L2")] = "libB"
mock.inext.all$site[which(mock.inext.all$site=="mockIM4p7L1")] = "libC"
mock.inext.all$site[which(mock.inext.all$site=="mockIM4p7L2")] = "libD"

## reset levels for plot 
mock.inext.all$Filt <- factor(mock.inext.all$Filt,levels = c("basic", "standard", "extra"))

## color palette
pal3 <- c('#9f9244', '#6c42b8', '#628a47', '#a9d190')

## plot function
theme_devon <- function () { 
  theme_bw(base_size=12, base_family="Avenir") %+replace% 
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill="transparent", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA)
    )
}

## split dataset by order (q=0 or q=2)
mock.inext.q0 <- mock.inext.all %>% filter(order==0)
mock.inext.q2 <- mock.inext.all %>% filter(order==2)

## split again into two plots to help visualize distinct axes for vsearch vs. deblur/dada2
dadbq0 <- ggplot(data = mock.inext.q0 %>% filter(Method != "vsearch"), aes(x=x, y=y, color=Filt, shape=Filt)) +
  geom_point(data = mock.inext.q0 %>% filter(method=="observed" & Method != "vsearch")) +
  geom_line(data = mock.inext.q0 %>% filter(method=="interpolated" & Method != "vsearch")) +
  geom_line(data = mock.inext.q0 %>% filter(method=="extrapolated" & Method != "vsearch"), linetype="dashed") +
  geom_ribbon(data = mock.inext.q0 %>% filter(Method != "vsearch"), mapping=aes(x=x, ymin=y.lwr, ymax=y.upr, fill=Filt), alpha=0.1) +
  facet_grid(Method ~ site, scales = "free") +
  scale_color_manual(values=pal3) +
  scale_x_continuous(labels=comma) +
  ylim(0,47) +
  labs(title="", x="", y="number of sequence variants", caption = "Shared yaxis scale for dada2 and deblur facets; alternate yaxis scale for vsearch. Note the x-axis varies freely to account for varying library sizes") +
  theme_devon() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "top", legend.title = element_blank())

vsq0 <- ggplot(data = mock.inext.q0 %>% filter(Method == "vsearch"), aes(x=x, y=y, color=Filt, shape=Filt)) +
  geom_point(data = mock.inext.q0 %>% filter(method=="observed" & Method == "vsearch")) +
  geom_line(data = mock.inext.q0 %>% filter(method=="interpolated" & Method == "vsearch")) +
  geom_line(data = mock.inext.q0 %>% filter(method=="extrapolated" & Method == "vsearch"), linetype="dashed") +
  geom_ribbon(data = mock.inext.q0 %>% filter(Method == "vsearch"), mapping=aes(x=x, ymin=y.lwr, ymax=y.upr, fill=Filt), alpha=0.1) +
  facet_grid(Method ~ site, scales = "free") +
  scale_x_continuous(labels=comma) +
  scale_color_manual(values=pal3) +
  ylim(0,1000) +
  #ylim(0,165) +
  labs(x="number of sequences", y="number of sequence variants") +
  theme_devon() +
  theme(axis.text.x = element_text(angle=22.5, hjust=1), legend.position = "none") +
  theme(strip.background.x = element_blank(), strip.text.x = element_blank())

##plot; save as 8_mock_specaccum_byFiltMethod_q0; export at 1000x800
q0 <- plot_grid(dadbq0, vsq0, ncol=1, rel_heights = c(2,1))
rm(dadbq0, vsq0)

##plot; save as 8_mock_specaccum_byFiltMethod_q2; export at 1000x800
ggplot(data = mock.inext.q2, aes(x=x, y=y, color=Filt, shape=Filt)) +
  geom_point(data = mock.inext.q2 %>% filter(method=="observed")) +
  geom_line(data = mock.inext.q2 %>% filter(method=="interpolated")) +
  geom_line(data = mock.inext.q2 %>% filter(method=="extrapolated"), linetype="dashed") +
  geom_ribbon(data = mock.inext.q2, mapping=aes(x=x, ymin=y.lwr, ymax=y.upr, fill=Filt), alpha=0.1) +
  facet_grid(Method ~ site, scales = "free_x") +
  scale_color_manual(values=pal3) +
  scale_x_continuous(labels=comma) +
  labs(title="", x="number of sequences", y="number of sequence variants") +
  theme_devon() + 
  theme(legend.position = "top", legend.title = element_blank(),
        axis.text.x = element_text(angle=22.5, hjust=1))

