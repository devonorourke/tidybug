# Plots used in paper

library(tidyverse)
library(reshape2)
library(vegan)
library(phyloseq)
library(stringr)

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
rm(df)

## alpha function applied to calculate diversity measures: observed OTUs, simpson, and shannon
alpha.function <- function(data, filter_exp, filter_exp2) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  Method <- tmp.df %>% distinct(Method)
  Filt <- tmp.df %>% distinct(Filt)
  Rare <- "unrarefied"
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  tmp.mat$SeqID <- NULL
  SeqID <- row.names(tmp.mat)
  Simpson <- diversity(tmp.mat, index="simpson")
  Shannon <- diversity(tmp.mat, index="shannon")
  OTUs <- specnumber(tmp.mat)
  data.frame(SeqID, OTUs, Simpson, Shannon, Method, Filt, Rare, row.names = NULL)
}


## calculate alphas for unrarefied data
dada2.basic <- alpha.function(fox.df, Method=="dada2", Filt=="basic")
dada2.standard <- alpha.function(fox.df, Method=="dada2", Filt=="standard")
dada2.extra <- alpha.function(fox.df, Method=="dada2", Filt=="extra")
deblur.basic <- alpha.function(fox.df, Method=="deblur", Filt=="basic")
deblur.standard <- alpha.function(fox.df, Method=="deblur", Filt=="standard")
deblur.extra <- alpha.function(fox.df, Method=="deblur", Filt=="extra")
vsearch.basic <- alpha.function(fox.df, Method=="vsearch", Filt=="basic")
vsearch.standard <- alpha.function(fox.df, Method=="vsearch", Filt=="standard")
vsearch.extra <- alpha.function(fox.df, Method=="vsearch", Filt=="extra")

## merge into single dataframe
all.fox.unrarefied <- rbind(dada2.basic, dada2.standard, dada2.extra, deblur.basic, deblur.standard, deblur.extra, vsearch.basic, vsearch.standard, vsearch.extra)
rm(dada2.basic, dada2.standard, dada2.extra, deblur.basic, deblur.standard, deblur.extra, vsearch.basic, vsearch.standard, vsearch.extra)

## repeat calculations, but rarefy data first 
## repeat, but this time rarefy the datasets.
## rarefying without replacement in Phyloseq
## alpha function applied to calculate diversity measures: observed OTUs, simpson, and shannon
alpha.function.rfy <- function(data, filter_exp, filter_exp2) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  Method <- tmp.df %>% distinct(Method)
  Filt <- tmp.df %>% distinct(Filt)
  Rare <- "rarefied"
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  rm(tmp.df)
  row.names(tmp.mat) <- tmp.mat$SeqID
  tmp.mat$SeqID <- NULL
  tmp.phy <- otu_table(tmp.mat, taxa_are_rows = TRUE)
  tmp.phy <- phyloseq(tmp.phy)
  set.seed(100)
  rare.phy <- rarefy_even_depth(tmp.phy, rngseed = FALSE, replace = FALSE, trimOTUs = TRUE)
  tmp.mat = as(otu_table(rare.phy), "matrix")
  rm(tmp.phy, rare.phy)
  SeqID <- row.names(tmp.mat)
  Simpson <- diversity(tmp.mat, index="simpson")
  Shannon <- diversity(tmp.mat, index="shannon")
  OTUs <- specnumber(tmp.mat)
  data.frame(SeqID, OTUs, Simpson, Shannon, Method, Filt, Rare, row.names = NULL)
}

r.dada2.basic <- alpha.function.rfy(fox.df, Method=="dada2", Filt=="basic")
r.dada2.standard <- alpha.function.rfy(fox.df, Method=="dada2", Filt=="standard")
r.dada2.extra <- alpha.function.rfy(fox.df, Method=="dada2", Filt=="extra")
r.deblur.basic <- alpha.function.rfy(fox.df, Method=="deblur", Filt=="basic")
r.deblur.standard <- alpha.function.rfy(fox.df, Method=="deblur", Filt=="standard")
r.deblur.extra <- alpha.function.rfy(fox.df, Method=="deblur", Filt=="extra")
r.vsearch.basic <- alpha.function.rfy(fox.df, Method=="vsearch", Filt=="basic")
r.vsearch.standard <- alpha.function.rfy(fox.df, Method=="vsearch", Filt=="standard")
r.vsearch.extra <- alpha.function.rfy(fox.df, Method=="vsearch", Filt=="extra")


all.fox.rarefied <- rbind(r.dada2.basic, r.dada2.standard, r.dada2.extra, r.deblur.basic, r.deblur.standard, r.deblur.extra, r.vsearch.basic, r.vsearch.standard, r.vsearch.extra)
rm(r.dada2.basic, r.dada2.standard, r.dada2.extra, r.deblur.basic, r.deblur.standard, r.deblur.extra, r.vsearch.basic, r.vsearch.standard, r.vsearch.extra)

## combine datasets 
all.fox.alpha <- rbind(all.fox.rarefied, all.fox.unrarefied)
rm(all.fox.rarefied, all.fox.unrarefied)

## add in a the WOY and MonthStart metadata
tmpmeta <- meta %>% select(SeqID, WOY, MonthStart)
all.fox.alpha <- merge(all.fox.alpha, tmpmeta)
rm(tmpmeta, meta)

## plotting only the months with sufficient data for this example (April,May,September,October)
SelectMonths=c("April", "May", "September", "October")
select.fox.alpha <- all.fox.alpha %>% filter(MonthStart %in% SelectMonths)

## add Labeler for plotting:
select.fox.alpha$Labeler <- paste(select.fox.alpha$Filt, select.fox.alpha$Rare, sep="-")

## set the levels
#select.fox.alpha$Rare <- factor(select.fox.alpha$Rare,levels = c("unrarefied", "rarefied"))
#select.fox.alpha$Filt <- factor(select.fox.alpha$Filt,levels = c("basic", "standard", "extra"))
select.fox.alpha$Method <- factor(select.fox.alpha$Method,levels = c("dada2", "deblur", "vsearch"))
select.fox.alpha$Labeler <- factor(select.fox.alpha$Labeler, 
                                   levels = c("basic-unrarefied", "basic-rarefied", "standard-unrarefied", "standard-rarefied", "extra-unrarefied", "extra-rarefied"))
select.fox.alpha$MonthStart <- factor(select.fox.alpha$MonthStart,levels = c("April", "May", "September", "October"))

## generate palette for 6 colors following plot5 color scheme but altering hue:
pal6 <- c('#ebdb8e', '#9f9244', '#c8b2e8', '#6c42b8', '#a9d190', '#628a47')
pal6 <- c('#9f9244', '#ebdb8e', '#6c42b8', '#c8b2e8', '#628a47', '#a9d190')

## plot differences in observed OTUs by $MonthStart
## save as 7_figure_FoxState_OTUs_byFilterMethodsandRarefy; export at 1000x1000
ggplot(select.fox.alpha, aes(x=MonthStart, y=OTUs, color=Labeler)) +
  #scale_y_continuous(trans = "log2") +
  ylim(0,157) +
  scale_color_manual(values = pal6) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.8, position=position_jitterdodge(jitter.width = 0.2)) +
  facet_grid(Method ~ .) +
  labs(title = "", x="", y="Number of unique sequences", color="",
       caption = "6 outliers > 150 OTUs not shown in plot - all vsearch+basic+unrarefied (range 157-223 OTUs / sample)\n Only select data collected at single site (Fox State Forest, Hillsboro NH) shown") +
  theme_devon() + theme(legend.position = "top", axis.text.x = element_text(angle=22.5, hjust=1))



## plot differences in observed Simpson Diversity by $MonthStart
ggplot(select.fox.alpha, aes(x=MonthStart, y=Simpson, color=Rare)) +
  #scale_y_continuous(trans = "log2") +
  scale_color_manual(values= c('#0571b0', '#92c5de')) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.8, position=position_jitterdodge()) +
  facet_grid(Filt ~ Method, scales = "free_y") +
  labs(title = "", x="", y="Simpsons's 1-D", color="") +
  theme_devon() + theme(legend.position = "top")


## plot differences in observed Shannon Diversity by $MonthStart
ggplot(select.fox.alpha, aes(x=MonthStart, y=Shannon, color=Rare)) +
  #scale_y_continuous(trans = "log2") +
  scale_color_manual(values= c('#0571b0', '#92c5de')) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.8, position=position_jitterdodge()) +
  facet_grid(Filt ~ Method, scales = "free_y") +
  labs(title = "", x="", y="Shannon's H", color="") +
  theme_devon() + theme(legend.position = "top")


## running Kruskal-Wallis rank sum test
## p <- capture.output(kruskal.test(AlphaTest ~ MonthStart, data=select.fox.alpha %>% filter(Method=="dada2" & Filt=="basic" & Rare=="unrarefied")))
p <- capture.output(kruskal.test(OTUs ~ MonthStart, data=select.fox.alpha %>% filter(Method=="dada2" & Filt=="basic" & Rare=="unrarefied")))
q <- str_extract(p[5], "p-value = .*")
r <- str_extract(q, "=.*")
s <- gsub("= ", "", r)

## function per batch
kruskal.function <- function(data, filter_exp, filter_exp2, filter_exp3) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  filter_exp_enq3 <- enquo(filter_exp3)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2) %>% filter(!!filter_exp_enq3)
  Method <- tmp.df %>% distinct(Method)
  Filt <- tmp.df %>% distinct(Filt)
  Rare <- tmp.df %>% distinct(Rare)
  p <- capture.output(kruskal.test(OTUs ~ MonthStart, data=tmp.df))
  q <- str_extract(p[5], "p-value = .*")
  r <- str_extract(q, "=.*")
  val <- gsub("= ", "", r)
  val <- round(as.numeric(val), 4)
  data.frame(Method, Filt, Rare, val)
}

dada2.basic.u <- kruskal.function(select.fox.alpha, Method=="dada2", Filt=="basic", Rare=="unrarefied")
dada2.standard.u <- kruskal.function(select.fox.alpha, Method=="dada2", Filt=="standard", Rare=="unrarefied")
dada2.extra.u <- kruskal.function(select.fox.alpha, Method=="dada2", Filt=="extra", Rare=="unrarefied")
deblur.basic.u <- kruskal.function(select.fox.alpha, Method=="deblur", Filt=="basic", Rare=="unrarefied")
deblur.standard.u <- kruskal.function(select.fox.alpha, Method=="deblur", Filt=="standard", Rare=="unrarefied")
deblur.extra.u <- kruskal.function(select.fox.alpha, Method=="deblur", Filt=="extra", Rare=="unrarefied")
vsearch.basic.u <- kruskal.function(select.fox.alpha, Method=="vsearch", Filt=="basic", Rare=="unrarefied")
vsearch.standard.u <- kruskal.function(select.fox.alpha, Method=="vsearch", Filt=="standard", Rare=="unrarefied")
vsearch.extra.u <- kruskal.function(select.fox.alpha, Method=="vsearch", Filt=="extra", Rare=="unrarefied")

dada2.basic.r <- kruskal.function(select.fox.alpha, Method=="dada2", Filt=="basic", Rare=="rarefied")
dada2.standard.r <- kruskal.function(select.fox.alpha, Method=="dada2", Filt=="standard", Rare=="rarefied")
dada2.extra.r <- kruskal.function(select.fox.alpha, Method=="dada2", Filt=="extra", Rare=="rarefied")
deblur.basic.r <- kruskal.function(select.fox.alpha, Method=="deblur", Filt=="basic", Rare=="rarefied")
deblur.standard.r <- kruskal.function(select.fox.alpha, Method=="deblur", Filt=="standard", Rare=="rarefied")
deblur.extra.r <- kruskal.function(select.fox.alpha, Method=="deblur", Filt=="extra", Rare=="rarefied")
vsearch.basic.r <- kruskal.function(select.fox.alpha, Method=="vsearch", Filt=="basic", Rare=="rarefied")
vsearch.standard.r <- kruskal.function(select.fox.alpha, Method=="vsearch", Filt=="standard", Rare=="rarefied")
vsearch.extra.r <- kruskal.function(select.fox.alpha, Method=="vsearch", Filt=="extra", Rare=="rarefied")

krustal.results.OTUs <- rbind(dada2.basic.u, dada2.basic.r, dada2.standard.u, dada2.standard.r, dada2.extra.u, dada2.extra.r,
                         deblur.basic.u, deblur.basic.r, deblur.standard.u, deblur.standard.r, deblur.extra.u, deblur.extra.r,
                         vsearch.basic.u, vsearch.basic.r, vsearch.standard.u, vsearch.standard.r, vsearch.extra.u, vsearch.extra.r)
krustal.results.OTUs$qval <- p.adjust(p = krustal.results.OTUs$val, method = "BH")

rm(dada2.basic.u, dada2.basic.r, dada2.standard.u, dada2.standard.r, dada2.extra.u, dada2.extra.r,
   deblur.basic.u, deblur.basic.r, deblur.standard.u, deblur.standard.r, deblur.extra.u, deblur.extra.r,
   vsearch.basic.u, vsearch.basic.r, vsearch.standard.u, vsearch.standard.r, vsearch.extra.u, vsearch.extra.r)

write.table(krustal.results.OTUs, file = "~/Repos/tidybug/data/text_tables/krustal.results.OTUs.csv", quote=FALSE, row.names = FALSE)


