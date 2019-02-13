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
mock <- df %>% filter(SampleType == "mock")
rm(df)
mock$SeqID <- as.character(mock$SeqID)
mock$SeqID[which(mock$SeqID=="mockIM4p4L1")] <- "libA"
mock$SeqID[which(mock$SeqID=="mockIM4p4L2")] <- "libB"
mock$SeqID[which(mock$SeqID=="mockIM4p7L1")] <- "libC"
mock$SeqID[which(mock$SeqID=="mockIM4p7L2")] <- "libD"

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
dada2.basic <- alpha.function(mock, Method=="dada2", Filt=="basic")
dada2.standard <- alpha.function(mock, Method=="dada2", Filt=="standard")
dada2.extra <- alpha.function(mock, Method=="dada2", Filt=="extra")
deblur.basic <- alpha.function(mock, Method=="deblur", Filt=="basic")
deblur.standard <- alpha.function(mock, Method=="deblur", Filt=="standard")
deblur.extra <- alpha.function(mock, Method=="deblur", Filt=="extra")
vsearch.basic <- alpha.function(mock, Method=="vsearch", Filt=="basic")
vsearch.standard <- alpha.function(mock, Method=="vsearch", Filt=="standard")
vsearch.extra <- alpha.function(mock, Method=="vsearch", Filt=="extra")

## merge into single dataframe
all.mock.unrarefied <- rbind(dada2.basic, dada2.standard, dada2.extra, deblur.basic, deblur.standard, deblur.extra, vsearch.basic, vsearch.standard, vsearch.extra)
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
  tmp.phy <- otu_table(tmp.mat, taxa_are_rows = FALSE)
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

r.dada2.basic <- alpha.function.rfy(mock, Method=="dada2", Filt=="basic")
r.dada2.standard <- alpha.function.rfy(mock, Method=="dada2", Filt=="standard")
r.dada2.extra <- alpha.function.rfy(mock, Method=="dada2", Filt=="extra")
r.deblur.basic <- alpha.function.rfy(mock, Method=="deblur", Filt=="basic")
r.deblur.standard <- alpha.function.rfy(mock, Method=="deblur", Filt=="standard")
r.deblur.extra <- alpha.function.rfy(mock, Method=="deblur", Filt=="extra")
r.vsearch.basic <- alpha.function.rfy(mock, Method=="vsearch", Filt=="basic")
r.vsearch.standard <- alpha.function.rfy(mock, Method=="vsearch", Filt=="standard")
r.vsearch.extra <- alpha.function.rfy(mock, Method=="vsearch", Filt=="extra")


all.mock.rarefied <- rbind(r.dada2.basic, r.dada2.standard, r.dada2.extra, r.deblur.basic, r.deblur.standard, r.deblur.extra, r.vsearch.basic, r.vsearch.standard, r.vsearch.extra)
rm(r.dada2.basic, r.dada2.standard, r.dada2.extra, r.deblur.basic, r.deblur.standard, r.deblur.extra, r.vsearch.basic, r.vsearch.standard, r.vsearch.extra)

## combine datasets 
all.mock.alpha <- rbind(all.mock.rarefied, all.mock.unrarefied)
rm(all.mock.rarefied, all.mock.unrarefied)

## add Labeler for plotting:
all.mock.alpha$Labeler <- paste(all.mock.alpha$Filt, all.mock.alpha$Rare, sep="-")

## set values to numeric:
all.mock.alpha$Simpson <- as.numeric(all.mock.alpha$Simpson)
all.mock.alpha$Shannon <- as.numeric(all.mock.alpha$Shannon)

## set the levels
all.mock.alpha$Method <- factor(all.mock.alpha$Method,levels = c("dada2", "deblur", "vsearch"))
all.mock.alpha$Labeler <- factor(all.mock.alpha$Labeler, 
                                   levels = c("basic-unrarefied", "basic-rarefied", "standard-unrarefied", "standard-rarefied", "extra-unrarefied", "extra-rarefied"))

## generate palette for 6 colors following plot5 color scheme but altering hue:
pal6 <- c('#9f9244', '#ebdb8e', '#6c42b8', '#c8b2e8', '#628a47', '#a9d190')


## save as 9_figure_mock_OTUs_byFilterMethodsandRarefy; export at 1000x1000
ggplot(all.mock.alpha, aes(x=SeqID, y=OTUs, fill=Labeler)) +
  scale_fill_manual(values = pal6) +
  geom_bar(stat="identity", position = "dodge", color="black") +
  #facet_grid(Method ~ ., scales = "free_y") +
  facet_grid(Method ~ .) +
  scale_y_continuous(trans="log2") +
  labs(title = "", x="", y="Number of unique sequences", fill="",
       caption="Yaxis scale varies by horizontal facet") +
  theme_devon() + theme(legend.position = "top", axis.text.x = element_text(angle=22.5, hjust=1))



## plot differences in observed Simpson Diversity
ggplot(all.mock.alpha, aes(x=SeqID, y=Simpson, fill=Labeler)) +
  scale_fill_manual(values = pal6) +
  geom_bar(stat="identity", position = "dodge", color="black") +
  facet_grid(Method ~ .) +
  scale_y_continuous(trans="log2") +
  labs(title = "", x="", y="Number of unique sequences", fill="",
       caption="Yaxis scale varies by horizontal facet") +
  theme_devon() + theme(legend.position = "top", axis.text.x = element_text(angle=22.5, hjust=1))


## plot differences in observed Shannon Diversity
ggplot(all.mock.alpha, aes(x=SeqID, y=Shannon, fill=Labeler)) +
  scale_fill_manual(values = pal6) +
  geom_bar(stat="identity", position = "dodge", color="black") +
  #facet_grid(Method ~ ., scales = "free_y") +
  facet_grid(Method ~ .) +
  scale_y_continuous(trans="log2") +
  labs(title = "", x="", y="Number of unique sequences", fill="",
       caption="Yaxis scale varies by horizontal facet") +
  theme_devon() + theme(legend.position = "top", axis.text.x = element_text(angle=22.5, hjust=1))


## running Kruskal-Wallis rank sum test
## p <- capture.output(kruskal.test(AlphaTest ~ MonthStart, data=all.mock.alpha %>% filter(Method=="dada2" & Filt=="basic" & Rare=="unrarefied")))
p <- capture.output(kruskal.test(OTUs ~ MonthStart, data=all.mock.alpha %>% filter(Method=="dada2" & Filt=="basic" & Rare=="unrarefied")))
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

dada2.basic.u <- kruskal.function(all.mock.alpha, Method=="dada2", Filt=="basic", Rare=="unrarefied")
dada2.standard.u <- kruskal.function(all.mock.alpha, Method=="dada2", Filt=="standard", Rare=="unrarefied")
dada2.extra.u <- kruskal.function(all.mock.alpha, Method=="dada2", Filt=="extra", Rare=="unrarefied")
deblur.basic.u <- kruskal.function(all.mock.alpha, Method=="deblur", Filt=="basic", Rare=="unrarefied")
deblur.standard.u <- kruskal.function(all.mock.alpha, Method=="deblur", Filt=="standard", Rare=="unrarefied")
deblur.extra.u <- kruskal.function(all.mock.alpha, Method=="deblur", Filt=="extra", Rare=="unrarefied")
vsearch.basic.u <- kruskal.function(all.mock.alpha, Method=="vsearch", Filt=="basic", Rare=="unrarefied")
vsearch.standard.u <- kruskal.function(all.mock.alpha, Method=="vsearch", Filt=="standard", Rare=="unrarefied")
vsearch.extra.u <- kruskal.function(all.mock.alpha, Method=="vsearch", Filt=="extra", Rare=="unrarefied")

dada2.basic.r <- kruskal.function(all.mock.alpha, Method=="dada2", Filt=="basic", Rare=="rarefied")
dada2.standard.r <- kruskal.function(all.mock.alpha, Method=="dada2", Filt=="standard", Rare=="rarefied")
dada2.extra.r <- kruskal.function(all.mock.alpha, Method=="dada2", Filt=="extra", Rare=="rarefied")
deblur.basic.r <- kruskal.function(all.mock.alpha, Method=="deblur", Filt=="basic", Rare=="rarefied")
deblur.standard.r <- kruskal.function(all.mock.alpha, Method=="deblur", Filt=="standard", Rare=="rarefied")
deblur.extra.r <- kruskal.function(all.mock.alpha, Method=="deblur", Filt=="extra", Rare=="rarefied")
vsearch.basic.r <- kruskal.function(all.mock.alpha, Method=="vsearch", Filt=="basic", Rare=="rarefied")
vsearch.standard.r <- kruskal.function(all.mock.alpha, Method=="vsearch", Filt=="standard", Rare=="rarefied")
vsearch.extra.r <- kruskal.function(all.mock.alpha, Method=="vsearch", Filt=="extra", Rare=="rarefied")

krustal.results.OTUs <- rbind(dada2.basic.u, dada2.basic.r, dada2.standard.u, dada2.standard.r, dada2.extra.u, dada2.extra.r,
                              deblur.basic.u, deblur.basic.r, deblur.standard.u, deblur.standard.r, deblur.extra.u, deblur.extra.r,
                              vsearch.basic.u, vsearch.basic.r, vsearch.standard.u, vsearch.standard.r, vsearch.extra.u, vsearch.extra.r)
krustal.results.OTUs$qval <- p.adjust(p = krustal.results.OTUs$val, method = "BH")

rm(dada2.basic.u, dada2.basic.r, dada2.standard.u, dada2.standard.r, dada2.extra.u, dada2.extra.r,
   deblur.basic.u, deblur.basic.r, deblur.standard.u, deblur.standard.r, deblur.extra.u, deblur.extra.r,
   vsearch.basic.u, vsearch.basic.r, vsearch.standard.u, vsearch.standard.r, vsearch.extra.u, vsearch.extra.r)

write.table(krustal.results.OTUs, file = "~/Repos/tidybug/data/text_tables/krustal.results.OTUs.csv", quote=FALSE, row.names = FALSE)


