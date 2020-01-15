library(tidyverse)
library(scales)
library(FSA)
library(twosamples)
library(grid)

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

## generate 3 color palette to distinguish between filtering pipelines:
pal3 <- c('#9f9244', '#6c42b8', '#628a47', '#a9d190')


################################################################################
## generate dataset for plot and analyses
## examining the per-ASV distribution of read depth and frequency of occurrence
################################################################################

## import data and select mock samples:
## retain only true samples for comparison and create dataframe for plots and stats
df <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz") %>% 
  filter(SampleType == "sample") %>% 
  mutate(log2reads = floor(log2(Reads))) %>% 
  group_by(Filt, Method) %>% 
  mutate(ReadSums=sum(Reads)) %>% 
  mutate(pReads=100*(Reads/ReadSums)) %>% 
  arrange(Reads) %>% 
  mutate(cd=cumsum(pReads))

################################################################################
## 2 plot: histogram
################################################################################
breakvals <- floor(min(df$log2reads)):ceiling(max(df$log2reads))

## first, make a histogram of ASV-read depth in log-scale
df2 <- df %>% 
  group_by(Filt, Method) %>% 
  count(log2bucket=log2reads, name = "Observations") %>% 
  mutate(nObservations=sum(Observations)) %>% 
  mutate(pObservations=100*(Observations/nObservations)) %>% 
  mutate(ReadScale=2^log2bucket) %>% 
  arrange(log2bucket) %>% 
  mutate(cd=cumsum(pObservations))

## set the levels for hist plot
df2$Filt <- factor(df2$Filt,levels = c("basic", "standard", "extra"))
df2$Method <- factor(df2$Method, levels = c("dada2", "deblur", "vsearch"))

## plot histogram
p <- ggplot(df2, aes(x=as.numeric(ReadScale),y=pObservations,fill=Method)) + 
  geom_bar(stat='identity') +
  facet_grid(Filt~Method) +
  scale_x_continuous(trans="log2", labels = comma) +
  labs(x="\nSequences per ASV", y="fraction of ASVs detected\n", fill="") +
  theme_devon() +
  scale_fill_manual(values=pal3) +
  scale_color_manual(values = c("firebrick", "black", "dodgerblue")) +
  theme(legend.position = "none",
        strip.text.x = element_text(size=13), strip.text.y = element_text(size=13),
        axis.text.x = element_text(angle=45, hjust = 1)) +
  geom_text(aes(x=ReadScale,y=pObservations,label=Vals,color=Filt),
            data=data.frame(
    Method=rep(c("dada2", "deblur", "vsearch"),3),
    Filt=c(rep("basic",3), rep("standard",3), rep("extra",3)),
    ReadScale=rep(262000,9),
    pObservations=rep(30,9),
    Vals=c("a", "b", "b", "a", "b", "b", "a", "a", "b")))
  
g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c("firebrick","dodgerblue","white")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid.draw(g)


################################################################################
## now run stats with this table to determine if distributions of Reads vary by denoiser (for each filtering method)
################################################################################

## read data are not normally distributed...
shapiro_function1 <- function(filtmethod){
  set.seed(10)  
  da_tmp <- sample(df %>% filter(Filt==filtmethod, Method=="dada2") %>% pull(Reads), 5000, replace = TRUE)
  db_tmp <- sample(df %>% filter(Filt==filtmethod, Method=="deblur") %>% pull(Reads), 5000, replace = TRUE)
  vs_tmp <- sample(df %>% filter(Filt==filtmethod, Method=="vsearch") %>% pull(Reads), 5000, replace = TRUE)
  print(shapiro.test(da_tmp))
  print(shapiro.test(db_tmp))
  print(shapiro.test(vs_tmp))
}

shapiro_function1("basic")    ## NOT normal for any dataset
shapiro_function1("standard")    ## NOT normal for any dataset
shapiro_function1("extra")    ## NOT normal for any dataset

## Tests to compare if empirical cumulative density functions vary between two samples 
## 3 tests explored - each one more stringent from local extremes to more global differences among entire distribution
## see: https://github.com/cdowd/twosamples for details about each test

## subsample dwn to fewer observations -- vsearch/extra has fewest with > 17,000 so we'll use 10,000
## note, this code takes ~20 min to execute per test... consider exporting/rewriting to cluster
## other ecdf functions availabe in package, but this Wasserstein test gives best balance of outlier detection across entire distribution
guano_ecdf_function <- function(filtmethod){
  dada2_tmp = df %>% filter(Filt==filtmethod, Method=="dada2") %>% pull(Reads) %>% sample(., 10000, replace=TRUE)
  deblur_tmp = df %>% filter(Filt==filtmethod, Method=="deblur") %>% pull(Reads) %>% sample(., 10000, replace=TRUE)
  vsearch_tmp =  df %>% filter(Filt==filtmethod, Method=="vsearch") %>% pull(Reads) %>% sample(., 10000, replace=TRUE)
  print(twosamples::wass_test(dada2_tmp, deblur_tmp, nboots = 500))
  print(twosamples::wass_test(dada2_tmp, vsearch_tmp, nboots = 500))
  print(twosamples::wass_test(deblur_tmp, vsearch_tmp, nboots = 500))
}

guano_ecdf_function("basic")
guano_ecdf_function("standard")
guano_ecdf_function("extra")

################################################################################
## unused code:
################################################################################

## are Observation data normal?
shapiro_function2 <- function(filtmethod){
  set.seed(10)  
  da_tmp <- sample(df2 %>% filter(Filt==filtmethod, Method=="dada2") %>% pull(Observations))
  db_tmp <- sample(df2 %>% filter(Filt==filtmethod, Method=="deblur") %>% pull(Observations))
  vs_tmp <- sample(df2 %>% filter(Filt==filtmethod, Method=="vsearch") %>% pull(Observations))
  print(shapiro.test(da_tmp))
  print(shapiro.test(db_tmp))
  print(shapiro.test(vs_tmp))
}

shapiro_function2("basic")    ## NOT normal for any dataset
shapiro_function2("standard")    ## NOT normal for any dataset
shapiro_function2("extra")    ## NOT normal for any dataset
## apply non-parametric tests given that data types are normally distributed

### what are the mean, median, and stdev of each read count distribution?
df %>% group_by(Method, Filt) %>% summarise(Mean=mean(Reads), Median=median(Reads), SD=sd(Reads))
### what about the observations for a distribution of 2^read? i.e. the number of times an ASV has 1, 2, 4, 8, 16, etc. reads?
df2 %>% group_by(Method, Filt) %>% summarise(Mean=mean(Observations), Median=median(Observations), SD=sd(Observations))

## plot cdf of Observations (not reads)
## this may be a confusing graphic because the statistic is on the distribution of READS
ggplot(df2, aes(x=as.numeric(ReadScale), y=cd, color=Method, fill=Method)) +
  geom_line(size=1) +
  geom_point(size=1.25) +
  facet_grid(Filt~Method) +
  scale_x_continuous(trans="log2", labels = comma, limits = c(1,2048)) +
  labs(x="\nSequences per ASV", y="fraction of library\n", fill="") +
  theme_devon() +
  scale_color_manual(values=pal3) +
  geom_hline(yintercept = 50, color="red", size=0.5, alpha=0.25) +
  theme(legend.position = "none",
        strip.text.x = element_text(size=13), strip.text.y = element_text(size=13),
        axis.text.x = element_text(angle=45, hjust = 1))

## run Kruskal-Wallis to test for rank ordered mean differences in # Reads between each Method (one per filter parameter)
## apply post hoc Dunn's Test for pairwise differences
guano_kw_function <- function(filtmethod){
  print(kruskal.test(Reads ~ as.factor(Method), data = df %>% filter(Filt == filtmethod)))
  dunnTest(Reads ~ as.factor(Method), data = df %>% filter(Filt == filtmethod), method = "bh")
}

guano_kw_function("basic")    ## kw sig diff; all pairwise comps different at p < 0.001
guano_kw_function("standard") ## kw sig diff; all pairwise comps different at p < 0.001
guano_kw_function("extra")    ## kw sig diff; all pairwise comps different at p < 0.001


library(formattable)
text_df <- data.frame(comparison=rep(c("dada2 : deblur", "dada2 : vsearch", "deblur : vsearch"),3),
                      difference=c(TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,FALSE,FALSE,TRUE))

p2 <- formattable(text_df, list(
  difference = formatter("span",
                         style = x ~ style(color = ifelse(x, "green", "red")),
                         x ~ icontext(ifelse(x, "ok", "remove"), ifelse(x, "Yes", "No"))
  )))
