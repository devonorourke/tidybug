library(tidyverse)
library(reshape2)
library(vegan)
library(phyloseq)

## theme for plot
theme_devon <- function () { 
  theme_bw(base_size=12, base_family="Avenir") %+replace% 
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill="transparent", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA)
    )
}

## ellipse function:
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}


## import data:
df <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz")
meta <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/metadata/meta_nomock.csv/")
meta <- meta %>% select(-StudyID, -SampleType, -Library)
df <- merge(df, meta)
rm(meta)
fox <- df %>% filter(Site=="FOX")
rm(df)
## only selected months
SelectMonths <- c("April", "May", "September", "October")
fox <- fox %>% filter(MonthStart %in% SelectMonths)
rm(SelectMonths)

## rarefied plots
rare.nmds.function <- function(data, filter_exp, filter_exp2, betatest, binaryval) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  Method <- tmp.df %>% distinct(Method)
  Filt <- tmp.df %>% distinct(Filt)
  BetaType <- binaryval
  tmp.meta <- tmp.df %>% distinct(SeqID, MonthStart, WOY)
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  SeqID <- row.names(tmp.mat)
  tmp.mat$SeqID <- NULL
  tmp.otu <- otu_table(tmp.mat, taxa_are_rows = FALSE)
  tmp.phy <- phyloseq(tmp.otu)
  set.seed(seed = 100)
  tmp.phy <- rarefy_even_depth(tmp.phy, sample.size = 5000, rngseed = TRUE, replace = FALSE, trimOTUs = TRUE)
  tmp.mat2 <- as(otu_table(tmp.phy), "matrix")
  tmp.betadiv <- vegdist(tmp.mat2, betatest, binary = binaryval)
  tmp.nmds <- metaMDS(tmp.betadiv, k=3,trymax = 250, autotransform = FALSE)
  nmds.df <- data.frame(tmp.nmds$points, Method, Filt, BetaType)
  nmds.df$SeqID <- row.names(nmds.df)
  row.names(nmds.df) <- NULL
  nmds.df <- merge(nmds.df, tmp.meta)
  nmds.df$MonthStart <- as.factor(nmds.df$MonthStart)
  nmds.hull <- nmds.df %>% group_by(MonthStart) %>% slice(chull(MDS1,MDS2))
  groups <- nmds.df$MonthStart
  plot(tmp.nmds)
  ord <- ordiellipse(tmp.nmds, groups=groups)
  nmds.ellipse <- data.frame()
  for(g in levels(nmds.df$MonthStart)) {
    nmds.ellipse <- rbind(nmds.ellipse, cbind(as.data.frame(with(nmds.df[nmds.df$MonthStart==g,],
                                                     veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale))) ,group=g))
  }
  colnames(nmds.ellipse) <- c("MDS1", "MDS2", "MonthStart")

  dada2.extra.unrarefied <- ggplot(nmds.df, aes(x=MDS1, y=MDS2, label=WOY, color=MonthStart)) +
    geom_text(size=3.5) +
    geom_polygon(data=nmds.hull, size=0.05, alpha=0.05, aes(color=MonthStart, fill=MonthStart)) +
    geom_path(data=nmds.ellipse, linetype="dashed", aes(x=MDS1, y=MDS2, color=MonthStart), inherit.aes = FALSE) +
    labs(x="mds1", y="mds2", title="", label="", color="") +
    xlim(-0.8, 0.8) +
    ylim(-0.8, 0.8) +
    theme_devon() +
    theme(legend.position = "none", axis.text = element_text(size = 7))
}

### ---- redo with alternative input... code the function 

dada2.basic.rarefied <- rare.nmds.function(fox, Method=="dada2", Filt=="basic", "bray", FALSE)
dada2.extra.rarefied <- rare.nmds.function(fox, Method=="dada2", Filt=="extra", "bray", FALSE)
vsearch.basic.rarefied <- rare.nmds.function(fox, Method=="vsearch", Filt=="basic", "bray", FALSE)
vsearch.extra.rarefied <- rare.nmds.function(fox, Method=="vsearch", Filt=="extra", "bray", FALSE)

plot_grid(dada2.basic.rarefied, dada2.extra.rarefied, vsearch.basic.rarefied, vsearch.extra.rarefied,
          ncol=2)
