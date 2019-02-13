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
  tmp.mat <- dcast(tmp.df, SeqID ~ Alias, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  SeqID <- row.names(tmp.mat)
  tmp.mat$SeqID <- NULL
  tmp.otu <- otu_table(tmp.mat, taxa_are_rows = FALSE)
  tmp.phy <- phyloseq(tmp.otu)
  set.seed(seed = 100)
  tmp.phy <- rarefy_even_depth(tmp.phy, sample.size = 5000, rngseed = TRUE, replace = FALSE, trimOTUs = TRUE)
  tmp.mat2 <- as(otu_table(tmp.phy), "matrix")
  tmp.mat2 <- as.data.frame(tmp.mat2)
  tmp.betadiv <- vegdist(tmp.mat2, betatest, binary = binaryval)
  tmp.pcoa <- capscale(formula = tmp.betadiv~1, distance = "bray")
  pcoa.df <- data.frame(scores(tmp.pcoa,display = "sites"), Method, Filt, BetaType)
  pcoa.df$SeqID <- row.names(pcoa.df)
  row.names(pcoa.df) <- NULL
  pcoa.df <- merge(pcoa.df, tmp.meta)
  #pcoa.df$MonthStart <- as.factor(pcoa.df$MonthStart)
  pcoa.df$WOY <- as.factor(pcoa.df$WOY)
  #pcoa.hull <- pcoa.df %>% group_by(MonthStart) %>% slice(chull(MDS1,MDS2))
  pcoa.hull <- pcoa.df %>% group_by(WOY) %>% slice(chull(MDS1,MDS2))
  #groups <- pcoa.df$MonthStart
  groups <- pcoa.df$WOY
  plot(tmp.pcoa)
  ord <- ordiellipse(tmp.pcoa, groups=groups)
  pcoa.ellipse <- data.frame()
  #for(g in levels(pcoa.df$MonthStart)) {
  #  pcoa.ellipse <- rbind(pcoa.ellipse, cbind(as.data.frame(with(pcoa.df[pcoa.df$MonthStart==g,],
  #                                                               veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale))) ,group=g))
  #}
  for(g in levels(pcoa.df$WOY)) {
      pcoa.ellipse <- rbind(pcoa.ellipse, cbind(as.data.frame(with(pcoa.df[pcoa.df$WOY==g,],
                                                                   veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale))) ,group=g))
    }
    
  screeplot(tmp.pcoa)
  eigvals <- data.frame(tmp.pcoa$CA$eig)
  colnames(eigvals) <- "eigenvalue"
  eigvals$MDSaxis <- row.names(eigvals)
  row.names(eigvals) <- NULL
  eigmds1 <- eigvals %>% filter(MDSaxis=="MDS1") %>% select(eigenvalue)
  eigmds1 <- round(eigmds1, 2)
  eigmds1 <- as.character(eigmds1$eigenvalue)
  eigmds2 <- eigvals %>% filter(MDSaxis=="MDS2") %>% select(eigenvalue)
  eigmds2 <- round(eigmds2, 2)
  eigmds2 <- as.character(eigmds2$eigenvalue)
  
  #colnames(pcoa.ellipse) <- c("MDS1", "MDS2", "MonthStart")
  colnames(pcoa.ellipse) <- c("MDS1", "MDS2", "WOY")
  
  #ggplot(pcoa.df, aes(x=MDS1, y=MDS2, label=WOY, color=MonthStart)) +
  ggplot(pcoa.df, aes(x=MDS1, y=MDS2, label=WOY, color=MonthStart)) +
    geom_text(size=3.5) +
    geom_polygon(data=pcoa.hull, size=0.05, alpha=0.05, aes(color=MonthStart, fill=MonthStart)) +
    geom_path(data=pcoa.ellipse, linetype="dashed", aes(x=MDS1, y=MDS2, color=MonthStart), inherit.aes = FALSE) +
    labs(x=paste0("mds1 [",eigmds1,"] %"),
         y=paste0("mds2 [",eigmds2,"] %"),
         title="", label="", color="") +
    #xlim(-0.8, 0.8) +
    #ylim(-0.8, 0.8) +
    theme_devon() +
    theme(legend.position = "none", axis.text = element_text(size = 7))
}


dada2.basic.rarefied <- rare.nmds.function(fox, Method=="dada2", Filt=="basic", "bray", FALSE)
dada2.extra.rarefied <- rare.nmds.function(fox, Method=="dada2", Filt=="extra", "bray", FALSE)
vsearch.basic.rarefied <- rare.nmds.function(fox, Method=="vsearch", Filt=="basic", "bray", FALSE)
vsearch.extra.rarefied <- rare.nmds.function(fox, Method=="vsearch", Filt=="extra", "bray", FALSE)

plot_grid(dada2.basic.rarefied, dada2.extra.rarefied, vsearch.basic.rarefied, vsearch.extra.rarefied,
          ncol=2)
