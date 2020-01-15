library(tidyverse)
library(reshape2)
library(vegan)

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
df <- read_csv("~/Repos/tidybug/data/text_tables/all.filtmethods.df.csv.gz")
mock <- df %>% filter(SampleType == "mock")
mock$Labeler <- paste(mock$Method, mock$Filt, sep="-")
mock$bigID <- paste(mock$SeqID, mock$Method, mock$Filt, sep="-")
rm(df)

## generate data to plot DISPERSIONS for each Method and Filt group with RAREFIED data:
rrarewithdrop <- function(x, sample) {
    rrarefy(x[rowSums(x) >= sample, , drop = FALSE], sample)
}

mock_betadisp.function <- function(betatest, binaryval) {
  tmp.meta <- mock %>% distinct(bigID, Labeler)
  groups <- as.character(tmp.meta$Labeler)
  tmp.mat <- dcast(mock, bigID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$bigID
  tmp.mat$bigID <- NULL
  tmp.mat2 <- data.frame(rrarewithdrop(tmp.mat, 5000))
  rm(tmp.mat)
  tmp.betadiv <- vegdist(tmp.mat2, betatest, binary = binaryval)
  tmp.disp <- betadisper(tmp.betadiv, groups, type = c("median"))
  rm(tmp.betadiv)
  tmp.disp_sco <- as.data.frame(scores(tmp.disp, display = "sites"))
  tmp.disp_sco$bigID <- row.names(tmp.disp_sco)
  row.names(tmp.disp_sco) <- NULL
  colnames(tmp.disp_sco) <- c("pc1end", "pc2end", "sample")
  tmp.disp_sco <- tmp.disp_sco %>% separate(., col="sample", into=c("SeqID", "Method", "Filt"), sep="-")
  tmp.disp_sco$Group <- paste(tmp.disp_sco$Method, tmp.disp_sco$Filt, sep = "-")
  tmp.disp_centroids <- as.data.frame(scores(tmp.disp, display="centroids"))
  tmp.disp_centroids$SeqID <- row.names(tmp.disp_centroids)
  row.names(tmp.disp_centroids) <- NULL
  colnames(tmp.disp_centroids) <- c("pc1start", "pc2start", "Group")
  tmp.disp.plot <- merge(tmp.disp_sco, tmp.disp_centroids)
  tmp.disp.plot$BetaTest <- betatest
  tmp.disp.plot
}


## gather data for plot
dice.betadisp_df <- mock_betadisp.function("bray", FALSE)
dice.betadisp_df$BetaTest <- gsub("bray", "dice", dice.betadisp_df$BetaTest)
bray.betadisp_df <- mock_betadisp.function("bray", TRUE)
mori.betadisp_df <- mock_betadisp.function("morisita", FALSE)
mock.betadisp <- rbind(dice.betadisp_df, bray.betadisp_df, mori.betadisp_df)

## set levels for plot
mock.betadisp$BetaTest <- factor(mock.betadisp$BetaTest, 
                                 levels = c("dice", "bray", "morisita"))

## colors and shapes for plot
pal3 <- c('#9f9244', '#6c42b8', '#628a47')

## plot dispersions; points are samples
## save as 14_mock_betadispersions; export at 1000x600
ggplot(mock.betadisp, aes(x=pc1end, y=pc2end, label=SeqID, color=Filt)) + 
  geom_point(data=mock.betadisp, aes(x=pc1end, y=pc2end), size=2, alpha=0.8) +
  scale_color_manual(values=pal3) + 
  scale_x_continuous(breaks = c(-0.2, 0, 0.2)) +
  geom_segment(data=mock.betadisp, aes(x=pc1end, y=pc2end, xend=pc1start, yend=pc2start), alpha=0.2) +
  facet_grid(BetaTest ~ Method + Filt) +
  labs(title="", color="", x="PCoA1", y="PCoA2", shape="", color="") +
  theme_devon() +
  theme(legend.position = "top")

## stats for dispersion
mock_betadisp_stats.function <- function(betatest, binaryval) {
  tmp.meta <- mock %>% distinct(bigID, Labeler)
  groups <- as.character(tmp.meta$Labeler)
  tmp.mat <- dcast(mock, bigID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$bigID
  tmp.mat$bigID <- NULL
  tmp.mat2 <- data.frame(rrarewithdrop(tmp.mat, 5000))
  rm(tmp.mat)
  tmp.betadiv <- vegdist(tmp.mat2, betatest, binary = binaryval)
  tmp.disp <- betadisper(tmp.betadiv, groups, type = c("median"))
}

## dispersion calculation; note we're calculating this two different ways...
## the first is following this conversation: https://stat.ethz.ch/pipermail/r-sig-ecology/2010-September/001524.html
## ... which is testing for dispersion among individual cells (one-way test where cells are treated as level of single factor)
## the second partitions the group into the main effects, but it's probably not worthwhile...
## all generate non-significant results regardless, so our differences in ADONIS are apparently not because of dipsersion

betadisp.dice <- mock_betadisp_stats.function("bray", TRUE)
anova(betadisp.dice)  ## anova on group treated as single factor
bdips.dice.dist <- data.frame(betadisp.dice$distances, betadisp.dice$group)
colnames(bdips.dice.dist) <- c("distances", "group")
bdips.dice.dist$sample <- row.names(bdips.dice.dist)
row.names(bdips.dice.dist) <- NULL
bdips.dice.dist <- bdips.dice.dist %>% separate(., col = "group", into = c("Method", "Filt"), sep = "-")
bdips.dice.anova <- aov(distances ~ Method * Filt, data = bdips.dice.dist)
summary(bdips.dice.anova) ## anova on group treated as separate factors

betadisp.bray <- mock_betadisp_stats.function("bray", FALSE)
anova(betadisp.bray)
bdips.bray.dist <- data.frame(betadisp.bray$distances, betadisp.bray$group)
colnames(bdips.bray.dist) <- c("distances", "group")
bdips.bray.dist$sample <- row.names(bdips.bray.dist)
row.names(bdips.bray.dist) <- NULL
bdips.bray.dist <- bdips.bray.dist %>% separate(., col = "group", into = c("Method", "Filt"), sep = "-")
bdips.bray.anova <- aov(distances ~ Method * Filt, data = bdips.bray.dist)
summary(bdips.bray.anova)


betadisp.mori <- mock_betadisp_stats.function("morisita", FALSE)
anova(betadisp.mori)
bdips.mori.dist <- data.frame(betadisp.mori$distances, betadisp.mori$group)
colnames(bdips.mori.dist) <- c("distances", "group")
bdips.mori.dist$sample <- row.names(bdips.mori.dist)
row.names(bdips.mori.dist) <- NULL
bdips.mori.dist <- bdips.mori.dist %>% separate(., col = "group", into = c("Method", "Filt"), sep = "-")
bdips.mori.anova <- aov(distances ~ Method * Filt, data = bdips.mori.dist)
summary(bdips.mori.anova)

