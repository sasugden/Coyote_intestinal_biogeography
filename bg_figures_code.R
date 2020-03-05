#### (done) FIGURE 1 - ALPHA DIVERSITY LINE GRAPHS ###############
# Ensure intestinal sites will plot in order.
desired_order <- c("duodenum","jejunum", "ileum", "caecum", "ascending.colon", "descending.colon", "feces")
alpha.summary$Type <- factor(alpha.summary$Type, levels=desired_order)

# Rename sites for plotting.
levels(alpha.summary$Type)[levels(alpha.summary$Type)=="duodenum"] <- "Duodenum"
levels(alpha.summary$Type)[levels(alpha.summary$Type)=="jejunum"] <- "Jejunum"
levels(alpha.summary$Type)[levels(alpha.summary$Type)=="ileum"] <- "Ileum"
levels(alpha.summary$Type)[levels(alpha.summary$Type)=="caecum"] <- "Caecum"
levels(alpha.summary$Type)[levels(alpha.summary$Type)=="ascending.colon"] <- "Ascending\nColon"
levels(alpha.summary$Type)[levels(alpha.summary$Type)=="descending.colon"] <- "Descending\nColon"
levels(alpha.summary$Type)[levels(alpha.summary$Type)=="feces"] <- "Feces"

plot1 <- ggplot(alpha.summary, aes(x=Type, y=mean.Observed)) +
  geom_errorbar(aes(ymin=mean.Observed-sd.Observed, ymax=mean.Observed+sd.Observed), width=.1) +
  geom_line(group=1) + geom_point() +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(color="black", size=9),
        axis.title=element_text(color="black", size=9),
        plot.tag=element_text(face="bold"),
        panel.spacing=unit(1, "lines"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position="none") +
  xlab(NULL) +
  ylab("ASV Richness\n") +
  labs(tag="a")

plot2 <- ggplot(alpha.summary, aes(x=Type, y=mean.Shannon)) +
  geom_errorbar(aes(ymin=mean.Shannon-sd.Shannon, ymax=mean.Shannon+sd.Shannon), width=.1) +
  geom_line(group=1) + geom_point() +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(color="black", size=9),
        axis.title=element_text(color="black", size=9),
        plot.tag=element_text(face="bold"),
        panel.spacing=unit(1, "lines"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position="none") +
  xlab(NULL) +
  ylab("Shannon index\n") +
  labs(tag="b")

plot3 <- ggplot(alpha.summary, aes(x=Type, y=cv.Observed)) +
  geom_line(group=1) + geom_point() +
  theme_bw() +
  scale_y_continuous(limits=c(0,65)) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, color="black", size=8),
        axis.text.y=element_text(color="black", size=9),
        axis.title=element_text(color="black", size=9),
        plot.tag=element_text(face="bold"),
        panel.spacing=unit(1, "lines"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position="none") +
  xlab(NULL) +
  ylab("ASV Richness\nCoefficient of variation") +
  labs(tag="c")

plot4 <- ggplot(alpha.summary, aes(x=Type, y=cv.Shannon)) +
  geom_line(group=1) + geom_point() +
  theme_bw() +
  scale_y_continuous(limits=c(0,40)) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, color="black", size=8),
        axis.text.y=element_text(color="black", size=8),
        axis.title=element_text(color="black", size=9),
        plot.tag=element_text(face="bold"),
        panel.spacing=unit(1, "lines"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position="none") +
  xlab(NULL) +
  ylab("Shannon index\nCoefficient of variation") +
  labs(tag="d")

plot1 <- ggplotGrob(plot1)
plot2 <- ggplotGrob(plot2)
plot3 <- ggplotGrob(plot3)
plot4 <- ggplotGrob(plot4)
maxWidth = grid::unit.pmax(plot1$widths[2:5], plot2$widths[2:5], plot3$widths[2:5], plot4$widths[2:5])
plot1$widths[2:5] <- as.list(maxWidth)
plot2$widths[2:5] <- as.list(maxWidth)
plot3$widths[2:5] <- as.list(maxWidth)
plot4$widths[2:5] <- as.list(maxWidth)

grid.arrange(gA, plot2, plot3, plot4, ncol=2, heights=c(0.9,1))

gs <- list(plot1, plot2, plot3, plot4)
Fig_2 <- arrangeGrob(grobs=gs, heights=c(0.9,1), ncol=2)

ggsave("Fig_1.jpg", plot = Fig_1, width=6.5, height=4, units="in", dpi=300)


#### (done) FIGURE 2 - PHYLUM RELATIVE ABUNDANCE ##############
# (a) Grouped by intestinal segment (small vs. large) 
# Agglomerate taxa to class level and merge by segment.
temp <- tax_glom(pseq.rarefied, taxrank="Class")
variable1 = as.character(get_variable(temp, "Segment"))
temp <- merge_samples(temp, "Segment")
temp.data <- as.data.frame(sample_data(temp))

# Replace data that gets overwritten in the merge.
temp.data$Segment <- rownames(temp.data)
sample_data(temp) <- temp.data

# Convert to relative abundance.
temp <- transform_sample_counts(temp, function(x) 100*x/sum(x))

# Rename all low-abundance or low-prevalence phyla to "Other"
psmelt.data <- psmelt(temp)
psmelt.data$Phylum <- as.character(psmelt.data$Phylum)
psmelt.data$Phylum[psmelt.data$Phylum %in% c("Acidobacteria",
                           "Deferribacteres", "Euryarchaeota","Fibrobacteres",
                           "Planctomycetes","Spirochaetes","Tenericutes","Verrucomicrobia")] <- "Other"

# Rename classes based on their phylum.
psmelt.data$Class <- as.character(psmelt.data$Class)
psmelt.data$Class[psmelt.data$Phylum=="Other"] <- "Taxa <1% Abundance"
psmelt.data$Class[psmelt.data$Phylum=="Proteobacteria"] <- "Proteobacteria"
psmelt.data$Class[psmelt.data$Class=="NA"] <- "Taxa <1% Abundance"
psmelt.data$Class[psmelt.data$Class %in% c("Bacteroidia", "Flavobacteriia", "Sphingobacteriia", "Cytophagia")] <- "Bacteroidetes"

# Convert Phylum and Class back to factors.
psmelt.data$Phylum <- as.factor(psmelt.data$Phylum)
psmelt.data$Class <- as.factor(psmelt.data$Class)

# Set plotting order for phyla.
desired_order <- c("Other","Actinobacteria","Proteobacteria","Fusobacteria","Bacteroidetes", "Firmicutes")
psmelt.data$Phylum <- factor(as.character(psmelt.data$Phylum), levels=desired_order)
psmelt.data <- psmelt.data[order(psmelt.data$Phylum, psmelt.data$Class),]

# Determine the number of colors needed.
colorCount <- length(unique(psmelt.data$Class))

# Set plotting order for classes.
desired_order <- c("Taxa <1% Abundance", "Actinobacteria",
                   "Proteobacteria", "Fusobacteriia", "Bacteroidetes",
                   "Bacilli", "Clostridia", "Erysipelotrichia", "Negativicutes")
psmelt.data$Class <- factor(as.character(psmelt.data$Class), levels=desired_order)

# Define colors for different taxa.
myColors <- c("#999999", # Other
              "#CC0000", # Actinobacteria
              "#0000FF", # Bacteroidetes
              "#99CC00", # Proteobacteria
              "#FF6699", # Fusobacteria
              "#FFCC00","#CC9900","#FFCC66","#996600") # Firmicutes
names(myColors) <- levels(psmelt.data$Class)

# Set plotting order for segments.
desired_order <- c("small.intestine","large.intestine")
psmelt.data$Segment <- factor(psmelt.data$Segment, levels=desired_order)

# Rename additional items for figure.
levels(psmelt.data$Segment)[levels(psmelt.data$Segment)=="small.intestine"] <- "Small\nIntestine"
levels(psmelt.data$Segment)[levels(psmelt.data$Segment)=="large.intestine"] <- "Large\nIntestine"
levels(psmelt.data$Class)[levels(psmelt.data$Class)=="Taxa <1% Abundance"] <- "Taxa <1% Abundance       "

Fig_2.legend <- cowplot::get_legend(ggplot(psmelt.data, aes(x=Segment, y=Abundance, fill=Class)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name="Class", values=myColors) +
  scale_y_continuous(expand=c(0,0), limits=c(0, 102)) +
  scale_x_discrete(labels=function(Segment) str_wrap(Segment, width=3)) +
  guides(fill=guide_legend(ncol=2)) +
  theme_bw() +
  theme(panel.spacing=unit(1, "lines"), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(color="black", size = 8, vjust=1),
        axis.title=element_text(size=9, color="black"),
        strip.background=element_blank(), 
        strip.text=element_text(size=9, color="black"),
        plot.tag=element_text(face="bold"),
        legend.title=element_blank(),
        legend.position="bottom",
        legend.text=element_text(size=8),
        legend.key.size=unit(0.5, "line"),
        legend.justification=0.5) +
  ylab("Mean relative abundance (%)") +
  xlab(NULL) +
  labs(tag="a"))

Fig_2a <- ggplot(psmelt.data, aes(x=Segment, y=Abundance, fill=Class)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name="Class", values=myColors) +
  scale_y_continuous(expand=c(0,0), limits=c(0, 102)) +
  scale_x_discrete(labels=function(Segment) str_wrap(Segment, width=3)) +
  guides(fill=guide_legend(ncol=2)) +
  theme_bw() +
  theme(panel.spacing=unit(1, "lines"), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(color="black", size = 8, vjust=0.5, hjust=1, angle=90),
        axis.text.y=element_text(color="black", size = 8, vjust=1),
        axis.title=element_text(size=9, color="black"),
        strip.background=element_blank(), 
        strip.text=element_text(size=9, color="black"),
        plot.tag=element_text(face="bold"),
        plot.margin=unit(c(5.5,5.5,17,5.5), "points"),
        legend.title=element_blank(),
        legend.position="none",
        legend.text=element_text(size=8),
        legend.key.size=unit(0.5, "line"),
        legend.justification=0.5) +
  ylab("Mean relative abundance (%)") +
  xlab(NULL) +
  labs(tag="a")

# Grouped by intestinal site - first agglomerate by class and group by site.
temp <- tax_glom(pseq.rarefied, taxrank="Class")
variable1 = as.character(get_variable(temp, "Type"))
temp <- merge_samples(temp, "Type")
temp.data <- as.data.frame(sample_data(temp))

# Replace data that gets overwritten in the merge.
temp.data$Type <- forcats::fct_inorder(rownames(temp.data))
sample_data(temp) <- temp.data

# Convert to relative abundance.
temp <- transform_sample_counts(temp, function(x) 100*x/sum(x))

# Melt the data.
psmelt.data <- psmelt(temp)

# Replace low-abundance phyla with "Other"
psmelt.data$Phylum <- as.character(psmelt.data$Phylum)
psmelt.data$Phylum[psmelt.data$Phylum %in% c("Acidobacteria",
                           "Deferribacteres", "Euryarchaeota","Fibrobacteres",
                           "Planctomycetes","Spirochaetes","Tenericutes","Verrucomicrobia")] <- "Other"

# Rename class category to the phylum level for plotting.
psmelt.data$Class <- as.character(psmelt.data$Class)
psmelt.data$Class[psmelt.data$Phylum=="Other"] <- "Taxa <1% Abundance"
psmelt.data$Class[psmelt.data$Class=="NA"] <- "Taxa <1% Abundance"
psmelt.data$Class[psmelt.data$Phylum=="Proteobacteria"] <- "Proteobacteria"
psmelt.data$Class[psmelt.data$Class %in% c("Bacteroidia", "Flavobacteriia", "Sphingobacteriia", "Cytophagia")] <- "Bacteroidetes"

# Convert variables back to factors.
psmelt.data$Phylum <- as.factor(psmelt.data$Phylum)
psmelt.data$Class <- as.factor(psmelt.data$Class)

# Set plotting order.
desired_order <- c("Other","Actinobacteria","Proteobacteria","Fusobacteria","Bacteroidetes", "Firmicutes")
psmelt.data$Phylum <- factor(as.character(psmelt.data$Phylum), levels=desired_order)
psmelt.data <- psmelt.data[order(psmelt.data$Phylum, psmelt.data$Class),]
colorCount <- length(unique(psmelt.data$Class))
desired_order <- c("Taxa <1% Abundance", "Actinobacteria",
                   "Proteobacteria", "Fusobacteriia", "Bacteroidetes",
                   "Bacilli", "Clostridia", "Erysipelotrichia", "Negativicutes")
psmelt.data$Class <- factor(as.character(psmelt.data$Class), levels=desired_order)

# Assign colors to taxa.
myColors <- c("#999999", # Other
              "#CC0000", # Actinobacteria
              "#0000FF", # Bacteroidetes
              "#99CC00", # Proteobacteria
              "#FF6699", # Fusobacteria
              "#FFCC00","#CC9900","#FFCC66","#996600") # Firmicutes
names(myColors) <- levels(psmelt.data$Class)

# Set plotting order.
desired_order <- c("duodenum","jejunum", "ileum", "caecum", "ascending.colon", "descending.colon", "feces")
psmelt.data$Type <- factor(psmelt.data$Type, levels=desired_order)

# Rename variables for plotting.
levels(psmelt.data$Type)[levels(psmelt.data$Type)=="duodenum"] <- "Duodenum"
levels(psmelt.data$Type)[levels(psmelt.data$Type)=="jejunum"] <- "Jejunum"
levels(psmelt.data$Type)[levels(psmelt.data$Type)=="ileum"] <- "Ileum"
levels(psmelt.data$Type)[levels(psmelt.data$Type)=="caecum"] <- "Caecum"
levels(psmelt.data$Type)[levels(psmelt.data$Type)=="ascending.colon"] <- "Ascending\nColon"
levels(psmelt.data$Type)[levels(psmelt.data$Type)=="descending.colon"] <- "Descending\nColon"
levels(psmelt.data$Type)[levels(psmelt.data$Type)=="feces"] <- "Feces"

# Create plot.
Fig_2b <- ggplot(psmelt.data, aes(x=Type, y=Abundance, fill=Class)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name="Class", values=myColors) +
  scale_y_continuous(expand=c(0,0), limits=c(0, 102)) +
  scale_x_discrete(labels=function(Type) str_wrap(Type, width=3)) +
  guides(fill=guide_legend(ncol=2)) +
  theme_bw() +
  theme(panel.spacing=unit(1, "lines"), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(color="black", size = 8, vjust=0.5, hjust=1, angle=90),
        axis.text.y=element_text(color="black", size = 8, vjust=1),
        axis.title=element_text(size=9, color="black"),
        strip.background=element_blank(), 
        strip.text=element_text(size=9, color="black"),
        plot.tag=element_text(face="bold"),
        legend.title=element_blank(),
        legend.position="none",
        legend.text=element_text(size=8),
        legend.key.size=unit(0.5, "line"),
        legend.justification=0) +
  ylab("Mean relative abundance (%)") +
  xlab(NULL) +
  labs(tag="b")

# Assemble and save plots.
layout <- rbind(c(1,2),
                c(3,3))

gs <- list(Fig_2a, Fig_2b, Fig_2.legend)
grid.arrange(grobs=gs, layout_matrix=layout, widths=c(0.35,1), heights=c(1,0.15))

Fig_2 <- arrangeGrob(grobs=gs, layout_matrix=layout, widths=c(0.35,1), heights=c(1,0.2))

ggsave("Fig_2.jpg", plot = Fig_2, width=6, height=4, units="in", dpi=300)

rm(temp, temp.data, psmelt.data, Fig_2a, Fig_2b, Fig_2.legend)

#### (done) FIGURE 3 - INTESTINAL BIOGEOGRAPHY ####
# Calculate summary data for heat map (families with mean abundance > 0.09%) #
# Convert to relative abundance and subset to only families present at >0.1% abundance.
types <- c("duodenum", "jejunum", "ileum", "caecum", "ascending.colon", "descending.colon", "feces")

plot.pseq <- tax_glom(pseq.rarefied, taxrank="Family")
plot.pseq <- transform_sample_counts(plot.pseq, function(x) 100*x/sum(x))
plot.pseq <- prune_taxa((taxa_sums(plot.pseq)/nsamples(plot.pseq)) > 0.09, plot.pseq)

# Subset to phyla of interest.
plot.pseq <- subset_taxa(plot.pseq, Phylum %in% c("Actinobacteria", "Proteobacteria", "Fusobacteria", "Bacteroidetes", "Firmicutes"))

# Calculate mean abundances in each intestinal site.
temp <- as.data.frame(otu_table(plot.pseq))
colnames(temp) <- as.data.frame(cbind(tax_table(plot.pseq)))$Family
temp$Type <- bg.metadata$Type
temp.duodenum <- subset(temp, Type=="duodenum")
temp.jejunum <- subset(temp, Type=="jejunum")
temp.ileum <- subset(temp, Type=="ileum")
temp.caecum <- subset(temp, Type=="caecum")
temp.asc.colon <- subset(temp, Type=="ascending.colon")
temp.des.colon <- subset(temp, Type=="descending.colon")
temp.feces <- subset(temp, Type=="feces")

means.duodenum <- colMeans(temp.duodenum[,1:ncol(temp.duodenum)-1])
means.jejunum <- colMeans(temp.jejunum[,1:ncol(temp.duodenum)-1])
means.ileum <- colMeans(temp.ileum[,1:ncol(temp.duodenum)-1])
means.caecum <- colMeans(temp.caecum[,1:ncol(temp.duodenum)-1])
means.asc.colon <- colMeans(temp.asc.colon[,1:ncol(temp.duodenum)-1])
means.des.colon <- colMeans(temp.des.colon[,1:ncol(temp.duodenum)-1])
means.feces <- colMeans(temp.feces[,1:ncol(temp.duodenum)-1])

plot.means <- as.data.frame(t(rbind(means.duodenum, means.jejunum, means.ileum, means.caecum, means.asc.colon, means.des.colon, means.feces)))

colnames(plot.means) <- c("duodenum", "jejunum", "ileum", "caecum", "ascending.colon", "descending.colon", "feces")
rm(temp, temp.duodenum, temp.jejunum, temp.ileum, temp.caecum, temp.asc.colon, temp.des.colon, temp.feces,
   means.duodenum, means.jejunum, means.ileum, means.caecum, means.asc.colon, means.des.colon, means.feces)

# Summarize data.
library(tidyverse)
heatmap.data <- plot.means %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
head(heatmap.data)
desired_order <- c("duodenum", "jejunum", "ileum", "caecum", "ascending.colon", "descending.colon", "feces")
heatmap.data$colname <- factor(as.character(heatmap.data$colname), levels=desired_order)
rm(plot.means)

# Convert any abundances >10% to 10%. This is to make it easier to visually distinguish shades of red in the final plot.
heatmap.data$value[heatmap.data$value > 10] <- 10

# Rename variables for plotting.
levels(heatmap.data$colname)[levels(heatmap.data$colname)=="duodenum"] <- "Duodenum"
levels(heatmap.data$colname)[levels(heatmap.data$colname)=="jejunum"] <- "Jejunum"
levels(heatmap.data$colname)[levels(heatmap.data$colname)=="ileum"] <- "Ileum"
levels(heatmap.data$colname)[levels(heatmap.data$colname)=="caecum"] <- "Caecum"
levels(heatmap.data$colname)[levels(heatmap.data$colname)=="ascending.colon"] <- "Asc. Colon"
levels(heatmap.data$colname)[levels(heatmap.data$colname)=="descending.colon"] <- "Des. Colon"
levels(heatmap.data$colname)[levels(heatmap.data$colname)=="feces"] <- "Feces"

# Order by phylum, then alphabetically by family.
temp.tax <- as.data.frame(cbind(tax_table(plot.pseq)))
temp.tax <- subset(temp.tax, Family %in% heatmap.data$rowname)
temp.tax <- arrange(temp.tax, Phylum, Family)
temp.tax <- mutate(temp.tax, Family = factor(Family, Family))
heatmap.data$rowname <- factor(as.character(heatmap.data$rowname), levels=levels(temp.tax$Family))

# Assign colors for each phylum.
colornames <- c("#CC0000", #Actinobacteria
                "magenta3", "magenta3", "magenta3", "magenta3",  # Bacteroidetes
                "#996600", "#996600", "#996600", "#996600", "#996600", "#996600", "#996600", "#996600", "#996600", "#996600", "#996600", "#996600", 
                "darkgreen", # Fusobacteria
                "#0000FF", "#0000FF", "#0000FF", "#0000FF", "#0000FF", "#0000FF", "#0000FF", "#0000FF", "#0000FF") # Proteobacteria
names(colornames) <- levels(heatmap.data$rowname)

# Extract summary data from mixed effects logistic regressions. #
logit.data <- coeff.logits.mix.eff[["Family"]]

# Obtain the same taxa that were used in the heat map, and sort by phylum, then alphabetically by family. 
logit.data$Family <- as.factor(logit.data$Family)
logit.data <- subset(logit.data, Family %in% temp.tax$Family)
logit.data$Family <- as.factor(logit.data$Family)

# Need to add two levels for families that were not used in logistic regression analysis.
levels(logit.data$Family) <- c(levels(logit.data$Family), "Veillonellaceae", "Desulfovibrionaceae")
logit.data$Family <- factor(as.character(logit.data$Family), levels=levels(temp.tax$Family))

# Add blank rows for these two families.
add.df <- rbind(c("Veillonellaceae", NA, NA, NA, NA, NA),
                c("Desulfovibrionaceae", NA, NA, NA, NA, NA))
colnames(add.df) <- colnames(logit.data)
logit.data <- rbind(logit.data, add.df)
logit.data$Family <- factor(logit.data$Family)
rm(add.df)

# Multiply coefficients and errors by -1 so that small intestines appear on the left (negative) and large intestinal indicators appear on the right (positive).
for(j in 2:6){
  logit.data[,j] <- as.numeric(logit.data[,j])
}
for(j in 2:4){
  logit.data[,j] <- -1*logit.data[,j]
}

# Replace values outside the range of the graph.
logit.data$coeff[logit.data$Family=="Vibrionaceae"] <- NA
logit.data$`2.5`[logit.data$Family=="Vibrionaceae"] <- NA
logit.data$`97.5`[logit.data$Family=="Vibrionaceae"] <- NA

logit.data$`2.5`[logit.data$Family=="Campylobacteraceae"] <- 24.9
logit.data$`2.5`[logit.data$Family=="Helicobacteraceae"] <- 24.9

# Replace family names with numbers (to avoid taking up too much space on the plot).
logit.data <- arrange(logit.data, Family)
levels(logit.data$Family) <- c(1:27)
logit.data$Family <- as.factor(logit.data$Family)

rm(temp.tax)

# Summarize data for ASV appearance bar chart.
# For the graph, "new to site" in the duodenum is all the taxa in the duodenum.
df.new.to.site$duodenum <- df.total.taxa$X1

# Calculate shared taxa with previous sites (total taxa - new taxa)
df.shared.taxa <- df.total.taxa - df.new.to.site

# Calculate means and clean data frame.
plot.data <- rbind(colMeans(df.shared.taxa, na.rm=TRUE), colMeans(df.new.to.site, na.rm=TRUE))
rownames(plot.data) <- c("shared", "new")
colnames(plot.data) <- c("duodenum", "jejunum", "ileum", "caecum", "ascending.colon", "descending.colon", "feces")
rm(df.shared.taxa, df.total.taxa, df.new.to.site, df.uniq.to.site)
plot.data <- as.data.frame(t(plot.data))
plot.data$Site <- rownames(plot.data)
plot.data <- melt(plot.data)

# Create plots.
Fig3.heatmap <- ggplot(heatmap.data, aes(x = colname, y = rowname, fill = value)) +
  geom_tile(color="black") +
  scale_fill_gradient(low="linen",
                      high="red",
                      limit=c(0.000001,10),
                      na.value="grey38",
                      space="Lab",
                      name="Relative Abundance (%)") +
  theme_bw() + 
  scale_y_discrete(limits = rev(levels(heatmap.data$rowname))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1, color="black"),
        axis.text.y = element_text(size=8, color=rev(colornames)),
        axis.title = element_text(color="black", size=9),
        axis.line = element_blank(),
        panel.spacing=unit(2, "lines"),
        panel.border = element_blank(),
        strip.text = element_text(color="black", size=9),
        strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title=element_text(color="black", size=9),
        legend.position="none",
        legend.text=element_text(color="black", size=8),
        legend.key.size = unit(1,"lines")) +
  xlab("Intestinal Location") +
  ylab("Family")

Fig3.bars <- ggplot(plot.data, aes(x=Site, y=value, fill=variable)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("darkgrey", "chartreuse3")) +
  labs(x=element_blank(), y="Number of ASVs") +
  theme_bw() +
  theme(legend.position="none",
        axis.text=element_text(color="black", size=8),
        axis.text.x=element_blank(),
        axis.title = element_text(color="black", size=9))

Fig3.logits <- ggplot(logit.data, aes(x=Family, y=coeff)) + 
  geom_point(color=colornames) +
  theme_bw() +
  geom_errorbar(aes(ymin=logit.data$`2.5`, ymax=logit.data$`97.5`, x=Family), position=position_dodge(0.9), width=0.2, color=rev(colornames)) +
  geom_hline(yintercept=0) +
  labs(x=element_blank(), y="Coefficient") +
  ylim(-25,NA) +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(logit.data$Family))) +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(color="black", size=8),
    axis.title = element_text(color="black", size=9))

# Ensure graphs have the same axes widths and heights so they align in the plot.
plots <- list(Fig3.bars, Fig3.heatmap, Fig3.logits)
grobs <- list()
widths <- list()
heights <- list()

for (i in 1:length(plots)){
  grobs[[i]] <- ggplotGrob(plots[[i]])
  widths[[i]] <- grobs[[i]]$widths[2:5]
  heights[[i]] <- grobs[[i]]$heights[1:12]
}

maxwidth <- do.call(grid::unit.pmax, widths)
maxheight <- do.call(grid::unit.pmax, heights)

# Assign the max width to the vertical grobs (bars and heatmap)
for (i in c(1,2)){
  grobs[[i]]$widths[2:5] <- as.list(maxwidth) }

# Assign the max height to the horizontal grobs (heatmap and logits plot)
for (i in c(2,3)){
  grobs[[i]]$heights[1:12] <- as.list(maxheight)
}

rm(maxwidth, maxheight, widths, heights, plots)

grid.arrange(grobs[[1]],
             blank,
             grobs[[2]],
             grobs[[3]],
             ncol=2, nrow=2, widths=c(0.6,0.4), heights=c(0.3,0.7))

Fig_3 <- arrangeGrob(grobs[[1]],
                     blank,
                     grobs[[2]],
                     grobs[[3]],
                     ncol=2, nrow=2, widths=c(0.6,0.4), heights=c(0.3,0.7))

ggsave(filename="~/Fig_3.jpg", plot=Fig_3, width=5.5, height=6.5, units="in", dpi=300)

# Extract figure legend for heat map.
heatplot.legend <- cowplot::get_legend(ggplot(heatmap.data, aes(x = rowname, y = colname, fill = value)) +
                                geom_tile(color="black") +
                                scale_fill_gradient(low="linen",
                                                    high="red",
                                                    limit=c(0.000001,10),
                                                    na.value="grey38",
                                                    space="Lab",
                                                    name="Relative\nAbundance (%)") +
                                theme_bw() + 
                                scale_y_discrete(limits = rev(levels(dt2$colname))) +
                                theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                                                 size = 8, hjust = 1, color=colornames),
                                      axis.text.y = element_text(color="black", size=8),
                                      axis.title = element_text(color="black", size=9),
                                      axis.line = element_blank(),
                                      panel.spacing=unit(2, "lines"),
                                      panel.border = element_blank(),
                                      strip.text = element_text(color="black", size=9),
                                      strip.background = element_blank(),
                                      panel.grid.major = element_blank(), 
                                      panel.grid.minor = element_blank(),
                                      legend.title=element_text(color="black", size=9),
                                      legend.position="right",
                                      legend.text=element_text(color="black", size=8),
                                      legend.key.size = unit(1,"lines")) +
                                xlab("Family") +
                                ylab("Intestinal Location"))

ggsave("~/heatplot.legend.jpg", heatplot.legend, width=1, height=2, dpi=300)

# Extract figure legend for bar graph.
bar.graph.legend <- plot.data
bar.graph.legend$variable <- factor(bar.graph.legend$variable, levels=c("new", "shared"))
levels(bar.graph.legend$variable)[levels(bar.graph.legend$variable)=="shared"] <- "Shared with proximal sites"
levels(bar.graph.legend$variable)[levels(bar.graph.legend$variable)=="new"] <- "New to intestine"
colnames(bar.graph.legend)[2] <- "ASV Classification"

bargraph.legend <- cowplot::get_legend(ggplot(bar.graph.legend, aes(x=Site, y=NumASVsALL, fill=`ASV classification`)) +
                                geom_bar(position="stack", stat="identity") +
                                scale_fill_manual(values=c("chartreuse3", "darkgrey")) +
                                labs(x=element_blank(), y="Number of ASVs") +
                                theme_bw() +
                                theme(legend.position="right",
                                      axis.text=element_text(color="black", size=8),
                                      axis.text.x=element_blank(),
                                      axis.title = element_text(color="black", size=9),
                                      legend.title=element_text(color="black", size=9),
                                      legend.text=element_text(color="black", size=8),
                                      legend.key.size = unit(1,"lines")))

ggsave("~/bargraph.legend.jpg", bargraph.legend, width=2, height=1, dpi=300)

# Clean workspace.
rm(Fig3.bars, Fig3.heatmap, Fig3.logits, heatmap.data, logit.data, plot.data, bar.graph.legend, heatmap.legend, grobs)


#### (done) FIGURE 4 - dbRDA STRUCTURING EFFECTS ##############
# (a) Complete dbRDA plot with segment error ellipses and hulls around individuals.
# Extract data from dbRDA model.
library(ggvegan)
dbrda.nofeces_fort <- fortify(dbrda.nofeces, display="sites")
dbrda.nofeces_fort <- cbind(bg.metadata.nf, dbrda.nofeces_fort)
desired_order <- c("small.intestine","large.intestine")
dbrda.nofeces_fort$Segment <- factor(dbrda.nofeces_fort$Segment, levels=desired_order)

# Calculate standard ellipse.
plot.new()
dbrda.nofeces.ellipse <- ordiellipse(dbrda.nofeces, bg.metadata.nf$Segment, display="sites", 
                                kind="sd", conf=0.75, label=T)
dbrda.nofeces.ell.data <- data.frame()
for(g in levels(bg.metadata.nf$Segment)){
  dbrda.nofeces.ell.data <- rbind(dbrda.nofeces.ell.data,
                             cbind(as.data.frame(with(bg.metadata.nf[dbrda.nofeces_fort$Segment==g,],
                                                      veganCovEllipse(dbrda.nofeces.ellipse[[g]]$cov,
                                                                      dbrda.nofeces.ellipse[[g]]$center,
                                                                      dbrda.nofeces.ellipse[[g]]$scale)))
                                   ,Type=g))}
rm(dbrda.nofeces.ellipse)

# Rename variables.
dbrda.nofeces_fort$Type <- as.character(dbrda.nofeces_fort$Type)
dbrda.nofeces_fort$Type[dbrda.nofeces_fort$Type=="duodenum"] <- "Duodenum"
dbrda.nofeces_fort$Type[dbrda.nofeces_fort$Type=="jejunum"] <- "Jejunum"
dbrda.nofeces_fort$Type[dbrda.nofeces_fort$Type=="ileum"] <- "Ileum"
dbrda.nofeces_fort$Type[dbrda.nofeces_fort$Type=="caecum"] <- "Caecum"
dbrda.nofeces_fort$Type[dbrda.nofeces_fort$Type=="ascending.colon"] <- "Ascending Colon"
dbrda.nofeces_fort$Type[dbrda.nofeces_fort$Type=="descending.colon"] <- "Descending Colon"
dbrda.nofeces_fort$Type <- factor(dbrda.nofeces_fort$Type, levels=c("Duodenum", "Jejunum", "Ileum", "Caecum", "Ascending Colon", "Descending Colon"))

# Create a synthetic variable that concatenates Individual + Segment, to create hulls for each individual within each intestinal segment.
bg.metadata.nf$SegInd <- paste0(bg.metadata.nf$Individual, bg.metadata.nf$Segment)
ind.hulls <- gg_ordiplot(dbrda.nofeces, groups=bg.metadata.nf$SegInd, hull=TRUE, ellipse=FALSE, label=TRUE, pt.size=2)

# Create first plot (ellipses for segment, hulls for individual).
dbrda.plot.a <- ggplot(dbrda.nofeces_fort, aes(x=CAP1, y=CAP2)) +
  geom_vline(xintercept = 0, color="grey") +
  geom_hline(yintercept = 0, color="grey") +
  
  # The individual hulls need to be scaled due to some discrepancy between the gg_ordiplot function and actual point values from dbrda.
  geom_path(data=ind.hulls[["df_hull"]], aes(x=2.18*x, y=2.75*y, linetype=Group), size=0.5, color="grey") +
  scale_linetype_manual("Individual", values=rep("solid", 20)) +
  
  geom_point(aes(color = Individual, shape=Type), size=2, fill="white") +
  
  # Manually add SI and LI labels at the segment centroids.
  geom_label(aes(x=0.70974, y=0.06806), color="black", label="LI") +
  geom_label(aes(x=-0.73421, y=-0.07041), color="black", label="SI") +
  new_scale("linetype") +
  
  # Add standard ellpse.
  geom_path(data=dbrda.nofeces.ell.data, aes(x=CAP1, y=CAP2, linetype=Type), color="black", size=0.5) +
  scale_linetype_manual(values=c("solid","solid"), guide=FALSE) +
  
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        axis.text=element_text(color="black", size = 8, vjust=1),
        axis.title=element_text(size=9),
        plot.title=element_text(size=10),
        strip.background=element_blank(), 
        plot.tag=element_text(face="bold"),
        strip.text=element_blank(),
        legend.title=element_blank(),
        legend.position="none") +
  scale_color_manual(values=dbrda.indv.pal(10), guide=FALSE) +
  scale_shape_manual(values=c(22,21,24,15,16,17)) +
  labs(x="dbRDA1 [20.3%]", y="dbRDA2 [12.6%]", tag="a")

# (b) Repeat dbRDA plot, with hulls distinguishing sites.
# Extract hull information.
site.hulls <- gg_ordiplot(dbrda.nofeces, groups=as.character(bg.metadata.nf$Type), hull=TRUE, ellipse=FALSE, label=TRUE, pt.size=2)

# Rename values for plotting.
site.hulls[["df_hull"]]$Group <- as.character(site.hulls[["df_hull"]]$Group)
site.hulls[["df_hull"]]$Group[site.hulls[["df_hull"]]$Group=="duodenum"] <- "Duodenum"
site.hulls[["df_hull"]]$Group[site.hulls[["df_hull"]]$Group=="jejunum"] <- "Jejunum"
site.hulls[["df_hull"]]$Group[site.hulls[["df_hull"]]$Group=="ileum"] <- "Ileum"
site.hulls[["df_hull"]]$Group[site.hulls[["df_hull"]]$Group=="caecum"] <- "Caecum"
site.hulls[["df_hull"]]$Group[site.hulls[["df_hull"]]$Group=="ascending.colon"] <- "Ascending Colon"
site.hulls[["df_hull"]]$Group[site.hulls[["df_hull"]]$Group=="descending.colon"] <- "Descending Colon"
site.hulls[["df_hull"]]$Group <- factor(site.hulls[["df_hull"]]$Group, levels=c("Duodenum", "Jejunum", "Ileum", "Caecum", "Ascending Colon", "Descending Colon"))

# Create plot.
dbrda.plot.b <- ggplot(dbrda.nofeces_fort, aes(x=CAP1, y=CAP2)) +
  geom_vline(xintercept = 0, color="grey") +
  geom_hline(yintercept = 0, color="grey") +
  geom_path(data=site.hulls[["df_hull"]], aes(x=2.18*x, y=2.75*y, color=Group), size=0.5) +
  geom_point(aes(color = Type, shape=Type), size=2, fill="white") +
  annotate("text", x=-1.45, y=1.05, hjust=0, color="#00A2BD", label="duodenum", size=3) +
  annotate("text", x=-0.70, y=0.25, hjust=0, color="#D544D1", label="jejunum", size=3) +
  annotate("text", x=-0.4, y=-0.6, hjust=0, color="#B97700", label="ileum", size=3) +
  annotate("text", x=0.65, y=1.4, hjust=0, color="#00A358", label="caecum", size=3) +
  annotate("text", x=0.05, y=1.7, hjust=0, color="#4C7FEE", label="asc.colon", size=3) +
  annotate("text", x=0.6, y=-0.85, hjust=0, color="#E14F77", label="des.colon", size=3) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        axis.text=element_text(color="black", size = 8, vjust=1),
        axis.title=element_text(size=9),
        plot.title=element_text(size=10),
        strip.background=element_blank(), 
        plot.tag=element_text(face="bold"),
        strip.text=element_blank(),
        legend.title=element_blank(),
        legend.position="none") +
  scale_color_manual(values=dbrda.indv.pal(6), guide=FALSE) +
  scale_shape_manual(values=c(22,21,24,15,16,17)) +
  labs(x="dbRDA1 [20.3%]", y="dbRDA2 [12.6%]", tag="b")

# Extract legend and save plot. #
dbrda.legend <- get_legend(ggplot(dbrda.nofeces_fort, aes(x=CAP1, y=CAP2)) +
                                   geom_point(aes(color = Individual, shape=Type), size=2) +
                                   theme(panel.border = element_blank(),
                                         panel.grid.major = element_blank(), 
                                         panel.grid.minor = element_blank(),
                                         panel.background = element_rect(colour = "black", size=0.5, fill=NA),
                                         axis.line = element_line(colour = "black", size=0.5),
                                         axis.text=element_text(color="black", size = 8, vjust=1),
                                         axis.title=element_text(size=9),
                                         plot.title=element_text(size=10),
                                         strip.background=element_blank(), 
                                         plot.tag=element_text(face="bold"),
                                         strip.text=element_blank(),
                                         legend.title=element_blank(),
                                         legend.position="bottom",
                                         legend.key=element_blank()) +
                                   scale_color_manual(values=dbrda.indv.pal(10), guide=FALSE) +
                                   scale_shape_manual(values=c(0,1,2,15,16,17)) +
                                   labs(x="dbRDA1 [20.3%]", y="dbRDA2 [12.6%]", tag="b"))

Fig_4 <- arrangeGrob(dbrda.plot.a, dbrda.plot.b, dbrda.legend,
                                  layout_matrix=rbind(c(1,2),c(3)), heights=c(0.8,0.2))

ggsave("~/Fig_4.jpg", Fig_4, width=6.5, height=4, units="in", dpi=300)
rm(dbrda.plot.a, dbrda.plot.b, dbrda.legend, ind.hulls, site.hulls)

#### (done) FIGURE 5 - FECAL EFFICACY PANEL ###################
# (a) Bray-Curtis distance from feces to each sample. 
feces.plot.a <- ggplot(dist.to.feces.BC.summary, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot(outlier.shape=NA) +
  scale_fill_manual(values=pal(6)) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_blank(),
        axis.title=element_text(size=9),
        strip.background=element_blank(), 
        strip.text=element_text(size=9),
        legend.title=element_blank(),
        plot.tag=element_text(face="bold"),
        legend.position="none") +
  labs(x="Site", y="Bray-Curtis distance to feces", tag="a")

# (b) Euler diagram showing ASV overlap (in percentages)
# These numbers were calculated in the analysis script using the function rowMeans(Venn.data)
vector <- c(16.1, 13.1, 1.7, 9.6, 29.9, 23.6, 5.9)
names(vector) <- c("SI", "SI&LI", "SI&F", "LI&F", "LI", "SI&LI&F", "F")

feces.plot.b <- plot(eulerr::euler(vector, shape="circle"), quantities=TRUE, fills=c("lightsteelblue1", "lightgoldenrod", "darkolivegreen1"),
             fill_alpha=0.4)
feces.plot.b <- ggplotify::as.grob(feces.plot.b)
rm(vector)

# (c) Line graph showing abundance and prevalence of fecal-undetected ASVs, by phylum.
# Rename factor levels for plotting.
feces.undetected.distr$Site <- as.character(feces.undetected.distr$Site)
feces.undetected.distr$Site[feces.undetected.distr$Site=="duo"] <- "Duodenum"
feces.undetected.distr$Site[feces.undetected.distr$Site=="jej"] <- "Jejunum"
feces.undetected.distr$Site[feces.undetected.distr$Site=="ile"] <- "Ileum"
feces.undetected.distr$Site[feces.undetected.distr$Site=="cae"] <- "Caecum"
feces.undetected.distr$Site[feces.undetected.distr$Site=="asc"] <- "Asc.\nColon"
feces.undetected.distr$Site[feces.undetected.distr$Site=="des"] <- "Des.\nColon"
feces.undetected.distr$Site <- as.factor(feces.undetected.distr$Site)
feces.undetected.distr$Site <- factor(feces.undetected.distr$Site, levels=c("Duodenum", "Jejunum", "Ileum", "Caecum", "Asc.\nColon", "Des.\nColon"))

# Create plot.
feces.plot.c <- ggplot(feces.undetected.distr, aes(x=Site, y=prevalence)) +
  geom_point(aes(size=abundance)) +
  geom_path(group=1) +
  theme_bw() +
  theme(legend.position="bottom", 
        legend.justification=c("top"),
        legend.title=element_text(size=9, color="black"),
        plot.tag=element_text(face="bold"),
        axis.line = element_line(colour = "black", size=0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1, color="black"),
        axis.text.y = element_text(size=8, color="black"),
        axis.title=element_text(size=9, color="black")) +
  #guides(colour=guide_legend(ncol=1, title.position="top")) +
  guides(size=guide_legend(ncol=2, title.position="top")) +
  labs(tag="c")

# (d) Abundance fold-changes from the SI and LI to fecal samples.
temp.melt <- melt(fold.ASV.seg.changes.summary)

# Subset to dominant phyla, and rename lower-abundance classes to their higher phylum.
temp.melt <- subset(temp.melt, Phylum %in% c("Actinobacteria", "Proteobacteria", "Fusobacteria",
                                             "Bacteroidetes", "Firmicutes"))
temp.melt$Class <- as.character(temp.melt$Class)
temp.melt$Class[temp.melt$Phylum=="Actinobacteria"] <- "Actinobacteria"
temp.melt$Class[temp.melt$Phylum=="Bacteroidetes"] <- "Bacteroidetes"
temp.melt$Class[temp.melt$Phylum=="Fusobacteria"] <- "Fusobacteria"
temp.melt$Class[temp.melt$Phylum=="Proteobacteria"] <- "Proteobacteria"
temp.melt <- subset(temp.melt, Class !="NA")
temp.melt$Class <- as.factor(temp.melt$Class)
temp.melt$Class <- factor(temp.melt$Class, levels=c("Actinobacteria", "Bacteroidetes", "Fusobacteria", "Proteobacteria",
                                                    "Bacilli", "Clostridia", "Erysipelotrichia", "Negativicutes"))

# Rename variables. 
temp.melt$variable <- as.character(temp.melt$variable)
temp.melt$variable[temp.melt$variable=="small.intestine"] <- "Small\nIntestine"
temp.melt$variable[temp.melt$variable=="large.intestine"] <- "Large\nIntestine"
temp.melt$variable <- as.factor(temp.melt$variable)
temp.melt$variable <- factor(temp.melt$variable, levels=c("Small\nIntestine", "Large\nIntestine"))

# Set fold changes of 0 to NA.
temp.melt[temp.melt==0] <- NA

# Plot.
feces.plot.d <- ggplot(temp.melt, aes(x=Class, y=value, fill=variable)) + 
  geom_boxplot(outlier.shape=NA) +
  scale_fill_manual(values=c("orangered1","deepskyblue3")) + 
  scale_y_continuous(trans="log10", breaks=c(0.01, 0.1, 1, 10, 100), labels=log_breaks(n=5, base=10)) + 
  geom_hline(yintercept=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1, color="black"),
        axis.text.y = element_text(size=8, color="black"),
        axis.title = element_text(size=9, color="black"),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.text = element_text(size=9, color="black"),
        legend.position="bottom",
        plot.tag = element_text(face="bold")) +
  guides(color=guide_legend(ncol=2)) +
  labs(x="Taxon", y="Log fold change vs. feces", tag="d")

# Ensure panels have the same heights and widths.
plots <- list(feces.plot.a,
              feces.plot.c, 
              feces.plot.d)

for (i in 1:length(plots)){
  grobs[[i]] <- ggplotGrob(plots[[i]])
  widths[[i]] <- grobs[[i]]$widths[2:5]
  heights[[i]] <- grobs[[i]]$heights[1:12]
}

maxwidth <- do.call(grid::unit.pmax, widths)
maxheight <- do.call(grid::unit.pmax, heights)

# Assign the max width to the vertical grobs (feces.plot.a and feces.plot.c)
for (i in c(1,2)){
  grobs[[i]]$widths[2:5] <- as.list(maxwidth)
}

# Assign the max height to the horizontal grob (feces.plot.c and feces.plot.d)
for (i in c(2,3)){
  grobs[[i]]$heights[1:12] <- as.list(maxheight)
}

grid.arrange(grobs[[1]],
             feces.plot.b,
             grobs[[2]],
             grobs[[3]],
             ncol=2, widths=c(0.5, 0.5), heights=c(0.4,0.6))
Fig_5 <- arrangeGrob(grobs[[1]],
                                  feces.plot.b,
                                  grobs[[2]],
                                  grobs[[3]],
                                  ncol=2, widths=c(0.5, 0.5), heights=c(0.4,0.6))

ggsave("~/Fig_5.jpg", Fig_5, width=6.5, height=8, units="in", dpi=300)

#### (done) FIGURE S1 - RAREFACTION CURVES ####
# Rename variable for plotting.
levels(rarefaction_curve_data_summary_verbose$Type)[levels(rarefaction_curve_data_summary_verbose$Type)=="duodenum"] <- "Duodenum"
levels(rarefaction_curve_data_summary_verbose$Type)[levels(rarefaction_curve_data_summary_verbose$Type)=="jejunum"] <- "Jejunum"
levels(rarefaction_curve_data_summary_verbose$Type)[levels(rarefaction_curve_data_summary_verbose$Type)=="ileum"] <- "Ileum"
levels(rarefaction_curve_data_summary_verbose$Type)[levels(rarefaction_curve_data_summary_verbose$Type)=="caecum"] <- "Caecum"
levels(rarefaction_curve_data_summary_verbose$Type)[levels(rarefaction_curve_data_summary_verbose$Type)=="ascending.colon"] <- "Ascending\nColon"
levels(rarefaction_curve_data_summary_verbose$Type)[levels(rarefaction_curve_data_summary_verbose$Type)=="descending.colon"] <- "Descending\nColon"
levels(rarefaction_curve_data_summary_verbose$Type)[levels(rarefaction_curve_data_summary_verbose$Type)=="feces"] <- "Feces"

# Plot rarefaction curves for species richness.
rare.plot.a <- ggplot(
  data = subset(rarefaction_curve_data_summary_verbose, Measure=="Observed"),
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    colour = Type,
    group = Sample)) +
  geom_line() +
  geom_pointrange(size=0.3) + 
  facet_wrap(facets = ~ Type, scales = 'free_y') +
  labs(x="Sequencing depth (number of reads)", y="ASV Richness", tag="a") +
  theme_bw() +
  theme(legend.position="none",
        panel.grid=element_blank(),
        axis.text.x=element_text(size=9, angle=45, color="black", hjust=1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title=element_text(size=10, color="black"),
        strip.background=element_blank(),
        plot.tag=element_text(face="bold"),
        strip.text=element_text(size=9, color="black"))

# Plot rarefaction curves for Shannon diversity. 
rare.plot.b <- ggplot(
  data = subset(rarefaction_curve_data_summary_verbose, Measure=="Shannon"),
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    colour = Type,
    group = Sample)) +
  geom_line() +
  geom_pointrange(size=0.3) + 
  facet_wrap(facets = ~ Type, scales = 'free_y') +
  labs(x="Sequencing depth (number of reads)", y="Shannon index", tag="b") +
  theme_bw() +
  theme(legend.position="none",
        panel.grid=element_blank(),
        axis.text.x=element_text(size=9, angle=45, color="black", hjust=1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title=element_text(size=10, color="black"),
        strip.background=element_blank(),
        plot.tag=element_text(face="bold"),
        strip.text=element_text(size=9, color="black"))

# Assemble and save files.
Fig_S1 <- arrangeGrob(grobs=list(rare.plot.a, rare.plot.b), ncol=2)
ggsave("~/Fig_S1.jpg", plot=Fig_S1, width=9, height=5.5, units="in", dpi=300)
rm(rare.plot.a, rare.plot.b)

#### (done) FIGURE S2 - MaAsLIN RESULTS ####
# Import the output data from Galaxy.
maaslin.graph.data <- read.csv("~/maaslin.graph.data.csv")

# Organize by Phylum, then by Family.
maaslin.graph.data <- arrange(maaslin.graph.data, Phylum, Feature)
maaslin.graph.data <- mutate(maaslin.graph.data, Feature = factor(Feature, Feature))

# Round adjusted p values to three decimal points.
maaslin.graph.data$Q.value <- round(maaslin.graph.data$Q.value, 3)
maaslin.graph.data$Q.value[maaslin.graph.data$Q.value < 0.001] <- "<0.001"

# Plot families against their MaAsLin coefficients.
Fig_S2 <- ggplot(maaslin.graph.data, aes(x=Feature, y=Coefficient)) + 
  geom_point(aes(color=Phylum), size=2) +
  geom_label(aes(x=Feature,  y=0.75, label=Q.value), size=3) +
  theme_bw() +
  geom_errorbar(aes(ymin=maaslin.graph.data$Min, ymax=maaslin.graph.data$Max, x=Feature, color=Phylum), position=position_dodge(0.9), width=0.2) +
  geom_hline(yintercept=0) +
  labs(x="Family", y="Coefficient") +
  scale_y_continuous(limits=c(NA,0.8), breaks=c(-0.25,0,0.25,0.5)) +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(maaslin.graph.data$Feature))) +
  scale_color_manual(values=c("#CC0000", "#e75480", "#996600", "darkgreen", "#0000FF")) +
  theme(
    axis.text.y = element_text(size=9, color=rev(c("#CC0000", "#e75480", "#e75480", "#996600","#996600","#996600","#996600","#996600","darkgreen", "#0000FF", "#0000FF", "#0000FF"))),
    axis.text.x = element_text(color="black", size=9),
    axis.title = element_text(color="black", size=10),
    legend.position="none")

# Save figure.
ggsave("~/Fig_S2.jpg", Fig_S2, width=5, height=3.5, units="in", dpi=300)



#### (done) FIGURE S3 - TAXON APPEARANCE AT EACH SITE ####
plot.data <- ids.genus.appear

# Add taxonomy information.
plot.data$Genus <- rownames(plot.data)
plot.data <- merge(plot.data, as.data.frame(cbind(tax_table(tax_glom(pseq.rarefied, taxrank="Genus")))), by="Genus")

# Identify families present in at least 8 individuals.
families <- subset(ids.family.appear, prev > 8)
families <- rownames(families)

# Subset genera table to only include prevalent families.
plot.data <- subset(plot.data, Family %in% families)

# Subset genera table to only include prevalent genera.
plot.data <- subset(plot.data, prev > 7)

# Clean data and melt.
plot.data$prev <- NULL
plot.data$Kingdom <- NULL
plot.data$Class <- NULL
plot.data$Order <- NULL

plot.data <- melt(plot.data)
plot.data$Genus <- factor(plot.data$Genus)

# Rename factor levels for plotting.
levels(plot.data$variable)[levels(plot.data$variable)=="duodenum"] <- "duo"
levels(plot.data$variable)[levels(plot.data$variable)=="jejunum"] <- "jej"
levels(plot.data$variable)[levels(plot.data$variable)=="ileum"] <- "ile"
levels(plot.data$variable)[levels(plot.data$variable)=="caecum"] <- "cae"
levels(plot.data$variable)[levels(plot.data$variable)=="asc.colon"] <- "asc"
levels(plot.data$variable)[levels(plot.data$variable)=="des.colon"] <- "des"
levels(plot.data$variable)[levels(plot.data$variable)=="feces"] <- "fec"
plot.data$Family <- factor(plot.data$Family)

# Assign linetype values for the plot.
plot.data$linetype <- rep(NA)
plot.data$linetype[plot.data$Family %in% c("Acidaminococcaceae", "Bacteroidaceae", "Coriobacteriaceae", "Enterobacteriaceae",
                                                             "Erysipelotrichaceae", "Helicobacteraceae", "Lactobacillaceae", "Prevotellaceae",
                                                             "Rubritaleaceae", "Succinivibrionaceae", "Sutterellaceae")] <- 1

plot.data$linetype[plot.data$Genus=="Fusobacterium"] <- 1
plot.data$linetype[plot.data$Genus=="Cetobacterium"] <- 2

plot.data$linetype[plot.data$Genus=="Blautia"] <- 1
plot.data$linetype[plot.data$Genus=="Clostridium_XlVa"] <- 2
plot.data$linetype[plot.data$Genus=="Clostridium_XlVb"] <- 3
plot.data$linetype[plot.data$Genus=="Dorea"] <- 4
plot.data$linetype[plot.data$Genus=="Fusicatenibacter"] <- 5
plot.data$linetype[plot.data$Genus=="Ruminococcus2"] <- 6

plot.data$linetype[plot.data$Genus=="Hathewaya"] <- 3
plot.data$linetype[plot.data$Genus=="Anaerobacter"] <- 2
plot.data$linetype[plot.data$Genus=="Clostridium_sensu_stricto"] <- 1

plot.data$linetype[plot.data$Genus=="Clostridium_XI"] <- 1
plot.data$linetype[plot.data$Genus=="Peptostreptococcus"] <- 2

plot.data$linetype[plot.data$Genus=="Faecalibacterium"] <- 1
plot.data$linetype[plot.data$Genus=="Hydrogenoanaerobacterium"] <- 2
plot.data$linetype[plot.data$Genus=="Butyricicoccus"] <- 3
plot.data$linetype <- as.factor(as.character(plot.data$linetype))

plot.data$value[is.na(plot.data$value)] <- 0

# Create plot
Fig_S3 <- ggplot(plot.data, aes(x=variable, y=value)) +
  geom_line(aes(group=Genus, linetype=linetype)) +
  facet_wrap(~Family, nrow=4, ncol=4) +
  scale_y_continuous(limits=c(0,9), breaks=c(0,2,4,6,8)) +
  theme_bw() +
  theme(
    axis.text.x=element_text(color="black", angle=90, size=9, hjust=1, vjust=0.5),
    axis.text.y=element_text(color="black", size=9),
    axis.title=element_text(color="black", size=11),
    legend.position="none",
    panel.grid.minor=element_blank(),
    strip.background=element_rect(fill="beige"),
    strip.text=element_text(color=c("black"))
  ) +
  labs(x="Site", y="Number of individuals")

ggsave("~/Fig_S3.jpg", Fig_S3, width=6, height=5, units="in", dpi=300)
rm(plot.data, families) 

#### (done) FIGURE S4 - UNIQUE ASVs AT EACH SITE ####
# (a) Boxplot showing number of unique ASVs by site.
Fig_S4a <- ggplot(df.uniq.to.site.prop, aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable)) +
  theme_bw() +
  scale_fill_manual(values=pal(7)) +
  theme(
    axis.text.x=element_text(color="black", angle=90, size=9, hjust=1, vjust=0.5),
    axis.text.y=element_text(color="black", size=9),
    axis.title=element_text(color="black", size=11),
    legend.position="none",
    plot.tag=element_text(face="bold")
  ) +
  labs(x="Site", y="Percent of ASVs unique to site", tag="a")

# (b) Faceted line graphs showing taxon affilitations.
# Add taxonomy information.
plot.data <- ids.genus.unique
plot.data$Genus <- rownames(plot.data)
plot.data <- merge(plot.data, as.data.frame(cbind(tax_table(tax_glom(pseq.rarefied, taxrank="Genus")))), by="Genus")

# Subset to taxa present in at least three individuals.
plot.data <- subset(plot.data, prev > 3)
plot.data$prev <- NULL
plot.data$Kingdom <- NULL
plot.data$Class <- NULL
plot.data$Order <- NULL

plot.data <- melt(plot.data)
plot.data$Genus <- factor(plot.data$Genus)

# Rename variables for plotting.
levels(plot.data$variable)[levels(plot.data$variable)=="duodenum"] <- "duo"
levels(plot.data$variable)[levels(plot.data$variable)=="jejunum"] <- "jej"
levels(plot.data$variable)[levels(plot.data$variable)=="ileum"] <- "ile"
levels(plot.data$variable)[levels(plot.data$variable)=="caecum"] <- "cae"
levels(plot.data$variable)[levels(plot.data$variable)=="asc.colon"] <- "asc"
levels(plot.data$variable)[levels(plot.data$variable)=="des.colon"] <- "des"
levels(plot.data$variable)[levels(plot.data$variable)=="feces"] <- "fec"
plot.data$Genus <- factor(plot.data$Genus)

# Remove two taxa from the plot. (This is a space-saving mechanism; both taxa are unlikely to have a functional role in the gut microbiome.
plot.data <- subset(plot.data, Genus != "Ralstonia" & Genus !="Parabacteroides")
plot.data$Genus <- factor(plot.data$Genus)
plot.data$value[is.na(plot.data$value)] <- 0

Fig_S4b <- ggplot(plot.data, aes(x=variable, y=value)) +
  geom_line(group=1) +
  facet_wrap(~Genus, ncol=4, nrow=3) +
  scale_y_continuous(limits=c(0,9), breaks=c(0,2,4,6,8)) +
  theme_bw() +
  theme(
    axis.text.x=element_text(color="black", angle=90, size=9, hjust=1, vjust=0.5),
    axis.text.y=element_text(color="black", size=9),
    axis.title=element_text(color="black", size=11),
    legend.position="none",
    panel.grid.minor=element_blank(),
    strip.background=element_rect(fill="beige"),
    strip.text=element_text(color=c("black")),
    plot.tag=element_text(face="bold")
  ) +
  labs(x="Site", y="Number of individuals", tag="b")

# Assemble figure. #
grid.arrange(Fig_S4a, Fig_S4b, ncol=2, layout_matrix=rbind(c(NA,2),c(1,2),c(1,2),c(NA,2)), widths=c(0.25,0.75))
Fig_S4 <- arrangeGrob(Fig_S4a, Fig_S4b, ncol=2, layout_matrix=rbind(c(NA,2),c(1,2),c(1,2),c(NA,2)), widths=c(0.25,0.75))
ggsave("~/Fig_S4.jpg", Fig_S4, width=6.5, height=4, units="in", dpi=300)
rm(plot.data)

#### (done) FIGURE S5 - dbRDA with UNIFRAC DISTANCE ###########
# Extract data from dbRDA model.
dbrda.wUF_fort <- fortify(dbrda.wUF, display="sites")
dbrda.wUF_fort <- cbind(bg.metadata.nf, dbrda.wUF_fort)
desired_order <- c("small.intestine","large.intestine")
dbrda.wUF_fort$Segment <- factor(dbrda.wUF_fort$Segment, levels=desired_order)

# Calculate standard ellipse.
plot.new()
dbrda.wUF.ellipse <- ordiellipse(dbrda.wUF, bg.metadata.nf$Segment, display="sites", 
                                kind="sd", conf=0.75, label=T)
dbrda.wUF.ell.data <- data.frame()
for(g in levels(bg.metadata.nf$Segment)){
  dbrda.wUF.ell.data <- rbind(dbrda.wUF.ell.data,
                             cbind(as.data.frame(with(bg.metadata.nf[dbrda.wUF_fort$Segment==g,],
                                                      veganCovEllipse(dbrda.wUF.ellipse[[g]]$cov,
                                                                      dbrda.wUF.ellipse[[g]]$center,
                                                                      dbrda.wUF.ellipse[[g]]$scale)))
                                   ,Type=g))}
rm(dbrda.wUF.ellipse)

# Rename variables.
dbrda.wUF_fort$Type <- as.character(dbrda.wUF_fort$Type)
dbrda.wUF_fort$Type[dbrda.wUF_fort$Type=="duodenum"] <- "Duodenum"
dbrda.wUF_fort$Type[dbrda.wUF_fort$Type=="jejunum"] <- "Jejunum"
dbrda.wUF_fort$Type[dbrda.wUF_fort$Type=="ileum"] <- "Ileum"
dbrda.wUF_fort$Type[dbrda.wUF_fort$Type=="caecum"] <- "Caecum"
dbrda.wUF_fort$Type[dbrda.wUF_fort$Type=="ascending.colon"] <- "Ascending Colon"
dbrda.wUF_fort$Type[dbrda.wUF_fort$Type=="descending.colon"] <- "Descending Colon"
dbrda.wUF_fort$Type <- factor(dbrda.wUF_fort$Type, levels=c("Duodenum", "Jejunum", "Ileum", "Caecum", "Ascending Colon", "Descending Colon"))

# Create first plot, with ellipses for segments. #
ind.hulls <- gg_ordiplot(dbrda.wUF, groups=bg.metadata.nf$SegInd, hull=TRUE, ellipse=FALSE, label=TRUE, pt.size=2)

dbrda.wUF.plot.a <- ggplot(dbrda.wUF_fort, aes(x=CAP1, y=CAP2)) +
  geom_vline(xintercept = 0, color="grey") +
  geom_hline(yintercept = 0, color="grey") +
  geom_path(data=ind.hulls[["df_hull"]], aes(x=1.96*x, y=3.48*y, linetype=Group), size=0.5, color="grey") +
  scale_linetype_manual("Individual", values=rep("solid", 20)) +
  geom_point(aes(color = Individual, shape=Type), size=2, fill="white") +
  geom_label(aes(x=0.50236, y=-2.399e-05), color="black", label="LI") +
  geom_label(aes(x=-0.51969, y=2.481e-05), color="black", label="SI") +
  new_scale("linetype") +
  geom_path(data=dbrda.wUF.ell.data, aes(x=CAP1, y=CAP2, linetype=Type), color="black", size=0.5) +
  scale_linetype_manual(values=c("solid","solid"), guide=FALSE) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        axis.text=element_text(color="black", size = 8, vjust=1),
        axis.title=element_text(size=9),
        plot.title=element_text(size=10),
        strip.background=element_blank(), 
        plot.tag=element_text(face="bold"),
        strip.text=element_blank(),
        legend.title=element_blank(),
        legend.position="none") +
  scale_color_manual(values=dbrda.indv.pal(10), guide=FALSE) +
  scale_shape_manual(values=c(22,21,24,15,16,17)) +
  labs(x="dbRDA1 [23.6%]", y="dbRDA2 [7.5%]", tag="a")

# Extract hull information for sites.
site.hulls <- gg_ordiplot(dbrda.wUF, groups=as.character(bg.metadata.nf$Type), hull=TRUE, ellipse=FALSE, label=TRUE, pt.size=2)
site.hulls[["df_hull"]]$Group <- as.character(site.hulls[["df_hull"]]$Group)
site.hulls[["df_hull"]]$Group[site.hulls[["df_hull"]]$Group=="duodenum"] <- "Duodenum"
site.hulls[["df_hull"]]$Group[site.hulls[["df_hull"]]$Group=="jejunum"] <- "Jejunum"
site.hulls[["df_hull"]]$Group[site.hulls[["df_hull"]]$Group=="ileum"] <- "Ileum"
site.hulls[["df_hull"]]$Group[site.hulls[["df_hull"]]$Group=="caecum"] <- "Caecum"
site.hulls[["df_hull"]]$Group[site.hulls[["df_hull"]]$Group=="ascending.colon"] <- "Ascending Colon"
site.hulls[["df_hull"]]$Group[site.hulls[["df_hull"]]$Group=="descending.colon"] <- "Descending Colon"
site.hulls[["df_hull"]]$Group <- factor(site.hulls[["df_hull"]]$Group, levels=c("Duodenum", "Jejunum", "Ileum", "Caecum", "Ascending Colon", "Descending Colon"))

# Create plot.
dbrda.wUF.plot.b <- ggplot(dbrda.wUF_fort, aes(x=CAP1, y=CAP2)) +
  geom_vline(xintercept = 0, color="grey") +
  geom_hline(yintercept = 0, color="grey") +
  geom_path(data=site.hulls[["df_hull"]], aes(x=1.965*x, y=3.5*y, color=Group), size=0.5) +
  geom_point(aes(color = Type, shape=Type), size=2, fill="white", alpha=1) +
  annotate("text", x=-1.1, y=-0.8, hjust=0, color="#00A2BD", label="duodenum", size=3) +
  annotate("text", x=-0.88, y=0.3, hjust=0, color="#D544D1", label="jejunum", size=3) +
  annotate("text", x=-0.55, y=1, hjust=0, color="#B97700", label="ileum", size=3) +
  annotate("text", x=0.6, y=2.48, hjust=0, color="#00A358", label="caecum", size=3) +
  annotate("text", x=0.1, y=1.2, hjust=0, color="#4C7FEE", label="asc.colon", size=3) +
  annotate("text", x=0.6, y=-1, hjust=0, color="#E14F77", label="des.colon", size=3) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        axis.text=element_text(color="black", size = 8, vjust=1),
        axis.title=element_text(size=9),
        plot.title=element_text(size=10),
        strip.background=element_blank(), 
        plot.tag=element_text(face="bold"),
        strip.text=element_blank(),
        legend.title=element_blank(),
        legend.position="none") +
  scale_color_manual(values=dbrda.indv.pal(6), guide=FALSE) +
  scale_shape_manual(values=c(22,21,24,15,16,17)) +
  labs(x="dbRDA1 [23.6%]", y="dbRDA2 [7.5%]", tag="b")

grid.arrange(dbrda.wUF.plot.a, dbrda.wUF.plot.b, dbrda.legend,
             layout_matrix=rbind(c(1,2),c(3)), heights=c(0.4,0.1))

Fig_S5 <- arrangeGrob(dbrda.wUF.plot.a, dbrda.wUF.plot.b, dbrda.legend,
                                      layout_matrix=rbind(c(1,2),c(3)), heights=c(0.4,0.1))

ggsave("~/Fig_S5.jpg", Fig_S5, width=6.5, height=3.5, units="in", dpi=300)
rm(dbrda.wUF.plot.a, dbrda.wUF.plot.b, dbrda.legend, site.hulls, ind.hulls)

#### (done) FIGURE S6 - RANDOM FOREST RESULTS ####
# Sort graph results by decreasing Gini scores.
rf.graph.data.seg <- rf.results[["Gini_scores"]][["Segment"]][order(-rf.results[["Gini_scores"]][["Segment"]]$MeanDecreaseGini),]
rf.graph.data.indv <- rf.results[["Gini_scores"]][["Individual"]][order(-rf.results[["Gini_scores"]][["Individual"]]$MeanDecreaseGini),]
rf.graph.data.type <- rf.results[["Gini_scores"]][["Type"]][order(-rf.results[["Gini_scores"]][["Type"]]$MeanDecreaseGini),]

# Insert ASVID as a row name.
rf.graph.data.seg$ASVID <- rownames(rf.graph.data.seg)
rf.graph.data.indv$ASVID <- rownames(rf.graph.data.indv)
rf.graph.data.type$ASVID <- rownames(rf.graph.data.type)

# Order plots by ASVID.
rf.graph.data.seg$ASVID <- forcats::fct_inorder(rf.graph.data.seg$ASVID)
rf.graph.data.indv$ASVID <- forcats::fct_inorder(rf.graph.data.indv$ASVID)
rf.graph.data.type$ASVID <- forcats::fct_inorder(rf.graph.data.type$ASVID)

# Determine which grouping (segment: small/large, individual: identity, site: intestinal site) each taxon is most abundant in.
rf.graph.data.seg$Max <- as.factor(colnames(rf.graph.data.seg[,8:9])[max.col(rf.graph.data.seg[,8:9], ties.method="first")])
rf.graph.data.indv$Max <- as.factor(colnames(rf.graph.data.indv[,8:17])[max.col(rf.graph.data.indv[,8:17], ties.method="first")])
rf.graph.data.type$Max <- as.factor(colnames(rf.graph.data.type[,8:14])[max.col(rf.graph.data.type[,8:14], ties.method="random")])
rf.graph.data.type$Max <- factor(rf.graph.data.type$Max, levels=c("duodenum", "jejunum", "ileum", "caecum", "ascending.colon", "descending.colon", "feces"))

# Plot for segment.
rf.plot.seg <- ggplot(rf.graph.data.seg, aes(x=MeanDecreaseGini, y=ASVID)) + 
  geom_point(aes(color=Max, size=Overall)) + 
  scale_y_discrete(limits=rev(levels(rf.graph.data.seg$ASVID)), labels=rev(rf.graph.data.seg$Genus)) +
  scale_color_manual(name="Segment", values=c("orangered1", "deepskyblue3")) +
  scale_size_continuous(name="Abundance (%)") +
  annotate("text", x=Inf, y=-Inf, hjust=1.1, vjust=-0.5, label="Model accuracy:\n88.06%", size=3) +
  theme_bw() +
  theme(legend.position=c("bottom"), 
        legend.justification=c("center"),
        legend.title=element_text(size=8, color="black"),
        legend.text=element_text(size=7, color="black"),
        plot.tag=element_text(face="bold"),
        axis.line = element_line(colour = "black", size=0.5),
        axis.text=element_text(size=7, color="black"),
        panel.grid.minor=element_blank(),
        axis.title=element_text(size=8, color="black")) +
  guides(colour=guide_legend(ncol=1, title.position="top", title.hjust=0.5, order=1)) +
  guides(size=guide_legend(ncol=1, title.position="top", title.hjust=0.5, order=0)) +
  labs(tag="a", y=NULL, x="Mean decrease Gini coefficient")

# Plot for individual.
rf.plot.indv <- ggplot(rf.graph.data.indv, aes(x=MeanDecreaseGini, y=ASVID)) + 
  geom_point(aes(color=Max, size=Overall)) + 
  scale_y_discrete(limits=rev(levels(rf.graph.data.indv$ASVID)), labels=rev(rf.graph.data.indv$Genus)) +
  scale_color_discrete(name="Individual", guide="none") +
  scale_size_continuous(name="Abundance (%)", breaks=c(1,3,5)) +
  annotate("text", x=Inf, y=-Inf, hjust=1.1, vjust=-0.5, label="Model accuracy:\n89.55%", size=3) +
  theme_bw() +
  theme(legend.position=c("bottom"), 
        legend.justification=c("center"),
        legend.title=element_text(size=8, color="black"),
        legend.text=element_text(size=7, color="black"),
        plot.tag=element_text(face="bold"),
        axis.line = element_line(colour = "black", size=0.5),
        axis.text=element_text(size=7, color="black"),
        panel.grid.minor=element_blank(),
        plot.title=element_text(size=8, color="black"),
        axis.title=element_text(size=8, color="black")) +
  #guides(colour=guide_legend(ncol=1, title.position="top", title.hjust=0.5, order=1)) +
  guides(size=guide_legend(ncol=1, title.position="top", title.hjust=0.5, order=0)) +
  labs(tag="b", y=NULL, x="Mean decrease Gini coefficient")

# Plot for intestinal site.
rf.plot.type <- ggplot(rf.graph.data.type, aes(x=MeanDecreaseGini, y=ASVID)) + 
  geom_point(aes(color=Max, size=Overall)) + 
  scale_y_discrete(limits=rev(levels(rf.graph.data.type$ASVID)), labels=rev(rf.graph.data.type$Genus)) +
  scale_color_discrete(name="Site") +
  scale_size_continuous(name="Abundance (%)") +
  annotate("text", x=Inf, y=-Inf, hjust=1.1, vjust=-0.5, label="Model accuracy:\n26.87%", size=3) +
  theme_bw() +
  theme(legend.position=c("bottom"), 
        legend.justification=c("center"),
        legend.title=element_text(size=8, color="black"),
        legend.text=element_text(size=7, color="black"),
        plot.tag=element_text(face="bold"),
        axis.line = element_line(colour = "black", size=0.5),
        axis.text=element_text(size=7, color="black"),
        panel.grid.minor=element_blank(),
        axis.title=element_text(size=8, color="black")) +
  guides(colour=guide_legend(ncol=2, title.position="top", title.hjust=0.5)) +
  guides(size=guide_legend(ncol=1, title.position="top", title.hjust=0.5)) +
  labs(tag="c", y=NULL, x="Mean decrease Gini coefficient")

# Make every plot the same height.
plots <- list(rf.plot.seg, rf.plot.indv, rf.plot.type)
grobs <- list()
heights <- list()

for (i in 1:length(plots)){
  grobs[[i]] <- ggplotGrob(plots[[i]])
  heights[[i]] <- grobs[[i]]$heights[1:12]
}

maxheight <- do.call(grid::unit.pmax, heights)

# Assign the max height the grobs.
for (i in 1:length(grobs)){
  grobs[[i]]$heights[1:12] <- as.list(maxheight)
}

# Arrange save graph.
layout_matrix <- rbind(c(1,1,2,2), c(NA,3,3,NA))
grid.arrange(grobs[[1]], grobs[[2]], grobs[[3]], nrow=1)
Fig_S6 <- arrangeGrob(grobs[[1]], grobs[[2]], grobs[[3]], nrow=1)
ggsave("~/Fig_S6.jpg", Fig_S6, width=9, height=5, units="in", dpi=300)

rm(rf.plot.seg, rf.plot.indv, rf.plot.type, rf.graph.data.indv, rf.graph.data.seg, rf.graph.data.type, grobs, plots, heights, maxheight, layout_matrix)

#### (done) FIGURE S7 - RICHNESS AND DIVERSITY, FACETED BY INDIVIDUAL ####
graph.data <- bg.metadata

# Rename individuals to say `Coyote #`
levels(graph.data$Individual)[levels(graph.data$Individual)=="34"] <- "Coyote 34"
levels(graph.data$Individual)[levels(graph.data$Individual)=="37"] <- "Coyote 37"
levels(graph.data$Individual)[levels(graph.data$Individual)=="38"] <- "Coyote 38"
levels(graph.data$Individual)[levels(graph.data$Individual)=="39"] <- "Coyote 39"
levels(graph.data$Individual)[levels(graph.data$Individual)=="40"] <- "Coyote 40"
levels(graph.data$Individual)[levels(graph.data$Individual)=="42"] <- "Coyote 42"
levels(graph.data$Individual)[levels(graph.data$Individual)=="46"] <- "Coyote 46"
levels(graph.data$Individual)[levels(graph.data$Individual)=="47"] <- "Coyote 47"
levels(graph.data$Individual)[levels(graph.data$Individual)=="63"] <- "Coyote 63"
levels(graph.data$Individual)[levels(graph.data$Individual)=="65"] <- "Coyote 65"

# Plot species richness per intestinal site, faceted by individual.
Fig_S7a <- ggplot(graph.data, aes(x=Type, y=Observed)) +
  geom_point() + geom_line(group=1) +
  facet_wrap(~Individual) +
  scale_x_discrete(labels=c('Duodenum','Jejunum','Ileum','Caecum','Asc. Colon',
                            'Des. Colon', 'Feces')) +
  theme_bw() +
  labs(x=NULL, y="ASV Richness", tag="a") +
  theme(panel.grid=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        axis.text.x=element_text(size=8, color="black", angle=90, hjust=1),
        axis.title=element_text(size=9, color="black"),
        strip.text=element_text(size=8, color="black"),
        strip.background=element_blank(),
        plot.tag=element_text(face="bold"))

# Plot Shannon diversity per intestinal site, faceted by individual.
Fig_S7b <- ggplot(graph.data, aes(x=Type, y=Shannon)) +
  geom_point() + geom_line(group=1) +
  facet_wrap(~Individual) +
  scale_x_discrete(labels=c('Duodenum','Jejunum','Ileum','Caecum','Asc. Colon',
                            'Des. Colon', 'Feces')) +
  theme_bw() +
  labs(x=NULL, tag="b", y="Shannon index") +
  theme(panel.grid=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        axis.text.x=element_text(size=8, color="black", angle=90, hjust=1),
        axis.title=element_text(size=9, color="black"),
        strip.text=element_text(size=8, color="black"),
        strip.background=element_blank(),
        plot.tag=element_text(face="bold"))

# Assemble figures and save.
grid.arrange(Fig_S7a, Fig_S7b, ncol=2)
gs <- list(Fig_S7a, Fig_S7b)

Fig_S7 <- arrangeGrob(grobs=gs, ncol=2)
ggsave("Fig_S7.jpg", plot=Fig_S7, width=9, height=4, dpi=300)
rm(graph.data)

#### (done) FIGURE S8 - TAXONOMIC BIAS in FECES-UNDETECTED TAXA ####
# Show which classes are more likely to be feces-undetected.
Fig_S8a <- ggplot(taxon.bias.stats, aes(x=variable, y=value)) + geom_boxplot(outlier.shape = NA, fill="grey") +
  theme_bw() +
  labs(x="Class", y="Percent of ASVs that were\nfecal-undetected", tag="a") +
  theme(axis.text.x=element_text(color="black", angle=45, size=9, hjust=1, vjust=1),
        axis.text.y=element_text(size=8, color="black"),
        axis.title=element_text(size=9, color="black"),
        plot.tag=element_text(color="black", face="bold"))

# Show the across-intestinal abundances of the two most undetected taxa (Alphaproteobacteria, Verrucomicrobia).
Fig_S8b <- ggplot(taxon.bias.Alpha.Verruco, aes(x=Type, y=value, fill=variable)) + 
  geom_boxplot() + 
  scale_y_continuous(trans="log10") + 
  facet_wrap(~variable) +
  theme_bw() +
  scale_fill_manual(values=c("darkorange1", "chartreuse3")) +
  theme(
    axis.text.x=element_text(color="black", angle=45, size=9, hjust=1, vjust=1),
    axis.text.y=element_text(color="black", size=9),
    axis.title=element_text(size=9, color="black"),
    legend.position="none",
    panel.grid.minor=element_blank(),
    strip.background=element_rect(fill="beige"),
    strip.text=element_text(color=c("black")),
    plot.tag=element_text(face="bold")) +
  labs(x="Site", y="log Abundance (%)", tag="b")

# Combine figures and save.
grid.arrange(Fig_S8a, Fig_S8b, nrow=2)
Fig_S8 <- arrangeGrob(Fig_S8a, Fig_S8b, nrow=2)
ggsave("~/Fig_S8.jpg", Fig_S8, width=4, height=5, units="in", dpi=300)

#### (done) FIGURE S9 - ASV ABUNDANCES SCATTER PLOT MATRIX ####
# This code inherits the object `means` from the analysis_code R workspace.
# Produce and save a scatter plot for each pairwise comparison of taxa means.
# Each scatter plot uses log-transformed axes. 
A <- ggplot(means, aes(duodenum, jejunum)) + geom_point(size=0.5) + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") + theme_bw() + labs(x=NULL,y=NULL) + theme(panel.grid=element_blank(), axis.text=element_blank())
B <- ggplot(means, aes(duodenum, ileum)) + geom_point(size=0.5) + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") + theme_bw() + labs(x=NULL,y=NULL) + theme(panel.grid=element_blank(), axis.text=element_blank())
D <- ggplot(means, aes(duodenum, caecum)) + geom_point(size=0.5) + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") + theme_bw() + labs(x=NULL,y=NULL) + theme(panel.grid=element_blank(), axis.text=element_blank())
G <- ggplot(means, aes(duodenum, asc.colon)) + geom_point(size=0.5) + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") + theme_bw() + labs(x=NULL,y=NULL) + theme(panel.grid=element_blank(), axis.text=element_blank())
K <- ggplot(means, aes(duodenum, des.colon)) + geom_point(size=0.5) + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") + theme_bw() + labs(x=NULL,y=NULL) + theme(panel.grid=element_blank(), axis.text=element_blank())
P <- ggplot(means, aes(duodenum, feces)) + geom_point(size=0.5) + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") + theme_bw() + labs(x=NULL,y=NULL) + theme(panel.grid=element_blank(), axis.text=element_blank())

C <- ggplot(means, aes(jejunum, ileum)) + geom_point(size=0.5) + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") + theme_bw() + labs(x=NULL,y=NULL) + theme(panel.grid=element_blank(), axis.text=element_blank())
E <- ggplot(means, aes(jejunum, caecum)) + geom_point(size=0.5) + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") + theme_bw() + labs(x=NULL,y=NULL) + theme(panel.grid=element_blank(), axis.text=element_blank())
H <- ggplot(means, aes(jejunum, asc.colon)) + geom_point(size=0.5) + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") + theme_bw() + labs(x=NULL,y=NULL) + theme(panel.grid=element_blank(), axis.text=element_blank())
L <- ggplot(means, aes(jejunum, des.colon)) + geom_point(size=0.5) + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") + theme_bw() + labs(x=NULL,y=NULL) + theme(panel.grid=element_blank(), axis.text=element_blank())
Q <- ggplot(means, aes(jejunum, feces)) + geom_point(size=0.5) + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") + theme_bw() + labs(x=NULL,y=NULL) + theme(panel.grid=element_blank(), axis.text=element_blank())

FX <- ggplot(means, aes(ileum, caecum)) + geom_point(size=0.5) + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") + theme_bw() + labs(x=NULL,y=NULL) + theme(panel.grid=element_blank(), axis.text=element_blank())
I <- ggplot(means, aes(ileum, asc.colon)) + geom_point(size=0.5) + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") + theme_bw() + labs(x=NULL,y=NULL) + theme(panel.grid=element_blank(), axis.text=element_blank())
M <- ggplot(means, aes(ileum, des.colon)) + geom_point(size=0.5) + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") + theme_bw() + labs(x=NULL,y=NULL) + theme(panel.grid=element_blank(), axis.text=element_blank())
R <- ggplot(means, aes(ileum, feces)) + geom_point(size=0.5) + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") + theme_bw() + labs(x=NULL,y=NULL) + theme(panel.grid=element_blank(), axis.text=element_blank())

J <- ggplot(means, aes(caecum, asc.colon)) + geom_point(size=0.5) + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") + theme_bw() + labs(x=NULL,y=NULL) + theme(panel.grid=element_blank(), axis.text=element_blank())
N <- ggplot(means, aes(caecum, des.colon)) + geom_point(size=0.5) + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") + theme_bw() + labs(x=NULL,y=NULL) + theme(panel.grid=element_blank(), axis.text=element_blank())
S <- ggplot(means, aes(caecum, feces)) + geom_point(size=0.5) + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") + theme_bw() + labs(x=NULL,y=NULL) + theme(panel.grid=element_blank(), axis.text=element_blank())

O <- ggplot(means, aes(asc.colon, des.colon)) + geom_point(size=0.5) + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") + theme_bw() + labs(x=NULL,y=NULL) + theme(panel.grid=element_blank(), axis.text=element_blank())
TX <- ggplot(means, aes(asc.colon, feces)) + geom_point(size=0.5) + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") + theme_bw() + labs(x=NULL,y=NULL) + theme(panel.grid=element_blank(), axis.text=element_blank())

U <- ggplot(means, aes(des.colon, feces)) + geom_point(size=0.5) + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") + theme_bw() + labs(x=NULL,y=NULL) + theme(panel.grid=element_blank(), axis.text=element_blank())

# Plot and save scatter plots to create a scatter plot matrix.  
grid.arrange(A, blank, blank, blank, blank, blank,
             B, C, blank, blank, blank, blank,
             D, E, FX, blank, blank, blank,
             G, H, I, J, blank, blank,
             K, L, M, N, O, blank,
             P, Q, R, S, TX, U,
             ncol=6, nrow=6)  

gs <- list(A, blank, blank, blank, blank, blank,
           B, C, blank, blank, blank, blank,
           D, E, FX, blank, blank, blank,
           G, H, I, J, blank, blank,
           K, L, M, N, O, blank,
           P, Q, R, S, TX, U)

Fig_S9 <- arrangeGrob(grobs=gs, ncol=6, nrow=6)

# Export figure. Correlation coefficients were manually added to the graphs using photo editing software.
ggsave("Fig_S9.jpg", plot = Fig_S9, width=6, height=6, units="in")
rm(gs, A, B, C, D, E, FX, G, H, I, J, K, L, M, N, O, P, Q, R, S, TX, U)