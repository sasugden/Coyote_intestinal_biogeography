# Define some variables that will be used later...
individuals <- levels(bg.metadata$Individual) # All individuals
ind.ids <- c("34", "38", "40", "42", "46", "47", "63", "65") # Individuals with fecal samples
desired_order_type <- c("duodenum", "jejunum", "ileum", "caecum", "ascending.colon", "descending.colon", "feces") # for sorting factors

#### SEQUENCE PRE-PROCESSING ##########
# Identify contaminants using decontam.
sample_data(pseq.original)$is.neg <- sample_data(pseq.original)$Type == "control" # Create logical identifying controls.
contamdf.prev <- isContaminant(pseq.original, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
contamdf.prev <- subset(contamdf.prev, contaminant==TRUE) # Subset to true contaminants
pseq.original <- subset_taxa(pseq.original, !(taxa_names(pseq.original) %in% as.character(rownames(contamdf.prev)))) # Remove.
sample_data(pseq.original)$is.neg <- NULL

# Remove negative controls from phyloseq object.
pseq.original <- subset_samples(pseq.original, Type !="control")

# Store taxa names as a 'refseq' object and rename taxa sequentially.
dna <- Biostrings::DNAStringSet(taxa_names(pseq.original))
names(dna) <- taxa_names(pseq.original)
pseq.original <- merge_phyloseq(pseq.original, dna)
rm(dna)

# Remove singletons prior to rarefaction.
pseq.rarefied <- prune_taxa(taxa_sums(pseq.original) > 1, pseq.original)

# Rarefy samples to the minimum remaining library size (8,279 reads). Average across 1,000 rarefactions.
rarefaction.average <- list()
for(i in 1:1000){
  temp.rarefy <- rarefy_even_depth(pseq.rarefied, 
                                   sample.size = min(sample_sums(pseq.rarefied)),
                                   replace = FALSE,
                                   trimOTUs = FALSE)
  rarefaction.average[[i]] <- as.data.frame(otu_table(temp.rarefy))
}

dfAvg <- Reduce("+", rarefaction.average)/length(rarefaction.average)
dfAvg <- round(dfAvg, 0)
dfAvg <- otu_table(dfAvg, taxa_are_rows=FALSE)

# Replace feature table in phyloseq object with new, rarefied feature table.
otu_table(pseq.rarefied) <- dfAvg
pseq.rarefied <- prune_taxa(taxa_sums(pseq.rarefied) > 0, pseq.rarefied)
rm(temp.rarefy, rarefaction.average, dfAvg)
taxa_names(pseq.rarefied) <- paste0("ASV", seq(ntaxa(pseq.rarefied)))

# Subset the original (unrarefied) phyloseq object to only include taxa that "survived" rarefaction.
# This phyloseq object will be the basis for all CLR-transformed analyses.
pseq.culled <- subset_taxa(pseq.original, taxa_names(pseq.original) %in% refseq(pseq.rarefied))
taxa_names(pseq.culled) <- taxa_names(pseq.rarefied)

# This `for` loop replaces all `NA` values in the taxonomy table with the placeholder `Uncl_(Higher Taxa)`.
# For example, a Firmicutes ASV not assigned to a class becomes "Uncl_Firmicutes"
# This is to avoid issues with duplicated taxon names (`NA`), inadvertently subsetting taxa, and to offer better taxon specificity.
pseq.rarefied.sub <- pseq.rarefied
for(i in 1:nrow(tax_table(pseq.rarefied.sub))){
  for(j in 2:6){
    if(is.na(tax_table(pseq.rarefied.sub)[i,j])==TRUE){
      if(substr(tax_table(pseq.rarefied.sub)[i,j-1], 1, 4)=="Uncl"){
        tax_table(pseq.rarefied.sub)[i,j] <- tax_table(pseq.rarefied.sub)[i,j-1]}
      else {
      tax_table(pseq.rarefied.sub)[i,j] <- paste0("Uncl_", tax_table(pseq.rarefied.sub)[i,j-1])}}
    }}

pseq.culled.sub <- pseq.culled
for(i in 1:nrow(tax_table(pseq.culled.sub))){
  for(j in 2:6){
    if(is.na(tax_table(pseq.culled.sub)[i,j])==TRUE){
      if(substr(tax_table(pseq.culled.sub)[i,j-1], 1, 4)=="Uncl"){
        tax_table(pseq.culled.sub)[i,j] <- tax_table(pseq.culled.sub)[i,j-1]}
      else {
        tax_table(pseq.culled.sub)[i,j] <- paste0("Uncl_", tax_table(pseq.culled.sub)[i,j-1])}}
  }}

#### CALCULATE RAREFACTION CURVES ##########
# Define rarefaction curve calculation function (obtained from https://github.com/joey711/phyloseq/issues/143)
calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

# Calculate rarefaction curves.
rarefaction_curve_data <- calculate_rarefaction_curves(pseq.rarefied, c('Observed', 'Shannon'), rep(c(1, 10, 100, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 8245), each = 10))

# Summarize alpha diversity and add sample data.
rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
rarefaction_curve_data_summary$Sample <- gsub("X", "", rarefaction_curve_data_summary$Sample)
rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(pseq.rarefied)), by.x = 'Sample', by.y = 'row.names')
rm(rarefaction_curve_data, rarefaction_curve_data_summary)

#### ALPHA DIVERSITY ##########
# Calculate alpha-diversity indices for all samples.
bg.metadata <- cbind(bg.metadata, estimate_richness(pseq.rarefied))

# Summarize the mean, SD, and coefficient of variance for species richness and Shannon diversity at each intestinal site.
alpha.summary <- bg.metadata %>%
  dplyr::group_by(Type) %>%
  dplyr::summarise(mean.Observed=mean(Observed),
            sd.Observed=sd(Observed),
            cv.Observed=100*sd(Observed)/mean(Observed),
            mean.Shannon=mean(Shannon),
            sd.Shannon=sd(Shannon),
            cv.Shannon=100*sd(Shannon)/mean(Shannon))

# Test for significant differences in species richness and Shannon diversity between segments.
alpha.segment <- subset(bg.metadata, Type !="feces") # Limit segment-based differences to exclude fecal samples.
alpha.segment <- alpha.segment %>% # Average richness and diversity for all three sites in each segment, per individual.
  dplyr::group_by(Individual, Segment) %>%
  dplyr::summarise(mean.Observed=mean(Observed),
                   sd.Observed=sd(Observed),
                   mean.Shannon=mean(Shannon),
                   sd.Shannon=sd(Shannon))

alpha.segment %>% # Summarize the mean and SD.
  dplyr::group_by(Segment) %>%
  dplyr::summarise(mean.Observed=mean(mean.Observed),
                   sd.Observed=0.5*sqrt(sum((sd.Observed)^2)),
                   mean.Shannon=mean(mean.Shannon),
                   sd.Shannon=0.5*sqrt(sum((sd.Shannon)^2)))

leveneTest(mean.Observed ~ Segment, testdata) # p=0.0386
t.test(mean.Observed ~ Segment, testdata, paired=TRUE, var.equal=FALSE) # t=4.3103, df=9, p=0.001961

leveneTest(mean.Shannon ~ Segment, testdata) # p=0.1925
t.test(mean.Shannon ~ Segment, testdata, paired=TRUE, var.equal=TRUE) # t=5.1238, df=9, p=0.000625

# Repeated measures ANOVA for significant differneces in species richness and Shannon diversity by intestinal site.
# Followed by Tukey's post hoc test. Implemented in a mixed-model framework to control for the effect of individual identity.
# Each segment is allowed to have a different slope within each individual.
model = glmer.nb(Observed ~ Type + (1|Individual),
            data=bg.metadata)
anova(model, type=3) # F=10.961, df=6, p<0.001
summary(glht(model, linfct = mcp(Type = "Tukey")), test = adjusted("bonferroni"))

model = lmer(Shannon ~ Type + (1|Individual),
            data=bg.metadata)
anova(model, type=3) # F=5.9857, df=6, p<0.001
summary(glht(model, linfct = mcp(Type = "Tukey")), test = adjusted("BH"))

#### DIFFERENTIAL ABUNDANCE AND MEAN/SD/MEDIAN ABUNDANCE AT PHYLUM, CLASS, and FAMILY LEVELS ####
# Calculate differential abundance at the phylum, class, and family levels. The `for` loop below performs the following functions for each level:
# (1) Calculate differential abundance between segments using ALDEx2.
# (2) Calculate differential abundance among sites using ALDEx2.
# (3) Summarize the mean, SD, and median relative abundance of each taxon at each segment and site.

# Define the taxonomic ranks we are studying.
taxranks <- c("Phylum", "Class", "Family")

# Create output lists where data can be stored. 
taxon.summary.TYPE <- list()
taxon.summary.SEGMENT <- list()

# Define a `for` loop for each taxonomic rank.
for(i in c("Phylum", "Class", "Family")){
  # Agglomerate taxa to the studied rank.
  temp.physeq <- tax_glom(pseq.rarefied.sub, taxrank=i, NArm=FALSE)
  
  # Extract feature and taxonomy tables.
  temp.otu.table <- as.data.frame(t(otu_table(temp.physeq)))
  temp.tax <- as.data.frame(cbind(tax_table(temp.physeq)))
  
  # Rename features with their taxon name.
  rownames(temp.otu.table) <- temp.tax[,i]
  
  # Define intestinal site as a covariate variable for differential abundance testing.
  covariates <- as.character(sample_data(pseq.rarefied)$Type)
  
  # Calculate differential abundance among sites using the ALDEx2 method. 
  sites <- aldex(temp.otu.table, covariates, mc.samples=128, test="kw", effect=FALSE, denom="all")
  
  # Calculate summary statistics for each taxon.
    # Transform to relative abundance.
    temp.transform <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
    
    # Extract OTU table and rename features with taxon names.
    temp.transform.otu <- as.data.frame(otu_table(temp.transform))
    colnames(temp.transform.otu) <- temp.tax[,i]
    
    # Append intestinal site.
    temp.transform.otu$Type <- sample_data(temp.transform)$Type
    
    # Define a data frame where calculation outputs can be stored. Number of rows is 7 intestinal sites * 3 statistics per site (mean, sd, median).
    summary.data <- data.frame(matrix(nrow=21)) 
    
    # For each taxon...
    for(j in 1:(ncol(temp.transform.otu)-1)){
      # Calculate mean relative abundance.
      means <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Type), mean)
      means$Group.1 <- NULL
      rownames(means) <- levels(temp.transform.otu$Type)
      rownames(means) <- paste("mean", rownames(means))

      # Calculate standard deviation of relative abundance.
      sds <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Type), sd)
      sds$Group.1 <- NULL
      rownames(sds) <- levels(temp.transform.otu$Type)
      rownames(sds) <- paste("stdev", rownames(sds))
          
      # Calculate median relative abundance.
      medians <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Type), median)
      medians$Group.1 <- NULL
      rownames(medians) <- levels(temp.transform.otu$Type)
      rownames(medians) <- paste("median", rownames(medians))
      
      # Combine mean, sd, and medians into one data frame.
      temp.summary <- rbind(means, sds, medians)
      colnames(temp.summary) <- colnames(temp.transform.otu)[j]
    
      # Append results for each taxon into one data frame.
      summary.data <- cbind(summary.data, temp.summary)
  }
  
    # Remove artefactual column.
    summary.data[,1] <- NULL
    
    # Add prevalence information to summary data. For this analysis, we consider prevalence the number of individuals that contain a given taxa, NOT the number of total samples.
    prev <- merge_samples(temp.physeq, "Individual")
    prev <- microbiome::prevalence(prev, detection=0, count=TRUE)
    
    # Add taxonomic information to summary data.
    summary.data <- as.data.frame(t(summary.data))
    summary.data <- cbind(temp.tax, prev, summary.data)
    rownames(summary.data) <- summary.data[,i]
    
    # Combine differential abundance results and summary data.
    taxon.summary.TYPE[[i]] <- transform(merge(sites, summary.data, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
    rm(means, sds, medians, sites, temp.transform, temp.transform.otu, summary.data, prev)
  
  # Repeat the above calculations, except using intestinal segment instead of intestinal site as the covariate.
  covariates <- as.character(sample_data(pseq.rarefied)$Segment)
  segments <- aldex(temp.otu.table, covariates, mc.samples=128, test="t", effect=FALSE, denom="all", paired=TRUE)
  
  temp.transform <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
  temp.transform.otu <- as.data.frame(otu_table(temp.transform))
  colnames(temp.transform.otu) <- temp.tax[,i]
  temp.transform.otu$Segment <- sample_data(temp.transform)$Segment
  
  summary.data <- data.frame(matrix(nrow=6)) 
    for(j in 1:(ncol(temp.transform.otu)-1)){
    means <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Segment), mean)
    means$Group.1 <- NULL
    rownames(means) <- levels(temp.transform.otu$Segment)
    rownames(means) <- paste("mean", rownames(means))
    
    sds <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Segment), sd)
    sds$Group.1 <- NULL
    rownames(sds) <- levels(temp.transform.otu$Segment)
    rownames(sds) <- paste("stdev", rownames(sds))
    
    medians <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Segment), median)
    medians$Group.1 <- NULL
    rownames(medians) <- levels(temp.transform.otu$Segment)
    rownames(medians) <- paste("median", rownames(medians))
    
    temp.summary <- rbind(means, sds, medians)
    colnames(temp.summary) <- colnames(temp.transform.otu)[j]
    
    summary.data <- cbind(summary.data, temp.summary)
  }
  
  summary.data[,1] <- NULL
  prev <- merge_samples(temp.physeq, "Individual")
  prev <- microbiome::prevalence(prev, detection=0, count=TRUE)
  summary.data <- as.data.frame(t(summary.data))
  summary.data <- cbind(temp.tax, prev, summary.data)
  rownames(summary.data) <- summary.data[,i]
  
  taxon.summary.SEGMENT[[i]] <- transform(merge(segments, summary.data, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
  
  rm(means, sds, medians, segments, temp.transform, temp.transform.otu, summary.data, prev,
     temp.physeq, temp.otu.table, temp.summary, temp.tax)
}

# What percent of the total sample is composed of the five most dominant phyla?
phylum <- tax_glom(pseq.rarefied, taxrank="Phylum")
phylum.summary <- colSums(otu_table(phylum))
names(phylum.summary) <- tax_table(phylum)[,"Phylum"]
phylum.summary <- 100*phylum.summary/sum(phylum.summary)
phylum.summary <- subset(phylum.summary, names(phylum.summary) %in% c("Firmicutes", "Fusobacteria", "Proteobacteria", "Bacteroidetes", "Actinobacteria"))
sum(phylum.summary) # =99.75298%
rm(phylum, phylum.summary)

#### MIXED EFFECTS TAXON REGRESSIONS ####
# Perform mixed effects logistic regressions for each bacterial taxa, with 'segment' (small or large) as the response and taxon abundance as the predictor.
# Individual is included as a random effect. 
taxranks <- c("Phylum", "Class", "Order", "Family", "Genus")

# Create a list to store the data.
coeff.logits.mix.eff <- list()

# Run a regression model at each taxonomic level.
library(lme4)
for(k in 1:length(taxranks)){ 
  # Agglomerate to the taxonomic rank of interest; CLR transform one object for regression, and use another for relative abundance filtering.
  temp.physeq <- tax_glom(pseq.culled.sub, taxrank=taxranks[[k]], NArm=FALSE)
  temp.physeq.clr <- microbiome::transform(temp.physeq, "clr")
  temp.physeq.ra <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
  temp.physeq.clr <- subset_taxa(temp.physeq.clr, taxranks[[k]] !="NA")
  temp.physeq.ra <- subset_taxa(temp.physeq.ra, taxranks[[k]] !="NA")
  
  # Extract OTU table from both objects and rename it with bacterial taxa.
  temp.vegan.otu.table = as(otu_table(temp.physeq.clr), "matrix")
  if(taxa_are_rows(temp.physeq.clr)){temp.vegan.otu.table <- t(temp.vegan.otu.table)}
  temp.vegan.otu.table = as.data.frame(temp.vegan.otu.table)
  colnames(temp.vegan.otu.table) <- as.data.frame(cbind(tax_table(temp.physeq.clr)))[,k+1]
  
  temp.vegan.ra.table = as(otu_table(temp.physeq.ra), "matrix")
  if(taxa_are_rows(temp.physeq.ra)){temp.vegan.ra.table <- t(temp.vegan.ra.table)}
  temp.vegan.ra.table = as.data.frame(temp.vegan.ra.table)
  colnames(temp.vegan.ra.table) <- as.data.frame(cbind(tax_table(temp.physeq.ra)))[,k+1]
  
  # Subset the OTU tables to taxa that are only present with at least 0.1% abundance.
  col_vector <- colMeans(temp.vegan.ra.table)
  temp.vegan.ra.table <- as.data.frame(t(temp.vegan.ra.table))
  temp.vegan.ra.table$Means <- col_vector
  temp.vegan.ra.table <- subset(temp.vegan.ra.table, Means > 0.1)
  temp.vegan.ra.table$Means <- NULL
  
  temp.vegan.otu.table <- as.data.frame(t(temp.vegan.otu.table))
  temp.vegan.otu.table <- subset(temp.vegan.otu.table, rownames(temp.vegan.otu.table) %in% rownames(temp.vegan.ra.table))
  temp.vegan.otu.table <- as.data.frame(t(temp.vegan.otu.table))
  
  # Add metadata to the OTU table.
  temp.vegan.otu.table$Segment <- bg.metadata$Segment
  temp.vegan.otu.table$Individual <- bg.metadata$Individual
  temp.vegan.otu.table$Type <- bg.metadata$Type
  temp.vegan.otu.table <- subset(temp.vegan.otu.table, Type !="feces")
  
  # Calculate regression for each taxon in the table.
  coefficients <- as.data.frame(matrix(ncol=2))
  
  for(i in 1:(ncol(temp.vegan.otu.table)-3)){
    logit <- lmer(temp.vegan.otu.table[,i] ~ Segment + (1+Segment|Individual),
                   data=temp.vegan.otu.table)
    
    # Use standard error to calculate 95% confidence intervals.
    se <- sqrt(diag(vcov(logit)))
    tab <- cbind(LL = fixef(logit) - 1.96 * se, Est = fixef(logit), UL = fixef(logit) + 1.96 * se)
    
    # Store the coefficient, 95% confidence intervals, and results from the model.
    coefficients[i,1] <- colnames(temp.vegan.otu.table)[i] # Taxon name
    coefficients[i,2] <- tab[2,1] # Lower limit of confidence interval
    coefficients[i,3] <- tab[2,2] # Coefficient
    coefficients[i,4] <- tab[2,3] # Upper limit of confidence interval
    coefficients[i,5] <- AIC(logit) # AIC of model
    coefficients[i,6] <- AIC(lmer(temp.vegan.otu.table[,i]~ 1 + (1+Segment|Individual),
                                  data=temp.vegan.otu.table))-AIC(logit) # Change in AIC (compared to null model)
  }
  
  # Name table and store in the `coeff.logits.mix.eff` object.
  colnames(coefficients) <- c(taxranks[[k]],"2.5","coeff","97.5", "AIC", "Delta_AIC")
  coeff.logits.mix.eff[[k]] <- coefficients
  rm(coefficients, tab, se, temp.vegan.otu.table, col_vector, temp.physeq, k, i, logit)
}
names(coeff.logits.mix.eff) <- taxranks
rm(taxranks)

#### OBTAIN FEATURE TABLE TO USE IN MaAsLIN ####
# Extract a feature table of family-level relative abundances.
maaslin.input <- microbiome::transform(tax_glom(pseq.culled.sub, taxrank="Genus"), "clr")
taxa_names(maaslin.input) <- tax_table(maaslin.input)[,"Genus"]
maaslin.input <- as.data.frame(otu_table(maaslin.input))
for(i in 1:ncol(maaslin.input)){maaslin.input[,i] <- scales::rescale(maaslin.input[,i], to=c(0,1))}

# Append sample metadata for `segment` and `individual` and rename columns appropriately.
maaslin.input <- cbind(bg.metadata$Segment, bg.metadata$Individual, maaslin.input)
colnames(maaslin.input) <- c("Segment", "Individual", colnames(maaslin.input[,c(3:ncol(maaslin.input))]))

# Transpose data frame and output this file as a tab-delimited file. 
maaslin.input <- as.data.frame(t(maaslin.input))
write.table(maaslin.input, file="~/maaslin.input.csv", sep="\t", quote=FALSE)
rm(maaslin.input)

#### EXTRACT ASV SEQUENCES FROM ABUNDANT TAXA FOR USE IN A BLAST SEARCH ####
# Subset phyloseq object to Helicobacter.
helicobacter <- subset_taxa(transform_sample_counts(pseq.rarefied, function(x) 100*x/sum(x)), Genus=="Helicobacter")

# Create a data frame with each ASV's mean abundance and DNA sequence.
helicobacter <- data.frame(rowMeans(as.data.frame(t(otu_table(helicobacter)))),
                           refseq(helicobacter))
colnames(helicobacter) <- c("mean_abund", "refseq")

# Create a FASTA file containing the sequences, named by ASV#.
fasta <- ShortRead(sread = DNAStringSet(helicobacter$refseq), id = BStringSet(rownames(helicobacter)))

# Export fasta file for a BLAST search and abundance information for later reference.
writeFasta(fasta, file = "~/helicobacter.fasta")
write.csv(helicobacter, "~/helicobacter.csv")
rm(helicobacter, fasta)

# Repeat for other taxa of interest.
clostridiaceae <- subset_taxa(transform_sample_counts(pseq.rarefied, function(x) 100*x/sum(x)), Family=="Clostridiaceae")
clostridiaceae <- data.frame(rowMeans(as.data.frame(t(otu_table(clostridiaceae)))),
                           refseq(clostridiaceae))
colnames(clostridiaceae) <- c("mean_abund", "refseq")
fasta <- ShortRead(sread = DNAStringSet(clostridiaceae$refseq), id = BStringSet(rownames(clostridiaceae)))
writeFasta(fasta, file = "~/clostridiaceae.fasta")
write.csv(clostridiaceae, "~/clostridiaceae.csv")
rm(clostridiaceae, fasta)

fusobacteria <- subset_taxa(transform_sample_counts(pseq.rarefied, function(x) 100*x/sum(x)), Phylum=="Fusobacteria")
fusobacteria <- data.frame(rowMeans(as.data.frame(t(otu_table(fusobacteria)))),
                             refseq(fusobacteria))
colnames(fusobacteria) <- c("mean_abund", "refseq")
fasta <- ShortRead(sread = DNAStringSet(fusobacteria$refseq), id = BStringSet(rownames(fusobacteria)))
writeFasta(fasta, file = "~/fusobacteria.fasta")
write.csv(fusobacteria, "~/fusobacteria.csv")
rm(fusobacteria, fasta)

#### 'NEW' and 'UNIQUE' TAXA AT EACH INTESTINAL SITE. ####
# This section produces six important data frames.
# df.new.to.site.prop = Number of ASVs that are first detected in each intestinal site, expressed as a percentage of the total ASVs. (Calculated separately for each individual).
# df.new.to.site.prop = Number of ASVs that are unique to each intestinal site, expressed as a percentage of the total ASVs. (Calculated separately for each individual).
# ids.family.appear = Names of the families that are first detected at each intestinal site. This table is stratified by intestinal site, with counts indicating number of individuals.
# ids.family.unique = Names of the families that aare unique to each each intestinal site. This table is stratified by intestinal site, with counts indicating number of individuals.
# ids.genus.appear = Names of the genera that are first detected at each intestinal site. This table is stratified by intestinal site, with counts indicating number of individuals.
# ids.genus.unique = Names of the genera that are unique to each intestinal site. This table is stratified by intestinal site, with counts indicating number of individuals.
# (a) Determine numbers of taxa at the ASV level ####
# Create data frames for storing data. Columns (7) are for intestinal sites; rows (10) are for individuals.
df.total.taxa <- data.frame(matrix(ncol=7, nrow=10))
df.new.to.site <- data.frame(matrix(ncol=7, nrow=10))
df.uniq.to.site <- data.frame(matrix(ncol=7, nrow=10))

# For each individual...
for(i in 1:length(individuals)){
  # Create a phyloseq object for that individual.
  temp <- subset_samples(pseq.rarefied, Individual==individuals[[i]])
  temp <- prune_taxa(taxa_sums(temp) > 0, temp)
  
  # Create independent phyloseq objects containing only the taxa present at each intestinal site.
  # "If" statement excludes individuals that are missing a particular sample (ex. no feces)
  if(i !=9){
  temp.duo <- subset_samples(temp, Type=="duodenum")
  temp.duo <- prune_taxa(taxa_sums(temp.duo) > 0, temp.duo) }
  temp.jej <- subset_samples(temp, Type=="jejunum")
  temp.jej <- prune_taxa(taxa_sums(temp.jej) > 0, temp.jej) 
  temp.ile <- subset_samples(temp, Type=="ileum")
  temp.ile <- prune_taxa(taxa_sums(temp.ile) > 0, temp.ile) 
  temp.cae <- subset_samples(temp, Type=="caecum")
  temp.cae <- prune_taxa(taxa_sums(temp.cae) > 0, temp.cae) 
  temp.asc <- subset_samples(temp, Type=="ascending.colon")
  temp.asc <- prune_taxa(taxa_sums(temp.asc) > 0, temp.asc) 
  temp.des <- subset_samples(temp, Type=="descending.colon")
  temp.des <- prune_taxa(taxa_sums(temp.des) > 0, temp.des) 
  if(i !=2 & i !=4){
  temp.fec <- subset_samples(temp, Type=="feces")
  temp.fec <- prune_taxa(taxa_sums(temp.fec) > 0, temp.fec)}
  
  # Create independent phyloseq objects based on a step-wise addition of sample sites.
  # That is, an object for "duodenum + jejunum", an object for "duodenum + jejunum + ileum", etc.
  temp.duo.jej <- subset_samples(temp, Type=="duodenum" | Type=="jejunum")
  temp.duo.jej <- prune_taxa(taxa_sums(temp.duo.jej) > 0, temp.duo.jej) 
  temp.duo.jej.ile <- subset_samples(temp, Type=="duodenum" | Type=="jejunum" | Type=="ileum")
  temp.duo.jej.ile <- prune_taxa(taxa_sums(temp.duo.jej.ile) > 0, temp.duo.jej.ile) 
  temp.duo.jej.ile.cae <- subset_samples(temp, Type=="duodenum" | Type=="jejunum" | Type=="ileum" | Type=="caecum")
  temp.duo.jej.ile.cae <- prune_taxa(taxa_sums(temp.duo.jej.ile.cae) > 0, temp.duo.jej.ile.cae) 
  temp.duo.jej.ile.cae.asc <- subset_samples(temp, Type=="duodenum" | Type=="jejunum" | Type=="ileum" | Type=="caecum" | Type=="ascending.colon")
  temp.duo.jej.ile.cae.asc <- prune_taxa(taxa_sums(temp.duo.jej.ile.cae.asc) > 0, temp.duo.jej.ile.cae.asc) 
  temp.duo.jej.ile.cae.asc.des <- subset_samples(temp, Type != "feces")
  temp.duo.jej.ile.cae.asc.des <- prune_taxa(taxa_sums(temp.duo.jej.ile.cae.asc.des) > 0, temp.duo.jej.ile.cae.asc.des) 

  # Create independent phyloseq objects for all taxa present outside of a given intestinal site.
  if(i !=9){
  temp.not.duo <- subset_samples(temp, Type !="duodenum")
  temp.not.duo <- prune_taxa(taxa_sums(temp.not.duo) > 0, temp.not.duo) }
  temp.not.jej <- subset_samples(temp, Type !="jejunum")
  temp.not.jej <- prune_taxa(taxa_sums(temp.not.jej) > 0, temp.not.jej)
  temp.not.ile <- subset_samples(temp, Type !="ileum")
  temp.not.ile <- prune_taxa(taxa_sums(temp.not.ile) > 0, temp.not.ile)
  temp.not.cae <- subset_samples(temp, Type !="caecum")
  temp.not.cae <- prune_taxa(taxa_sums(temp.not.cae) > 0, temp.not.cae)
  temp.not.asc <- subset_samples(temp, Type !="ascending.colon")
  temp.not.asc <- prune_taxa(taxa_sums(temp.not.asc) > 0, temp.not.asc)
  temp.not.des <- subset_samples(temp, Type !="descending.colon")
  temp.not.des <- prune_taxa(taxa_sums(temp.not.des) > 0, temp.not.des)
  temp.not.fec <- subset_samples(temp, Type !="feces")
  temp.not.fec <- prune_taxa(taxa_sums(temp.not.fec) > 0, temp.not.fec)
  
  # Populate a data frame giving the total number of ASVs present at each intestinal site, in this one individual.
  if(i !=9){df.total.taxa[i,1] <- ntaxa(temp.duo)}
  df.total.taxa[i,2] <- ntaxa(temp.jej)
  df.total.taxa[i,3] <- ntaxa(temp.ile)
  df.total.taxa[i,4] <- ntaxa(temp.cae)
  df.total.taxa[i,5] <- ntaxa(temp.asc)
  df.total.taxa[i,6] <- ntaxa(temp.des)
  if(i !=2 & i !=4){df.total.taxa[i,7] <- ntaxa(temp.fec)}
  
  # Populate a data frame giving the number of ASVs that appear for the first time in each intestinal site (e.g. not detected upstream)
  if(i !=9){df.new.to.site[i,2] <- length(setdiff(taxa_names(temp.jej), taxa_names(temp.duo)))}
  df.new.to.site[i,3] <- length(setdiff(taxa_names(temp.ile), taxa_names(temp.duo.jej)))
  df.new.to.site[i,4] <- length(setdiff(taxa_names(temp.cae), taxa_names(temp.duo.jej.ile)))
  df.new.to.site[i,5] <- length(setdiff(taxa_names(temp.asc), taxa_names(temp.duo.jej.ile.cae)))
  df.new.to.site[i,6] <- length(setdiff(taxa_names(temp.des), taxa_names(temp.duo.jej.ile.cae.asc)))
  if(i !=2 & i !=4){df.new.to.site[i,7] <- length(setdiff(taxa_names(temp.fec), taxa_names(temp.duo.jej.ile.cae.asc.des)))}
  
  # Populate a data frame giving the number of ASVs that are unique to each intestinal site (e.g. not present elsewhere in the intestine)
  if(i !=9){df.uniq.to.site[i,1] <- length(setdiff(taxa_names(temp.duo), taxa_names(temp.not.duo)))}
  df.uniq.to.site[i,2] <- length(setdiff(taxa_names(temp.jej), taxa_names(temp.not.jej)))
  df.uniq.to.site[i,3] <- length(setdiff(taxa_names(temp.ile), taxa_names(temp.not.ile)))
  df.uniq.to.site[i,4] <- length(setdiff(taxa_names(temp.cae), taxa_names(temp.not.cae)))
  df.uniq.to.site[i,5] <- length(setdiff(taxa_names(temp.asc), taxa_names(temp.not.asc)))
  df.uniq.to.site[i,6] <- length(setdiff(taxa_names(temp.des), taxa_names(temp.not.des)))
  if(i !=2 & i !=4){df.uniq.to.site[i,7] <- length(setdiff(taxa_names(temp.fec), taxa_names(temp.not.fec)))}
  
  rm(temp, temp.duo, temp.jej, temp.ile, temp.cae, temp.asc, temp.des, temp.fec,
     temp.duo.jej, temp.duo.jej.ile, temp.duo.jej.ile.cae, temp.duo.jej.ile.cae.asc, temp.duo.jej.ile.cae.asc.des,
     temp.not.duo, temp.not.jej, temp.not.ile, temp.not.cae, temp.not.asc, temp.not.des, temp.not.fec)
}

colnames(df.new.to.site) <- levels(bg.metadata$Type)
colnames(df.uniq.to.site) <- levels(bg.metadata$Type)

# (b) Determine names of taxa at the genus and family levels ####
pseq.glom <- tax_glom(pseq.rarefied, taxrank="Genus")
taxa_names(pseq.glom) <- as.data.frame(cbind(tax_table(pseq.glom)))$Genus

ids.genus.new.to.site <- list()
ids.genus.uniq.to.site <- list()

for(i in 1:length(individuals)){
  # Create a phyloseq object for that individual.
  temp <- subset_samples(pseq.glom, Individual==individuals[[i]])
  temp <- prune_taxa(taxa_sums(temp) > 0, temp)
  
  # Create independent phyloseq objects containing only the taxa present at each intestinal site.
  # "If" statement excludes individuals that are missing a particular sample (ex. no feces)
  if(i !=9){
    temp.duo <- subset_samples(temp, Type=="duodenum")
    temp.duo <- prune_taxa(taxa_sums(temp.duo) > 0, temp.duo) }
  temp.jej <- subset_samples(temp, Type=="jejunum")
  temp.jej <- prune_taxa(taxa_sums(temp.jej) > 0, temp.jej) 
  temp.ile <- subset_samples(temp, Type=="ileum")
  temp.ile <- prune_taxa(taxa_sums(temp.ile) > 0, temp.ile) 
  temp.cae <- subset_samples(temp, Type=="caecum")
  temp.cae <- prune_taxa(taxa_sums(temp.cae) > 0, temp.cae) 
  temp.asc <- subset_samples(temp, Type=="ascending.colon")
  temp.asc <- prune_taxa(taxa_sums(temp.asc) > 0, temp.asc) 
  temp.des <- subset_samples(temp, Type=="descending.colon")
  temp.des <- prune_taxa(taxa_sums(temp.des) > 0, temp.des) 
  if(i !=2 & i !=4){
    temp.fec <- subset_samples(temp, Type=="feces")
    temp.fec <- prune_taxa(taxa_sums(temp.fec) > 0, temp.fec)}
  
  # Create independent phyloseq objects based on a step-wise addition of sample sites.
  # That is, an object for "duodenum + jejunum", an object for "duodenum + jejunum + ileum", etc.
  temp.duo.jej <- subset_samples(temp, Type=="duodenum" | Type=="jejunum")
  temp.duo.jej <- prune_taxa(taxa_sums(temp.duo.jej) > 0, temp.duo.jej) 
  temp.duo.jej.ile <- subset_samples(temp, Type=="duodenum" | Type=="jejunum" | Type=="ileum")
  temp.duo.jej.ile <- prune_taxa(taxa_sums(temp.duo.jej.ile) > 0, temp.duo.jej.ile) 
  temp.duo.jej.ile.cae <- subset_samples(temp, Type=="duodenum" | Type=="jejunum" | Type=="ileum" | Type=="caecum")
  temp.duo.jej.ile.cae <- prune_taxa(taxa_sums(temp.duo.jej.ile.cae) > 0, temp.duo.jej.ile.cae) 
  temp.duo.jej.ile.cae.asc <- subset_samples(temp, Type=="duodenum" | Type=="jejunum" | Type=="ileum" | Type=="caecum" | Type=="ascending.colon")
  temp.duo.jej.ile.cae.asc <- prune_taxa(taxa_sums(temp.duo.jej.ile.cae.asc) > 0, temp.duo.jej.ile.cae.asc) 
  temp.duo.jej.ile.cae.asc.des <- subset_samples(temp, Type != "feces")
  temp.duo.jej.ile.cae.asc.des <- prune_taxa(taxa_sums(temp.duo.jej.ile.cae.asc.des) > 0, temp.duo.jej.ile.cae.asc.des) 
  
  # Create independent phyloseq objects for all taxa present outside of a given intestinal site.
  if(i !=9){
    temp.not.duo <- subset_samples(temp, Type !="duodenum")
    temp.not.duo <- prune_taxa(taxa_sums(temp.not.duo) > 0, temp.not.duo) }
  temp.not.jej <- subset_samples(temp, Type !="jejunum")
  temp.not.jej <- prune_taxa(taxa_sums(temp.not.jej) > 0, temp.not.jej)
  temp.not.ile <- subset_samples(temp, Type !="ileum")
  temp.not.ile <- prune_taxa(taxa_sums(temp.not.ile) > 0, temp.not.ile)
  temp.not.cae <- subset_samples(temp, Type !="caecum")
  temp.not.cae <- prune_taxa(taxa_sums(temp.not.cae) > 0, temp.not.cae)
  temp.not.asc <- subset_samples(temp, Type !="ascending.colon")
  temp.not.asc <- prune_taxa(taxa_sums(temp.not.asc) > 0, temp.not.asc)
  temp.not.des <- subset_samples(temp, Type !="descending.colon")
  temp.not.des <- prune_taxa(taxa_sums(temp.not.des) > 0, temp.not.des)
  temp.not.fec <- subset_samples(temp, Type !="feces")
  temp.not.fec <- prune_taxa(taxa_sums(temp.not.fec) > 0, temp.not.fec)

  # Populate a data frame giving the name of the genera that are appearing for the first time.
  ids.genus.new.to.site[[i]] <- data.frame()
  if(i !=9){ids.genus.new.to.site[[i]] <- cbind.fill(ids.genus.new.to.site[[i]], taxa_names(temp.duo), fill=NA)}
  if(i==9){ids.genus.new.to.site[[i]] <- cbind.fill(ids.genus.new.to.site[[i]], rep(NA, 6), fill=NA)}
  if(i !=9 & i !=10){ids.genus.new.to.site[[i]] <- cbind.fill(ids.genus.new.to.site[[i]], setdiff(taxa_names(temp.jej), taxa_names(temp.duo)), fill=NA)}
  if(i==9){ids.genus.new.to.site[[i]] <- cbind.fill(ids.genus.new.to.site[[i]], rep(NA, 4), fill=NA)}
  if(i==10){ids.genus.new.to.site[[i]] <- cbind.fill(ids.genus.new.to.site[[i]], rep(NA, 4), fill=NA)} # This is here because `setdiff(taxa_names(temp.jej), taxa_names(temp.duo))` has length of 0 for the tenth individual.
  ids.genus.new.to.site[[i]] <- cbind.fill(ids.genus.new.to.site[[i]], setdiff(taxa_names(temp.ile), taxa_names(temp.duo.jej)), fill=NA)
  ids.genus.new.to.site[[i]] <- cbind.fill(ids.genus.new.to.site[[i]], setdiff(taxa_names(temp.cae), taxa_names(temp.duo.jej.ile)), fill=NA)
  ids.genus.new.to.site[[i]] <- cbind.fill(ids.genus.new.to.site[[i]], setdiff(taxa_names(temp.asc), taxa_names(temp.duo.jej.ile.cae)), fill=NA)
  ids.genus.new.to.site[[i]] <- cbind.fill(ids.genus.new.to.site[[i]], setdiff(taxa_names(temp.des), taxa_names(temp.duo.jej.ile.cae.asc)), fill=NA)
  if(i !=2 & i !=4){ids.genus.new.to.site[[i]] <- cbind.fill(ids.genus.new.to.site[[i]], setdiff(taxa_names(temp.fec), taxa_names(temp.duo.jej.ile.cae.asc.des)), fill=NA)}
  if(i==2 | i==4){ids.genus.new.to.site[[i]] <- cbind.fill(ids.genus.new.to.site[[i]], rep(NA))}
  ids.genus.new.to.site[[i]]$init <- NULL
  colnames(ids.genus.new.to.site[[i]]) <- c("duodenum", "jejunum", "ileum", "caecum", "asc.colon", "des.colon", "feces")
  
  # Populate a data frame giving the name of genera that are unique to each intestinal site.
  ids.genus.uniq.to.site[[i]] <- data.frame()
  if(i !=9){ids.genus.uniq.to.site[[i]] <- cbind.fill(ids.genus.uniq.to.site[[i]], setdiff(taxa_names(temp.duo), taxa_names(temp.not.duo)), fill=NA)}
  if(i==9){ids.genus.uniq.to.site[[i]] <- cbind.fill(ids.genus.uniq.to.site[[i]], rep(NA, 4), fill=NA)}
  ids.genus.uniq.to.site[[i]] <- cbind.fill(ids.genus.uniq.to.site[[i]], setdiff(taxa_names(temp.jej), taxa_names(temp.not.jej)), fill=NA)
  ids.genus.uniq.to.site[[i]] <- cbind.fill(ids.genus.uniq.to.site[[i]], setdiff(taxa_names(temp.ile), taxa_names(temp.not.ile)), fill=NA)
  ids.genus.uniq.to.site[[i]] <- cbind.fill(ids.genus.uniq.to.site[[i]], setdiff(taxa_names(temp.cae), taxa_names(temp.not.cae)), fill=NA)
  ids.genus.uniq.to.site[[i]] <- cbind.fill(ids.genus.uniq.to.site[[i]], setdiff(taxa_names(temp.asc), taxa_names(temp.not.asc)), fill=NA)
  ids.genus.uniq.to.site[[i]] <- cbind.fill(ids.genus.uniq.to.site[[i]], setdiff(taxa_names(temp.des), taxa_names(temp.not.des)), fill=NA)
  if(i !=2 & i !=4){ids.genus.uniq.to.site[[i]] <- cbind.fill(ids.genus.uniq.to.site[[i]], setdiff(taxa_names(temp.fec), taxa_names(temp.not.fec)), fill=NA)}
  if(i==2 | i==4){ids.genus.uniq.to.site[[i]] <- cbind.fill(ids.genus.uniq.to.site[[i]], rep(NA))}
  ids.genus.uniq.to.site[[i]]$init <- NULL
  colnames(ids.genus.uniq.to.site[[i]]) <- c("duodenum", "jejunum", "ileum", "caecum", "asc.colon", "des.colon", "feces")
  
  rm(temp, temp.duo, temp.jej, temp.ile, temp.cae, temp.asc, temp.des, temp.fec,
     temp.duo.jej, temp.duo.jej.ile, temp.duo.jej.ile.cae, temp.duo.jej.ile.cae.asc, temp.duo.jej.ile.cae.asc.des,
     temp.not.duo, temp.not.jej, temp.not.ile, temp.not.cae, temp.not.asc, temp.not.des, temp.not.fec)
}

pseq.glom <- tax_glom(pseq.rarefied, taxrank="Family")
taxa_names(pseq.glom) <- as.data.frame(cbind(tax_table(pseq.glom)))$Family

ids.family.new.to.site <- list()
ids.family.uniq.to.site <- list()

for(i in 1:length(individuals)){
  # Create a phyloseq object for that individual.
  temp <- subset_samples(pseq.glom, Individual==individuals[[i]])
  temp <- prune_taxa(taxa_sums(temp) > 0, temp)
  
  # Create independent phyloseq objects containing only the taxa present at each intestinal site.
  # "If" statement excludes individuals that are missing a particular sample (ex. no feces)
  if(i !=9){
    temp.duo <- subset_samples(temp, Type=="duodenum")
    temp.duo <- prune_taxa(taxa_sums(temp.duo) > 0, temp.duo) }
  temp.jej <- subset_samples(temp, Type=="jejunum")
  temp.jej <- prune_taxa(taxa_sums(temp.jej) > 0, temp.jej) 
  temp.ile <- subset_samples(temp, Type=="ileum")
  temp.ile <- prune_taxa(taxa_sums(temp.ile) > 0, temp.ile) 
  temp.cae <- subset_samples(temp, Type=="caecum")
  temp.cae <- prune_taxa(taxa_sums(temp.cae) > 0, temp.cae) 
  temp.asc <- subset_samples(temp, Type=="ascending.colon")
  temp.asc <- prune_taxa(taxa_sums(temp.asc) > 0, temp.asc) 
  temp.des <- subset_samples(temp, Type=="descending.colon")
  temp.des <- prune_taxa(taxa_sums(temp.des) > 0, temp.des) 
  if(i !=2 & i !=4){
    temp.fec <- subset_samples(temp, Type=="feces")
    temp.fec <- prune_taxa(taxa_sums(temp.fec) > 0, temp.fec)}
  
  # Create independent phyloseq objects based on a step-wise addition of sample sites.
  # That is, an object for "duodenum + jejunum", an object for "duodenum + jejunum + ileum", etc.
  temp.duo.jej <- subset_samples(temp, Type=="duodenum" | Type=="jejunum")
  temp.duo.jej <- prune_taxa(taxa_sums(temp.duo.jej) > 0, temp.duo.jej) 
  temp.duo.jej.ile <- subset_samples(temp, Type=="duodenum" | Type=="jejunum" | Type=="ileum")
  temp.duo.jej.ile <- prune_taxa(taxa_sums(temp.duo.jej.ile) > 0, temp.duo.jej.ile) 
  temp.duo.jej.ile.cae <- subset_samples(temp, Type=="duodenum" | Type=="jejunum" | Type=="ileum" | Type=="caecum")
  temp.duo.jej.ile.cae <- prune_taxa(taxa_sums(temp.duo.jej.ile.cae) > 0, temp.duo.jej.ile.cae) 
  temp.duo.jej.ile.cae.asc <- subset_samples(temp, Type=="duodenum" | Type=="jejunum" | Type=="ileum" | Type=="caecum" | Type=="ascending.colon")
  temp.duo.jej.ile.cae.asc <- prune_taxa(taxa_sums(temp.duo.jej.ile.cae.asc) > 0, temp.duo.jej.ile.cae.asc) 
  temp.duo.jej.ile.cae.asc.des <- subset_samples(temp, Type != "feces")
  temp.duo.jej.ile.cae.asc.des <- prune_taxa(taxa_sums(temp.duo.jej.ile.cae.asc.des) > 0, temp.duo.jej.ile.cae.asc.des) 
  
  # Create independent phyloseq objects for all taxa present outside of a given intestinal site.
  if(i !=9){
    temp.not.duo <- subset_samples(temp, Type !="duodenum")
    temp.not.duo <- prune_taxa(taxa_sums(temp.not.duo) > 0, temp.not.duo) }
  temp.not.jej <- subset_samples(temp, Type !="jejunum")
  temp.not.jej <- prune_taxa(taxa_sums(temp.not.jej) > 0, temp.not.jej)
  temp.not.ile <- subset_samples(temp, Type !="ileum")
  temp.not.ile <- prune_taxa(taxa_sums(temp.not.ile) > 0, temp.not.ile)
  temp.not.cae <- subset_samples(temp, Type !="caecum")
  temp.not.cae <- prune_taxa(taxa_sums(temp.not.cae) > 0, temp.not.cae)
  temp.not.asc <- subset_samples(temp, Type !="ascending.colon")
  temp.not.asc <- prune_taxa(taxa_sums(temp.not.asc) > 0, temp.not.asc)
  temp.not.des <- subset_samples(temp, Type !="descending.colon")
  temp.not.des <- prune_taxa(taxa_sums(temp.not.des) > 0, temp.not.des)
  temp.not.fec <- subset_samples(temp, Type !="feces")
  temp.not.fec <- prune_taxa(taxa_sums(temp.not.fec) > 0, temp.not.fec)
  
  # Populate a data frame giving the name of the families that are appearing for the first time.
  ids.family.new.to.site[[i]] <- data.frame()
  if(i !=9){ids.family.new.to.site[[i]] <- cbind.fill(ids.family.new.to.site[[i]], taxa_names(temp.duo), fill=NA)}
  if(i==9){ids.family.new.to.site[[i]] <- cbind.fill(ids.family.new.to.site[[i]], rep(NA, 6), fill=NA)}
  if(i !=9){ids.family.new.to.site[[i]] <- cbind.fill(ids.family.new.to.site[[i]], setdiff(taxa_names(temp.jej), taxa_names(temp.duo)), fill=NA)}
  if(i==9){ids.family.new.to.site[[i]] <- cbind.fill(ids.family.new.to.site[[i]], rep(NA, 4), fill=NA)}
  ids.family.new.to.site[[i]] <- cbind.fill(ids.family.new.to.site[[i]], setdiff(taxa_names(temp.ile), taxa_names(temp.duo.jej)), fill=NA)
  ids.family.new.to.site[[i]] <- cbind.fill(ids.family.new.to.site[[i]], setdiff(taxa_names(temp.cae), taxa_names(temp.duo.jej.ile)), fill=NA)
  ids.family.new.to.site[[i]] <- cbind.fill(ids.family.new.to.site[[i]], setdiff(taxa_names(temp.asc), taxa_names(temp.duo.jej.ile.cae)), fill=NA)
  ids.family.new.to.site[[i]] <- cbind.fill(ids.family.new.to.site[[i]], setdiff(taxa_names(temp.des), taxa_names(temp.duo.jej.ile.cae.asc)), fill=NA)
  if(i !=2 & i !=4){ids.family.new.to.site[[i]] <- cbind.fill(ids.family.new.to.site[[i]], setdiff(taxa_names(temp.fec), taxa_names(temp.duo.jej.ile.cae.asc.des)), fill=NA)}
  if(i==2 | i==4){ids.family.new.to.site[[i]] <- cbind.fill(ids.family.new.to.site[[i]], rep(NA))}
  ids.family.new.to.site[[i]]$init <- NULL
  colnames(ids.family.new.to.site[[i]]) <- c("duodenum", "jejunum", "ileum", "caecum", "asc.colon", "des.colon", "feces")
  
  # Populate a data frame giving the name of families that are unique to each intestinal site.
  ids.family.uniq.to.site[[i]] <- data.frame()
  if(i !=9){ids.family.uniq.to.site[[i]] <- cbind.fill(ids.family.uniq.to.site[[i]], setdiff(taxa_names(temp.duo), taxa_names(temp.not.duo)), fill=NA)}
  if(i==9){ids.family.uniq.to.site[[i]] <- cbind.fill(ids.family.uniq.to.site[[i]], rep(NA, 4), fill=NA)}
  ids.family.uniq.to.site[[i]] <- cbind.fill(ids.family.uniq.to.site[[i]], setdiff(taxa_names(temp.jej), taxa_names(temp.not.jej)), fill=NA)
  ids.family.uniq.to.site[[i]] <- cbind.fill(ids.family.uniq.to.site[[i]], setdiff(taxa_names(temp.ile), taxa_names(temp.not.ile)), fill=NA)
  ids.family.uniq.to.site[[i]] <- cbind.fill(ids.family.uniq.to.site[[i]], setdiff(taxa_names(temp.cae), taxa_names(temp.not.cae)), fill=NA)
  ids.family.uniq.to.site[[i]] <- cbind.fill(ids.family.uniq.to.site[[i]], setdiff(taxa_names(temp.asc), taxa_names(temp.not.asc)), fill=NA)
  ids.family.uniq.to.site[[i]] <- cbind.fill(ids.family.uniq.to.site[[i]], setdiff(taxa_names(temp.des), taxa_names(temp.not.des)), fill=NA)
  if(i !=2 & i !=4){ids.family.uniq.to.site[[i]] <- cbind.fill(ids.family.uniq.to.site[[i]], setdiff(taxa_names(temp.fec), taxa_names(temp.not.fec)), fill=NA)}
  if(i==2 | i==4){ids.family.uniq.to.site[[i]] <- cbind.fill(ids.family.uniq.to.site[[i]], rep(NA))}
  ids.family.uniq.to.site[[i]]$init <- NULL
  colnames(ids.family.uniq.to.site[[i]]) <- c("duodenum", "jejunum", "ileum", "caecum", "asc.colon", "des.colon", "feces")
  
  rm(temp, temp.duo, temp.jej, temp.ile, temp.cae, temp.asc, temp.des, temp.fec,
     temp.duo.jej, temp.duo.jej.ile, temp.duo.jej.ile.cae, temp.duo.jej.ile.cae.asc, temp.duo.jej.ile.cae.asc.des,
     temp.not.duo, temp.not.jej, temp.not.ile, temp.not.cae, temp.not.asc, temp.not.des, temp.not.fec)
}

rm(pseq.glom)

# (c) Calculate summary statistics #####
# Convert absolute ASV numbers to percentages of total ASVs at each site.
df.new.to.site.prop <- 100*df.new.to.site/df.total.taxa
df.uniq.to.site.prop <- 100*df.uniq.to.site/df.total.taxa

# Prepare data for mixed models.
df.new.to.site.prop.model <- df.new.to.site.prop
df.new.to.site.prop.model$Individual <- levels(bg.metadata$Individual)
df.new.to.site.prop.model$duodenum <- NULL
df.new.to.site.prop.model <- melt(df.new.to.site.prop.model)

df.uniq.to.site.prop.model <- df.uniq.to.site.prop
df.uniq.to.site.prop.model$Individual <- levels(bg.metadata$Individual)
df.uniq.to.site.prop.model <- melt(df.uniq.to.site.prop.model)

# Melt data into long format.
df.new.to.site.prop <- melt(df.new.to.site.prop)
df.new.to.site.prop <- df.new.to.site.prop[complete.cases(df.new.to.site.prop),]
df.uniq.to.site.prop <- melt(df.uniq.to.site.prop)
df.uniq.to.site.prop <- df.uniq.to.site.prop[complete.cases(df.uniq.to.site.prop),]

# Calculate summary statistics for `new to site` ASVs and `unique to site` ASVs.
# Use both a standard ANOVA and a mixed model approach (to control for repeated measures in individuals).
leveneTest(value ~ variable, df.new.to.site.prop) # p=0.1357
summary(aov(value ~ variable, df.new.to.site.prop)) # F=13.44, df=5, p<0.001
model <- lmer(value ~ variable + (1|Individual),
              data=df.new.to.site.prop.model)
anova(model) # F=13.44, df=5, p<0.001
summary(glht(model, linfct = mcp(variable = "Tukey")), test = adjusted("BH"))

leveneTest(value ~ variable, df.uniq.to.site.prop) # p=0.4671
summary(aov(value ~ variable, df.uniq.to.site.prop)) # F=3.351, df=6, p=0.00649
model <- lmer(value ~ variable + (1|Individual),
              data=df.uniq.to.site.prop.model)
anova(model) # F=4.5638, df=6, p<0.001 
summary(glht(model, linfct = mcp(variable = "Tukey")), test = adjusted("BH"))

# (d) Determine patterns for which taxa appear most frequently at each intestinal site. ####
# First, row-bind all data frames from each individual. This occurs at the family level first.
ids.family.new.to.site.cat <- rbind(ids.family.new.to.site[[1]], ids.family.new.to.site[[2]], ids.family.new.to.site[[3]], ids.family.new.to.site[[4]], ids.family.new.to.site[[5]], 
                                    ids.family.new.to.site[[6]], ids.family.new.to.site[[7]], ids.family.new.to.site[[8]], ids.family.new.to.site[[9]], ids.family.new.to.site[[10]])

# For each intestinal site, count the number of times each bacterial family appears there.
# Store this information as a named number string (a1, a2, a3, a4, a5, a6, a7)
# This section of code will throw error messages ("Factor contains implicit NA" and "Setting row names on a tibble is deprecated") that do not affect the results.
duo <- ids.family.new.to.site.cat %>% dplyr::count(duodenum)
duo$duodenum <- forcats::fct_explicit_na(duo$duodenum, na_level = "Blank")
rownames(duo) <- duo$duodenum
duo$duodenum <- NULL

jej <- ids.family.new.to.site.cat %>% dplyr::count(jejunum)
jej$jejunum <- forcats::fct_explicit_na(jej$jejunum, na_level = "Blank")
rownames(jej) <- jej$jejunum
jej$jejunum <- NULL

ile <- ids.family.new.to.site.cat %>% dplyr::count(ileum)
ile$ileum <- forcats::fct_explicit_na(ile$ileum, na_level = "Blank")
rownames(ile) <- ile$ileum
ile$ileum <- NULL

cae <- ids.family.new.to.site.cat %>% dplyr::count(caecum)
cae$caecum <- forcats::fct_explicit_na(cae$caecum, na_level = "Blank")
rownames(cae) <- cae$caecum
cae$caecum <- NULL

asc <- ids.family.new.to.site.cat %>% dplyr::count(asc.colon)
asc$asc.colon <- forcats::fct_explicit_na(asc$asc.colon, na_level = "Blank")
rownames(asc) <- asc$asc.colon
asc$asc.colon <- NULL

des <- ids.family.new.to.site.cat %>% dplyr::count(des.colon)
des$des.colon <- forcats::fct_explicit_na(des$des.colon, na_level = "Blank")
rownames(des) <- des$des.colon
des$des.colon <- NULL

fec <- ids.family.new.to.site.cat %>% dplyr::count(feces)
fec$feces <- forcats::fct_explicit_na(fec$feces, na_level = "Blank")
rownames(fec) <- fec$feces
fec$feces <- NULL

rm(ids.family.new.to.site.cat, ids.family.new.to.site)

# Merge all of these together into one data frame.
ids.family.appear <- transform(merge(duo, jej, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
ids.family.appear <- transform(merge(ids.family.appear, ile, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
ids.family.appear <- transform(merge(ids.family.appear, cae, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
ids.family.appear <- transform(merge(ids.family.appear, asc, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
ids.family.appear <- transform(merge(ids.family.appear, des, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
ids.family.appear <- transform(merge(ids.family.appear, fec, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
colnames(ids.family.appear) <- c("duodenum", "jejunum", "ileum", "caecum", "asc.colon", "des.colon", "feces")
rm(duo, jej, ile, cae, asc, des, fec)

# Calculate the prevalence of each of these families (number of individuals that contain this family).
ids.family.appear$prev <- rowSums(ids.family.appear, na.rm=TRUE)
ids.family.appear <- subset(ids.family.appear, rownames(ids.family.appear) !="Blank")

# Repeat all of the above code at the genus level.
ids.genus.new.to.site.cat <- rbind(ids.genus.new.to.site[[1]], ids.genus.new.to.site[[2]], ids.genus.new.to.site[[3]], ids.genus.new.to.site[[4]], ids.genus.new.to.site[[5]], 
                                    ids.genus.new.to.site[[6]], ids.genus.new.to.site[[7]], ids.genus.new.to.site[[8]], ids.genus.new.to.site[[9]], ids.genus.new.to.site[[10]])

duo <- ids.genus.new.to.site.cat %>% dplyr::count(duodenum)
duo$duodenum <- forcats::fct_explicit_na(duo$duodenum, na_level = "Blank")
rownames(duo) <- duo$duodenum
duo$duodenum <- NULL

jej <- ids.genus.new.to.site.cat %>% dplyr::count(jejunum)
jej$jejunum <- forcats::fct_explicit_na(jej$jejunum, na_level = "Blank")
rownames(jej) <- jej$jejunum
jej$jejunum <- NULL

ile <- ids.genus.new.to.site.cat %>% dplyr::count(ileum)
ile$ileum <- forcats::fct_explicit_na(ile$ileum, na_level = "Blank")
rownames(ile) <- ile$ileum
ile$ileum <- NULL

cae <- ids.genus.new.to.site.cat %>% dplyr::count(caecum)
cae$caecum <- forcats::fct_explicit_na(cae$caecum, na_level = "Blank")
rownames(cae) <- cae$caecum
cae$caecum <- NULL

asc <- ids.genus.new.to.site.cat %>% dplyr::count(asc.colon)
asc$asc.colon <- forcats::fct_explicit_na(asc$asc.colon, na_level = "Blank")
rownames(asc) <- asc$asc.colon
asc$asc.colon <- NULL

des <- ids.genus.new.to.site.cat %>% dplyr::count(des.colon)
des$des.colon <- forcats::fct_explicit_na(des$des.colon, na_level = "Blank")
rownames(des) <- des$des.colon
des$des.colon <- NULL

fec <- ids.genus.new.to.site.cat %>% dplyr::count(feces)
fec$feces <- forcats::fct_explicit_na(fec$feces, na_level = "Blank")
rownames(fec) <- fec$feces
fec$feces <- NULL

rm(ids.genus.new.to.site.cat, ids.genus.new.to.site)

ids.genus.appear <- transform(merge(duo, jej, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
ids.genus.appear <- transform(merge(ids.genus.appear, ile, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
ids.genus.appear <- transform(merge(ids.genus.appear, cae, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
ids.genus.appear <- transform(merge(ids.genus.appear, asc, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
ids.genus.appear <- transform(merge(ids.genus.appear, des, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
ids.genus.appear <- transform(merge(ids.genus.appear, fec, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
colnames(ids.genus.appear) <- c("duodenum", "jejunum", "ileum", "caecum", "asc.colon", "des.colon", "feces")
rm(duo, jej, ile, cae, asc, des, fec)

ids.genus.appear$prev <- rowSums(ids.genus.appear, na.rm=TRUE)
ids.genus.appear <- subset(ids.genus.appear, rownames(ids.genus.appear) !="Blank")

# (e) Determine patterns for which taxa are unique most frequently at each intestinal site. ####
ids.family.uniq.to.site.cat <- rbind(ids.family.uniq.to.site[[1]], ids.family.uniq.to.site[[2]], ids.family.uniq.to.site[[3]], ids.family.uniq.to.site[[4]], ids.family.uniq.to.site[[5]], 
                                    ids.family.uniq.to.site[[6]], ids.family.uniq.to.site[[7]], ids.family.uniq.to.site[[8]], ids.family.uniq.to.site[[9]], ids.family.uniq.to.site[[10]])

duo <- ids.family.uniq.to.site.cat %>% dplyr::count(duodenum)
duo$duodenum <- forcats::fct_explicit_na(duo$duodenum, na_level = "Blank")
rownames(duo) <- duo$duodenum
duo$duodenum <- NULL

jej <- ids.family.uniq.to.site.cat %>% dplyr::count(jejunum)
jej$jejunum <- forcats::fct_explicit_na(jej$jejunum, na_level = "Blank")
rownames(jej) <- jej$jejunum
jej$jejunum <- NULL

ile <- ids.family.uniq.to.site.cat %>% dplyr::count(ileum)
ile$ileum <- forcats::fct_explicit_na(ile$ileum, na_level = "Blank")
rownames(ile) <- ile$ileum
ile$ileum <- NULL

cae <- ids.family.uniq.to.site.cat %>% dplyr::count(caecum)
cae$caecum <- forcats::fct_explicit_na(cae$caecum, na_level = "Blank")
rownames(cae) <- cae$caecum
cae$caecum <- NULL

asc <- ids.family.uniq.to.site.cat %>% dplyr::count(asc.colon)
asc$asc.colon <- forcats::fct_explicit_na(asc$asc.colon, na_level = "Blank")
rownames(asc) <- asc$asc.colon
asc$asc.colon <- NULL

des <- ids.family.uniq.to.site.cat %>% dplyr::count(des.colon)
des$des.colon <- forcats::fct_explicit_na(des$des.colon, na_level = "Blank")
rownames(des) <- des$des.colon
des$des.colon <- NULL

fec <- ids.family.uniq.to.site.cat %>% dplyr::count(feces)
fec$feces <- forcats::fct_explicit_na(fec$feces, na_level = "Blank")
rownames(fec) <- fec$feces
fec$feces <- NULL

rm(ids.family.uniq.to.site.cat, ids.family.uniq.to.site)

ids.family.unique <- transform(merge(duo, jej, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
ids.family.unique <- transform(merge(ids.family.unique, ile, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
ids.family.unique <- transform(merge(ids.family.unique, cae, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
ids.family.unique <- transform(merge(ids.family.unique, asc, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
ids.family.unique <- transform(merge(ids.family.unique, des, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
ids.family.unique <- transform(merge(ids.family.unique, fec, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
colnames(ids.family.unique) <- c("duodenum", "jejunum", "ileum", "caecum", "asc.colon", "des.colon", "feces")
rm(duo, jej, ile, cae, asc, des, fec)

# Calculate the prevalence of each of these families (number of individuals that contain this family).
ids.family.unique$prev <- rowSums(ids.family.unique, na.rm=TRUE)
ids.family.unique <- subset(ids.family.unique, rownames(ids.family.unique) !="Blank")

ids.genus.uniq.to.site.cat <- rbind(ids.genus.uniq.to.site[[1]], ids.genus.uniq.to.site[[2]], ids.genus.uniq.to.site[[3]], ids.genus.uniq.to.site[[4]], ids.genus.uniq.to.site[[5]], 
                                   ids.genus.uniq.to.site[[6]], ids.genus.uniq.to.site[[7]], ids.genus.uniq.to.site[[8]], ids.genus.uniq.to.site[[9]], ids.genus.uniq.to.site[[10]])

duo <- ids.genus.uniq.to.site.cat %>% dplyr::count(duodenum)
duo$duodenum <- forcats::fct_explicit_na(duo$duodenum, na_level = "Blank")
rownames(duo) <- duo$duodenum
duo$duodenum <- NULL

jej <- ids.genus.uniq.to.site.cat %>% dplyr::count(jejunum)
jej$jejunum <- forcats::fct_explicit_na(jej$jejunum, na_level = "Blank")
rownames(jej) <- jej$jejunum
jej$jejunum <- NULL

ile <- ids.genus.uniq.to.site.cat %>% dplyr::count(ileum)
ile$ileum <- forcats::fct_explicit_na(ile$ileum, na_level = "Blank")
rownames(ile) <- ile$ileum
ile$ileum <- NULL

cae <- ids.genus.uniq.to.site.cat %>% dplyr::count(caecum)
cae$caecum <- forcats::fct_explicit_na(cae$caecum, na_level = "Blank")
rownames(cae) <- cae$caecum
cae$caecum <- NULL

asc <- ids.genus.uniq.to.site.cat %>% dplyr::count(asc.colon)
asc$asc.colon <- forcats::fct_explicit_na(asc$asc.colon, na_level = "Blank")
rownames(asc) <- asc$asc.colon
asc$asc.colon <- NULL

des <- ids.genus.uniq.to.site.cat %>% dplyr::count(des.colon)
des$des.colon <- forcats::fct_explicit_na(des$des.colon, na_level = "Blank")
rownames(des) <- des$des.colon
des$des.colon <- NULL

fec <- ids.genus.uniq.to.site.cat %>% dplyr::count(feces)
fec$feces <- forcats::fct_explicit_na(fec$feces, na_level = "Blank")
rownames(fec) <- fec$feces
fec$feces <- NULL

rm(ids.genus.uniq.to.site.cat, ids.genus.uniq.to.site)

ids.genus.unique <- transform(merge(duo, jej, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
ids.genus.unique <- transform(merge(ids.genus.unique, ile, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
ids.genus.unique <- transform(merge(ids.genus.unique, cae, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
ids.genus.unique <- transform(merge(ids.genus.unique, asc, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
ids.genus.unique <- transform(merge(ids.genus.unique, des, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
ids.genus.unique <- transform(merge(ids.genus.unique, fec, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
colnames(ids.genus.unique) <- c("duodenum", "jejunum", "ileum", "caecum", "asc.colon", "des.colon", "feces")
rm(duo, jej, ile, cae, asc, des, fec)

ids.genus.unique$prev <- rowSums(ids.genus.unique, na.rm=TRUE)
ids.genus.unique <- subset(ids.genus.unique, rownames(ids.genus.unique) !="Blank")



#### dbRDA MODELS AND PERMANOVA ####
# Create a phyloseq object without fecal samples, both for rarefied data and for CLR-transformed data.
pseq.rarefied.nf <- subset_samples(pseq.rarefied, Type !="feces")
pseq.rarefied.nf <- prune_taxa(taxa_sums(pseq.rarefied.nf) > 0, pseq.rarefied.nf)
bg.metadata.nf <- subset(bg.metadata, Type !="feces")

pseq.culled.nf <- subset_samples(pseq.culled, Type !="feces")
pseq.culled.nf <- prune_taxa(taxa_sums(pseq.culled.nf) > 0, pseq.culled.nf)
pseq.clr <- microbiome::transform(pseq.culled, "clr")
pseq.clr.nf <- microbiome::transform(pseq.culled.nf, "clr")

# Create distance matrices (Bray-Curtis, weighted/unweighted UniFrac, Jaccard, and Aitchison), for feature tables including feces.
dist.BC <- vegdist(as.data.frame(otu_table(pseq.rarefied)), method="bray")
dist.wUF <- UniFrac(pseq.rarefied, weighted=TRUE, normalized=TRUE)
dist.uwUF <- UniFrac(pseq.rarefied, weighted=FALSE, normalized=TRUE)
dist.AIT <- vegdist(as.data.frame(otu_table(pseq.clr)), method="euclidean")

otu <- as.data.frame(otu_table(pseq.rarefied))
otu[otu != 0] <- 1
dist.JAC <- vegdist(otu, method="jaccard")
rm(otu)

# Create distance matrices (Bray-Curtis, weighted/unweighted UniFrac, Jaccard, and Aitchison), for feature tables excluding feces.
dist.BC.nofeces <- vegdist(as.data.frame(otu_table(pseq.rarefied.nf)), method="bray")
dist.wUF.nofeces <- UniFrac(pseq.rarefied.nf, weighted=TRUE, normalized=TRUE)
dist.uwUF.nofeces <- UniFrac(pseq.rarefied.nf, weighted=FALSE, normalized=TRUE)
dist.AIT.nofeces <- vegdist(as.data.frame(otu_table(pseq.clr.nf)), method="euclidean")

otu <- as.data.frame(otu_table(pseq.rarefied.nf))
otu[otu != 0] <- 1
dist.JAC.nofeces <- vegdist(otu, method="jaccard")
rm(otu)

# Create a dbRDA model based on these variables and evaluate predictive power.
# Because the Aitchison distance appears in the manuscript, that is the distance matrix used here.
# All of the other distance matrices were tested and could be substituted in the code below.
dbrda.nofeces <- capscale(dist.AIT.nofeces ~ Segment/Type+Individual, data = subset(bg.metadata, Type !="feces"))
RsquareAdj(dbrda.nofeces)

# Test the constraints for significance.
anova(dbrda.nofeces, by="term")

# Evaluate the axis loadings based on ANOVA against the dbRDA axes.
# The percentage explained by each variable as calculated as the sum of squares for that variable,
# divided by  the total sum of sqquares for the model.
dbrda.nofeces_fort <- fortify(dbrda.nofeces, display="sites")
dbrda.nofeces_fort <- cbind(bg.metadata.nf, dbrda.nofeces_fort)

anova(lm(CAP1 ~ Individual+Segment/Type, dbrda.nofeces_fort))
anova(lm(CAP1 ~ Segment/Type+Individual, dbrda.nofeces_fort))
anova(lm(CAP2 ~ Individual+Segment/Type, dbrda.nofeces_fort))
anova(lm(CAP2 ~ Segment/Type+Individual, dbrda.nofeces_fort))
anova(lm(CAP3 ~ Individual+Segment/Type, dbrda.nofeces_fort))
anova(lm(CAP3 ~ Segment/Type+Individual, dbrda.nofeces_fort))

# Define pairwise PERMANOVA function.
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni'){
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
}

# Test clustering effects using PERMANOVA.
adonis(dist.AIT ~ Segment/Type+Individual, bg.metadata, permutations=1000)

temp.otu.table <- as.data.frame(otu_table(pseq.clr.nf)) # Without feces
pairwise.adonis(temp.otu.table, bg.metadata.nf$Type, p.adjust.m = "BH", sim.method="euclidean")

# Compare the effects of including or excluding feces in the analysis.
# Run a dbRDA model including feces and test variance explained by the model and each axis.
dbrda.all <- capscale(dist.AIT ~ Individual+Segment/Type, bg.metadata)
RsquareAdj(dbrda.all)
anova(dbrda.all, by="term")
dbrda.all_fort <- fortify(dbrda.all, display="sites")
dbrda.all_fort <- cbind(bg.metadata, dbrda.all_fort)

anova(lm(CAP1 ~ Individual+Segment/Type, dbrda.all_fort))
anova(lm(CAP1 ~ Segment/Type+Individual, dbrda.all_fort))
anova(lm(CAP2 ~ Individual+Segment/Type, dbrda.all_fort))
anova(lm(CAP2 ~ Segment/Type+Individual, dbrda.all_fort))
anova(lm(CAP3 ~ Individual+Segment/Type, dbrda.all_fort))
anova(lm(CAP3 ~ Segment/Type+Individual, dbrda.all_fort))

# Compare the effects of classifying the ileum as the large intestine or the caecum as the small intestine.
# Note that for this analysis, intestinal site cannot be included in the CCA models because it masks the puposefully misclassified segments.
# Classify the ileum as the large intestine.
bg.metadata.nf.il <- bg.metadata.nf
bg.metadata.nf.il$Segment[bg.metadata.nf.il$Type=="ileum"] <- "large.intestine"

# Classify the caecum as the small intestine.
bg.metadata.nf.ca <- bg.metadata.nf
bg.metadata.nf.ca$Segment[bg.metadata.nf.ca$Type=="caecum"] <- "small.intestine"

# Evaluate the model for misclassified ileum.
RsquareAdj(capscale(dist.AIT.nofeces ~ Individual+Segment, data = bg.metadata.nf.il)) # Rsq = 50.90%, Adj-Rsq = 40.67%
RsquareAdj(capscale(dist.AIT.nofeces ~ Individual+Segment, data = bg.metadata.nf.ca)) # Rsq = 45.07%, Adj-Rsq = 33.63%
RsquareAdj(capscale(dist.AIT.nofeces ~ Individual+Segment, data = bg.metadata.nf)) # Rsq = 51.77%, Adj-Rsq = 41.7%

rm(dbrda.ca, dbrda.il, bg.metadata.nf.ca, bg.metadata.nf.il)

# Repeat the dbRDA analysis within each intestinal segment. Start by subsetting phyloseq objects and calculating new distance matrices.
pseq.rarefied.small <- subset_samples(pseq.rarefied, Segment=="small.intestine")
pseq.rarefied.small <- prune_taxa(taxa_sums(pseq.rarefied.small) > 0, pseq.rarefied.small)
pseq.rarefied.large <- subset_samples(pseq.rarefied, Segment=="large.intestine")
pseq.rarefied.large <- prune_taxa(taxa_sums(pseq.rarefied.large) > 0, pseq.rarefied.large)
pseq.culled.small <- subset_samples(pseq.culled, Segment=="small.intestine")
pseq.culled.small <- prune_taxa(taxa_sums(pseq.culled.small) > 0, pseq.culled.small)
pseq.culled.large <- subset_samples(pseq.culled, Segment=="large.intestine")
pseq.culled.large <- prune_taxa(taxa_sums(pseq.culled.large) > 0, pseq.culled.large)
pseq.culled.large.nf <- subset_samples(pseq.culled.nf, Segment=="large.intestine")
pseq.culled.large.nf <- prune_taxa(taxa_sums(pseq.culled.large.nf) > 0, pseq.culled.large.nf)

bg.metadata.small <- subset(bg.metadata, Segment=="small.intestine")
bg.metadata.large <- subset(bg.metadata, Segment=="large.intestine")
bg.metadata.large.nf <- subset(bg.metadata.nf, Segment=="large.intestine" & Type !="feces")

dist.BC.small <- vegdist(as.data.frame(otu_table(pseq.rarefied.small)), distance="bray")
dist.BC.large <- vegdist(as.data.frame(otu_table(pseq.rarefied.large)), distance="bray")

dist.wUF.small <- UniFrac(pseq.rarefied.small, weighted=TRUE, normalized=TRUE)
dist.wUF.large <- UniFrac(pseq.rarefied.large, weighted=TRUE, normalized=TRUE)

dist.AIT.small <- vegdist(as.data.frame(otu_table(microbiome::transform(pseq.culled.small, "clr"))), method="euclidean")
dist.AIT.large <- vegdist(as.data.frame(otu_table(microbiome::transform(pseq.culled.large, "clr"))), method="euclidean")
dist.AIT.large.nf <- vegdist(as.data.frame(otu_table(microbiome::transform(pseq.culled.large.nf, "clr"))), method="euclidean")

dist.uwUF.small <- UniFrac(pseq.rarefied.small, weighted=FALSE, normalized=TRUE)
dist.uwUF.large <- UniFrac(pseq.rarefied.large, weighted=FALSE, normalized=TRUE)

otu <- as.data.frame(otu_table(pseq.rarefied.small))
otu[otu != 0] <- 1
dist.JAC.small <- vegdist(otu, method="jaccard")

otu <- as.data.frame(otu_table(pseq.rarefied.large))
otu[otu != 0] <- 1
dist.JAC.large <- vegdist(otu, method="jaccard")
rm(otu)

# Perform the model selection procedure and analysis for the small intestine.
dbrda.small <- capscale(dist.AIT.small ~ Individual+Type, data = bg.metadata.small)
anova(dbrda.small, by="term")
RsquareAdj(dbrda.small)
adonis(dist.AIT.small ~ Type+Individual, bg.metadata.small, permutations=1000)

dbrda.small_fort <- fortify(dbrda.small, display="sites")
dbrda.small_fort <- cbind(bg.metadata.small, dbrda.small_fort)
anova(lm(CAP1 ~ Individual+Type, dbrda.small_fort))
anova(lm(CAP1 ~ Type+Individual, dbrda.small_fort))
anova(lm(CAP2 ~ Individual+Type, dbrda.small_fort))
anova(lm(CAP2 ~ Type+Individual, dbrda.small_fort))

# Perform the model selection procedure and analysis for the large intestine.
dbrda.large <- capscale(dist.AIT.large ~ Type+Individual, data = bg.metadata.large) # With feces
anova(dbrda.large, by="term")
RsquareAdj(dbrda.large) 
adonis(dist.AIT.large ~ Type+Individual, bg.metadata.large, permutations=1000)

dbrda.large_fort <- fortify(dbrda.large, display="sites")
dbrda.large_fort <- cbind(bg.metadata.large, dbrda.large_fort)
anova(lm(CAP1 ~ Individual+Type, dbrda.large_fort))
anova(lm(CAP1 ~ Type+Individual, dbrda.large_fort))
anova(lm(CAP2 ~ Individual+Type, dbrda.large_fort))
anova(lm(CAP2 ~ Type+Individual, dbrda.large_fort))

dbrda.large <- capscale(dist.AIT.large.nf ~ Type+Individual, data = bg.metadata.large.nf) # Without feces
anova(dbrda.large, by="term")
RsquareAdj(dbrda.large) 
adonis(dist.AIT.large ~ Type+Individual, bg.metadata.large, permutations=1000)

dbrda.large_fort <- fortify(dbrda.large, display="sites")
dbrda.large_fort <- cbind(bg.metadata.large, dbrda.large_fort)
anova(lm(CAP1 ~ Individual+Type, dbrda.large_fort))
anova(lm(CAP1 ~ Type+Individual, dbrda.large_fort))
anova(lm(CAP2 ~ Individual+Type, dbrda.large_fort))
anova(lm(CAP2 ~ Type+Individual, dbrda.large_fort))

rm(pseq.rarefied.large, pseq.rarefied.small, pseq.rarefied.nf,
   pseq.culled.large, pseq.culled.small, pseq.culled.nf)

#### MULTIVARIATE DISPERSION #####
bg.metadata.nf <- subset(bg.metadata, Type !="feces")

# Calculate the dispersion around centroids for each categorical variable.
ind <- betadisper(dist.AIT.nofeces, bg.metadata.nf$Individual)
seg <- betadisper(dist.AIT.nofeces, bg.metadata.nf$Segment)
type <- betadisper(dist.AIT.nofeces, bg.metadata.nf$Type)

# Test for significant differences.
permutest(ind, pairwise=TRUE, permutations=1000) # F=0.610, df=9, p=0.7872
permutest(seg, pairwise=TRUE, permutations=1000) # F=17.447, df=1, p=0.001
permutest(type, pairwise=TRUE, permutations=1000) # F=8.2796, df=5, p=0.001

# Add distances to centroids as new data in the metadata object.
bg.metadata.nf$ind_dist <- ind[["distances"]]
bg.metadata.nf$seg_dist <- seg[["distances"]]
bg.metadata.nf$type_dist <- type[["distances"]]

# Calculate mean dispersions for each group.
mean(bg.metadata.nf$ind_dist) # 30.22029
mean(bg.metadata.nf$seg_dist) # 35.87815
mean(bg.metadata.nf$type_dist) # 34.45207

# Test for significant differences in dispersion around different groups - significantly more dispersion around 'segment' centroids.
df1 <- data.frame(Group = rep("Seg", 59),
                  Distance = bg.metadata.nf$seg_dist)
df2 <- data.frame(Group = rep("Indv", 59),
                  Distance = bg.metadata.nf$ind_dist)
df3 <- data.frame(Group = rep("Type", 59),
                  Distance = bg.metadata.nf$type_dist)
df <- rbind(df1, df2, df3)

leveneTest(Distance ~ Group, df)
oneway.test(Distance ~ Group, df, var.equal = TRUE) # F=5.1641, df=2, p=0.006627
TukeyHSD(aov(df$Distance ~ df$Group), p.adjust="BH")

rm(ind, seg, type, df1, df2, df3, df)

# Repeat this analysis within the small and large intestine separately.
# Calculate the dispersion around centroids for each categorical variable.
ind_small <- betadisper(dist.AIT.small, bg.metadata.small$Individual)
type_small <- betadisper(dist.AIT.small, bg.metadata.small$Type)
ind_large <- betadisper(dist.AIT.large.nf, bg.metadata.large.nf$Individual)
type_large <- betadisper(dist.AIT.large.nf, bg.metadata.large.nf$Type)

# Test for significant differences.
permutest(ind_small, pairwise=TRUE, permutations=1000) # F=0.5911, df=9, p=0.7792
permutest(type_small, pairwise=TRUE, permutations=1000) # F=7.7856, df=2, p=0.004995
permutest(ind_large, pairwise=TRUE, permutations=1000) # F=0.8652, df=9, p=0.5684
permutest(type_large, pairwise=TRUE, permutations=1000) # F=0.361, df=2, p=0.7153

# Append distances to centroid.
bg.metadata.small$ind_dist <- ind_small[["distances"]]
bg.metadata.small$type_dist <- type_small[["distances"]]
bg.metadata.large.nf$ind_dist <- ind_large[["distances"]]
bg.metadata.large.nf$type_dist <- type_large[["distances"]]

# Calculate mean dispersions around each group.
mean(bg.metadata.small$ind_dist) # 21.4749
mean(bg.metadata.small$type_dist) # 27.92973
mean(bg.metadata.large.nf$ind_dist) # 20.23386
mean(bg.metadata.large.nf$type_dist) # 39.70994

# Test for significant differences in dispersions around different groups, within the small and large intestine.
df1 <- data.frame(Group = rep("Indv", 29),
                  Distance = bg.metadata.small$ind_dist)
df2 <- data.frame(Group = rep("Site", 29),
                  Distance = bg.metadata.small$type_dist)
df <- rbind(df1, df2)

leveneTest(Distance ~ Group, df)
t.test(Distance ~ Group, df, var.equal=TRUE) # t = -2.4886, df = 56, p-value = 0.01583

df1 <- data.frame(Group = rep("Indv", 30),
                  Distance = bg.metadata.large.nf$ind_dist)
df2 <- data.frame(Group = rep("Site", 30),
                  Distance = bg.metadata.large.nf$type_dist)
df <- rbind(df1, df2)

leveneTest(Distance ~ Group, df)
t.test(Distance ~ Group, df, var.equal=TRUE) # t = -9.3946, df = 58, p-value = < 0.001

rm(df, df1, df2, bg.metadata.large, bg.metadata.small, bg.metadata.nf, 
   ind_large, ind_small, type_large, type_small,
   dist.BC.large, dist.BC.small, dist.BC, dist.BC.nofeces,
   dist.wUF.large, dist.wUF.small, dist.wUF, dist.wUF.nofeces,
   dist.uwUF.large, dist.uwUF.small, dist.uwUF, dist.uwUF.nofeces,
   dist.AIT.large, dist.AIT.small, dist.AIT, dist.AIT.nofeces, dist.AIT.large.nf,
   dist.JAC.large, dist.JAC.small, dist.JAC, dist.JAC.nofeces)

#### RANDOM FOREST MODELS ####
# Create a place to store model outputs.
rf.results <- list()
rf.results[["Models"]] <- list()
rf.results[["Gini_scores"]] <- list()

# This code is adapted from https://rpubs.com/michberr/randomforestmicrobe
# Repeat the model 3 times (once each for segment, individual identity, and intestinal site)
for(i in c("Individual","Segment","Type")){
  # Filter taxa with a relative abundance below 0.01%
  prunescale = 0.0001
  minlib <- min(sample_sums(pseq.rarefied))
  tax.mean <- taxa_sums(pseq.rarefied)/nsamples(pseq.rarefied)
  sites.prune <- prune_taxa(tax.mean > prunescale*minlib, pseq.rarefied)
  
  sites.prune <- subset_taxa(pseq.culled.sub, taxa_names(pseq.culled.sub) %in% taxa_names(sites.prune))
  sites.prune <- microbiome::transform(sites.prune, "clr")
  
  # Define the model predictors as ASV abundances.
  predictors <- otu_table(sites.prune)
  
  # Define the model response as the variable of interest (segment, individual identity, or intestinal site).
  response <- bg.metadata[,i]

  # Combine the predictors and response into a single data frame.
  rf.data <- data.frame(response, predictors)
  
  # Perform random forest model.
  set.seed(2)
  rf.raw.results <- randomForest(response~., data = rf.data, ntree = 1000)

  # Extract the 40 most important ASVs, based on their mean decrease in the Gini coefficient.
  rf.importance <- randomForest::importance(rf.raw.results)
  rf.importance <- data.frame(predictors = rownames(rf.importance), rf.importance)
  rf.importance <- arrange(rf.importance, desc(MeanDecreaseGini))
  rf.importance$predictors <- factor(rf.importance$predictors,
                                      levels = rf.importance$predictors)
  rf.importance <- rf.importance[1:40, ]
  rownames(rf.importance) <- rf.importance$predictors
  rf.importance$predictors <- NULL
  
  # Add taxonomy information to the predictors.
  temp <- tax_table(sites.prune)
  temp <- subset(temp, rownames(temp) %in% rownames(rf.importance))
  rf.importance <- merge(rf.importance, temp, by="row.names", all=TRUE)

  # Organize by Gini score.  
  rf.importance <- arrange(rf.importance, desc(MeanDecreaseGini))

  # Add information about taxon abundances for the given predictor.
  taxa_names(pseq.nsrare.glom) <- taxa_names(pseq.glom)
  temp.pseq <- transform_sample_counts(pseq.nsrare.glom, function(x) 100*x/sum(x))
  temp.pseq <- subset_taxa(temp.pseq, taxa_names(temp.pseq) %in% rf.importance$Row.names)
  temp.pseq.otu <- as.data.frame(otu_table(temp.pseq))
  
  # Calculate taxon means based on each predictor, and the overall mean in the sample.
  temp.pseq.otu[,i] <- bg.metadata[,i]
  temp.pseq.tab <- temp.pseq.otu %>% group_by(temp.pseq.otu[,i]) %>% summarise_all(list(mean=mean))
  temp.pseq.tab2 <- temp.pseq.otu %>% summarise_all(list(mean=mean))
  temp.pseq.tab <- as.data.frame(temp.pseq.tab)
  temp.pseq.tab2 <- as.data.frame(temp.pseq.tab2)
  rownames(temp.pseq.tab) <- temp.pseq.tab[,1]
  rownames(temp.pseq.tab2) <- "Overall"
  
  temp.pseq.tab[,c(1,ncol(temp.pseq.tab))] <- NULL
  temp.pseq.tab2[,c(ncol(temp.pseq.tab2))] <- NULL
  
  colnames(temp.pseq.tab) <- taxa_names(temp.pseq)
  colnames(temp.pseq.tab2) <- taxa_names(temp.pseq)
  
  temp.pseq.tab <- as.data.frame(t(rbind(temp.pseq.tab, temp.pseq.tab2)))

  rownames(rf.importance) <- rf.importance$Row.names
  rf.importance <- merge(rf.importance, temp.pseq.tab, by="row.names", all=TRUE)
  rownames(rf.importance) <- rf.importance$Row.names
  rf.importance$Row.names <- NULL
  rf.importance$Row.names <- NULL
  rm(temp.pseq, temp.pseq.otu, temp.pseq.tab, temp.pseq.tab2)
  
  # Store data and clean workspace.
  rf.results[["Models"]][[i]] <- rf.raw.results
  rf.results[["Gini_scores"]][[i]] <- rf.importance
  rm(prunescale, minlib, tax.mean, sites.prune, predictors, response, rf.data, temp, rf.raw.results, rf.importance)
  }

# Evaluate confusion matrices for each model.
print(rf.results[["Models"]][["Individual"]])
print(rf.results[["Models"]][["Segment"]])
print(rf.results[["Models"]][["Type"]])

#### ALPHA DIVERSITY: VARIANCE DUE TO INDIVIDUAL VS. SITE ####
library(MuMIn)
# Create regression models predicting ASV richness based on intestinal location, with and without a random effect.
fm1.0 <- glm.nb(Observed ~ 1,
            data=bg.metadata)

fm1.1 <- glm.nb(Observed~Type,
             data=bg.metadata)

fm1.2 <- glmer.nb(Observed~Type + (1|Individual),
              data=bg.metadata)

# Ensure that intestinal site has a significant effect on species richness
AIC(fm1.0) # 727.3853
AIC(fm1.1) # 700.4605
AIC(fm1.2) # 698.1156

anova.negbin(fm1.0, fm1.1) # Species richness - null model vs. intestinal site
# Model    theta  Resid. df    2 x log-lik.   Test    df LR stat.      Pr(Chi)
#     1 3.965372        66       -723.3853                                   
#  Type 7.070310        60       -684.4605 1 vs 2     6 38.92484 7.404534e-07

# Estimate the proportion of variance explained by the fixed (R2m) and fixed+random (R2c) effects
r.squaredGLMM(fm1.2)
#                 R2m       R2c
# delta     0.4566037 0.5602257
# lognormal 0.4682755 0.5745463
# trigamma  0.4441504 0.5449462

(0.5449462-0.4441504)*100 # 10.08%

# 44.4% of variance explained by intestinal location / 10.08% explained by individual.

# Repeat for Shannon diversity.
fm2.0 <- lm(Shannon~1, data=bg.metadata)
fm2.1 <- lm(Shannon~Type, data=bg.metadata)
fm2.2 <- lmer(Shannon~Type+(1|Individual), data=bg.metadata)

anova(fm2.0, fm2.1) # Shannon diversity - null model vs. intestinal site
#   Res.Df    RSS Df Sum of Sq      F    Pr(>F)    
# 1     66 40.931                                  
# 2     60 21.226  6    19.705 9.2836 3.428e-07 ***

AIC(fm2.0) # 161.1198
AIC(fm2.1) # 129.123
AIC(fm2.2) # 145.4818

# Estimate the proportion of variance explained by the fixed (R2m) and fixed+random (R2c) effects
r.squaredGLMM(fm2.2)
#          R2m       R2c
# [1,] 0.45636 0.4720273

(0.4720273-0.45636)*100 # 1.56673%

rm(fm1.0, fm1.1, fm1.2, fm2.0, fm2.1, fm2.2)

#### DISTANCE FROM FECAL COMMUNITIES TO OTHER INTESTINAL SITES ####
# Obtain a table of the Bray-Curtis and weighted UniFrac distances from feces to each intestinal site.
# Create a place to store the output tables.
dist.to.feces.BC <- list()
dist.to.feces.wUF <- list()
dist.to.feces.AIT <- list()

# For each individual with a fecal sample...
for(i in 1:length(ind.ids)){
  # Subset phyloseq object to that individual.
  temp <- subset_samples(pseq.rarefied, Individual==ind.ids[[i]])
  temp <- prune_taxa(taxa_sums(temp) > 0, temp)
  
  temp2 <- subset_samples(pseq.culled, Individual==ind.ids[[i]])
  temp2 <- prune_taxa(taxa_sums(temp2) > 0, temp2)
  temp2 <- microbiome::transform(temp2, "clr")
  
  # Calculate Bray-Curtis and weighted UniFrac distances for the phyloseq object.
  tempBC <- as.matrix(vegdist(as.data.frame(otu_table(temp)), method = "bray"))
  tempwUF <- as.matrix(UniFrac(temp), weighted=TRUE, normalized=TRUE)
  tempAIT <- as.matrix(vegdist(as.data.frame(otu_table(temp2)), method = "euclidean"))
  
  # Name distance matrices by intestinal sites.
  colnames(tempBC) <- as.data.frame(cbind(sample_data(temp)))$Type
  rownames(tempBC) <- as.data.frame(cbind(sample_data(temp)))$Type
  tempBC[tempBC==0] <- NA
  
  colnames(tempwUF) <- as.data.frame(cbind(sample_data(temp)))$Type
  rownames(tempwUF) <- as.data.frame(cbind(sample_data(temp)))$Type
  tempwUF[tempwUF==0] <- NA
  
  colnames(tempAIT) <- as.data.frame(cbind(sample_data(temp2)))$Type
  rownames(tempAIT) <- as.data.frame(cbind(sample_data(temp2)))$Type
  tempAIT[tempAIT==0] <- NA
  
  # Store distance matrices as data frames.
  dist.to.feces.BC[[i]] <- as.data.frame(tempBC)
  dist.to.feces.wUF[[i]] <- as.data.frame(tempwUF)
  dist.to.feces.AIT[[i]] <- as.data.frame(tempAIT)
  
  rm(temp, temp2, tempBC, tempwUF, tempAIT)
}

# Create a data frame containing ONLY distances from each site to the feces, for each individual.
# Because one individual lacks a duodenum sample, this code is somewhat complicated.
dist.to.feces.BC.summary <- cbind(dist.to.feces.BC[[1]]$feces, dist.to.feces.BC[[2]]$feces, dist.to.feces.BC[[3]]$feces, dist.to.feces.BC[[4]]$feces,
                           dist.to.feces.BC[[5]]$feces, dist.to.feces.BC[[6]]$feces, dist.to.feces.BC[[8]]$feces)
rownames(dist.to.feces.BC.summary) <- rownames(dist.to.feces.BC[[1]])
dist.to.feces.BC.summary <- transform(merge(dist.to.feces.BC.summary, data.frame(row.names=rownames(dist.to.feces.BC[[7]]), dist.to.feces.BC[[7]]$feces), by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
colnames(dist.to.feces.BC.summary) <- c("1", "2", "3", "4", "5", "6", "7", "8")
rm(dist.to.feces.BC)

dist.to.feces.wUF.summary <- cbind(dist.to.feces.wUF[[1]]$feces, dist.to.feces.wUF[[2]]$feces, dist.to.feces.wUF[[3]]$feces, dist.to.feces.wUF[[4]]$feces,
                                  dist.to.feces.wUF[[5]]$feces, dist.to.feces.wUF[[6]]$feces, dist.to.feces.wUF[[8]]$feces)
rownames(dist.to.feces.wUF.summary) <- rownames(dist.to.feces.wUF[[1]])
dist.to.feces.wUF.summary <- transform(merge(dist.to.feces.wUF.summary, data.frame(row.names=rownames(dist.to.feces.wUF[[7]]), dist.to.feces.wUF[[7]]$feces), by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
colnames(dist.to.feces.wUF.summary) <- c("1", "2", "3", "4", "5", "6", "7", "8")
rm(dist.to.feces.wUF)

dist.to.feces.AIT.summary <- cbind(dist.to.feces.AIT[[1]]$feces, dist.to.feces.AIT[[2]]$feces, dist.to.feces.AIT[[3]]$feces, dist.to.feces.AIT[[4]]$feces,
                                   dist.to.feces.AIT[[5]]$feces, dist.to.feces.AIT[[6]]$feces, dist.to.feces.AIT[[8]]$feces)
rownames(dist.to.feces.AIT.summary) <- rownames(dist.to.feces.AIT[[1]])
dist.to.feces.AIT.summary <- transform(merge(dist.to.feces.AIT.summary, data.frame(row.names=rownames(dist.to.feces.AIT[[7]]), dist.to.feces.AIT[[7]]$feces), by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
colnames(dist.to.feces.AIT.summary) <- c("1", "2", "3", "4", "5", "6", "7", "8")
rm(dist.to.feces.AIT)

# Melt and organize data frames for statistical tests.
dist.to.feces.BC.summary <- as.data.frame(t(dist.to.feces.BC.summary))
dist.to.feces.BC.summary <- melt(dist.to.feces.BC.summary)
dist.to.feces.BC.summary <- dist.to.feces.BC.summary[complete.cases(dist.to.feces.BC.summary),]

dist.to.feces.wUF.summary <- as.data.frame(t(dist.to.feces.wUF.summary))
dist.to.feces.wUF.summary <- melt(dist.to.feces.wUF.summary)
dist.to.feces.wUF.summary <- dist.to.feces.wUF.summary[complete.cases(dist.to.feces.wUF.summary),]

dist.to.feces.AIT.summary <- as.data.frame(t(dist.to.feces.AIT.summary))
dist.to.feces.AIT.summary <- melt(dist.to.feces.AIT.summary)
dist.to.feces.AIT.summary <- dist.to.feces.AIT.summary[complete.cases(dist.to.feces.AIT.summary),]

dist.to.feces.BC.summary$variable <- factor(as.character(dist.to.feces.BC.summary$variable), levels=desired_order_type)
dist.to.feces.wUF.summary$variable <- factor(as.character(dist.to.feces.wUF.summary$variable), levels=desired_order_type)
dist.to.feces.AIT.summary$variable <- factor(as.character(dist.to.feces.AIT.summary$variable), levels=desired_order_type)

# Perform statistical tests (ANOVA with Tukey's post hoc test).
leveneTest(dist.to.feces.BC.summary$value, dist.to.feces.BC.summary$variable) # p = 0.2928
summary(aov(dist.to.feces.BC.summary$value ~ dist.to.feces.BC.summary$variable)) # F=6.183, df=5, p<0.001
TukeyHSD(aov(dist.to.feces.BC.summary$value ~ dist.to.feces.BC.summary$variable), p.adjust="BH")

leveneTest(dist.to.feces.wUF.summary$value, dist.to.feces.wUF.summary$variable) # p = 0.9254
summary(aov(dist.to.feces.wUF.summary$value ~ dist.to.feces.wUF.summary$variable)) # F=5.531, df=5, p<0.001
TukeyHSD(aov(dist.to.feces.wUF.summary$value ~ dist.to.feces.wUF.summary$variable), p.adjust="BH")

leveneTest(dist.to.feces.AIT.summary$value, dist.to.feces.AIT.summary$variable) # p = 0.9569
summary(aov(dist.to.feces.AIT.summary$value ~ dist.to.feces.AIT.summary$variable)) # F=2.41, df=5, p=0.0526
TukeyHSD(aov(dist.to.feces.AIT.summary$value ~ dist.to.feces.AIT.summary$variable), p.adjust="BH")

#### ASV DISTRIBUTION AMONG SMALL INTESTINE, LARGE INTESTINE, AND FECES ####
# Subset only to ASVs present with at least 0.001% relative abundance.
temp.transform <- transform_sample_counts(pseq.rarefied, function(x) 100*x/sum(x))
temp.transform <- prune_taxa((taxa_sums(temp.transform)/nsamples(temp.transform)) > 0.001, temp.transform)

# Reclassify fecal samples from "large intestine" to "feces."
temp <- sample_data(temp.transform)
temp$Segment <- as.character(temp$Segment)
temp$Segment[temp$Type == "feces"] <- "feces"
temp$Segment <- as.factor(temp$Segment)
sample_data(temp.transform) <- temp
rm(temp)

# Create a list for data storage.
Venn.data <- data.frame(nrow=8, ncol=10)

# For each individual with fecal samples...
for(i in 1:length(ind.ids)){
  # Subset phyloseq object to that individual.
  grand.merge.temp <- subset_samples(temp.transform, Individual==ind.ids[[i]])
  
  # Merge samples by segment (producing 3 samples per individual: small intestine, large intestine, and feces).
  grand.merge.temp <- merge_samples(grand.merge.temp, "Segment")
  grand.merge.temp <- prune_taxa(taxa_sums(grand.merge.temp) > 0, grand.merge.temp)
  
  # Subset to only the small intestine.
  grand.merge.temp.SI <- subset_samples(grand.merge.temp, sample_names(grand.merge.temp)=="small.intestine")
  grand.merge.temp.SI <- prune_taxa(taxa_sums(grand.merge.temp.SI) > 0, grand.merge.temp.SI)
  
  # Subset to only the large intestine.
  grand.merge.temp.LI <- subset_samples(grand.merge.temp, sample_names(grand.merge.temp)=="large.intestine")
  grand.merge.temp.LI <- prune_taxa(taxa_sums(grand.merge.temp.LI) > 0, grand.merge.temp.LI)
  
  # Subset to only feces.
  grand.merge.temp.F <- subset_samples(grand.merge.temp, sample_names(grand.merge.temp)=="feces")
  grand.merge.temp.F <- prune_taxa(taxa_sums(grand.merge.temp.F) > 0, grand.merge.temp.F)
  
  # Create a data frame that counts number of ASVs unique to SI, unique to LI, unique to F, shared by SI/LI, shared by SI/F, shared by LI/F, and shared among all three sites.
  Venn.data[1,i] <- length(setdiff(setdiff(taxa_names(grand.merge.temp.SI), taxa_names(grand.merge.temp.LI)), taxa_names(grand.merge.temp.F))) # SI -- only
  Venn.data[2,i] <- length(setdiff(setdiff(taxa_names(grand.merge.temp.LI), taxa_names(grand.merge.temp.SI)), taxa_names(grand.merge.temp.F))) # LI -- only
  Venn.data[3,i] <- length(setdiff(setdiff(taxa_names(grand.merge.temp.F), taxa_names(grand.merge.temp.SI)), taxa_names(grand.merge.temp.LI))) # F -- only
  Venn.data[4,i] <- length(setdiff(intersect(taxa_names(grand.merge.temp.SI), taxa_names(grand.merge.temp.LI)), taxa_names(grand.merge.temp.F))) # SI & LI, but not F
  Venn.data[5,i] <- length(setdiff(intersect(taxa_names(grand.merge.temp.SI), taxa_names(grand.merge.temp.F)), taxa_names(grand.merge.temp.LI))) # SI & F, but not LI
  Venn.data[6,i] <- length(setdiff(intersect(taxa_names(grand.merge.temp.LI), taxa_names(grand.merge.temp.F)), taxa_names(grand.merge.temp.SI))) # LI & F, but not SI
  Venn.data[7,i] <- length(intersect(intersect(taxa_names(grand.merge.temp.LI), taxa_names(grand.merge.temp.F)), taxa_names(grand.merge.temp.SI))) # SI, LI, and F
  Venn.data[8,i] <- ntaxa(grand.merge.temp) # Total number of taxa present in this individual.
  
  rm(grand.merge.temp.F, grand.merge.temp.SI, grand.merge.temp.LI, grand.merge.temp)
}

# Name rows based on their Venn diagram category.
rownames(Venn.data) <- c("SI", "LI", "F", "SI+LI", "SI+F", "LI+F", "ALL", "TOTAL")
colnames(Venn.data) <- ind.ids

# Sum the total number of ASVs not in feces.
Venn.data <- as.data.frame(t(Venn.data))
Venn.data$NOT.F <- Venn.data$SI+Venn.data$LI+Venn.data$`SI+LI`
Venn.data$SI.NOT.F <- Venn.data$SI+Venn.data$`SI+LI`
Venn.data$LI.NOT.F <- Venn.data$LI+Venn.data$`SI+LI`
Venn.data <- as.data.frame(t(Venn.data))

# Divide ASV counts in each Venn diagram category by total number of ASVs in the individual, to get percentages of ASVs.
for(i in 1:ncol(Venn.data)){
  for(j in c(1:7,9:11)){
  Venn.data[j,i] <- 100*Venn.data[j,i]/Venn.data[8,i] }}

# Calculate means and standard deviations to obtain approximate ASV distribution throughout the intestine.
rowMeans(Venn.data)
#        SI         LI          F      SI+LI       SI+F       LI+F        ALL      TOTAL 
# 16.098999  30.219091   5.547450  13.003483   1.687217   9.682151  23.761608 287.500000 

rowSds(as.matrix(Venn.data))
#        SI        LI        F      SI+LI      SI+F      LI+F       ALL     TOTAL 
#  7.939141 16.259965  2.220451  6.756601  1.103028  7.866859 11.527924 44.493980

#### PREVALENCE AND ABUNDANCE OF FECAL-UNDETECTED ASVs #####
# This section create a summary of fecal-undetected taxa for each individual.
# First, determine the names of fecal-undetected ASVs in each individual. 
feces.undetected.names <- data.frame(matrix(ncol=0, nrow=0))

# For each individual with fecal samples...
for(i in 1:length(ind.ids)){
  # Subset the rarefied phyloseq object to that single individual.
  grandt <- subset_samples(pseq.rarefied, Individual==ind.ids[[i]])
  grandt <- prune_taxa(taxa_sums(grandt) > 0, grandt)
  
  # Create a phyloseq object with only the feces.
  grandf <- subset_samples(grandt, Type=="feces")
  grandf <- prune_taxa(taxa_sums(grandf) > 0, grandf)
  
  # Create a phyloseq object with only intestinal sites.
  grandnf <- subset_samples(grandt, Type != "feces")
  grandnf <- prune_taxa(taxa_sums(grandnf) > 0, grandnf)
  
  # Create a single data frame with taxa names in the intestinal sites and taxa names in the fecal sites.
  temp.list <- cbind.fill(rownames(as.data.frame(tax_table(grandnf))),
                          rownames(as.data.frame(tax_table(grandf))),
                          fill=NA)
  colnames(temp.list) <- c("int","feces")
  
  # Create a list of taxa present in intestinal sites but not in feces (i.e., fecal-undetected taxa).
  feces.undetected <- setdiff(temp.list$int, temp.list$feces)
  
  # Append to data storage.
  feces.undetected.names <- cbind.fill(feces.undetected.names, feces.undetected, fill=NA)
  rm(grandt, grandf, grandnf, temp.list, feces.undetected)
}

feces.undetected.names$init <- NULL

# Summarize the abundances and identities of fecal-undetected ASVs for each individual.
summary.ASV.indv <- list()

# For each individual with a fecal sample...
for(i in 1:length(ind.ids)){
  # Obtain phyloseq object for the given individual, containing only ASVs present in the intestines.
  temp.physeq <- transform_sample_counts(pseq.rarefied, function(x) 100*x/sum(x))
  temp.physeq <- subset_samples(temp.physeq, Individual==ind.ids[[i]])
  temp.physeq <- subset_samples(temp.physeq, Type != "feces")
  
  # Subset to only include fecal-undetected taxa.
  temp.physeq <- subset_taxa(temp.physeq, taxa_names(temp.physeq) %in% feces.undetected.names[,i])

  # Extract feature table and taxonomy table. Name the feature table based on intestinal site.
  temp.summary <- as.data.frame(t(otu_table(temp.physeq)))
  temp.tax.table <- as.data.frame(cbind(tax_table(temp.physeq)))
  colnames(temp.summary) <- sample_data(temp.physeq)$Type
  
  # Count the number of intestinal sites where this fecal-undetected taxa is present.
  temp.summary$prev <- rowSums(temp.summary != 0)
  
  # Merge taxonomy table and abundance information together.
  temp.summary <- transform(merge(temp.tax.table, temp.summary, by="row.names", all=TRUE), row.names=Row.names, Row.names=NULL)
  
  # Store information on fecal-undetected taxa.
  summary.ASV.indv[[i]] <- temp.summary
  rm(temp.summary, temp.physeq, temp.tax.table)
}

names(summary.ASV.indv) <- ind.ids

# Examine the distribution of fecal-undetected ASVs, based on prevalence and abundance at each intestinal site. #

# Make a new copy of summary.ASV.indv before manipulating the data.
summary.undetected <- summary.ASV.indv

# Define data frames for information on prevalence ("count") and abundance ("abund")
abund.summary <- data.frame(matrix(ncol=6, nrow=0))
count.summary <- data.frame(matrix(ncol=6, nrow=0))
total.summary <- data.frame(matrix(ncol=6, nrow=0))

# For each individual with fecal samples...
for(i in 1:length(ind.ids)){
  
  # Convert all zero-abundance values to NA, so that abundances of 0 don't affect the mean.
  summary.undetected[[i]][summary.undetected[[i]]==0] <- NA
  
  if(i !=7){ # This tag is here to account for the individual without a duodenal sample.
    # Calculate the mean abundance of fecal-undetected ASVs at each sampling site.
  abund <- summary.undetected[[i]] %>%
    summarise(mean.duo=mean(duodenum, na.rm=TRUE),
              mean.jej=mean(jejunum, na.rm=TRUE),
              mean.ile=mean(ileum, na.rm=TRUE),
              mean.cae=mean(caecum, na.rm=TRUE),
              mean.asc=mean(ascending.colon, na.rm=TRUE),
              mean.des=mean(descending.colon, na.rm=TRUE))
  
  # Calculate the total abundance of fecal-undetected ASVs at each sampling site.
  total <- summary.undetected[[i]] %>%
    summarise(sum.duo=sum(duodenum, na.rm=TRUE),
              sum.jej=sum(jejunum, na.rm=TRUE),
              sum.ile=sum(ileum, na.rm=TRUE),
              sum.cae=sum(caecum, na.rm=TRUE),
              sum.asc=sum(ascending.colon, na.rm=TRUE),
              sum.des=sum(descending.colon, na.rm=TRUE))
  
  # Calculate the number of fecal-undetected ASVs at each sampling site.
  count <- summary.undetected[[i]] %>%
    summarise(count.duo=sum(!is.na(duodenum)),
              count.jej=sum(!is.na(jejunum)),
              count.ile=sum(!is.na(ileum)),
              count.cae=sum(!is.na(caecum)),
              count.asc=sum(!is.na(ascending.colon)),
              count.des=sum(!is.na(descending.colon)))}
  
  # Repeat for the individual without a duodenal sample.
  if(i==7){
  abund <- summary.undetected[[i]] %>%
    summarise(
      mean.jej=mean(jejunum, na.rm=TRUE),
      mean.ile=mean(ileum, na.rm=TRUE),
      mean.cae=mean(caecum, na.rm=TRUE),
      mean.asc=mean(ascending.colon, na.rm=TRUE),
      mean.des=mean(descending.colon, na.rm=TRUE))
    abund <- cbind(rep(NA), abund)
  
    total <- summary.undetected[[i]] %>%
      summarise(sum.jej=sum(jejunum, na.rm=TRUE),
                sum.ile=sum(ileum, na.rm=TRUE),
                sum.cae=sum(caecum, na.rm=TRUE),
                sum.asc=sum(ascending.colon, na.rm=TRUE),
                sum.des=sum(descending.colon, na.rm=TRUE))
    total <- cbind(rep(NA), total)
    
  count <- summary.undetected[[i]] %>%
    summarise(
      count.jej=sum(!is.na(jejunum)),
      count.ile=sum(!is.na(ileum)),
      count.cae=sum(!is.na(caecum)),
      count.asc=sum(!is.na(ascending.colon)),
      count.des=sum(!is.na(descending.colon)))
  count <- cbind(rep(NA), count)}
  
  # Convert taxon counts to prevalence, expressed as a percent of the total number of fecal-undetected ASVs.
  for(j in 1:ncol(count)){
    count[,j] <- 100*count[,j]/nrow(summary.undetected[[i]])
  }
  
  # Rename columns for clarity.
  colnames(abund) <- c("mean.duo", "mean. jej", "mean.ile", "mean.cae", "mean.asc", "mean.des")
  colnames(count) <- c("count.duo", "count.jej", "count.ile", "count.cae", "count.asc", "count.des")
  colnames(total) <- c("total.duo", "total.jej", "total.ile", "total.cae", "total.asc", "total.des")
  
  # Append data to storage data frames.
  abund.summary <- rbind(abund.summary, abund)
  count.summary <- rbind(count.summary, count)
  total.summary <- rbind(total.summary, total)
  rm(abund, count, total) 
}

# Name columns based on intestinal site.
colnames(abund.summary) <- c("duo", "jej", "ile", "cae", "asc", "des")
colnames(count.summary) <- c("duo", "jej", "ile", "cae", "asc", "des")
colnames(total.summary) <- c("duo", "jej", "ile", "cae", "asc", "des")

# Create a data frame that shows, for each intestinal site, 1) prevalence of fecal-undetected ASVs; 2) mean relative abundance of fecal-undetected ASVs;
# and 3) mean total abundance of fecal-undetected ASVs (i.e. proportion of the community that is fecal-undetected).
feces.undetected.distr <- data.frame(
  prevalence = colMeans(count.summary, na.rm=TRUE),
  prevalence.sd = colSds(as.matrix(count.summary), na.rm=TRUE),
  abundance = colMeans(abund.summary, na.rm=TRUE),
  abundance.sd = colSds(as.matrix(abund.summary), na.rm=TRUE),
  total = colMeans(subset(total.summary, rownames(total.summary) !=4), na.rm=TRUE), ## NOTE AN OUTLIER SAMPLE IS BEING REMOVED
  total.sd = colSds(as.matrix(subset(total.summary, rownames(total.summary) !=4)), na.rm=TRUE)) ## NOTE AN OUTLIER SAMPLE IS BEING REMOVED
feces.undetected.distr$Site <- factor(rownames(feces.undetected.distr), levels=c("duo", "jej", "ile", "cae", "asc", "des"))
rm(abund.summary, count.summary, summary.undetected, total.summary)

# Calculate the total abundances of fecal-undetected taxa in each intestinal segment.
mean(feces.undetected.distr[c(1:3),"total"]) # 15.69%
sd(feces.undetected.distr[c(1:3),"total"]) # 7.006%
mean(feces.undetected.distr[c(4:6),"total"]) # 6.533%
sd(feces.undetected.distr[c(4:6),"total"]) # 5.1019%



#### TAXON BIAS IN FECAL-UNDETECTED ASVs #####
# Determine if any taxonomic groups are preferentially excluded from feces.
# Create a place to store the data.
taxon.bias.results <- list()

# For each individual with fecal samples...
for(i in 1:length(summary.ASV.indv)){
  # Replace unassigned classes with a new factor ("Unclassified")
  summary.ASV.indv[[i]]$Class <- forcats::fct_explicit_na(summary.ASV.indv[[i]]$Class, na_level = "Unclassified")
  
  # Subset a phyloseq object to exclude fecal samples and extract the taxonomy table.
  temp.pseq <- subset_samples(pseq.rarefied, Type !="feces" & Individual==ind.ids[[i]])
  temp.pseq <- prune_taxa(taxa_sums(temp.pseq) > 0, temp.pseq)
  temp.tax <- as.data.frame(cbind(tax_table(temp.pseq)))
  
  # Replace unassigned classes in the taxonomy table with the  same ("Unclassified") factor.
  temp.tax$Class <- forcats::fct_explicit_na(temp.tax$Class, na_level = "Unclassified")
  
  # Count the total number of ASVs present from each Class.
  total <- temp.tax %>%
    group_by(Class) %>%
    tally()
  colnames(total) <- c("Class", "TotalASVs")
  
  # Count the number of fecal-undetected ASVs from each Class.
  undetected <- summary.ASV.indv[[i]] %>%
    group_by(Class) %>%
    tally()
  colnames(undetected) <- c("Class", "FecesUndetected")
  
  # Merge data frames based on Class.
  summary <- merge(total, undetected, by="Class", all.x=TRUE)
  
  # Express the number of fecal-undetected ASVs as a percentage of the total number of ASVs.
  summary$Proportion <- 100*summary$FecesUndetected/summary$TotalASVs
  summary$TotalASVs <- NULL
  summary$FecesUndetected <- NULL
  
  # Store data.
  taxon.bias.results[[i]] <- summary
  rm(temp.pseq, temp.tax, total, undetected, summary)
}

# Merge results into a single data frame. This expresses, for each individual, the percentage of ASVs from each Class that are fecal-undetected.
taxon.bias <- merge(taxon.bias.results[[1]], taxon.bias.results[[2]], by="Class", all=TRUE)
taxon.bias <- merge(taxon.bias, taxon.bias.results[[3]], by="Class", all=TRUE)
taxon.bias <- merge(taxon.bias, taxon.bias.results[[4]], by="Class", all=TRUE)
taxon.bias <- merge(taxon.bias, taxon.bias.results[[5]], by="Class", all=TRUE)
taxon.bias <- merge(taxon.bias, taxon.bias.results[[6]], by="Class", all=TRUE)
taxon.bias <- merge(taxon.bias, taxon.bias.results[[7]], by="Class", all=TRUE)
taxon.bias <- merge(taxon.bias, taxon.bias.results[[8]], by="Class", all=TRUE)

# Clean data.
rownames(taxon.bias) <- taxon.bias$Class
taxon.bias$Class <- NULL
colnames(taxon.bias) <- ind.ids

# Calculate mean percentage of undetected ASVs from each class and the prevalence, or number of individuals that contain that Class.
taxon.bias$means = rowMeans(taxon.bias, na.rm=TRUE)
taxon.bias$prev = rowSums(taxon.bias[,1:8] != 0, na.rm=TRUE)
rm(taxon.bias.results)

# For statistical tests, subset to taxa that are present in at least 6 individuals.
taxon.bias.stats <- taxon.bias
taxon.bias.stats <- subset(taxon.bias.stats, prev > 6)

# Remove unnecessary data columns and melt.
taxon.bias.stats$prev <- NULL
taxon.bias.stats$means <- NULL
taxon.bias.stats <- melt(as.data.frame(t(taxon.bias.stats)))

# Determine if there is any significant taxonomic bias in fecal-undetected taxa.
summary(aov(value ~ variable, taxon.bias.stats)) # F=6.045, df=11, p<0.001
TukeyHSD(aov(value ~ variable, taxon.bias.stats), p.adjust="BH")

# Note from the data that Alphaproteobacteria and Verrucomicrobia are highly undetected; evaluate their abundance at each intestinal site.
taxon.bias.Alpha.Verruco <- transform_sample_counts(pseq.rarefied, function(x) 100*x/sum(x))
taxon.bias.Alpha.Verruco <- subset_taxa(taxon.bias.Alpha.Verruco, Class=="Alphaproteobacteria" | Phylum=="Verrucomicrobia")
tax_table(taxon.bias.Alpha.Verruco)[,"Class"][tax_table(taxon.bias.Alpha.Verruco)[,"Phylum"]=="Verrucomicrobia"] <- "Verrucomicrobia"
taxon.bias.Alpha.Verruco <- tax_glom(taxon.bias.Alpha.Verruco, taxrank="Class", NArm=FALSE)
taxon.bias.Alpha.Verruco <- as.data.frame(otu_table(taxon.bias.Alpha.Verruco))
taxon.bias.Alpha.Verruco$Type <- bg.metadata$Type
colnames(taxon.bias.Alpha.Verruco) <- c("Alphaproteobacteria", "Verrucomicrobia", "Type")
taxon.bias.Alpha.Verruco <- melt(taxon.bias.Alpha.Verruco)

#### FECAL-ONLY TAXA ####
feces.unique <- list()
# For each individual with fecal samples...
for(i in 1:length(ind.ids)){
  # Subset to fecal samples.
  temp.feces <- subset_samples(pseq.rarefied, Type=="feces" & Individual==ind.ids[[i]])
  temp.feces <- prune_taxa(taxa_sums(temp.feces) > 0, temp.feces)
  
  # Subset to intestinal samples.
  temp.int <- subset_samples(pseq.rarefied, Type!="feces" & Individual==ind.ids[[i]])
  temp.int <- prune_taxa(taxa_sums(temp.int) > 0, temp.int)
  
  # Identify fecal-unique taxa.
  feces.only.names <- setdiff(taxa_names(temp.feces), taxa_names(temp.int))
  
  # Subset a phyloseq object to fecal-unique taxa.
  temp.pseq <- transform_sample_counts(pseq.rarefied, function(x) 100*x/sum(x))
  temp.pseq <- subset_samples(temp.pseq, Individual==ind.ids[[i]] & Type=="feces")
  temp.pseq <- subset_taxa(temp.pseq, taxa_names(temp.pseq) %in% feces.only.names)
  
  # Append the taxonomy table and OTU table. Rename unassigned Classes to "Unclassified."
  feces.only.summary <- cbind(as.data.frame(cbind(tax_table(temp.pseq))), as.data.frame(t(otu_table(temp.pseq))))
  feces.only.summary$Class <- forcats::fct_explicit_na(feces.only.summary$Class, na_level = "Unclassified")
  
  # Store data.
  feces.unique[[i]] <- feces.only.summary
  feces.unique[[i]][,1:6] <- NULL
  rm(temp.feces, temp.int, feces.only.names, temp.pseq, feces.only.summary)
}

# Merge individual summaries into a single data frame.
feces.unique.summary <- transform(merge(feces.unique[[1]], feces.unique[[2]], by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
feces.unique.summary <- transform(merge(feces.unique.summary, feces.unique[[3]], by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
feces.unique.summary <- transform(merge(feces.unique.summary, feces.unique[[4]], by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
feces.unique.summary <- transform(merge(feces.unique.summary, feces.unique[[5]], by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
feces.unique.summary <- transform(merge(feces.unique.summary, feces.unique[[6]], by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
feces.unique.summary <- transform(merge(feces.unique.summary, feces.unique[[7]], by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
feces.unique.summary <- transform(merge(feces.unique.summary, feces.unique[[8]], by=0, all=TRUE), row.names=Row.names, Row.names=NULL)

# Add taxonomy information.
feces.unique.summary <- transform(merge(feces.unique.summary, as.data.frame(cbind(tax_table(pseq.rarefied))), by=0, all=FALSE), row.names=Row.names, Row.names=NULL)

# Label unassigned classes as "Unclassified."
feces.unique.summary$Class <- forcats::fct_explicit_na(feces.unique.summary$Class, na_level = "Unclassified")
rm(feces.unique)

# Summarize total abundance of fecal-unique classes at each sampling site.
feces.unique.class <- feces.unique.summary %>% dplyr::group_by(Class) %>%
  dplyr::summarise(sum34=sum(X34F, na.rm=TRUE),
            sum38=sum(X38F, na.rm=TRUE),
            sum40=sum(X40F, na.rm=TRUE),
            sum42=sum(X42F, na.rm=TRUE),
            sum46=sum(X46F, na.rm=TRUE),
            sum47=sum(X47F, na.rm=TRUE),
            sum63=sum(X63F, na.rm=TRUE),
            sum65=sum(X65F, na.rm=TRUE))

# Summarize total abundances of fecal-unique taxa at each site.
feces.unique.class <- rbind(feces.unique.class, c(NA, colSums(feces.unique.class[,2:9])))

# Convert everything back to numeric because the last command undid that...
feces.unique.class$sum34 <- as.numeric(feces.unique.class$sum34)
feces.unique.class$sum38 <- as.numeric(feces.unique.class$sum38)
feces.unique.class$sum40 <- as.numeric(feces.unique.class$sum40)
feces.unique.class$sum42 <- as.numeric(feces.unique.class$sum42)
feces.unique.class$sum46 <- as.numeric(feces.unique.class$sum46)
feces.unique.class$sum47 <- as.numeric(feces.unique.class$sum47)
feces.unique.class$sum63 <- as.numeric(feces.unique.class$sum63)
feces.unique.class$sum65 <- as.numeric(feces.unique.class$sum65)

# Calculate the mean and SD for class-level abundances of fecal-unique taxa.
feces.unique.class$mean <- rowMeans(feces.unique.class[,2:9])
feces.unique.class$sd <- rowSds(as.matrix(feces.unique.class[,2:9]))








#### FECAL ASV ABUNDANCES: FOLD CHANGE AGAINST OTHER INTESTINAL SITES #####
# Summarize average fold-changes between ASV relative abundance in each intestinal segment and feces, per individual.
# To exclude 'feces' from segment based analysis, create a new level for the 'segment' variable.
temp.physeq <- pseq.rarefied
temp <- as.data.frame(cbind(sample_data(temp.physeq)))
levels(temp$Segment) <- c("large.intestine", "small.intestine", "feces")
temp$Segment[temp$Type=="feces"] <- "feces"
sample_data(temp.physeq) <- temp
rm(temp)

# Merge samples based on both individual identity and intestinal segment.
variable1 = as.character(get_variable(temp.physeq, "Individual"))
variable2 = as.character(get_variable(temp.physeq, "Segment"))
sample_data(temp.physeq)$MergeVariable <- mapply(paste0, variable1, variable2, 
                                                        collapse = "_")
temp.physeq <- merge_samples(temp.physeq, "MergeVariable")

# Rename variables that got replaced with "NA" in the merge.
temp <- as.data.frame(cbind(sample_data(temp.physeq)))
temp$Individual <- substring(rownames(temp), 1, 2)
temp$Individual <- as.factor(temp$Individual)
temp$Segment <- substring(rownames(temp), 3, 20)
temp$Segment <- as.factor(temp$Segment)
sample_data(temp.physeq) <- temp
rm(temp, variable1, variable2)

# Convert to relative abundance.
temp.physeq <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))

# Define data storage.
fold.seg.change <- list() # For fold changes
fold.seg.orig.abund <- list() # For raw abundances

# For each individual with fecal samples...
for(i in 1:length(ind.ids)){
  # Subset to that individual
  temp <- subset_samples(temp.physeq, Individual==ind.ids[[i]])
  
  # Extract taxonomy and OTU tables.
  temp.tax <- as.data.frame(cbind(tax_table(temp)))
  temp.otu <- as.data.frame(t(otu_table(temp)))
  
  # Name OTU table by segment.
  colnames(temp.otu) <- sample_data(temp)$Segment
  
  # Create a new data frame and calculate abundance changes (small intestine / feces) and (large intestine / feces)
  df <- data.frame(matrix(ncol=2, nrow=length(rownames(temp.otu))))
  df[,1] <- (temp.otu$small.intestine)/(temp.otu$feces)
  df[,2] <- (temp.otu$large.intestine)/(temp.otu$feces)
  colnames(df) <- c("SI.F", "LI.F")
  
  # Append taxonomy information to OTU table.
  temp.table <- cbind(temp.tax, temp.otu)
  
  # Append taxonomy information to fold-change table.
  temp.table2 <- cbind(temp.tax, df)
  
  # Store data.
  fold.seg.change[[i]] <- temp.table2
  fold.seg.orig.abund[[i]] <- temp.table
  
  rm(temp, temp.tax, temp.otu, df, temp.table, temp.table2)
}

# Create a place to store data - now calculating averages for each intestinal segment.
fold.ASV.seg.averages <- list()

# Name both intestinal segments.
typesnf <- c("small.intestine", "large.intestine")

# For each segment...
for(i in 1:length(typesnf)){
  # Create a data frame with the number of ASVs in the sample.
  df <- data.frame(matrix(ncol=0, nrow=1109))
  
  # Extract taxonomy table.
  temp.tax <- as.data.frame(cbind(tax_table(pseq.rarefied)))
  
  # For each individual with fecal samples...
  for(j in 1:length(ind.ids)){
    # Extract the information about fold change for each segment. 
    df <- cbind(df, fold.seg.change[[j]][,(i+6)])  
  }
  
  # Label this data frame (with taxa as rows and each individual as a column).
  colnames(df) <- ind.ids
  
  # Replace infinite values (fold change divided by 0) with NA.
  df <- as.matrix(df)
  df[!is.finite(df)] <- NA
  df <- as.data.frame(df)
  
  # Add taxonomy information.
  df <- cbind(temp.tax, df)
  rownames(df) <- taxa_names(pseq.rarefied)
  fold.ASV.seg.averages[[i]] <- df
  rm(df)
}
names(fold.ASV.seg.averages) <- typesnf

# Calculate mean fold changes across all individuals.
for(i in 1:length(fold.ASV.seg.averages)){
  fold.ASV.seg.averages[[i]]$Mean <- rowMeans(fold.ASV.seg.averages[[i]][,c(7:14)], na.rm=TRUE)
}

# Create a data table that combines average from small intestinal fold changes and large intestinal fold changes.
fold.ASV.seg.changes.summary <- data.frame(matrix(ncol=0, nrow=1109))

# This puts small and large intestinal measurements into the same data frame.
for(i in 1:length(fold.ASV.seg.averages)){
  fold.ASV.seg.changes.summary <- cbind(fold.ASV.seg.changes.summary, fold.ASV.seg.averages[[i]]$Mean)
}

# Remove infinite values.
fold.ASV.seg.changes.summary <- as.matrix(fold.ASV.seg.changes.summary)
fold.ASV.seg.changes.summary[!is.finite(fold.ASV.seg.changes.summary)] <- NA
fold.ASV.seg.changes.summary <- as.data.frame(fold.ASV.seg.changes.summary)

# Clean data.
colnames(fold.ASV.seg.changes.summary) <- typesnf
rownames(fold.ASV.seg.changes.summary) <- taxa_names(pseq.rarefied)

# Add taxonomy information.
fold.ASV.seg.changes.summary <- cbind(as.data.frame(cbind(tax_table(pseq.rarefied))),
                                      fold.ASV.seg.changes.summary)

rm(fold.ASV.seg.changes, fold.seg.change)

#### ASV ABUNDANCE CORRELATIONS AMONG INTESTINAL SITES ####
# Remove ASVs with fewer than 15 counts, convert to relative abundance, and extract OTU table.
temp.pseq <- transform_sample_counts(pseq.rarefied, function(x) 100*x/sum(x))
temp.pseq <- prune_taxa((taxa_sums(temp.pseq)/nsamples(temp.pseq)) > 0.001, temp.pseq)

temp.otu <- as.data.frame(otu_table(temp.pseq))
temp.otu$Type <- bg.metadata$Type

# Create individual OTU tables for each intestinal segment, without any additional metadata.
temp.duodenum <- subset(temp.otu, Type=="duodenum")
temp.jejunum <- subset(temp.otu, Type=="jejunum")
temp.ileum <- subset(temp.otu, Type=="ileum")
temp.caecum <- subset(temp.otu, Type=="caecum")
temp.asc.colon <- subset(temp.otu, Type=="ascending.colon")
temp.des.colon <- subset(temp.otu, Type=="descending.colon")
temp.feces <- subset(temp.otu, Type=="feces")

# Calculate mean abundance for each ASV in each intestinal site.
means.duodenum <- colMeans(temp.duodenum[,1:length(temp.duodenum)-1])
means.jejunum <- colMeans(temp.jejunum[,1:length(temp.jejunum)-1])
means.ileum <- colMeans(temp.ileum[,1:length(temp.ileum)-1])
means.caecum <- colMeans(temp.caecum[,1:length(temp.caecum)-1])
means.asc.colon <- colMeans(temp.asc.colon[,1:length(temp.asc.colon)-1])
means.des.colon <- colMeans(temp.des.colon[,1:length(temp.des.colon)-1])
means.feces <- colMeans(temp.feces[,1:length(temp.feces)-1])

# Compile data into a single table.
means <- as.data.frame(rbind(means.duodenum, means.jejunum, means.ileum, means.caecum, means.asc.colon, means.des.colon, means.feces))
rownames(means) <- c("duodenum","jejunum", "ileum", "caecum", "asc.colon", "des.colon", "feces")
means <- as.data.frame(t(means))

# Replace zeroes with NA - we are only correlating if taxa are present in both sites.
means[means==0] <- NA

# Calculate Spearman correlations between taxon abundances. (These are column-column correlations).
means.corr <- Hmisc::rcorr(as.matrix(means), type="spearman")

# Add taxonomy information.
means <- transform(merge(as.data.frame(tax_table(temp.pseq)), means, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
rm(temp.pseq, temp.otu, temp.asc.colon, temp.caecum, temp.des.colon, temp.duodenum, temp.feces, temp.ileum,
   temp.jejunum, means.asc.colon, means.caecum, means.des.colon, means.duodenum, means.feces, means.ileum, means.jejunum)
