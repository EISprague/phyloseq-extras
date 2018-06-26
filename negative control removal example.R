library(phyloseq)
data(GlobalPatterns)

###Sorting OTUs into negatives and samples----
samps = subset_samples(physeq, Depth != "control")
negs = subset_samples(physeq, Depth == "control")

#Generating list of negative OTUs
#OTU.all = sort(taxa_names(physeq)) ###All taxa
#OTU.negs = sort(taxa_names(prune_taxa(taxa_sums(negs) > 0, negs))) ###Taxa in both negs and samps
#OTU.samp.only = sort(OTU.all[!(OTU.all %in% OTU.negs)]) ###Taxa in samps only
###Ratios of all OTUs found in samples AND in controls for a given extraction day;
### higher than 1 means there were more sequences in the sample than its corresponding control so you should probably keep that OTU
neg.ratios = neg_ratios_together(samps, negs)

less.than.OTUs = less_than_taxa(neg.ratios, 1) #2nd argument determines ratio cutoff. For example, 1 will find OTUs with a consistent 1:1 ratio or lower, while 0.5 will find OTUs that consistently appeared up to 2x as often in the negative controls compared to the samples

###Removing OTUs identified as contaminants
samps.final = prune_taxa(!(taxa_names(samps) %in% less.than.OTUs), samps) %>%
  prune_taxa(taxa_sums(.) > 0, .)

rm(samps, negs, OTU.negs, OTU.all, OTU.samp.only, neg.ratios, less.than.OTUs) #these objects are no longer necessary
