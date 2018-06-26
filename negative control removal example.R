library(phyloseq)
library(magrittr)
data(GlobalPatterns)

###Sorting OTUs into negatives and samples----
samps = subset_samples(physeq, Depth != "control")
negs = subset_samples(physeq, Depth == "control")
sample_data(physeq)$Control_type[sample_data(physeq)$Depth == "control"] = "negative"
sample_data(physeq)$Control_type[is.na(sample_data(physeq)$Control_type)] = "sample"
#Generating list of negative OTUs
#OTU.all = sort(taxa_names(physeq)) ###All taxa
#OTU.negs = sort(taxa_names(prune_taxa(taxa_sums(negs) > 0, negs))) ###Taxa in both negs and samps
#OTU.samp.only = sort(OTU.all[!(OTU.all %in% OTU.negs)]) ###Taxa in samps only
###Ratios of all OTUs found in samples AND in controls for a given extraction day;
### higher than 1 means there were more sequences in the sample than its corresponding control so you should probably keep that OTU
#neg.ratios.testing = neg_OTU_ratios_by_date(samps, negs)
#neg.ratios = neg_ratios_matrix(samps, negs)

less.than.OTUs = contaminants_by_date_ratio(samps, negs, 1, "Extraction.date") #2nd argument determines ratio cutoff. For example, 1 will find OTUs with a consistent 1:1 ratio or lower, while 0.5 will find OTUs that consistently appeared up to 2x as often in the negative controls compared to the samples

###Removing OTUs identified as contaminants
samps.final = prune_taxa(!(taxa_names(samps) %in% less.than.OTUs), samps) %>%
  prune_taxa(taxa_sums(.) > 0, .)

samps.final = remove_contaminants(physeq, 1, "negative", "Extraction.date")

rm(samps, negs, OTU.negs, OTU.all, OTU.samp.only, neg.ratios, less.than.OTUs) #these objects are no longer necessary
