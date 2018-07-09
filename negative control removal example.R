source("negative control functions.R")
library(phyloseq)
library(magrittr)
data(GlobalPatterns)

#First set up a column in the sample_data() slot that says whether a given sample is a negative control or a sample. The column is called "Control_type"
sample_data(GlobalPatterns)$Control_type[sample_data(GlobalPatterns)$SampleType == "Mock"] = "negative"
sample_data(GlobalPatterns)$Control_type[is.na(sample_data(GlobalPatterns)$Control_type)] = "sample"

#Then set up a column called "Extr_date" that assigns a unique date to each negative control, and the same date for all of the corresponding samples. As in this example, it doesn't have to be an actual date, just an ID works too.
#GlobalPatterns has 3 mock communities, which for demo purposes I am treating as negatives controls, one for each of three extraction batches.
sample_data(GlobalPatterns)$Extr_date[sample_data(GlobalPatterns)$Control_type == "negative"] = c("A", "B", "C")
#Randomly assigning extraction batch ID to samples (in my project we randomly grouped samples for extraction)
set.seed(711)
sample_data(GlobalPatterns)$Extr_date[sample_data(GlobalPatterns)$Control_type == "sample"] = sample(c("A", "B", "C"), size = 23, replace = T)

#remove_contaminants(physeq, ratio, negcontrol, date_colname) is the easiest function to use. See "negative control functions.R" for explanation of each argument and assumptions about the physeq object.
cleaned = remove_contaminants(GlobalPatterns, 1, "negative", "Extr_date") #returns phyloseq object with negative controls and contaminant OTUs removed in one step (using the following 3 functions). If there were OTUs that were present ONLY in the negative controls, they will still be included in the output, but should be easy to remove with prune_taxa() since all their abundance values will be zero at this point.

#If you want to use one of the other functions instead for exploratory purposes, you'll have to break your phyloseq object into samps consisting of only samples, and negs consisting only of negative controls
samps = subset_samples(GlobalPatterns, Control_type != "negative")
negs = subset_samples(GlobalPatterns, Control_type == "negative")
contaminants = contaminants_by_date_ratio(samps, negs, 1, "Extr_date") #returns vector of OTU IDs determined to be contaminants by ratio
ratios.matrix = neg_ratios_matrix(samps, negs, "Extr_date") #returns single matrix of ratios of sample sequence counts to negative control counts
ratios.list = neg_OTU_ratios_by_date(samps, negs, "Extr_date") #returns list of ratio matrices, one element for each extraction date
