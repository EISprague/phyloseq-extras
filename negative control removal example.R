source("negative control functions.R")
library(phyloseq)
library(magrittr)
data(GlobalPatterns)



sample_data(GlobalPatterns)$Control_type[sample_data(GlobalPatterns)$SampleType == "Mock"] = "negative"
sample_data(GlobalPatterns)$Control_type[is.na(sample_data(GlobalPatterns)$Control_type)] = "sample"

sample_data(GlobalPatterns)$Extr_date[sample_data(GlobalPatterns)$Control_type == "negative"] = c("A", "B", "C")
set.seed(711)
sample_data(GlobalPatterns)$Extr_date[sample_data(GlobalPatterns)$Control_type == "sample"] = sample(c("A", "B", "C"), size = 23, replace = T)

cleaned = remove_contaminants(GlobalPatterns, 1, "negative", "Extr_date") #returns phyloseq object with contaminant OTUs removed in one step (using the following 3 functions)

#If you want to use one of the other functions instead for exploratory purposes, you'll have to break your phyloseq object into one consisting of only samples, and one consistent only of negative controls
samps = subset_samples(GlobalPatterns, Control_type != "negative")
negs = subset_samples(GlobalPatterns, Control_type == "negative")
contaminants = contaminants_by_date_ratio(samps, negs, 1, "Extr_date") #returns vector of contaminant OTUs as determined by ratio
ratios.matrix = neg_ratios_matrix(samps, negs, "Extr_date") #returns single matrix of ratios of sample sequence counts to negative control counts
ratios.list = neg_OTU_ratios_by_date(samps, negs, "Extr_date") #returns list of ratio matrices, one element for each extraction date


