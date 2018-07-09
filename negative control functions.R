remove_contaminants = function(physeq, ratio, negcontrol, date_colname) {
  #Assumes there is a column "Control_type" in sample_data(physeq) to identify what type of control each sample is (or isn't)
  #Assumes all samples are either negative controls or samples (positive controls are treated like samples)
  #Assumes there is a column "X.SampleID" with sample names
  #physeq is the phyloseq object
  #ratio is the cutoff ratio of sample sequences : negative control sequences. So if you want OTUs to be considered contaminants if they are equally or more abundant in the negative control than in the corresponding samples, use ratio = 1. If you want them to be considered contaminants if there is 1 sequence in the sample for every 2 or more sequences in the negative control, use ratio = 1/2.
  #negcontrol is a string that tells the function which value in the Control_type column means a given sample is a negative control
  #date_colname is a string telling the function which column in sample_data(physeq) will group the samples with negative controls by date
  #The ratio has to be consistent for every sample. So for a given OTU, if there are 5 samples where the ratio is less than the ratio cutoff given by the user (and therefore the OTU is likely a contaminant), but a sixth sample where the ratio is higher, the OTU will NOT be included in the final output and therefore will not be considered a contaminant.
  require(phyloseq)
  orig.sample.data = as(sample_data(physeq), "data.frame") #Turns sample_data(physeq) into data frame
  samps = physeq
  sample_data(samps) = sample_data(orig.sample.data[orig.sample.data$Control_type != negcontrol, ]) #creates phyloseq object with all non-negative controls
  negs = physeq
  sample_data(negs) = sample_data(orig.sample.data[orig.sample.data$Control_type == negcontrol, ]) #creates phyloseq object with only negative controls
  
  contaminants = contaminants_by_date_ratio(samps, negs, ratio, date_colname) #makes vector of OTUs identified as contaminants
  out = prune_taxa(!(taxa_names(samps) %in% contaminants), samps) #removes contaminant OTUs from samps object
  return(out) #returns a phyloseq object with contaminant OTUs removed and no negative controls
}

contaminants_by_date_ratio = function(samps, negs, ratio, date_colname) {
  #Returns vector of OTU names of contaminants as determined by the ratio input.
  ratio.mat = neg_ratios_matrix(samps, negs, date_colname) #makes single matrix of OTUs (rows) and samples (columns), and the ratio of sample sequences : negative control sequences. If the ratio is zero, it means that OTU did NOT appear in BOTH the sample and its corresponding negative control
  out = rownames(ratio.mat[apply(X = ratio.mat, MARGIN = 1, 
                                 FUN = function(y) { all(y <= ratio) } ), 
                           ]) #trims matrix so OTUs are only returned/identified as contaminants if ALL of their ratios are less than the ratio argument. In other words, for every pair of sample and negative control, a given OTU has to CONSISTENTLY be more common in the negative control than the sample (or it had a ratio of zero). A single instance of the OTU having more sequences in the sample than the corresponding negative control will prevent it from being ID'd as a contaminant.
  return(out)
}

neg_ratios_matrix = function(samps, negs, date_colname) {
  #Calls neg_OTU_ratios_by_date, combines list into a single matrix to make everything else easier
  ###Assumes the user has left the sample name column in sample_data(physeq) as "X.SampleID"
  input = neg_OTU_ratios_by_date(samps, negs, date_colname) #gets the list from the previous function
  neg.OTUs = lapply(input, function(x) { unlist(dimnames(x)[1]) } ) #makes a list of vectors; each vector contains all the OTUs found in both the negative control and samples for an extraction date
  all.neg.OTUs = sort(unique(unlist(neg.OTUs))) #makes a single vector containing all unique OTU IDs that were found to occur in both a sample and its corresponding negative control at least once
  new.data = matrix(nrow = length(all.neg.OTUs), ncol = length(sample_data(samps)$X.SampleID),
                    dimnames = list(all.neg.OTUs, sample_data(samps)$X.SampleID)) #makes an empty matrix with one row for each OTU and one column for each sample (negative controls don't need their own columns anymore)
  
  for (i in 1:length(input)) { #starts a loop with one cycle for each extraction date
    rnames = rownames(input[[i]]) #generates vector of OTU IDs from matrix for current cycle's date
    cnames = colnames(input[[i]]) #generates vector of sample IDs from matrix for current cycle's date
    new.data[rnames, cnames] = input[[i]]
  }
  
  new.data[is.na(new.data)] = 0 #replaces any remaining NA values with 0
  new.data = new.data[rowSums(new.data) > 0, ] #double check that all OTUs have a non-zero ratio for at least one sample
  new.data = signif(new.data, digits = 3) #more than 3 sig figs seemed annoying and unnecessary
  return(new.data)
}

neg_OTU_ratios_by_date = function(samps, negs, date_colname) {
  #Outputs list of OTU ratio matrices of samples to controls
  #Each element of output list is an extraction date, each column of each element is a sample, each row is an OTU
  ###If there's an extraction date with only one sample, this won't work
  dates = unique(unlist(samps@sam_data@.Data[samps@sam_data@names == date_colname])) #character vector of extraction dates
  ratio.list = setNames(vector("list", length(dates)), dates) #empty list to fill, with dates as names of elements
  for(i in 1:length(ratio.list)) { #loop that makes each element of tab.ratio.comb into a matrix
    negs.date = negs@otu_table@.Data[, unlist(negs@sam_data@.Data[negs@sam_data@names == date_colname]) == dates[i]] #subsets the negs object to the negative control used on the current loop's date
    negs.OTUs = negs.date[negs.date > 0] #removes all OTUs that didn't appear on current date
    negs.OTUs = negs.OTUs[order(names(negs.OTUs))] #alphabetizes the OTUs
    samps.date = as.matrix(samps@otu_table@.Data[, unlist(samps@sam_data@.Data[samps@sam_data@names == date_colname]) == dates[i]]) #subsets the samps object to the samples extracted on the current date
    samps.date = samps.date[rownames(samps.date) %in% names(negs.OTUs), ] #subsets to OTUs that also appeared in that date's negative control
    samps.date = samps.date[order(names(negs.OTUs)), ] #puts OTUs in same order as negs.OTUs object
    date.ratio = samps.date/negs.OTUs #ratio calculation; ratios > 1 for a given OTU & sample mean that there were more sequences in the sample than the negative control for that extraction date
    date.ratio = date.ratio[rowSums(date.ratio) > 0, ] #removes OTUs in the matrix that did not generate a ratio
    ratio.list[[ dates[i] ]] = date.ratio #adds matrix from previous step to list as element, and restarts the loop
  }
  
  return(ratio.list) #You now have a list with a matrix element for each extraction date.
}
