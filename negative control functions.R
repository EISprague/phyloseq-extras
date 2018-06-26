by_date = "Extraction.date"
i = 1
neg_OTU_ratios_by_date = function(samps, negs, by_date = "Extraction.date") {
  #Outputs list of OTU ratio matrices of samples to controls
  #Each element of output list is an extraction date, each column of each element is a sample, each row is an OTU
  ###If there's an extraction date with only one sample, this won't work
  dates = levels(unlist(samps@sam_data@.Data[samps@sam_data@names == by_date])) #character vector of extraction dates
  ratio.list = setNames(vector("list", length(dates)), dates) #empty list to fill, with dates as names of elements
  for(i in 1:length(ratio.list)) { #loop that makes each element of tab.ratio.comb into a matrix
    negs.date = negs@otu_table@.Data[, unlist(negs@sam_data@.Data[negs@sam_data@names == by_date]) == dates[i]] #subsets the negs object to the negative control used on the current loop's date
    negs.OTUs = negs.date[negs.date > 0] #removes all OTUs that didn't appear on current date
    negs.OTUs = negs.OTUs[order(names(negs.OTUs))] #alphabetizes the OTUs
    samps.date = as.matrix(samps@otu_table@.Data[, unlist(samps@sam_data@.Data[samps@sam_data@names == by_date]) == dates[i]]) #subsets the samps object to the samples extracted on the current date
    samps.date = samps.date[rownames(samps.date) %in% names(negs.OTUs), ] #subsets to OTUs that also appeared in that date's negative control
    samps.date = samps.date[order(names(negs.OTUs)), ] #puts OTUs in same order as negs.OTUs object
    date.ratio = samps.date/negs.OTUs #ratio calculation; ratios > 1 for a given OTU & sample mean that there were more sequences in the sample than the negative control for that extraction date
    date.ratio = date.ratio[rowSums(date.ratio) > 0, ] #removes OTUs in the matrix that did not generate a ratio
    ratio.list[[ dates[i] ]] = date.ratio #adds matrix from previous step to list as element, and restarts the loop
  }
  
  return(ratio.list) #You now have a list with a matrix element for each extraction date.
}
rm(by_date, i, dates, ratio.list, negs.date, negs.OTUs, samps.date, date.ratio)

neg_ratios_matrix = function(samps, negs) {
  #Calls neg_OTU_ratios_by_date, combines list into a single matrix to make everything else easier
  ###Assumes the user has left the sample name column in their sample_data(physeq) as "X.SampleID"
  input = neg_OTU_ratios_by_date(samps, negs) #gets the list from the previous function
  neg.OTUs = lapply(input, function(x) { unlist(dimnames(x)[1]) } ) #makes a list of vectors; each vector contains all the OTUs found in both the negative controls and samples for an extraction date
  all.neg.OTUs = sort(unique(unlist(neg.OTUs))) #makes a single vector containing all unique OTU IDs that were found to occur in both a sample and its corresponding negative control at least once
  new.data = matrix(nrow = length(all.neg.OTUs), ncol = length(sample_data(samps)$X.SampleID),
                    dimnames = list(all.neg.OTUs, sample_data(samps)$X.SampleID)) #makes an empty matrix with one row for each OTU and one column for each sample (negative controls don't need their own columns anymore)
  for (i in 1:length(input)) { #starts a loop with one cycle for each extraction date
    rnames = rownames(input[[i]]) #generates vector of OTU IDs from matrix for current cycle's date
    cnames = colnames(input[[i]]) #generates vector of sample IDs from matrix for current cycle's date
    #The following two loops might not actually be necessary, but this is how I wrote it several months ago when I was still figuring stuff out
    for (j in 1:nrow(input[[i]])) { #starts a sub-loop, one cycle for each OTU ID from current loop's date
      for (k in 1:ncol(input[[i]])) { #starts a sub-sub-loop, one cycle for each sample ID from current loop's date
        new.data[rnames[j], cnames[k]] = input[[i]][j, k] #Fills cells in matrix with ratio values from the original list
      }
    }
  }
  new.data[is.na(new.data)] = 0 #replaces any remaining NA values with 0
  new.data = new.data[rowSums(new.data) > 0, ] #double check that all OTUs have a ratio for at least one sample
  new.data = signif(new.data, digits = 3) #more than 3 sig figs seemed annoying and unnecessary
  return(new.data)
}

less_than_taxa = function(rat.mat, ratios) {
  #Returns list of OTU name vectors, one vector per number in ratio argument; rat.mat is matrix of neg ratios
  # produced by neg_ratios_together; only finds OTUs that are always <= in given ratio 
  out = vector("list", length = length(ratios))
  names(out) = paste("OTU", ratios, "or.ls", sep = ".")
  for(i in 1:length(out)) {
    out[[i]] = rownames(rat.mat[apply(X = rat.mat, MARGIN = 1, FUN = function(y) { (all(y <= ratios[i])) }),])
  }
  
  if (length(out) > 1) {
    return(out) }
  
  if (length(out) == 1) {
    return(unname(unlist(out))) }
}
