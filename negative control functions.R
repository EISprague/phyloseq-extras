neg_OTU_ratios_groups = function(samps, negs) {
  #Outputs list of OTU ratios of samples to controls, each element of list is an extraction date
  #Extraction date column in sample data needs to be called "Extraction.date"
  ###If there's an extraction date with only one sample, this won't work
  tab.ratio.comb = vector("list", length(levels(unlist(samps@sam_data@.Data[samps@sam_data@names == "Extraction.date"])))) #empty list to fill
  dates = levels(unlist(samps@sam_data@.Data[samps@sam_data@names == "Extraction.date"])) #vector of extraction dates
  names(tab.ratio.comb) = dates #each element of the list has a date for a name
  for(i in 1:length(dates)) { #loop that creates a matrix for each date
    negs.date = as.matrix(negs@otu_table@.Data[, (unlist(negs@sam_data@.Data[negs@sam_data@names == "Extraction.date"]) == dates[i])]) #subsets the negs object to the negative control used on the current cycle's date
    tab.negs = negs.date[rowSums(negs.date) > 0, ] #removes all OTUs that didn't appear on current date
    tab.negs = tab.negs[order(names(tab.negs))] #alphabetizes the OTUs
    samps.date = as.matrix(samps@otu_table@.Data[, unlist(samps@sam_data@.Data[samps@sam_data@names == "Extraction.date"]) == dates[i]]) #subsets the samps object to the samples extracted on the current date
    samps.date = samps.date[rownames(samps.date) %in% names(tab.negs), ] #subsets to OTUs that also appeared in that date's negative control
    samps.date = samps.date[order(rownames(samps.date)), ] #orders the matrix by OTU ID
    tab.ratio = samps.date/tab.negs #ratio calculation
    tab.ratio = tab.ratio[rowSums(tab.ratio) > 0, ] #makes double-y sure there are no OTUs left in the matrix that did not generate a ratio
    tab.ratio.comb[[dates[i]]] = tab.ratio #adds matrix from previous step to list as element, and restarts the loop
  }
  out = tab.ratio.comb
  return(out) #You now have a list with a matrix element for each extraction date. It gets fed into the next function.
}

neg_ratios_together = function(samps, negs) { #Calls neg_OTU_ratios_groups, combines list into a single matrix to make everything else easier
  tab = neg_OTU_ratios_groups(samps, negs) #gets the list from the previous function
  neg.OTUs = lapply(tab, function(x) { dimnames(x)[1] } ) #makes a list of lists, each list corresponds to one extraction date and contains all the OTUs found in both the negative controls and samples for that date
  neg.OTUs = lapply(neg.OTUs, '[[', 1) #changes the sublists to vectors
  all.neg.OTUs = sort(unique(unlist(neg.OTUs))) #makes a single vector containing all unique OTU IDs that were found to occur in both a sample and its corresponding negative control at least once
  new.data = matrix(nrow = length(all.neg.OTUs), ncol = length(sample_data(samps)$X.SampleID),
                    dimnames = list(all.neg.OTUs, (sample_data(samps)$X.SampleID))) #makes an empty matrix with one row for each OTU and one column for each actual sample (no negative controls)
  for (i in 1:length(tab)) { #starts a loop with one cycle for each extraction date
    rnames = rownames(tab[[i]]) #generates vector of OTU IDs from matrix for current cycle's date
    cnames = colnames(tab[[i]]) #generates vector of sample IDs from matrix for current cycle's date
    #The following two loops might not actually be necessary, but this is how I wrote it several months ago when I was still figuring stuff out
    for (j in 1:nrow(tab[[i]])) { #starts a sub-loop, one cycle for each OTU ID from current loop's date
      for (k in 1:ncol(tab[[i]])) { #starts a sub-sub-loop, one cycle for each sample ID from current loop's date
        new.data[rnames[j], cnames[k]] = tab[[i]][j, k] #Fills cells in matrix with values from the original list
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
