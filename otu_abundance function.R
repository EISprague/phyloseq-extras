otu_abundance = function(physeq, rank, cutoff) { 
  ##physeq is the phyloseq object; it should already be converted to relative abundance if desired
  #rank is which taxonomy rank you want to glom and examine such as "Phylum"
  #cutoff determines which taxa keep their names and which are renamed "Other"
  #If cutoff < 1, all taxa with an overall relative abundance greater than cutoff will keep their names. ex: if cutoff = 0.01, all taxa that make up 1% or more of all sequences across all samples will keep their names. All other taxa will be labeled "Other"
  #If cutoff > 1, the most abundant taxa will keep their names. ex: if cutoff = 12, the 12 most abundant taxa across all samples will keep their names, and all other taxa will be renamed "Other". If you want to decide the number of colors to use ahead of time, this is the route to choose
  #If cutoff = 0, no taxa will be renamed as "Other" and a useable data frame is returned
  require(phyloseq)
  require(magrittr)
  require(plyr)
  require(dplyr)
  
  physeq.glom = tax_glom(physeq, taxrank = rank) #Merges taxa at the given rank
  total = taxa_sums(physeq.glom) #sums abundance for each phylum across all samples; named with OTU tags, not actual phyla names at this step
  overall = data.frame(OTU = names(total), raw.count = as.numeric(total))
  overall$rel.abund = overall$raw.count/sum(overall$raw.count) #calculates the relative abundance of each taxon across samples
  physeq.melt = psmelt(physeq.glom) #turns phyloseq object into dataframe
  physeq.melt[, rank] = as.character(physeq.melt[, rank]) #makes rank into character vector instead of factor
  merged = merge(physeq.melt, overall, by = "OTU") #adds the overall calculations to the dataframe
  
  if (cutoff < 1) { #this statement is used for percentage-based cutoffs
    merged[merged$rel.abund < cutoff, rank] = "Other" #labels taxa that didn't make the cutoff as "Other"
  }
  
  else { #this statement is used for rankings-based cutoffs
    top = data.frame(rel.abund = merged$rel.abund, rank = merged[rank]) %>% #orders taxa by abundance
      distinct() %>% arrange(desc(rel.abund)) %>% top_n(cutoff, rel.abund)
    merged[!(merged[, rank] %in% top[, rank]), rank] = "Other" #labels taxa that didn't make the cutoff as "Other"
    
  }
  
  all_ranks = colnames(physeq.glom@tax_table@.Data)[1:which(colnames(physeq.glom@tax_table@.Data) == rank)] #used later
  colnames(merged)[which(colnames(merged) == rank)] = "Glommed_rank" #renames the rank column for the following steps
  others.combined = aggregate(cbind(Abundance, rel.abund, raw.count) ~ X.SampleID + Glommed_rank, #sums abundance columns for each combination of taxon and sample
                              FUN = sum, data = merged)
  merged$Abundance = NULL #deleting columns unecessary to the end user (contain numbers for overall abundance instead of by sample)
  merged$raw.count = NULL
  merged$rel.abund = NULL
  merged$OTU = NULL
  merged = distinct(merged[, !(colnames(merged) %in% all_ranks)]) #removes repeated rows
  
  out = merge(others.combined, merged, by = c("X.SampleID", "Glommed_rank"), all = F) #merges sample_data from original physeq input with newly calculated abundance data and renamed taxa
  colnames(out)[which(colnames(out) == "Glommed_rank")] = rank #renames column with rank again
  out$raw.count = NULL #unecessary to end user
  
  return(distinct(out))
}
