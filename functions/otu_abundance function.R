otu_abundance = function(physeq, rank = NULL, cutoff, keep_ranks = T) { 
  #physeq is the phyloseq object
  #rank is the taxonomic rank you want to examine such as "Phylum"; this will "agglomerate" the data to the given rank. Use rank = NULL if your physeq object is already the output of the tax_glom or tip_glom functions to the taxonomic rank you wish to examine. If you use rank = NULL and your input has not been agglomerated (ie, still has individual OTUs or sequence variants), the individual OTUs/SVs will be ranked and returned.
  #cutoff (if < 1) is the minimum percentage of the relative abundance of the taxa you wish to visualize. If > 1, it is the x most abundant number of taxa you wish to visualize. All taxa with either less than the minimum percent or not in the top x most abundant will be grouped together and returned as "Other". If = 0, it returns all taxa.
  #keep_ranks: if TRUE, all taxonomic levels above the glommed level will be returned in the final data frame. If FALSE they will be omitted
  require(phyloseq)
  require(magrittr)
  require(plyr)
  require(dplyr)
  
  if (is.null(rank)) {#If you use this option and haven't previously agglomerated your input physeq object, the final column name in your taxonomy table will be treated as the rank for all subsequent steps
    physeq.new = physeq
    rank = colnames(physeq.new@tax_table@.Data)[ncol(physeq.new@tax_table@.Data)]
  } else {
    physeq.new = tax_glom(physeq, taxrank = rank) #Merges taxa with the same taxonomy at the given rank
  }
  nranks = ncol(physeq.new@tax_table@.Data[, 1:which(colnames(physeq.new@tax_table@.Data) == rank)])
  first_rank = colnames(physeq.new@tax_table@.Data)[1]
  total = taxa_sums(physeq.new) #sums raw abundance for each phylum across all samples; named with OTU/SV tags, not taxonomic rank names at this point
  overall = data.frame(OTU = names(total), raw.count = as.numeric(total))
  overall$rel.abund = overall$raw.count/sum(overall$raw.count) #calculates the relative abundance across all samples of each taxon to determine if it is above or below the cutoff
  physeq.melt = psmelt(physeq.new) #turns phyloseq object into dataframe
  #physeq.melt = physeq.melt[, -4]
  rank_end = which(colnames(physeq.melt) == rank)
  rank_start = which(colnames(physeq.melt) == first_rank)
  physeq.melt[, rank_start:rank_end] = apply(physeq.melt[, rank_start:rank_end], 2, as.character)
  #colnames(physeq.melt)[1] = "OTU"
  merged = merge(physeq.melt, overall, by = "OTU") #adds the cutoff calculations to the dataframe
  
  if (cutoff < 1) { 
    merged[merged$rel.abund < cutoff, rank] = "Other" #labels taxa that didn't make the cutoff as "Other"
  } else {
    top = data.frame(rel.abund = merged$rel.abund, rank = merged[rank]) %>%
      distinct() %>%
      arrange(desc(rel.abund)) %>%
      top_n(cutoff, rel.abund)
    merged[!(merged[, rank] %in% top[, rank]), rank] = "Other"
  }
  
  #colnames(merged)[which(colnames(merged) == rank)] = "Glommed_rank"
  others.combined = aggregate(cbind(Abundance, rel.abund, raw.count) ~ Sample + OTU,
                              FUN = sum, data = merged)
  merged$Abundance = NULL
  merged$raw.count = NULL
  merged$rel.abund = NULL
  if (keep_ranks) {
    merged = distinct(merged)
  } else {
    all_ranks = colnames(physeq.new@tax_table@.Data)[1:which(colnames(physeq.new@tax_table@.Data) == rank)]
    merged = distinct(merged[, !(colnames(merged) %in% all_ranks)])
    nranks = 1
  }
  
  out = merge(others.combined, merged, by = c("Sample", "OTU"), all = F)
  out$raw.count = NULL
  out$OTU = NULL
  out$rel.abund = NULL
  rank_end = which(colnames(out) == rank)
  rank_start = which(colnames(out) == first_rank)
  out[out[rank] == "Other", rank_start:rank_end] = "Other"
  #Now sum "Other" for each individual sample
  out2 = aggregate(Abundance ~ ., data = out, sum)
  out2 = distinct(out2)
  return(out2)
}
