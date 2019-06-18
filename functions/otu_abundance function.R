function(physeq, rank = NULL, cutoff, keep_ranks = T) { 
  #physeq is the phyloseq object
  #rank is the taxonomic rank you want to examine such as "Phylum". Use rank = NULL if your physeq object is already the output of the tax_glom or tip_glom functions to the level you wish to examine. Function assumes you wish to view taxa that have been agglomerated to at least the Species or other lowest taxonomic level, not individual OTUs or sequence variants.
  #cutoff (if < 1) is the minimum percentage of the relative abundance of the taxa you wish to visualize. If > 1, it is the x most abundant number of taxa you wish to visualize. All taxa with either less than the minimum percent or not in the top x most abundant will be grouped together and returned as "Other". If = 0, it returns all taxa.
  #keep_ranks: if TRUE, all taxonomic levels above the glommed level will be returned in the final data frame. If FALSE they will be omitted
  require(phyloseq)
  require(magrittr)
  require(plyr)
  require(dplyr)
  
  if (is.null(rank)) {
    physeq.new = physeq
    rank = colnames(physeq.new@tax_table@.Data)[ncol(physeq.new@tax_table@.Data)]
  } else {
    physeq.new = tax_glom(physeq, taxrank = rank) #Merges taxa with the same taxonomy at the given rank
  }
  total = taxa_sums(physeq.new) #sums raw abundance for each phylum across all samples; named with OTU tags, not actual phyla names
  overall = data.frame(OTU = names(total), raw.count = as.numeric(total))
  overall$rel.abund = overall$raw.count/sum(overall$raw.count) #calculates the overall relative abundance of each taxon to determine if it is above or below the cutoff
  physeq.melt = psmelt(physeq.new) #turns phyloseq object into dataframe
  physeq.melt = physeq.melt[, -4]
  physeq.melt[, rank] = as.character(physeq.melt[, rank])
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
  }
  
  out = merge(others.combined, merged, by = c("Sample", "OTU"), all = F)
  out$raw.count = NULL
  out$OTU = NULL
  
  return(distinct(out))
}
