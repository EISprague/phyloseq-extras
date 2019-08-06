taxa_abundance = function(physeq, rank = NULL, cutoff) { 

  require(phyloseq)
  require(magrittr)
  require(plyr)
  require(dplyr)
  
  if (is.null(rank)) {
    physeq.new = physeq
    rank = colnames(tax_table(physeq.new))[ncol(tax_table(physeq.new))]
  } else {
    physeq.new = tax_glom(physeq, taxrank = rank)
  }
  nranks = which(colnames(tax_table(physeq.new)) == rank)
  first_rank = colnames(tax_table(physeq.new))[1]
  total_across = taxa_sums(physeq.new)
  across_samps = data.frame(OTU = names(total_across), raw.count = as.numeric(total_across))
  across_samps$rel.abund = across_samps$raw.count/sum(across_samps$raw.count)
  
  physeq.melt = psmelt(physeq.new)
  rank_end = which(colnames(physeq.melt) == rank)
  rank_start = which(colnames(physeq.melt) == first_rank)
  physeq.melt[, rank_start:rank_end] = apply(physeq.melt[, rank_start:rank_end], 2, as.character)
  
  merged = merge(physeq.melt, across_samps, by = "OTU")
  rank_end = which(colnames(merged) == rank)
  rank_start = which(colnames(merged) == first_rank)
  
  if (cutoff < 1) { 
    merged[merged$rel.abund < cutoff, rank] = "Other"
  } else {
    top = data.frame(rel.abund = merged$rel.abund, rank = merged[rank]) %>%
      distinct() %>%
      arrange(desc(rel.abund)) %>%
      top_n(cutoff, rel.abund)
    merged[!(merged[, rank] %in% top[, rank]), rank] = "Other"
  }
  merged[merged[rank] == "Other", rank_start:rank_end] = "Other"
  
  merged$raw.count = NULL
  merged$rel.abund = NULL
  merged$OTU = NULL
  out = aggregate(Abundance ~ ., data = merged, sum)
  return(out)

}
