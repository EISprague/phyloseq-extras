get_taxonomy = function(OTUlist, physeq) {
  tax = physeq@tax_table@.Data
  out = tax[rownames(tax) %in% OTUlist, ]
  return(out)
}
