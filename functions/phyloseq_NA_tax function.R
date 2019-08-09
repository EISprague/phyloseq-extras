phyloseq_NA_tax = function(physeq) {
  ###This function assumes that the taxonomy table of physeq has exactly seven columns with the classic Linnaean taxonomic level names as the column names.
  require(phyloseq)
  out = physeq
  
  tax_table(out)[is.na(tax_table(out)[, "Kingdom"]), ] = "Unassigned_kingdom"
  tax_table(out)[is.na(tax_table(out)[, "Phylum"]), 2:7] = 
    paste0("Kingdom_", tax_table(out)[is.na(tax_table(out)[, "Phylum"]), "Kingdom"])
  tax_table(out)[is.na(tax_table(out)[, "Class"]), 3:7] = 
    paste0("Phylum_", tax_table(out)[is.na(tax_table(out)[, "Class"]), "Phylum"])
  tax_table(out)[is.na(tax_table(out)[, "Order"]), 4:7] = 
    paste0("Class_", tax_table(out)[is.na(tax_table(out)[, "Order"]), "Class"])
  tax_table(out)[is.na(tax_table(out)[, "Family"]), 5:7] = 
    paste0("Order_", tax_table(out)[is.na(tax_table(out)[, "Family"]), "Order"])
  tax_table(out)[is.na(tax_table(out)[, "Genus"]), 6:7] = 
    paste0("Family_", tax_table(out)[is.na(tax_table(out)[, "Genus"]), "Family"])
  tax_table(out)[is.na(tax_table(out)[, "Species"]), 7] = 
    paste0("Genus_", tax_table(out)[is.na(tax_table(out)[, "Species"]), "Genus"])

  return(out)
}
