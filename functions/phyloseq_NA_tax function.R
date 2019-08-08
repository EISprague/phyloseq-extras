phyloseq_NA_tax = function(physeq) {
  
  out = physeq
  
  for(row in 1:nrow(tax_table(out))) {
    if(is.na(tax_table(out))[row, "Kingdom"]) {
      tax_table(out)[row, 1:7] = "Unassigned_kingdom"
    } else if(is.na(tax_table(out)[row, "Phylum"])) {
      tax_table(out)[row, 2:7] = paste0("Kingdom_", tax_table(out)[row, "Kingdom"])
    } else if(is.na(tax_table(out)[row, "Class"])) {
      tax_table(out)[row, 3:7] = paste0("Phylum_", tax_table(out)[row, "Phylum"])
    } else if(is.na(tax_table(out)[row, "Order"])) {
      tax_table(out)[row, 4:7] = paste0("Class_", tax_table(out)[row, "Class"])
    } else if(is.na(tax_table(out)[row, "Family"])) {
      tax_table(out)[row, 5:7] = paste0("Order_", tax_table(out)[row, "Order"])
    } else if(is.na(tax_table(out)[row, "Genus"])) {
      tax_table(out)[row, 6:7] = paste0("Family_", tax_table(out)[row, "Family"])
    } else if(is.na(tax_table(out)[row, "Species"])) {
      tax_table(out)[row, 7] = paste0("Genus_", tax_table(out)[row, "Genus"])
    }
  }
  return(out)
}
