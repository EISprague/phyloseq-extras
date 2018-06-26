#Your working directory should already have the rep_set.tre, mapping file, and the json version of your .biom file, along with the phyloseq extras.R file

source("phyloseq extras.R")
library(phyloseq)
library(magrittr)

tree = read_tree("rep_set.tre")
otu_table = import_biom("otu_table_mc2_w_tax_no_pynast_failures_ANC_gone_json.biom") #
colnames(tax_table(otu_table)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

###Replacing taxonomy NAs to prevent removal of OTUs because some phyloseq functions are picky
otu_table@tax_table@.Data[is.na(otu_table@tax_table@.Data[, "Kingdom"]), "Kingdom"] = "Unassigned_domain"
otu_table@tax_table@.Data[is.na(otu_table@tax_table@.Data[, "Phylum"]), "Phylum"] = "p__Unassigned"
otu_table@tax_table@.Data[is.na(otu_table@tax_table@.Data[, "Class"]), "Class"] = "c__Unassigned"
otu_table@tax_table@.Data[is.na(otu_table@tax_table@.Data[, "Order"]), "Order"] = "o__Unassigned"
otu_table@tax_table@.Data[is.na(otu_table@tax_table@.Data[, "Family"]), "Family"] = "f__Unassigned"
otu_table@tax_table@.Data[is.na(otu_table@tax_table@.Data[, "Genus"]), "Genus"] = "g__Unassigned"
otu_table@tax_table@.Data[is.na(otu_table@tax_table@.Data[, "Species"]), "Species"] = "s__Unassigned"
otu_table@tax_table@.Data[which(otu_table@tax_table@.Data[, "Kingdom"] == "k__"), "Kingdom"] = "Unassigned_domain"
otu_table@tax_table@.Data[which(otu_table@tax_table@.Data[, "Kingdom"] == "Unassigned"), "Kingdom"] = "Unassigned_domain"
otu_table@tax_table@.Data[which(otu_table@tax_table@.Data[, "Phylum"] == "p__"), "Phylum"] = "p__Unassigned"
otu_table@tax_table@.Data[which(otu_table@tax_table@.Data[, "Class"] == "c__"), "Class"] = "c__Unassigned"
otu_table@tax_table@.Data[which(otu_table@tax_table@.Data[, "Order"] == "o__"), "Order"] = "o__Unassigned"
otu_table@tax_table@.Data[which(otu_table@tax_table@.Data[, "Family"] == "f__"), "Family"] = "f__Unassigned"
otu_table@tax_table@.Data[which(otu_table@tax_table@.Data[, "Genus"] == "g__"), "Genus"] = "g__Unassigned"
otu_table@tax_table@.Data[which(otu_table@tax_table@.Data[, "Species"] == "s__"), "Species"] = "s__Unassigned"

###Use this part to check that the otu_table and tre imports match perfectly if you're feeling suspicious
#noverlap.t = sort(tree$tip.label[(which(!(tree$tip.label %in% rownames(otu_table@otu_table@.Data))))])
#noverlap.o = sort(rownames(otu_table@otu_table@.Data)[(which(!(rownames(otu_table@otu_table@.Data) %in% tree$tip.label)))])
#noverlap.samps = sort(rownames(samps.final@otu_table@.Data)[(which(!(rownames(samps.final@otu_table@.Data) %in% tree$tip.label)))])

#rm(noverlap.o,noverlap.t) ###Remove once import is confirmed

map = import_qiime_sample_data("lakedata_map.txt")
map$InputFastaFileName = NULL #No longer necessary
physeq.mc = merge_phyloseq(otu_table, map, tree) %>% #version of the phyloseq object with mitochondria and chloroplasts included
  subset_samples(X.SampleID != "May.N.1") %>% #Samples that are probably contaminated
  subset_samples(X.SampleID != "Feb.W.4") %>%
  subset_samples(X.SampleID != "Nov.W.0")
physeq = subset_taxa(physeq.mc, Class != "c__Chloroplast") %>% #Removing chloroplasts, mitochondria, and genera that are probably contaminants
  subset_taxa(Family != "f__mitochondria") %>%
  subset_taxa(Genus != "g__Staphylococcus") %>%
  subset_taxa(Genus != "g__Streptococcus") %>%
  subset_taxa(Genus != "g__Micrococcus") %>%
  subset_taxa(Genus != "g__Propionibacterium") %>%
  subset_taxa(Genus != "g__Ralstonia") %>%
  subset_taxa(Kingdom != "Unassigned_domain") %>%
  prune_taxa(taxa_sums(.) > 0, .)
