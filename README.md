Functions and scripts that make microbiome data analysis using phyloseq easier. Compatible with phyloseq version 1.24.0  


### Functions to help with use of phyloseq (functions folder)  
**negative control functions.R**: Contains four functions; the first removes OTUs identified as contaminants in a single step and outputs the new phyloseq object. If the user wishes to examine the abundance and/or identity of the contaminant OTUs before removing them, use one of the other functions instead. The function assumes there is a column named "Control_type" in the sample_data() slot of the physeq input which designates whether or not samples are controls. See annotations under "remove_contaminants" function for other assumptions about each argument (in particular the assumption about the ratio argument).

**taxa_abundance function.R**: Outputs dataframe similar to the plot_bar() function that comes with phyloseq (using the justDF = T option). However, this function gives the additional option to set a cutoff for the most abundant taxa in the given taxonomic rank (e.g. "Phylum"); any taxa that are below this cutoff will be relabeled as "Other" and grouped together as the same color in the bar graph. This makes it easier for the user to control how the graph looks and how many colors are needed for the stacked bars. This function assumes the user has already converted abundance to relative abundance if desired.  

**phyloseq_NA_tax function.R**: Replaces NAs in phyloseq object taxonomy table with the name of lowest identified taxonomic level. Useful to prep for taxa_abundance function.

**get_taxonomy function.R**: Simple function to access full taxonomy from tax_table for a vector of OTU IDs.  

### Example graphs and workflows that use the above functions (demos folder)  
**negative control removal example.R**: A workflow example of my functions for removing OTUs identified as potential contaminants found in negative controls (see "negative control functions.R" above). Uses the GlobalPatterns data set, and for demonstration purposes treats the mock community samples like negative controls.  

**taxa_abundance_demo.Rmd**: A workflow example of how to use the taxa_abundance and phyloseq_NA_tax functions. See the HTML output here: http://htmlpreview.github.io/?https://github.com/EISprague/phyloseq-extras/blob/master/demos/taxa_abundance_demo.html
