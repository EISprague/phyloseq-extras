Functions and scripts that I wrote during my M.S. that made microbiome data analysis easier. Compatible with phyloseq version 1.24.0


**GlobalPatterns raw abundance graph script.R**: A workflow example of my otu_abundance() function (see "functions/otu_abundance function.R" below) using the GlobalPatterns data set that comes with the phyloseq package. This particular example uses the raw abundance data.

**GlobalPatterns raw graphed.png**: Graphical output of "GlobalPatterns raw abundance graph script.R"

**GlobalPatterns relative abundance graph script.R**: Another workflow example of my otu_abundance() function, but in this verison the data was converted to relative abundance first.

**GlobalPatterns relative graphed.png**: Graphical output of "GlobalPatterns relative abundance graph script.R"

**negative control functions.R**: Contains four functions; the first removes OTUs identified as contaminants in a single step and outputs the new phyloseq object. If the user wishes to examine the abundance and/or identity of the contaminant OTUs before removing them, use one of the other functions instead. The function assumes there is a column named "Control_type" in the sample_data() slot of the physeq input which designates whether or not samples are controls. See annotations under "remove_contaminants" function for other assumptions about each argument (in particular the assumption about the ratio argument).

**negative control removal example.R**: A workflow example of my functions for removing OTUs identified as potential contaminants found in negative controls (see "negative control functions.R" above). Uses the GlobalPatterns data set, and for demonstration purposes treats the mock community samples like negative controls.

**otu_abundance function.R**: Outputs dataframe similar to the plot_bar() function that comes with phyloseq (using the justDF = T option). However, this function gives the additional option to set a cutoff for the most abundant taxa in the given taxonomic rank (e.g. "Phylum"); any taxa that are below this cutoff will be relabeled as "Other" and grouped together as the same color in the bar graph. This makes it easier for the user to control how the graph looks and how many colors are needed for the stacked bars. This function assumes the user has already converted abundance to relative abundance if desired.
