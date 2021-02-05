#_______________________________________________________________________
#
# Pipeline to analyze metaproteomics of MICROBIAL GUT
# Data are generated from High-throughput Liquid Chromatography 
# Mass Spectrometry instruments 
# 
#_______________________________________________________________________

# ------------- Set working directory --------------


setwd("~/metaprotr")


# ------------- Install the required packages -------------


# Before to start, install the following packages :
# ad4
# dendextend
# ggforce
# ggrepel
# reshape2
# stringr
# tidyverse

# ex. 
# install.packages("ade4")
# install.packages("dendextend")
# install.packages("ggforce")
# install.packages("ggrepel")
# install.packages("reshape2")
# install.packages("stringr")
# install.packages("tidyverse")


# ------------- Load required libraries -------------


library("metaprotr")
library("ade4")
library("dendextend")
library("dplyr")
library("ggforce")
library("ggrepel")
library("reshape2")
library("stringr")
library("tidyverse")

# The following warnings could appear, you can continue the analysis
# ── Conflicts ────────────────────────────────────────────────── tidyverse_conflicts() ──
# x dplyr::filter() masks stats::filter()
# x dplyr::lag()    masks stats::lag()


# ------------- Load Mass Spectrometry data and database with taxonomical levels -------------


# DB containing the taxonomic classification for each protein
# The Integrated non-redundant Gene Catalog 9.9 with taxonomic annotation can be downloaded from
# https://zenodo.org/record/3997093#.X0UYI6Zb_mE

meta99_full_taxo <- read.csv2("full_taxonomy_MetaHIT99.tsv", header = T, sep = "\t")

# Characters indicating the localisation (including directories) of files
# check that the items from the column 'msrunfile' (in metadata) are identical to the samples names in 'peptides_file'

proteins_file <- "data_directory/protein_list_test_condor_22janv.txt"
peptides_file <- "data_directory/peptide_counting_test_condor_22janv.txt"
metadata_file <- "data_directory/metadata.csv"

# Combine the files containing proteins names, peptides abundance in spectral counts and metadata

metaproteome <- load_protspeps(proteins_file, peptides_file, metadata_file)
dim(metaproteome$peptides_proteins) # verification


# ------------- Integrate taxonomic database into metaproteome_object -------------


metaproteome <- add_taxonomy(metaproteome, meta99_full_taxo) 
head(metaproteome$peptides_proteins) # verification


# ------------- Get spectral counts from specific peptides  -------------


# Spectral counts can be arranged by peptides, subgroups or groups
# SC_peps <- getsc_specific(metaproteome, "sc_specific_peptides") 
# SC_groups <- getsc_specific(metaproteome, "sc_groups")

# For this pipepline the analysis will be carried out from subgroups, also referred as "metaproteins"

SC_subgroups <- getsc_specific(metaproteome, "sc_subgroups") 


# ------------- Inspection of spectral data  -------------


# 1. The function below allows to observe the number of elements (subgroups)
# shared by the samples in the experiment

inspect_sample_elements(SC_subgroups)

# 2. Portrays the abundance intensities of the items (peptides, subgroups, 
# groups or taxonomic elements) in the conditions or samples
# To define the second argument of the function, you can inspect the 
# names of columns from metadata

colnames(SC_subgroups$metadata)
# [1] "SC_name"  "msrunfile"  "SampleID"  "Order_injection"  "Condition"

plot_intensities(SC_subgroups, "SC_name", "Spectral counts of subgroups")
plot_intensities(SC_subgroups, "Condition", "Spectral counts of subgroups")

# 3. The commands below show the clustering of samples after a Principal Components Analysis (PCA). 
# This allows to observe the overall grouping of the samples. 
# To define the second argument you can inspect the names of columns from metadata

colnames(SC_subgroups$metadata)
# [1] "SC_name"  "msrunfile"  "SampleID"  "Order_injection"  "Condition" 

plot_pca(SC_subgroups, "Condition", c(1, 2))
plot_pca(SC_subgroups, "Condition", c(2, 3))
# In this example, we observe that the 3 first components explain the major part of the variance


# ------------- Non supervised Clustering -------------


# To visualize the similarities between the different samples, an unsupervised clustering analysis can be informative
# First, you need to set a the name of ONE column from metadata
# The colors of the dendogram will be set in function of the different levels in the slected column  

str(SC_subgroups$metadata)

# 'data.frame':	11 obs. of  4 variables:
# $ SC_name        : chr  
# $ msrunfile      : chr  
# $ SampleID       : chr  
# $ Condition      : chr  

# The function will ask for the colors for a given sample or condition
# Provide the colors in the suitable format:
# ex. red, yellow, darkorange, blue
# ex. #440154FF, #FDE725FF, #404788FF, #55C667FF, #B8DE29FF

plot_dendocluster(SC_subgroups, "Condition", "graph_dendogram_by_groups")


# ------------- Reformating  -------------


# After data inspection, if some samples, conditions or elements (peptide, 
# proteins, subgroups, groups, etc.) should be removed there are some tools
#for reformating data

?remove_element
?select_element
?filter_text

# For instance, this metaproteome analysis requires to separate 
# human and microbial proteins.

SC_sgp_human <- filter_text(SC_subgroups, "Description", "HUMAN", "keep")
SC_sgp_microbial <- filter_text(SC_subgroups, "Description", "HUMAN", "discard")


# ------------- Comparison between two samples or conditions -------------


# The following lines allow to compare globally two conditions or two samples 

plot_intensities_ratio(SC_sgp_microbial, "Condition", c("CTRL", "S"))

# To verify that replicates are similar in a given condition

plot_intensities_ratio(SC_sgp_microbial, "SampleID", c("Q1_batch1_prot", "Q2_batch1_prot"))


# ------------- Venn diagrams  -------------

# Option 1

venn_Q1_Q2 <- plot_venn(SC_sgp_microbial, "SampleID", c("Q1_batch1_prot", "Q2_batch1_prot"))

# Export of the elements of each part of the venn diagram on the working directory
export_vennlists(venn_Q1_Q2)

# Option 2

# Venn plots can also be performed by comparing three conditions

SC_sgp_microbial$metadata$SampleID
# [1] "STD_env"            "CTRL"               "Q1+FW1_batch1_prot" "Q2+FW2_batch1_prot"
# [5] "Q3+FW3_batch1_prot" "Q1_batch1_prot"     "Q2_batch1_prot"     "Q3_batch1_prot"    
# [9] "FW1_batch1_prot"    "FW2_batch1_prot"    "FW3_batch1_prot" 

# Set the three conditions
list_conditions <- c("Q1_batch1_prot", "Q2_batch1_prot", "Q3_batch1_prot")
venn_Q1_Q2_Q3 <- plot_venn(SC_sgp_microbial, "SampleID", list_conditions)

# Export of "venn_object" lists in a subdirectory
dir.create("list_venn_Q1_Q2_Q3")
export_vennlists(venn_Q1_Q2_Q3, "list_venn_Q1_Q2_Q3")


# ------------- Information about the taxonomy levels in all samples of an "spectral_count_object" -------------


# ~7 min depending on the experiment size
plot_fulltaxo(SC_sgp_microbial) 


# ------------- Calculates spectral counts (MICROBIAL origin for instance) per taxonomic level -------------


# All in ~ 16 mins
SCsgp_superkingdom <- crumble_taxonomy(SC_sgp_microbial, "superkingdom")
SCsgp_phylum <- crumble_taxonomy(SC_sgp_microbial, "phylum")
SCsgp_class <- crumble_taxonomy(SC_sgp_microbial, "class")
SCsgp_order <- crumble_taxonomy(SC_sgp_microbial, "order")
SCsgp_family <- crumble_taxonomy(SC_sgp_microbial, "family")
SCsgp_genus <- crumble_taxonomy(SC_sgp_microbial, "genus")
SCsgp_species <- crumble_taxonomy(SC_sgp_microbial, "species")


# ------------- Plot stacked bars of spectral_count_objects with defined taxonomy  ------------- 


# The lines below represent the distribution of spectral counts per taxonomic level (ex. SCsgp_superkingdom, SCsgp_phylum, SCsgp_class, etc.)

plot_stackedtaxo(SCsgp_family, "SC_name", "percent")

# Third argument is optional and can be adjusted to improve visual 
# representaion by increasing/decreasing the elements that pass the "filter_percent" argument

plot_stackedtaxo(SCsgp_family, "SC_name", "percent", 2)

# Other columns from metadata can also be selected (ex. conditions, order of injection, etc.)

plot_stackedtaxo(SCsgp_family, "Condition", "percent", 2)

# Bar plot can be also represented as spectral counts by setting the second argument as "numbers"

plot_stackedtaxo(SCsgp_family, "SampleID", "numbers", 2)


# ------------- Plot pie a chart of a sample or condition from spectral_count_objects with defined taxonomy  ------------- 


plot_pietaxo(SCsgp_genus, "Condition", "CTRL")

# Third argument is optional and can be adjusted to improve visual representaion
# by increasing/decreasing the elements that pass the "filter_percent" argument

plot_pietaxo(SCsgp_species, "Condition", "STD", 0.8)


# ------------- Principal Components Analysis (PCA) from taxonomic entities -------------


plot_pca(SCsgp_species, "Condition", c(1, 2))
plot_pca(SCsgp_genus, "Condition", c(1, 2))
plot_pca(SCsgp_family, "Condition", c(1, 2))


# ------------- Differential taxonomic elements (over / under represented) ------------- 


# To display the most divergent elements between two conditions or samples 
# The differential elements are selected by a log2(condition / reference) > 3
# this is equivalent to a 8 fold-change between two conditions

identify_differences(SCsgp_species, "Condition", c("CTRL", "S"))

identify_differences(SCsgp_genus, "Condition", c("CTRL", "S"))

identify_differences(SCsgp_family, "SampleID", c("Q1_batch1_prot", "Q2_batch1_prot"))

# If required, the third argument can be adjusted to identify elements with lower log2(fold change)

identify_differences(SCsgp_genus, "Condition", c("CTRL", "S"), 1.5)


# ------------- Functional Annotation from KEGG database -------------


# Import the KEGG annotations and the MetaHit99 database
# The IGC database with KEGG functionnal annotation can be downloaded from
# https://zenodo.org/record/3997093#.X0UYI6Zb_mE

kegg_db <- read.csv2("KEGG89_IGC_hs99.table", header = T, sep = "\t")

# The following objects were loaded at the beginning of the pipeline, they are required for functional annotation

meta99_full_taxo <- read.csv2("full_taxonomy_MetaHIT99.tsv", header = T, sep = "\t")
proteins_file <- "data_directory/protein_list_test_condor_22janv.txt"
peptides_file <- "data_directory/peptide_counting_test_condor_22janv.txt"
metadata_file <- "data_directory/metadata.csv"
metaproteome <- load_protspeps(proteins_file, peptides_file, metadata_file)

# Assignment of functional annotation per taxonomic level

SCsgp_species_annot <- add_kegg(
  SCsgp_species, 
  kegg_db, meta99_full_taxo, metaproteome, proteins_file, peptides_file, 
  text_to_filter = "HUMAN"
)

SCsgp_genus_annot <- add_kegg(
  SCsgp_genus, 
  kegg_db, meta99_full_taxo, metaproteome, proteins_file, peptides_file, 
  text_to_filter = "HUMAN"
)

SCsgp_family_annot <- add_kegg(
  SCsgp_family, 
  kegg_db, meta99_full_taxo, metaproteome, proteins_file, peptides_file, 
  text_to_filter = "HUMAN"
)

SCsgp_phylum_annot <- add_kegg(
  SCsgp_phylum, 
  kegg_db, meta99_full_taxo, metaproteome, proteins_file, peptides_file, 
  text_to_filter = "HUMAN"
)

SCsgp_superkingdom_annot <- add_kegg(
  SCsgp_superkingdom, 
  kegg_db, meta99_full_taxo, metaproteome, proteins_file, peptides_file, 
  text_to_filter = "HUMAN"
)


# ------------- Export of KO terms to draw metabolic pathways in iPATH3 (https://pathways.embl.de/) -------------


export_ipath3(SCsgp_phylum_annot, "all", "SampleID", "CTRL", "#840AA3")

taxonomic_levels <- c("Bacteroidetes", "Actinobacteria", "Proteobacteria")
export_ipath3(SCsgp_phylum_annot, "selection", "Condition", "S", "#28c1df", taxonomic_levels)

taxonomic_levels <- c("Eubacteriaceae")
export_ipath3(SCsgp_family_annot, "selection", "Condition", "S", "#097737", taxonomic_levels)


# ------------- Export spectral counts expressed in percentages, considering the sample, for multi-omics integration -------------


export_robject(SCsgp_species_annot, "spectral_percent", "RDATA")

export_robject(SCsgp_family_annot, "spectral_percent", "RDATA")

export_robject(SCsgp_phylum_annot, "spectral_percent", "RDATA")

