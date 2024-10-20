# Microbiome_NCWS_Bread_analysis
Gut microbiome analysis in non-celiac wheat-sensitive patients and effects of different bread types

# Data comes from the in vitro study, where the effect of different bread types on fecal microbiome composition in healthy (HC) and Non-Celiac Wheat-Sensitive Patients (NCWS) was investigated.
# Brief description of the data:
- 138 Samples
- Taxonomy data: 1830 taxa with ranks from kingdom to species
- Many sample data variables:
	-Patient types: Healthy control individuals, Non-Celiac Wheat-Sensitive Patients (NCWS).
	-Time points: T0 (at zero time), T5 (after 5 hours of adding grane and ferment)
	-Grain types: spelt, wheat, emmer, and medium.
	-Ferment types: yeast, sourdough, and medium.
	-Age
	-Sex

# The data is in a phyloseq object from the phyloseq package in R containing:
- Amplicon Sequence Variants (ASVs), in so called otu_table
- metadata in sample_data
- taxonomy data in tax_table
- phylogenetic tree in phy_tree

# Project aims: 

- Finding if there is an association between gut microbiome and the Non-Celiac Wheat-Sensitive(NCWS) patients
- Finding if there is an association between the gut microbiome and different bread types in NCWS patients

# Analysis was done using R studio (version 2024.09.0+375) and R (version 4.4.1)

# From the repository (Microbiome_NCWS_Bread_analysis) do the following:
# Download the file of the phyloseq object of the data named (WoW_project_finPS.rds) 
# Download the file of R.script to run the analysis named (Scient._program_project.R) 
# Put both files together in the same folder of your working directory.
# Run the R.script using the R studio (version 2024.09.0+375)
