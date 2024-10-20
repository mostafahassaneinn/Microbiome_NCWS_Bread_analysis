#=====================================================================================================#
# MSB1015_Scient_program_project.R                                                                    #
#	Gut microbiome analysis in non-celiac wheat-sensitive patients and effects of different bread types	#															                                       		  #
# Date: Oct 20, 2024											                                                            #
# Author: Mostafa Hassanein, Maastricht University                                                    #
#=====================================================================================================#

# Part 1) Data cleaning
## a) Removing duplication and missing values (NA) from the sample data
## b) Forming one variable of bread type from grane and ferment
## c) Check if Sample filtering is needed
## d) Inspecting/Fixing the phyloseq object
## e) Aggregation, filtration, imputation, transformation, and mean centring (will be done per each aim)

# Part 2) Aims:

## Aim 1) Investigating the potential association between the gut microbiome and Non-Celiac Wheat-Sensitive (NCWS) patients
## 1.a) Aggregation
## 1.b) Filtration, imputation, trasnformation, and mean centring
## 1.b.1) Filtration, imputation, transformation and mean centring for beta diversity
## 1.b.2) Filtration, imputation, transformation and mean centring for differential abundance 
## 1.c) Alpha diversity
## 1.d) Visualisation using Composition barplot
## 1.e) Statistics with alpha diversity
## 1.f) Dissimilarity and ordination
## 1.g) Statistical analysis using PERMANOVA
## 1.h) Differential Abundance testing

## Aim 2) Investigating the potential association between the gut microbiome and different bread types in Non-Celiac Wheat-Sensitive (NCWS) patients.

## 2.a) Aggregation
## 1.b) Filtration, imputation, trasnformation, and mean centring
## 1.b.1) Filtration, imputation, transformation and mean centring for beta diversity
## 1.b.2) Filtration, imputation, transformation and mean centring for differential abundance
## 2.c) Alpha diversity
## 2.d) Visualisation using Composition barplot
## 2.e) Statistics with alpha diversity
## 2.f) Dissimilarity and ordination
## 2.g) Statistical analysis using PERMANOVA
## 2.h) Differential Abundance testing


set.seed(1) # for reproducible stochastic processes

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"),update=FALSE)

# Packages necessary to run this script
package_list <- c("base",
                  "utils",
                  "stats",
                  "methods",
                  "grDevices", 
                  "datasets",
                  "dplyr",
                  "ggplot2",
                  "patchwork",
                  "corncob",
                  "graphics",
                  "tibble",
                  "rstudioapi",
                  "phyloseq",
                  "microbiome",
                  "ComplexHeatmap",
                  "ggtext",
                  "ggraph",
                  "DT",
                  "ggside")

# Loading each package in the list
for(package in package_list){
  if (!require(package,character.only = TRUE)) {
    # If the package is not yet installed, install it
    install.packages(package)
  }
  library(package,character.only = TRUE)
} 

# Check if the following package is already installed
"microViz" %in% rownames(installed.packages())

# if TRUE --> Comment the following four lines (no need to re-install)
# install.packages(
#   "microViz",
#   repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
# )

# Loading the package
library(microViz)

# Setting up your directory

wd <- getActiveDocumentContext()$path
setwd(dirname(wd))
########################################################################################################
# Read the RDS file:
ps_data <- readRDS("WoW_project_finPS.rds")

# Explore the phyloseq object:
ps_data
# This illustrate that the object has 1830 taxa, 138 samples, 10 variables

# Part of OTU_table to see how it lookslike:
otu_get(ps_data, taxa = 1:3, samples = 1:3)
# The table contain counts not proportions

# More exploration of the phyloseq object:
sample_names(ps_data) %>% head()
taxa_names(ps_data) %>% head()
sample_variables(ps_data)
rank_names(ps_data)
# This gives an idea of the naming of the samples, taxa, variables, ranks.
#########################################################################################################
## 1) Data Cleaning:

### 1.a) Removing duplication and missing values (NA) from the sample data in the phyloseq: _______________

# Explore the sample data:
samdat_tbl(ps_data)  
# There are 138 samples and 10 variables

# Remove duplicated/repeated samples(identifying "duplicates" via the interaction of multiple columns):
ps_non_dup <- ps_dedupe(
  ps = ps_data, verbose = TRUE,
  vars = c("Sample_Type", "time_point", "grane", "ferment", "patient_type","patient", "age", "seks", "DNA_con")
)
# No duplicates are detected

#Remove samples with missing values in sample_data: 
ps_non_dup <- ps_drop_incomplete(ps_non_dup, vars = NA, verbose = "max")
# Just in case there are different missing values in sample data, in this case the samples have not missing values

# Changing one variable name:
ps_non_dup <- ps_non_dup %>%
  ps_mutate(sex = seks)

# Extract sample metadata
sample_metadata <- sample_data(ps_non_dup)

# 1.b) Forming one variable of bread type from grane and ferment and removing not needed variables (pool, sample_type)
sample_metadata$bread_type <- paste(sample_metadata$ferment, sample_metadata$grane, sep = "_")
sample_metadata <- sample_metadata[,-c(1, 2, 5, 6, 10)]
ps_non_dup@sam_data <- sample_data(sample_metadata)

# 1.c) Check if Sample filtering is needed:

ps_non_dup %>%
  ps_mutate(reads = sample_sums(ps_non_dup)) %>%
  samdat_tbl() %>%
  ggplot(aes(x = reads)) +
  geom_freqpoly(bins = 100) +
  geom_rug(alpha = 0.5) +
  labs(x = "Number of classified reads", y = "Number of samples") +
  theme_bw()
# Most of the reads for the samples are around the average (normally distributed)
# No need to remove any sample according to this criteria

### 1.d) Inspecting/Fixing the phyloseq object: ____________________________________________________________________
phyloseq_validate(ps_non_dup)
# This function returns "NA" it indicates that there are inconsistencies or missing values within the phyloseq object
# This illustrates that many taxa have "NA" in different ranks.

# Fixing the phyloseq:
ps_data_fixed<- ps_non_dup %>%
  tax_fix(
    min_length = 4,
    unknowns = c("s__NA", "g__NA", "f__NA", "o__NA", "c__NA", "p__NA"),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified"
    )
# IT is important to fix this problem before aggregation to correct for short values, like “g__NA”, etc.
# By using annon_unique: This will rename the detected "NA" with the name of the higher rank in its taxonomy

# Inspecting the fixed phyloseq object:
phyloseq_validate(ps_data_fixed)
# No problems are detected

################################################################################################
################################################################################################
################################################################################################

### Aim 1) Investigating the potential association between the gut microbiome and Non-Celiac Wheat-Sensitive (NCWS) patients
# (Use only T0) (Use HC vs NCWS) **********************************************************

### 1.a) Aggregation: 

# Before filtration: Aggregation without filteration, only for visualisation and determining the threshold
# Important to apply filtration per group
ps_data_agg_no_fil <- ps_data_fixed %>% 
  ps_filter(time_point == "T0")  %>%
  tax_agg(rank = "Family")

### 1.b) Filtarion: 

# To remove rare taxa:
# Low prevalence if taxon only detected in a small number of samples in dataset.
# Low abundance when relatively few reads assigned to that taxon.

# Table of summary statistics for the unique taxa in in HC and NCWS groups:

ps_data_agg_HC <- ps_data_agg_no_fil %>%
  ps_filter(patient_type == "HC")

ps_data_agg_NCWS <- ps_data_agg_no_fil %>%
  ps_filter(patient_type == "NCWS")

Taxa_Stats_HC <- tibble(
  taxon = taxa_names(ps_data_agg_HC@tax_table),
  prevalence = microbiome::prevalence(t(ps_data_agg_HC@otu_table)),
  total_abundance = taxa_sums(ps_data_agg_HC@otu_table)
)

Taxa_Stats_NCWS <- tibble(
  taxon = taxa_names(ps_data_agg_NCWS@tax_table),
  prevalence = microbiome::prevalence(t(ps_data_agg_NCWS@otu_table)),
  total_abundance = taxa_sums(ps_data_agg_NCWS@otu_table)
)

# Visulaisation using the prevalence-abundance plot:
Taxa_Stats_HC %>%
  ggplot(aes(x = total_abundance, y = prevalence)) +
  geom_vline(xintercept = 10000, color = "red", linetype = "dotted") +
  geom_hline(yintercept = 3 / 100, color = "red", linetype = "dotted") +
  ggtitle("Taxa statistics in HC") +
  geom_point(alpha = 0.5) +
  geom_rug(alpha = 0.1) +
  scale_x_log10(labels = scales::label_number(), name = "Total Abundance") +
  scale_y_log10(
    labels = scales::label_percent(), breaks = scales::breaks_log(n = 9),
    name = "Prevalence (%)",
    sec.axis = sec_axis(
      trans = ~ . * nsamples(ps_data_agg_HC), breaks = scales::breaks_log(n = 9),
      name = "Prevalence (N samples)"
    )
  )  +
  theme_bw()

Taxa_Stats_NCWS %>%
  ggplot(aes(x = total_abundance, y = prevalence)) +
  geom_vline(xintercept = 10000, color = "red", linetype = "dotted") +
  geom_hline(yintercept = 3 / 100, color = "red", linetype = "dotted") +
  ggtitle("Taxa statistics in NCWS") +
  geom_point(alpha = 0.5) +
  geom_rug(alpha = 0.1) +
  scale_x_log10(labels = scales::label_number(), name = "Total Abundance") +
  scale_y_log10(
    labels = scales::label_percent(), breaks = scales::breaks_log(n = 9),
    name = "Prevalence (%)",
    sec.axis = sec_axis(
      trans = ~ . * nsamples(ps_data_agg_NCWS), breaks = scales::breaks_log(n = 9),
      name = "Prevalence (N samples)"
    )
  )  +
  theme_bw()
# not interested in unique taxa that occur in fewer than 3% of samples
# and they have to have at least 10000 reads in total across all samples.

# For clarification:
# Prevalence 100% means that this taxon presents in all samples
# Total abundance means the count of this taxon in all samples

#Depends on what analysis method:
# alpha diversity = No need for filtration
# beta diversity = need filtering
# differential abundance testing = more stringent filtering 

### 1.b.1) Filtration, Transformation and mean centring for beta-diversity

ps_data_filt_HC <- ps_data_agg_no_fil %>%
  ps_filter(patient_type == "HC")%>% 
  tax_filter(min_prevalence = 3 / 100, prev_detection_threshold = 1, min_total_abundance = 10000, verbose = FALSE, tax_level = "Family")

ps_data_filt_NCWS <- ps_data_agg_no_fil %>%
  ps_filter(patient_type == "NCWS")%>% 
  tax_filter(min_prevalence = 3 / 100, prev_detection_threshold = 1,min_total_abundance = 10000,verbose = FALSE, tax_level = "Family")
ps_filtered_T0 <- merge_phyloseq(ps_data_filt_HC, ps_data_filt_NCWS)

# This removes many taxa that contain many zeros. However, this threshold removes many taxa that are of low abundance
# and can be distinctive. In this project, I am only focus on the main key taxa that of high abundant.
# in the future projects, I can only remove most of the taxa with the zero count before using a less strict cuttoff if
# it is important to see the differences of the low abundant taxa.

# Imputation and Transfomration:
ps_trans_clr_T0 <- ps_filtered_T0 %>%
  tax_agg("Family") %>%
  tax_transform("clr", pseudocount =1) 

# Aggregation on family level due to undetermined genus level in more than 30% of the taxa
# Imputation using a pseudocount =1, to replace zeros in the table count, however this can introduce bias specifically in 
# the low abundant taxa. This is not a problem here I already don't focus on the low abundant ones.
# clr: centered log-ratio (CLR) transformation to the taxonomic abundances at the family level. 
#This normalize the data and account for the compositional nature of microbiome data. 

### Mean centering scaling over clr transformed data:
# Extract otu_table:
otu_table <- otu_table(ps_trans_clr_T0)
otu_table_df <- as.data.frame(otu_table)
sample_metadata_T0 <- sample_metadata[sample_metadata$time_point == "T0", ]

patient <- sample_metadata_T0$patient
otu_table_df <- cbind(patient, otu_table_df)
otu_table_df <- rownames_to_column(otu_table_df, var = "row_name")

centered_otu_table <- otu_table_df %>%
  group_by(patient) %>%
  mutate_at(vars(starts_with("asv")), ~ . - mean(.)) %>%
  ungroup() %>%
  as.data.frame()
centered_otu_table <- column_to_rownames(centered_otu_table, var = "row_name")
centered_otu_table[1] <- NULL

ps_trans_clr_T0@otu_table <- otu_table(centered_otu_table, taxa_are_rows = FALSE)
ps_trans_clr_T0

### 1.b.2) Filtration, Transformation and mean centring for differential abundance testing 

# Again, repeat the previous filtration steps, but with more strict threshold for differential abundance:
ps_data_filt_HC_DA <- ps_data_agg_no_fil %>%
  ps_filter(patient_type == "HC")%>% 
            tax_filter(min_prevalence = 10 / 100, prev_detection_threshold = 1, min_total_abundance = 10000, verbose = FALSE,
            undetected = 0, use_counts = TRUE, tax_level = "Family")

ps_data_filt_NCWS_DA <- ps_data_agg_no_fil %>%
  ps_filter(patient_type == "NCWS")%>% 
  tax_filter(min_prevalence = 10 / 100, prev_detection_threshold = 1, min_total_abundance = 10000, verbose = FALSE,
            undetected = 0, use_counts = TRUE, tax_level = "Family")

ps_filtered_DA_T0 <- merge_phyloseq(ps_data_filt_HC_DA, ps_data_filt_NCWS_DA)

# Transformation:
ps_trans_comp_T0 <- ps_filtered_DA_T0%>%
  tax_agg("Family") %>%
  tax_transform(trans = "compositional", rank = "Family", keep_counts = TRUE)
# for various statistical, biological, and practical reasons, let's strictly filter taxa

# Mean centering over compositional-transformed data
# Extract otu_table:
otu_table <- otu_table(ps_trans_comp_T0)
otu_table_df <- as.data.frame(otu_table)
sample_metadata_T0 <- sample_metadata[sample_metadata$time_point == "T0", ]

patient <- sample_metadata_T0$patient
otu_table_df <- cbind(patient, otu_table_df)
otu_table_df <- rownames_to_column(otu_table_df, var = "row_name")

centered_otu_table <- otu_table_df %>%
  group_by(patient) %>%
  mutate_at(vars(starts_with("asv")), ~ . - mean(.)) %>%
  ungroup() %>%
  as.data.frame()
centered_otu_table <- column_to_rownames(centered_otu_table, var = "row_name")
centered_otu_table[1] <- NULL

ps_trans_comp_T0@otu_table <- otu_table(centered_otu_table, taxa_are_rows = FALSE)
ps_trans_comp_T0

# 1.c) Alpha- Diversity: ((NO FILTRATION NEEDED))

#WHY?
#insights into the richness and evenness of microbial communities within individual samples.

# Notes:
# The library sizes can dominate the biology in determining the result of the diversity analysis
# Therefore, it is likely to observe higher numbers of different taxa in the sample with more microbial reads

#HOW?
# Checking the noise of readcount
# It is valid to compare richness across samples, as the readcount variation is only random noise, 
# as can be seen from the following:

ps_data_fixed %>%
  ps_filter(time_point == "T0") %>%
  ps_calc_diversity(index = "shannon", rank = "Family", varname = "Shannon_Family") %>% 
  ps_mutate(Effective_Shannon_Family = exp(Shannon_Family)) %>%
  samdat_tbl() %>%
  ggplot(aes(DNA_con, Effective_Shannon_Family)) +
  geom_point(alpha = 0.4, size = 2.5) +
  theme_bw(14)

# 1.d) Visualisation using Composition barplot

ps_data_fixed %>% 
  ps_filter(time_point == "T0") %>%
  comp_barplot(
    tax_level = "Family",
    label = "patient",
    n_taxa = 10, # give 10 taxa unique colours
    taxon_renamer = function(x) stringr::str_replace_all(x, c("o__"="","f__" ="")), # this will remove underscores
    other_name = "Other families", # to set a custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, 
    bar_outline_colour = "grey5" 
  )  +
  facet_wrap(vars(patient_type), scales = "free") +
  coord_flip()

# Insights:
# Some families are more enriched in NCWS such as (Prevotellaceae and Rikenellaceae)

# 1.e) Statistics with alpha diversity

ps_data_fixed %>%
  ps_filter(time_point == "T0") %>%
  ps_calc_diversity(index = "shannon", rank = "Family", varname = "Shannon_Family") %>% 
  ps_mutate(Effective_Shannon_Family = exp(Shannon_Family)) %>%
  samdat_tbl() %>%
  ggplot(aes(y = patient_type, x = Effective_Shannon_Family)) +
  geom_boxplot(width = 0.3) +
  geom_point(position = position_jitter(height = 0.2), alpha = 0.5) +
  theme_bw()
# It seems that NCWS group has less diversity, let's do a statistical test
# Effective Shannon index: combine richness and evenness. Therefore, it is better than richness alone.

# Statistical test:
ps_data_fixed %>%
  ps_filter(time_point == "T0") %>%
  ps_calc_diversity(index = "shannon", rank = "Family", varname = "Shannon_Family") %>% 
  ps_mutate(Effective_Shannon_Family = exp(Shannon_Family)) %>%
  samdat_tbl() %>%
  wilcox.test(formula = Effective_Shannon_Family ~ patient_type, data = .)
# Not significant different, wilcox is  non-parametric test

# Linear model:
ps_data_fixed %>%
  ps_filter(time_point == "T0") %>%
  ps_calc_diversity(index = "shannon", rank = "Family", varname = "Shannon_Family") %>% 
  ps_mutate(Effective_Shannon_Family = exp(Shannon_Family)) %>%
  samdat_tbl() %>%
  lm(formula = Effective_Shannon_Family ~ patient_type + age +sex, data = .) %>%
  summary()
# This shows that age and sex are significant predictors but patient_type is NOT!

# 1.f) Dissimilarity and ordination: ______________________________________________________________________________

## similar samples will be close to each other 

## 1.f.1) PCA:

ps_trans_clr_T0 %>%
  ord_calc("PCA") %>%
  ord_plot(color = "patient_type", size = 2) +
  scale_colour_brewer(palette = "Dark2")
## Separation of the 2 groups based on PC1
ps_trans_clr_T0 %>%
  ord_calc("PCA") %>%
  ord_plot(axes = c(3, 4), color = "patient_type", size = 2) +
  scale_colour_brewer(palette = "Dark2")

## 1.f.2) PCA Biplot:

pca_bi <- ps_trans_clr_T0 %>%
  ord_calc("PCA") %>%
  ord_plot(axes = c(1, 2),
    colour = "patient_type", alpha = 0.7, size = 2,
    plot_taxa = 1:30, tax_vec_length = 0.5,
    taxon_renamer = function(x) stringr::str_replace_all(x, c("o__"="","f__" ="")), # remove underscores
    center = TRUE,
    tax_lab_style = tax_lab_style(
      type = "label", max_angle = 90, fontface = "bold",
      alpha = 0.7, size = 3.2
    )
  ) 
pca_bi
## PCA Biplot shows colouring samples by group and adding taxon loading arrows to visualise which taxa generally differ across samples 

## 1.f.3) Non-metric Multidimensional Scaling (NMDS):
ps_trans_clr_T0 %>%
#  tax_transform("identity", rank = "Family") %>% # don't transform!
  dist_calc("euclidean") %>%
  ord_calc("NMDS") %>%
  ord_plot(color = "patient_type", shape = "sex", size = 2) +
  scale_colour_brewer(palette = "Dark2", aesthetics = c("fill", "colour"), name = "Patient type") +
  theme_bw() +
  ggside::geom_xsidedensity(aes(fill = patient_type), alpha = 0.5, show.legend = FALSE) +
  ggside::geom_ysidedensity(aes(fill = patient_type), alpha = 0.5, show.legend = FALSE) +
  ggside::theme_ggside_void()
## NMDS is a dimensionality reduction technique that aims to preserve the rank order of distances between samples in a lower-dimensional space.
## Using euclidean distance over clr transformed data
# No Outliers are visualised


# 1.g) Statistical analysis using PERMANOVA:

# WHY? What variables is the overall microbial community variation associated with?

ps_trans_clr_T0 %>%
  dist_calc(dist = "euclidean") %>%
  dist_permanova(variables = c("patient_type", "sex", "age"), n_perms = 9999, seed = 123, 
                 complete_cases = TRUE, verbose = "max") %>%
  perm_get()
# A good approachthat can work on distance matrices and don't assume normality becasue it depends on permutating the samples
# however, other methods such as ASCA could be a good future approach.
# P-value is significant means that there is difference in microbiota composition (statistically)
# Therefore, there is an association between microbiome and being NCWS


# 1.h) Differential abundance testing

#Why? To test statistically which taxa have a higher relative abundance

ps_data_filtered_DAT <- ps_trans_comp_T0 %>%
  tax_prepend_ranks()
ps_data_treeStats <- ps_data_filtered_DAT %>%
  # run all the statistical models
  taxatree_models(
    ranks = c("Phylum", "Class", "Order","Family"),
    trans = "log2", trans_args = list(zero_replace = "halfmin"),
    variables = "patient_type", type = lm # modelling function
  ) %>%
  # extract stats from the models
  taxatree_models2stats(.keep_models = TRUE) %>%
  # adjust the p values for multiple testing, within each rank
  taxatree_stats_p_adjust(method = "fdr", grouping = "rank")

# Using half-min imputation that depends on using the half for the minimum value per group to impute the zero values
# Then, log2 transformation on the imputed compositional  data
taxatree_stats_get(ps_data_treeStats)

treePlotsSimple <- ps_data_treeStats %>%
  # specify which taxa will get labeled (adds a "label" variable to the stats tibble)
  taxatree_label(p.adj.fdr.rank < 0.01, rank %in% c("Phylum", "Class", "Order","Family")) %>%
  # make the plots (1 per predictor variable, so a list of 1 in this example)
  taxatree_plots(
    sig_stat = "p.adj.fdr.rank", sig_threshold = 0.01,
    drop_ranks = FALSE # drop_ranks = TRUE has a bug in version 0.10.4 :(
  )

treePlotsSimple %>% str(max.level = 1)

treePlotsSimple$patient_type %>%
  # add labels to the plot, only for the taxa indicated earlier
  taxatree_plot_labels(
    taxon_renamer = function(x) {
      gsub(pattern = "[pg]: ", replacement = "", x = x) %>%
        stringr::str_replace_all(c("o__" = "", "f__" = "", "k__" ="", "p__" = "", "c__" = ""))
    }
  ) +
  coord_fixed(expand = TRUE, clip = "off") + # allow scale expansion
  scale_x_continuous(expand = expansion(mult = 0.2)) # make space for labels

# At one taxa using linear regression:

ps_DAT_Prevot <- ps_trans_comp_T0 %>%
  tax_filter(min_prevalence = 10/100, min_total_abundance = 10000, use_counts = TRUE, tax_level = "Family")
# strictly filter taxa


Prevotellaceae_lm <- ps_DAT_Prevot %>%
  tax_model(
    type = "lm", rank = "Family", taxa = "f__Prevotellaceae",
    trans = "log2", trans_args = list(zero_replace = "halfmin"),
    variables = c("patient_type", "sex", "age"),
    return_psx = FALSE
  )

summary(Prevotellaceae_lm$f__Prevotellaceae)
##############################################################################################################
# Findings:
### P-value in PERMANOVA is significant means that there is difference in microbiota composition (statistically)
### Therefore, there is an association between microbiome and being NCWS
## From DAT:
### Many taxa have a higher relative abundance in NCWS
##############################################################################################################
##############################################################################################################
##############################################################################################################


# Aim 2) Finding an association between the gut microbiome and different bread types in NCWS patients

# (Use only NCWS) **********************************************************

## 2.a) Aggregation

# Before filtration: Aggregation without filteration, only for visualisation and determining the threshold

# Selecting only NCWS group:
ps2_data_agg_no_fil <- ps_data_fixed %>% 
  ps_filter(patient_type == "NCWS")  %>%
  tax_agg(rank = "Family")


### 2.b) Filtarion: 
### 2.b.1) Filtration, Transformation and mean centring for beta-diveristy

bread_types <- c("sourdough_wheat", "yeast_spelt", "sourdough_spelt", "yeast_emmer", "sourdough_emmer", "yeast_wheat", "medium_medium")

ps2_list <- list()
for (bread_type in bread_types) {
  ps2_list[[bread_type]] <- ps2_data_agg_no_fil %>%
    ps_filter(bread_type == bread_type) %>%
    tax_filter(min_prevalence = 3 / 100, min_total_abundance = 10000, prev_detection_threshold = 1, verbose = FALSE, tax_level = "Family")
}

ps2_filtered_NCWS <- merge_phyloseq(ps2_list[[1]], ps2_list[[2]], ps2_list[[3]], ps2_list[[4]], ps2_list[[5]],
                                    ps2_list[[6]], ps2_list[[7]])

# Transfomration:
ps2_trans_clr_NCWS <- ps2_filtered_NCWS %>%
  tax_agg("Family") %>%
  tax_transform("clr")

# Aggregation on family level due to undetermined genus level in more than 30% of the taxa
# clr: centered log-ratio (CLR) transformation to the taxonomic abundances at the family level. 
#This normalize the data and account for the compositional nature of microbiome data. 

### Mean centering scaling over clr transformed data:
# Extract otu_table:
otu2_table <- otu_table(ps2_trans_clr_NCWS)
otu2_table_df <- as.data.frame(otu2_table)
sample_metadata_NCWS <- sample_metadata[sample_metadata$patient_type == "NCWS", ]

patient2 <- sample_metadata_NCWS$patient
otu2_table_df <- cbind(patient2, otu2_table_df)
otu2_table_df <- rownames_to_column(otu2_table_df, var = "row_name")

centered_otu2_table <- otu2_table_df %>%
  group_by(patient2) %>%
  mutate_at(vars(starts_with("asv")), ~ . - mean(.)) %>%
  ungroup() %>%
  as.data.frame()
centered_otu2_table <- column_to_rownames(centered_otu2_table, var = "row_name")
centered_otu2_table[1] <- NULL

ps2_trans_clr_NCWS@otu_table <- otu_table(centered_otu2_table, taxa_are_rows = FALSE)
ps2_trans_clr_NCWS

### 2.b.2) Filtration, Transformation and mean centring for differential abundance 

# Again, repeat the previous filtration steps, but with more strict threshold for differential abundance:

bread_types <- c("sourdough_wheat", "yeast_spelt", "sourdough_spelt", "yeast_emmer", "sourdough_emmer", "yeast_wheat", "medium_medium")

ps2_list_DA <- list()
for (bread_type in bread_types) {
  ps2_list_DA[[bread_type]] <- ps2_data_agg_no_fil %>%
    ps_filter(bread_type == bread_type) %>%
    tax_filter(min_prevalence = 10 / 100, prev_detection_threshold = 1, verbose = FALSE,
               min_total_abundance = 10000, undetected = 0, use_counts = TRUE, tax_level = "Family")
}

ps2_filtered_DA <- merge_phyloseq(ps2_list_DA[[1]], ps2_list_DA[[2]], ps2_list_DA[[3]], ps2_list_DA[[4]], ps2_list_DA[[5]],
                                  ps2_list_DA[[6]], ps2_list_DA[[7]])

# Transformation:
ps2_trans_comp_NCWS <- ps2_filtered_DA%>%
  tax_agg("Family") %>%
  tax_transform(trans = "compositional", rank = "Family", keep_counts = TRUE)
# for various statistical, biological, and practical reasons, let's strictly filter taxa

# Mean centering over compositional-transformed data
# Extract otu_table:
otu2_table <- otu_table(ps2_trans_comp_NCWS)
otu2_table_df <- as.data.frame(otu2_table)
sample_metadata_NCWS <- sample_metadata[sample_metadata$patient_type == "NCWS", ]

patient2 <- sample_metadata_NCWS$patient
otu2_table_df <- cbind(patient2, otu2_table_df)
otu2_table_df <- rownames_to_column(otu2_table_df, var = "row_name")

centered_otu2_table <- otu2_table_df %>%
  group_by(patient2) %>%
  mutate_at(vars(starts_with("asv")), ~ . - mean(.)) %>%
  ungroup() %>%
  as.data.frame()
centered_otu2_table <- column_to_rownames(centered_otu2_table, var = "row_name")
centered_otu2_table[1] <- NULL

ps2_trans_comp_NCWS@otu_table <- otu_table(centered_otu2_table, taxa_are_rows = FALSE)
ps2_trans_comp_NCWS


# 2.c) Alpha- Diversity: ((NO FILTRATION NEEDED))

# Checking the noise of readcount
#It is usually valid to compare richness across groups of samples:

ps_data_fixed %>%
  ps_filter(patient_type == "NCWS") %>%
  ps_calc_diversity(index = "shannon", rank = "Family", varname = "Shannon_Family") %>% 
  ps_mutate(Effective_Shannon_Family = exp(Shannon_Family)) %>%
  samdat_tbl() %>%
  ggplot(aes(DNA_con, Effective_Shannon_Family)) +
  geom_point(alpha = 0.4, size = 2.5) +
  theme_bw(14)
# This shows that readcount variation is only random noise, 

# 2.d) Visualisation using Composition barplot

ps_data_fixed %>% 
  ps_filter(patient_type == "NCWS") %>%
  comp_barplot(
    tax_level = "Family",
    label = "time_point", 
    n_taxa = 10, # this gives more taxa unique colours
    taxon_renamer = function(x) stringr::str_replace_all(x, c("o__"="","f__" ="")), # to remove underscores
    other_name = "Other families", # this is to set custom name for the "other" category more than the 10 taxa
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # to reduce the bar width to 70% of one row
    bar_outline_colour = "grey5"
  )  +
  facet_wrap(vars(bread_type), scales = "free") +
  coord_flip()

## No differences can be detected in richness or distribution of taxa

# 2.e) Statistics with alpha diversity

ps_data_fixed %>%
  ps_filter(patient_type == "NCWS") %>%
  ps_calc_diversity(index = "shannon", rank = "Family", varname = "Shannon_Family") %>% 
  ps_mutate(Effective_Shannon_Family = exp(Shannon_Family)) %>%
  samdat_tbl() %>%
  ggplot(aes(y = bread_type, x = Effective_Shannon_Family)) +
  geom_boxplot(width = 0.3) +
  geom_point(position = position_jitter(height = 0.2), alpha = 0.5) +
  theme_bw()
# It seems that there is no statistical differences in richness and diversity bet. different bread types

# Linear model:
ps_data_fixed %>%
  ps_filter(patient_type == "NCWS") %>%
  ps_calc_diversity(index = "shannon", rank = "Family", varname = "Shannon_Family") %>% 
  ps_mutate(Effective_Shannon_Family = exp(Shannon_Family)) %>%
  samdat_tbl() %>%
  lm(formula = Effective_Shannon_Family ~ bread_type + age +sex + time_point, data = .) %>%
  summary()
# This linear model including the covariates shows that bread_type is not significant

# 2.f) Dissimilarity and ordination: ______________________________________________________________________________

## similar samples will be close to each other 

## 2.f.1) PCA:

ps2_trans_clr_NCWS %>%
  ord_calc("PCA") %>%
  ord_plot(axes = c(1, 2), color = "bread_type", shape = "time_point", size = 2) +
  scale_colour_brewer(palette = "Dark2")
## PC1 and PC2 are not informative!
ps2_trans_clr_NCWS %>%
  ord_calc("PCA") %>%
  ord_plot(axes = c(3, 4), color = "bread_type", shape = "time_point", size = 2) +
  scale_colour_brewer(palette = "Dark2")
# PC3 !! informative for Time_point not bread_type
ps2_trans_clr_NCWS %>%
  ord_calc("PCA") %>%
  ord_plot(axes = c(3, 4), color = "bread_type", size = 2, shape = "time_point") + theme_classic(12) +
  coord_fixed(0.7) +
  stat_ellipse(aes(shape = time_point)) +
  scale_color_brewer(palette = "Set1") + 
  theme_bw() +
  ggside::geom_xsidedensity(aes(fill = time_point), alpha = 0.5, show.legend = FALSE) +
  ggside::geom_ysidedensity(aes(fill = time_point), alpha = 0.5, show.legend = FALSE) +
  ggside::theme_ggside_void()

## 2.f.2) Non-metric Multidimensional Scaling (NMDS):
ps2_trans_clr_NCWS %>%
  dist_calc("euclidean") %>%
  ord_calc("NMDS") %>%
  ord_plot(color = "bread_type", shape = "time_point",size = 2) +
  scale_colour_brewer(palette = "Dark2", aesthetics = c("colour"), name = "bread_type") +
  theme_bw() +
  ggside::geom_xsidedensity(aes(fill = time_point), alpha = 0.5, show.legend = FALSE) +
  ggside::geom_ysidedensity(aes(fill = time_point), alpha = 0.5, show.legend = FALSE) +
  ggside::theme_ggside_void()
## NMDS is a dimensionality reduction technique that aims to preserve the rank order of distances between samples in a lower-dimensional space.
# No Outliers are visualised

# 2.g) Statistical analysis using PERMANOVA:

# WHY? What variables is the overall microbial community variation associated with?

ps2_trans_clr_NCWS %>%
  dist_calc(dist = "euclidean") %>%
  dist_permanova(variables = c("bread_type", "sex", "age","time_point"), n_perms = 9999, seed = 123, 
                 complete_cases = TRUE, verbose = "max") %>%
  perm_get()
# P-value is not significant means that there is no difference in microbiota composition (statistically)
# Therefore, there is NO association between microbiome and different bread types

# 2.h) Differential abundance testing

#Why? To test statistically which taxa have a higher relative abundance

ps2_data_filtered_DAT <- ps2_trans_comp_NCWS %>%
  tax_prepend_ranks()
ps2_data_treeStats <- ps2_data_filtered_DAT %>%
  # run all the statistical models
  taxatree_models(
    ranks = c("Family"),
    trans = "log2", trans_args = list(zero_replace = "halfmin"),
    variables = "bread_type", type = lm # modelling function
  ) %>%
  # extract stats from the models
  taxatree_models2stats(.keep_models = TRUE) %>%
  # adjust the p values for multiple testing, within each rank
  taxatree_stats_p_adjust(method = "fdr", grouping = "rank")

## Imputed and transformed count table

taxatree_stats_get(ps2_data_treeStats)

treePlotsSimple2 <- ps2_data_treeStats %>%
  # specify which taxa will get labeled (adds a "label" variable to the stats tibble)
  taxatree_label(p.adj.fdr.rank < 0.05, rank %in% c("Family")) %>%
  # make the plots (1 per predictor variable, so a list of 1 in this example)
  taxatree_plots(
    sig_stat = "p.adj.fdr.rank", sig_threshold = 0.05,
    drop_ranks = FALSE # drop_ranks = TRUE has a bug in version 0.10.4 :(
  )

treePlotsSimple2 %>% str(max.level = 1)

p1 <- treePlotsSimple2$bread_typesourdough_emmer %>%
  # add labels to the plot, only for the taxa indicated earlier
  taxatree_plot_labels(
    taxon_renamer = function(x) {
      gsub(pattern = "[pg]: ", replacement = "", x = x) %>%
        stringr::str_replace_all(c("o__" = "", "f__" = "", "k__" ="", "p__" = "", "c__" = ""))
    }
  ) +
  coord_fixed(expand = TRUE, clip = "off") + # allow scale expansion
  scale_x_continuous(expand = expansion(mult = 0.2)) + theme(legend.position = "none")

p2 <- treePlotsSimple2$bread_typesourdough_spelt %>%
  # add labels to the plot, only for the taxa indicated earlier
  taxatree_plot_labels(
    taxon_renamer = function(x) {
      gsub(pattern = "[pg]: ", replacement = "", x = x) %>%
        stringr::str_replace_all(c("o__" = "", "f__" = "", "k__" ="", "p__" = "", "c__" = ""))
    }
  ) +
  coord_fixed(expand = TRUE, clip = "off") + # allow scale expansion
  scale_x_continuous(expand = expansion(mult = 0.2)) + theme(legend.position = "none")

p3 <- treePlotsSimple2$bread_typesourdough_wheat %>%
  # add labels to the plot, only for the taxa indicated earlier
  taxatree_plot_labels(
    taxon_renamer = function(x) {
      gsub(pattern = "[pg]: ", replacement = "", x = x) %>%
        stringr::str_replace_all(c("o__" = "", "f__" = "", "k__" ="", "p__" = "", "c__" = ""))
    }
  ) +
  coord_fixed(expand = TRUE, clip = "off") + # allow scale expansion
  scale_x_continuous(expand = expansion(mult = 0.2)) 

p4 <- treePlotsSimple2$bread_typeyeast_emmer %>%
  # add labels to the plot, only for the taxa indicated earlier
  taxatree_plot_labels(
    taxon_renamer = function(x) {
      gsub(pattern = "[pg]: ", replacement = "", x = x) %>%
        stringr::str_replace_all(c("o__" = "", "f__" = "", "k__" ="", "p__" = "", "c__" = ""))
    }
  ) +
  coord_fixed(expand = TRUE, clip = "off") + # allow scale expansion
  scale_x_continuous(expand = expansion(mult = 0.2)) + theme(legend.position = "none")

p5 <- treePlotsSimple2$bread_typeyeast_spelt %>%
  # add labels to the plot, only for the taxa indicated earlier
  taxatree_plot_labels(
    taxon_renamer = function(x) {
      gsub(pattern = "[pg]: ", replacement = "", x = x) %>%
        stringr::str_replace_all(c("o__" = "", "f__" = "", "k__" ="", "p__" = "", "c__" = ""))
    }
  ) +
  coord_fixed(expand = TRUE, clip = "off") + # allow scale expansion
  scale_x_continuous(expand = expansion(mult = 0.2)) + theme(legend.position = "none")

p6 <- treePlotsSimple2$bread_typeyeast_wheat %>%
  # add labels to the plot, only for the taxa indicated earlier
  taxatree_plot_labels(
    taxon_renamer = function(x) {
      gsub(pattern = "[pg]: ", replacement = "", x = x) %>%
        stringr::str_replace_all(c("o__" = "", "f__" = "", "k__" ="", "p__" = "", "c__" = ""))
    }
  ) +
  coord_fixed(expand = TRUE, clip = "off") + # allow scale expansion
  scale_x_continuous(expand = expansion(mult = 0.2)) + theme(legend.position = "none")

#Combine all plots using patchwork
patchwork::wrap_plots(p1,p2,p3,p4,p5,p6, ncol = 3, guides = "collect")

##############################################################################################################
# Findings of Aim 2:
### P-value in PERMANOVA is not significant means that there is no difference in microbiota composition (statistically) with different bread-types
### Therefore, there is no association between microbiome and bread types
## From DAT:
### No differential abundance
##############################################################################################################
##############################################################################################################
##############################################################################################################
