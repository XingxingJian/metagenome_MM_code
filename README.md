# metagenome_MM_code
"Alterations of gut microbiome accelerate multiple myeloma progression by increasing the relative abundances of nitrogen recycling bacteria"

All code using in the manuscript, Alterations of gut microbiome accelerate multiple myeloma progression by increasing the relative abundances of nitrogen recycling bacteria, is presented here, including two Shell script in Linux and several R script in Windows.
1.	run_kraken_standard: this script was used for taxa classification using Kraken.
2.	run_metacv: this script was used for function annotation using MetaCV.
3.	Code_1_species_analysis.R: this script was used to calculate alpha-diversity within the samples and beta- diversity between the two groups, as well as visualization.
4.	Code_2_phylum_genus.R: this script was used to analyze statistically taxa classification at the phylum and genus level, as well as visualization.
5.	Code_3_diff_species.R: this script was used to identify significantly differential species using R package DESeq2, as well as visualization.
6.	Code_4_diff_ko.R: this script was used to identify significantly differential KOs using R package DESeq2, and calculate beta- diversity between the two groups.
7.	Code_5_metabolite.R: this script was used to screen significant metabolites which was p_value <0.05 (two-tailed T-test) and VIP >1 (SIMCA), as well as visualization.
8.	Code_6_interaction.R: this script was used to calculate the Spearman’s correlation between the differential species and differential metabolites, as well as visualization.
9.	Code_7_microbiotaViz.R: this script was used to visualize the Cladogram of the significant species using R package microbiomeViz.
10.	Code_RDA_analysis.R: this script was used to perform redundancy analysis (RDA), which presented the variation of the relative abundance of several significant strains in FMT_HC, PBS, and FMT_MM mice with time (week 0, week 2, week 4, week 6).
11.	sample_group: all samples using in this study were divided into two groups, HC and MM.
