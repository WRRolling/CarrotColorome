This GitHub repository contains the analyses completed as part of a project to complete genomic analysis of carotenoid accumulation in carrots. Analyses include phenotypic data organization,  filtering, outlier removal, and presentation. There are analyses related to genotypic data variant calling, filtering, and formatting. Genomic analysis included population structure, genome-wide association analysis (GWAS), genomic prediction (GP), and producing associated figures. 

A full written description of this project are found in the publication: (TBD) 
Images of carrot roots and a file genotypic data file from polyRAD is available at [CarrotOmics.org/.](https://carrotomics.org/). 
Software used in these analysis include: (not an all software listed)
https://github.com/lvclark/polyRAD


# Colorome Paper
Analyses used in publication: ""

1) The "Analysis_1_HPLCDat_CA19" folder includes analyses and outputs related to QC-filtered HPLC data from the CA_2019_732PI trial. Data from these analyses are presented in Supplemental Table 2 and Supplemental Figure 3.
2) The "Analysis_2_FigsHPLCdat_CA19" folder contains analyses and outputs related to HPLC data that was further filtered based on appropriate HPLC data for each color category.
   Test of normality were completed and Outliers were removed here. 
   The counts of each color score were calculated for supplemental 4.
3) The folder "Analysis_3_HPLCDat_WI18" contains analyses and output related to QC-filter HPLC data from the WI_2018_432 trial.
     Data from these analyses are included in Supplemental Table 2 and Supplemental Figure 1. 
4) The folder "Analysis_4_FigsHPLCdat_WI18" contains analyses and output-related contains analyses and output-related HPLC data was further filtered checking for appropriate HPLC data for each color category.
5) The folder "Analysis_5_DescribePhenotypeData_WI18" Removes outliers from the WI_2018_432 trial, and calculates the mean for each carotenoid. Mean separation for Supplemental Figure 2a,b is included here as well as the boxplot in the figure. These analyses include the correlation test found in Table 1. A phenotype file for  Genomic Prediction and GWA analyses was created during these analyses. 
6)  The folder: "Analysis_6_DescribePhenotypeData_CA19" contains analyses and output related to mean separation between colors, and the boxplot found in Supplemental Figure 4 and Supplemental Table 6. The correlation and carotenoid means for Table 1 were calculated with this analysis. This analysis calculated the core color carotenoid comparisons found in Supplemental Figure 2. A comparison of HPLC data (i.e., mean separation and correlation) found in supplemental table 7 was calculated here). 
7) The folder: "Analysis_7_CA19.VCF" contains a text file describing how TASSEL was used to identify variants and SNP data was filtered.
8) The folder: "Analysis_8_CombineandFilterVCF" contains analysis to combine genotypic data from two trial (Resequencing data dn GBS data).
9) The folder: "Analysis_9_LD_PopStructure" contains files to complete population structure analysis with admxture and fastStructure and calculate LD decay.
10) The folder: Analysis_10_polyRAD" contains files related to variant calling and filtering in polyRAD
11) The folder: "Analysis_11_GWA" contains files related to completing GWAS.
12) The folder: "Analysis_12_HaploviewLDblocks" contains the files used in Haploview.
13) The folder: "Analysis_13_GenomicPredictionHPLC" contains example files of how BGLR was used to complete Genomic Prediction
14) The folder: "Analysis_14_GPFigures" contains R scripts used to create figures from the GP figures.
15) The folder: "Analysis_15_CarotenoidsEnvironment" contains R analysis used to complete calculations related to intra-plot, intra-field, and location differences in order to estimate the minimum effective sample to calculate an average carotenoid concentration.  
