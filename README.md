This GitHub repository houses the analyses conducted for a project investigating the genomic basis of carotenoid accumulation in carrots. The work encompassed phenotypic data management (organization, filtering, outlier removal, visualization), genotypic data processing (variant calling, filtering, formatting), and population-level analyses (population structure, genome-wide association study, genomic prediction, and associated visualizations).

--------------------------------------------------------------------------
A full written description of this project can be found in the following publication: (TBD) 

Images of carrot roots and a file containing genotypic data from polyRAD are available on [CarrotOmics](https://carrotomics.org/). 

----------------------------------------------------------------------------
Software used in these analyses includes: 
  - [BGLR](https://github.com/gdlc/BGLR-R/) 
  - [fastStructure](https://rajanil.github.io/fastStructure/)
  - [GAPIT](https://github.com/jiabowang/GAPIT)
  - [Haploview](https://www.broadinstitute.org/haploview/haploview)
  - [Plink](https://www.cog-genomics.org/plink/)
  - [polyRAD](https://github.com/lvclark/polyRAD/) 
  - [qqman](https://github.com/stephenturner/qqman)
  - [simpleM](https://simplem.sourceforge.net/)
  - [TASSEL](https://github.com/maize-genetics/tassel-6-source)
  - [Tidyverse](https://www.tidyverse.org/) 
  - [VCFtools](https://github.com/vcftools)
-------------------------------------------
Folder Descriptions:

1) The "Analysis_1_HPLCDat_CA19" folder includes analyses and outputs related to QC-filtered HPLC data from the CA_2019_730PI trial. A description of the germplasm and the analyzed HPLC data can be found in Supplemental Table 2.
2) The "Analysis_2_FigsHPLCdat_CA19" folder contains analyses and outputs related to HPLC data analysis and presentation. Results of this can be found in Supplemental Tables 2,3,4. 
3) The "Analysis_3_HPLCDat_WI18" folder includes analyses and outputs related to QC-filtered HPLC data from the WI_2019_605PI trial. A description of the germplasm and the analyzed HPLC data can be found in Supplemental Table 1.
4) The "Analysis_4_FigsHPLCdat_WI18" folder contains analyses and outputs related to HPLC data analysis and presentation. Results of this can be found in Supplemental Tables 1,3,4. 
5) The "Analysis_5_DescribePhenotypeData_WI18" folder contains analyses that remove outliers from the WI_2018_605PI trial and calculate mean carotenoid levels. Mean separation for Supplemental Figure S2, S4 and the corresponding boxplot are included. Additionally, the correlation test results presented in Table 1 were conducted here. A phenotype file, utilized for subsequent Genomic Prediction and GWAS analyses, was generated during this process
6)  The "Analysis_6_DescribePhenotypeData_CA19" folder contains analyses that remove outliers from the CA_2019_730PI trial and calculate mean carotenoid levels. Mean separation for Supplemental Figure S3, S4 and the corresponding boxplot are included. Additionally, the correlation test results presented in Table 1 were conducted here. A phenotype file, utilized for subsequent Genomic Prediction and GWAS analyses, was generated during this process
7) The "Analysis_7_CA19.VCF" folder includes a sample batch file demonstrating how TASSEL was used to identify variants in the CA_2019_730PI trial. The same TASSEL pipeline was employed for the WI_2018_605PI trial.
8) The "Analysis_8_CombineandFilterVCF" folder contains analysis used to combine genotypic data from WI_2018_605PI and CA_2019_730PI trials and filter the genotypic data. There are some custom functions referred to in this file including an in-house  perl script use to wrangle and filter vcf files. 
9) The "Analysis_9_LD_PopStructure" folder contains files to complete population structure analysis with admixture, fastStructure, and calculate LD decay.
10) The Analysis_10_polyRAD" folder contains files related to variant calling and filtering in polyRAD. A genotypic file resulting from this analysis was used in GWAS and genomic prediction. Genotypic data is available at CarrotOmics.org/. 
11) The folder "Analysis_11_GWA" contains files related to completing GWAS in GAPIT. Analyses include GWAS, significance threshold calculation with simpleM and manhattan plot creation with qqman. Results from these analyses are presented in Figure 1, Figure 2,  Figure S5, Figure S6, Table S6- S10
12) The folder: "Analysis_12_HaploviewLDblocks" contains the files used in Haploview. Haploview is a GUI to visualize LD available through the Broad institute. The input format is unique, so this folder contains formatted files used to identify LD blocks to define loci and positional candidate genes TableS7 and S8. 
13) The folder: "Analysis_13_GenomicPredictionHPLC" contains example files of how BGLR was used to complete Genomic Prediction. This folder contains many scripts used to iterate through leave-one-out cross-validation. Ensure that the multi-threaded job has results from all processing units, the accuracy (correlation) was calculated, the best model was identified. The same type of analyses, sans Leave-one-out cross-validation was applied to complete across trial genomic prediction with results presented in Table S11
14) The folder: "Analysis_14_GPFigures" contains R scripts used to create figures from the GP figures. Results from Genomic Prediction are presented in Figure 3, 4 and Supplemental Figures 6,7.  
15) The folder: "Analysis_15_CarotenoidsEnvironment" contains R analysis used to complete calculations related to intra-plot, intra-field, and location differences in order to estimate the minimum effective sample to calculate an average carotenoid concentration. Results from these analyses are presented in Figure 5 and Supplemental Tables S12, and S13.    
