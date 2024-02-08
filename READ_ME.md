# Function files

### `Functions.R`
Contains numerous functions used in the other scripts for data analysis and data visualization.

### `Function_GSEA_plot.R`
Contains a modified function to display GSEA plot with annotations. 

### `Dataset_pre_processing.R`
Script to filter and normalised Wang et al. datset. Similar scripts can be adpted for each dataset for examples for Xiao et al. dataset.

# Script files

### `Main_1a_b_Supp_4_GSEA_Conservation.R`
GSEA analysis using both conserved and non-conserved targets or conserved targets only, threshold on targets expression >  4 RPKM in log2 scale and TCS < -0.1 with Pearson correlations and log2 scale for miRNA and mRNA expressions, for Wang et al. dataset. Similar script can be used for Xiao et al. dataset. 

Outputs:  
Table of results from GSEA analysis for all miRNAs selected for analysis.  
Plots of GSEA results log10-pvalue as a function of enrichment score for each conditions.  
Plot of GSEA results +/-log10 (p-values) of one condition as a function of the other condition.  

### `Main_1c_d_Supp_8_15_GSEA_Efficacy.R`
GSEA analysis using both conserved and non-conserved targets, threshold on targets expression > 4 RPKM in log2 scale and TCS < -0.5, -0.4, -0.3, -0.2,-0.1 or  <= 0 with Pearson correlations and log2 scale for miRNA and mRNA expressions, for Wang et al. dataset. Similar script can be used for Xiao et al. dataset.  

Outputs:  
Table of results from GSEA analysis for all miRNAs selected for analysis.  
Plots of GSEA results log10-pvalue as a function of enrichment score.  
Plots of standard deviation of miRNA expression in log2 scale as a function of mean miRNA expression in log2 scale, colorized according to GSEA results.


### `Supp_1_KS_distribution_correlations.R` 
Kolmogorov-Smirnov test, Cummulative distribution function and density distribution fonction for miR-92a-3p, miR-26a-5p, miR-125a-5p and let-7i-5p. 
KS test is performed for targets with mean(log2(RPKM expression)) > 4 and with either all non targets genes for H0 or non targets genes after selection on expression using the same threshold as for targets. Box and whiskers plot of the correlations are build  for the 4 miRNAs for non targets genes only and using -7, 0, 4 or 8 as threshold on expression to select genes.

Output:  
Plots of cummulative distribution function for miR-92a-3p, miR-26a-5p, miR-125a-5p and let-7i-5p, with KS test
Box and whiskers plots of correlations for miR-92a-3p, miR-26a-5p, miR-125a-5p and let-7i-5p.


### `Supp_2a_b_P-values_simulations.R`
GSEA analysis and Kolmogorov-Smirnov test on simulated genes and sets of targets for different number of used targets.

Outputs:  
Plot with density plot of distributions of simulated genes and targets.  
Plot of log10 (pvalues) from GSEA analysis and KS test as a function of number of targets used.

 
### `Supp_2c_d_P-values_Simulations_with_different_shift.R`
GSEA analysis and Kolmogorov-Smirnov test on simulated genes and sets of targets for different number of used targets. The greater the number of targets used, the less the gap between the distribution of simulated targets compared to the distribution of simulated genes (H0)

Outputs:  
Plot with density plot of distributions of simulated genes and targets.  
Plot of log10 (pvalues) from GSEA analysis and KS test as a function of number of targets used.


### `Supp_3_Barplot_number_targets.R`
Barplot of number of conserved and non-conserved targets for each miRNAs expressed in the dataset. 
Targets are selected using two options:
- Both conserved and non-conserved targets, threshold on targets expression > -7 RPKM in log2 scale and TCS <= 0
- Both conserved and non-conserved targets, threshold on targets expression > 4 RPKM in log2 scale and TCS <= 0

Outputs:  
Plots with barplot of number of targets for each miRNAs selected for analysis. 


### `Supp_5_Box_plot_and_GSEA_plot.R`
Box and whiskers plot of correlations for miR-92a-3p for non targets, conserved targets, non conserved targets and both conserved and non conserved targets. 
GSEA plot for miR-92a-3p with conserved targets only, non-conserved targets only and both conserved and non-conserved targets. Similar script can be used for Xiao et al. dataset and let-7i-5p.

Outputs: 
Box and whiskers plot of correlations.
GSEA plots with different type of targets used. 


### `Supp_6_14_GSEA_Expression.R`
GSEA analysis using both conserved and non-conserved targets, threshold on targets expression > -7, -4, 0, 4, 6 or 8 RPKM in log2 scale and TCS < -0.1 with Pearson correlations and log2 scale for miRNA and mRNA expressions,for Wang et al. dataset. Similar script can be used for Xiao et al. dataset.  

Outputs:  
Table of results from GSEA analysis for all miRNAs selected for analysis.  
Plots of GSEA results log10-pvalue as a function of enrichment score.  


### `Supp_7_Simulation_Correlation_miRNA_targets.R`
Simulation of the repression of a miRNA on its mRNA targets for different expression levels of the mRNA, with noise restricted to sampling one. miRNA expression follows a normal distribution with mean =1 and sd =0.1. miRNA repression efficacy was set between 0.2 and 0.002. We have set the mRNA level of expression before repression to obtain an average expression after miRNA repression between -7 and 8 (log2 (RPKM expression)). For each miRNA efficacy and each mRNA mean, we ran 1000 replicates. After repression of mRNAs, we applied a Poisson law on mRNA counts to model the uncertainty due to sampling, and then computed Pearson correlations and corresponding p-values.

Outputs:   
Plot of correlations as a function of mean mRNA expression for different miRNA repression efficacy.  
Plot of a targeted mRNA expression as a function of miRNA expression as examples of correlations.   


### `Supp_9a_b_GSEA_H0.R` 
GSEA analysis using both conserved and non-conserved targets, threshold on targets expression > 4 RPKM in log2 scale and TCS <= 0. With either all non targets genes used as H0 (Null hypothesis) or genes selected if their mean log2 expressio is above 4 RPKM for H0, for Wang et al. dataset. Similar script can be used for Xiao et al. dataset.  

Output:  
Plot of GSEA +/-log10 (p-values) of one condition as a function of the other condition.  


### `Supp_10a_c_Histogram_jaccard_index.R`
Intersections between the set of predicted targets of each miRNAs with the set of core targets of hsa-miR-92a-3p are compute and compare with random intersection expected. Targets are selected using different options: 
- Both conserved and non-conserved targets, threshold on targets expression > 4 RPKM in log2 scale and TCS <= 0 
- Conserved targets, threshold on targets expression > 4 RPKM in log2 scale and TCS < -0.1 
- Conserved targets, threshold on targets expression > -7 RPKM in log2 scale and TCS < -0.1 

Outputs:   
Plot with histogram of percentage of common targets between miRNAs and core targets of hsa-miR-92a-3p and histogram of random percentage of intersection expected (H0). 


### `Supp_10b_GSEA_Negative_Control.R`
GSEA analysis with correlations computed using the expression of hsa-miR-92a-3p in every cases, the targets selection is made for every miRNAs using different options and with Pearson correlations and log2 scale for miRNA and mRNA expression. 
- Conserved targets only, threshold on targets expression > -7 RPKM in log2 scale and TCS < -0.1
- Both conserved and non-conserved targets, threshold on targets expression > 4 RPKM in log2 scale and TCS <=0

Outputs:   
Table of results from GSEA analysis for all miRNAs selected for analysis.  
Plots of GSEA results log10-pvalue as a function of enrichment score.


### `Supp_10d_GSEA_simulation_Jaccard_core.R`
GSEA analysis on simulated sets of 588 targets which have a Jaccard index between 0 and 0.07 with miR-92a-3p core targets. For each Jaccard index, the simulated sets of targets were generated 100 times and the p-value of GSEA analysis were plot as boxplots.

Outputs:  
Boxplot of GSEA log10 (pvalues) obtained for miR-92a-3p correlations with simulated sets of 570 targets with varying Jaccard index with the
miR-92a-3p target core. 


### `Supp_11a_b_GSEA_shuffle_miR.R`
GSEA analysis with correlations computed using the expression of hsa-miR-92a-3p in every cases, the targets selection is made for 100 artificials miRNAs obtain by shuffling miR-92a-3p sequence or seed sequence. Correlations computed with Pearson correlations and log2 scale for miRNA and mRNA expression. 

Outputs:   
Table of results from GSEA analysis for artificial miRNAs.  
Plots of GSEA results log10-pvalue as a function of enrichment score for miRNA obtain by shuffling miR-92a-3p sequence and for miRNA obtain by shuffling miR-92a-3p seed sequence.


### `Supp_11c_d_Targetome_and_UTR_length.R`
Targetome and 3'UTR length relationship analysis. 
Analysis of 3'UTR length of genes. Genes are split into two groups, one with genes less expressed in cells in which miR-92a-3p is highly expressed. The other one being genes more expressed in cells in which miR-92a-3p is weakly expressed.   

Outputs:   
Plot of number of miR trageting each genes as a function of 3'UTR length of these genes.
Cummulative distribution of 3'UTR length of genes represented with two groups. One group of gene have a positive log fold change between the cells in which miR-92a-3p is highly expressed and the cells in which miR-92a-3p is less expressed. The other group having a negative log flod change. 


### `Supp_12_GSEA_KS_Wilcoxon_test.R`
GSEA analysis, Kolmogorov-Smirnov (KS) and Wilcoxon tests using different options on targets selection with Pearson correlations and log2 scale for miRNA and mRNA expression: 
- Both conserved and non-conserved targets, threshold on targets expression > 4 RPKM in log2 scale and TCS <= 0 
- Conserved targets, threshold on targets expression > 4 RPKM in log2 scale and TCS < -0.1 

Outputs:   
Table of results from GSEA analysis, KS and Wilcoxon tests for all miRNAs selected for analysis.  
Plots of GSEA, KS and Wilcoxon tests results log10-pvalue as a function of enrichment score or D (compute with a sign).  
Plot of GSEA and KS results +/-log10 (p-values) of one condition as a function of the other condition. 
Plot of GSEA and Wilcoxon results +/-log10 (p-values) of one condition as a function of the other condition.  


### `Supp_13_GSEA_miR-92a-3p_correlation.R`
Use GSEA results from Main_1a_Supp_4_GSEA_Conservation.R in the conditions with both conserved and non conserved targets, thrshold on targets expression > 4 and TCS <= 0 and compute the correlations of expression between all miRNAs and miR-92a-3p. similar script can be used for Xiao et al. dataset with let-7i-5p.  

Outputs:
Plot of GSEA results by correlation between the expression of miRNAs and miR-92a-3p. 


### `Supp_16_GSEA_Scale.R`
GSEA analysis using both conserved and non-conserved targets, threshold on targets expression > 4 RPKM in log2 scale and TCS <= 0 with Pearson correlations. Using either log2 scale for miRNA and mRNA expressions or linear scale. 

Outputs:  
Table of results from GSEA analysis for all miRNAs selected for analysis.  
Plots of GSEA results log10-pvalue as a function of enrichment score.  
Plot of GSEA results +/-log10 (p-values) of one condition as a function of the other condition.    
Plots of standard deviation of miRNA expression in log2 scale as a function of mean miRNA expression in log2 scale, colorized according to GSEA results.


### `Supp_17_GSEA_Correlation.R`
GSEA analysis using both conserved and non-conserved targets, threshold on targets expression > 4 RPKM in log2 scale and TCS <= 0. Using either Pearson coefficient for computing correlations between miRNA expression and mRNA expression or Spearman coefficient and log2 scale for miRNA and mRNA expression. 

Outputs:  
Table of results from GSEA analysis for all miRNAs selected for analysis.  
Plots of GSEA results log10-pvalue as a function of enrichment score.  
Plot of GSEA results +/-log10 (p-values) of one condition as a function of the other condition.  
Plots of standard deviation of miRNA expression in log2 scale as a function of mean miRNA expression in log2 scale, colorized according to GSEA results.
