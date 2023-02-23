# Function files

### `Functions.R`
Contains numerous functions used in the other scripts for data analysis and data visualization.

### `Function_GSEA_plot.R`
Contains a modified function to display GSEA plot with annotations. 

# Script files

### `Main_1a_Supp_4_GSEA_Conservation.R`
GSEA analysis using both conserved and non-conserved targets or conserved targets only, threshold on targets expression >  4 RPKM in log2 scale and TCS < -0.1 with Pearson correlations and log2 scale for miRNA and mRNA expressions. 

Outputs: 
Table of results from GSEA analysis for all miRNAs selected for analysis.
Plots of GSEA results log10-pvalue as a function of enrichment score for each conditions
Plot of GSEA results +/-log10 (p-values) of one condition as a function of the other condition.  

### `Main_1b_Supp_6_GSEA_Efficacy.R`
GSEA analysis using both conserved and non-conserved targets, threshold on targets expression > 4 RPKM in log2 scale and TCS < -0.5, -0.4, -0.3, -0.2,-0.1 or  <= 0 with Pearson correlations and log2 scale for miRNA and mRNA expressions. 

Outputs:
Table of results from GSEA analysis for all miRNAs selected for analysis.
Plots of GSEA results log10-pvalue as a function of enrichment score.
Plots of standard deviation of miRNA expression in log2 scale as a function of mean miRNA expression in log2 scale, colorized according to GSEA results.

### `Supp_1_GSEA_miR-26a-5p.R` 
GSEA analysis of hsa-miR-26a-5p using different conditions on targets selection with Pearson correlation and log2 scale for miRNA and mRNA expression. 

Output: GSEA plots from GSEA analysis of hsa-miR-26a-5p using different conditions of targets selection. 

### `Supp_2a_b_P-values_simulations.R`
GSEA analysis and Kolmogorov-Smirnov test on simulated genes and sets of targets for different number of used targets.

Outputs:
Plot with density plot of distributions of simulated genes and targets
Plot of log10 (pvalues) from GSEA analysis and KS test as a function of number of targets used
 
### `Supp_2c_d_P-values_Simulations_with_different_shift.R`
GSEA analysis and Kolmogorov-Smirnov test on simulated genes and sets of targets for different number of used targets. The greater the number of targets used, the less the gap between the distribution of simulated targets compared to the distribution of simulated genes (H0)

Outputs:
Plot with density plot of distributions of simulated genes and targets
Plot of log10 (pvalues) from GSEA analysis and KS test as a function of number of targets used


### `Supp_3_Barplot_number_targets.R`
Barplot of number of conserved and non-conserved targets for each miRNAs expressed in the dataset. 
Targets are selected using two options:
- Both conserved and non-conserved targets, threshold on targets expression > -7 RPKM in log2 scale and TCS <= 0
- Both conserved and non-conserved targets, threshold on targets expression > 4 RPKM in log2 scale and TCS <= 0

Outputs:
Plots with barplot of number of targets for each miRNAs selected for analysis. 

### `Supp_5_GSEA_Expression.R`
GSEA analysis using both conserved and non-conserved targets, threshold on targets expression > -7, -4, 0, 4, 6 or 8 RPKM in log2 scale and TCS < -0.1 with Pearson correlations and log2 scale for miRNA and mRNA expressions. 

Outputs:
Table of results from GSEA analysis for all miRNAs selected for analysis.
Plots of GSEA results log10-pvalue as a function of enrichment score.
Plots of standard deviation of miRNA expression in log2 scale as a function of mean miRNA expression in log2 scale, colorized according to GSEA results.


### `Supp_7_Simulation_Correlation_miRNA_targets.R`
Simulation of the repression of a miRNA on its mRNA targets for different expression levels of the mRNA, with noise restricted to sampling one. miRNA expression follows a normal distribution with mean =1 and sd =0.1. miRNA repression efficacy was set between 0.2 and 0.002. We have set the mRNA level of expression before repression to obtain an average expression after miRNA repression between -7 and 8 (log2 (RPKM expression)). For each miRNA efficacy and each mRNA mean, we ran 1000 replicates. After repression of mRNAs, we applied a Poisson law on mRNA counts to model the uncertainty due to sampling, and then computed Pearson correlations and corresponding p-values.

Outputs: 
Plot of correlations as a function of mean mRNA expression for different miRNA repression efficacy.
Plot of a targeted mRNA expression as a function of miRNA expression as examples of correlations.   

### `Supp_8a_c_Histogram_jaccard_index.R`
Intersections between the set of predicted targets of each miRNAs with the set of core targets of hsa-miR-92a-3p are compute and compare with random intersection expected. Targets are selected using different options: 
- Both conserved and non-conserved targets, threshold on targets expression > 4 RPKM in log2 scale and TCS <= 0 
- Conserved targets, threshold on targets expression > 4 RPKM in log2 scale and TCS < -0.1 
- Conserved targets, threshold on targets expression > -7 RPKM in log2 scale and TCS < -0.1 

Outputs: 
Plot with histogram of percentage of common targets between miRNAs and core targets of hsa-miR-92a-3p and histogram of random percentage of intersection expected (H0). 

### `Supp_8b_GSEA_Negative_Control.R`
GSEA analysis with correlations computed using the expression of hsa-miR-92a-3p in every cases, the targets selection is made for every miRNAs using different options and with Pearson correlations and log2 scale for miRNA and mRNA expression. 
- Conserved targets only, threshold on targets expression > -7 RPKM in log2 scale and TCS < -0.1
- Both conserved and non-conserved targets, threshold on targets expression > 4 RPKM in log2 scale and TCS <=0

Outputs: 
Table of results from GSEA analysis for all miRNAs selected for analysis.
Plots of GSEA results log10-pvalue as a function of enrichment score.

### `Supp_8d_GSEA_simulation_Jaccard_core.R`
GSEA analysis on simulated sets of 588 targets which have a Jaccard index between 0 and 0.07 with miR-92a-3p core targets. For each Jaccard index, the simulated sets of targets were generated 100 times and the p-value of GSEA analysis were plot as boxplots.

Outputs:
Boxplot of GSEA log10 (pvalues) obtained for miR-92a-3p correlations with simulated sets of 588 targets with varying Jaccard index with the
miR-92a-3p target core. 

### `Supp_9_GSEA_and_KS_test.R`
GSEA analysis and Kolmogorov-Smirnov (KS) tests using different options on targets selection with Pearson correlations and log2 scale for miRNA and mRNA expression: 
- Both conserved and non-conserved targets, threshold on targets expression > 4 RPKM in log2 scale and TCS <= 0 
- Conserved targets, threshold on targets expression > 4 RPKM in log2 scale and TCS < -0.1 

Outputs: 
Table of results from GSEA analysis and KS tests for all miRNAs selected for analysis.
Plots of GSEA and KS test results log10-pvalue as a function of enrichment score or D (compute with a sign).
Plot of GSEA and KS results +/-log10 (p-values) of one condition as a function of the other condition.  

### `Supp_10_GSEA_Scale.R`
GSEA analysis using both conserved and non-conserved targets, threshold on targets expression > 4 RPKM in log2 scale and TCS <= 0 with Pearson correlations. Using either log2 scale for miRNA and mRNA expressions or linear scale. 

Outputs:
Table of results from GSEA analysis for all miRNAs selected for analysis.
Plots of GSEA results log10-pvalue as a function of enrichment score.
Plot of GSEA results +/-log10 (p-values) of one condition as a function of the other condition.  
Plots of standard deviation of miRNA expression in log2 scale as a function of mean miRNA expression in log2 scale, colorized according to GSEA results.

### `Supp_11_GSEA_Correlation.R`
GSEA analysis using both conserved and non-conserved targets, threshold on targets expression > 4 RPKM in log2 scale and TCS <= 0. Using either Pearson coefficient for computing correlations between miRNA expression and mRNA expression or Spearman coefficient and log2 scale for miRNA and mRNA expression. 

Outputs:
Table of results from GSEA analysis for all miRNAs selected for analysis.
Plots of GSEA results log10-pvalue as a function of enrichment score.
Plot of GSEA results +/-log10 (p-values) of one condition as a function of the other condition.
Plots of standard deviation of miRNA expression in log2 scale as a function of mean miRNA expression in log2 scale, colorized according to GSEA results.
