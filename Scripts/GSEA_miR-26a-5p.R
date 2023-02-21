###########################################################################
## Script for Matter Arising 
## GSEA plot for hsa-miR-26a-5p
## Parameters : Conservation = conserved or both, 
## Expression > -7 or 4 and TCS < -0.1
###########################################################################

## Importation 
source("Functions.R")# functions importation

## Data importation -------------------------------------------------------
## TargetScan v7.1 data importation 
if (!exists('TS'))
  TS <- TargetScan_importation()
TargetScan <- TS[[1]] ; families <- TS[[2]]

## Single-cell data importation
data_imported <- import_SCdata()
data_RNA_19 <- data_imported[[1]] ; data_miRNA_19 <- data_imported[[2]] 

## Sort miRNAs by mean expression 
mean_miRNAs <- apply(data_miRNA_19,1,mean)
all_miRNAs <- names(sort(mean_miRNAs[mean_miRNAs > -13],decreasing = TRUE))
list_miRNA <- all_miRNAs [-which(all_miRNAs %in% 
                                   c("hsa-miR-183-5p","hsa-miR-140-3p"))]

###########################################################################
## GSEA analysis on hsa-miR-26a-5p ----------------------------------------
###########################################################################


## GSEA analysis ----------------------------------------------------------
## Plot GSEA for miR-26-a-5p in article conditions
file1 <- "R.results/Supp_1a_GSEA_miR-26a-5p_article.pdf"
pdf (file1, width = 14.07, height = 8)

apply_GSEA (miRNA = 'hsa-miR-26a-5p', conservation ='conserved',
            selection = 'TCS', threshold = -0.1, thr_exp= -7, 
            display = TRUE)

dev.off()




## Plot GSEA for miR-2-a-5p in article conditions with highly
file2 <- "R.results/Supp_1b_GSEA_miR-26a-5p_article_highly.pdf"
pdf (file2, width = 14.07, height = 8)

apply_GSEA (miRNA = 'hsa-miR-26a-5p', conservation ='conserved',
            selection = 'TCS', threshold = -0.1, thr_exp= 4, 
            display = TRUE)

dev.off()




## Plot GSEA for miR-26-a-5p in both conditions with highly
file3 <- "R.results/Supp_1c_GSEA_miR-26a-5p_both_highly.pdf"
pdf (file3, width = 14.07, height = 8)

apply_GSEA (miRNA = 'hsa-miR-26a-5p', conservation ='both',
            selection = 'TCS', threshold = -0.1, thr_exp= 4, 
            display = TRUE)

dev.off()






