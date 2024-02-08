###########################################################################
## Script for Matter Arising
## GSEA analysis on hsa-miR-92a-3p
## GSEA plot for conserved targets, non conserved targets
## or both conserved and non conserved targets
## Box plot of correlations
###########################################################################


## Importation ------------------------------------------------------------
source("Functions.R")# functions importation

## Data importation -------------------------------------------------------
## TargetScan v7.1 data importation 
if(!exists('TS'))
  TS <- TargetScan_importation()
TargetScan <- TS[[1]] ; families <- TS[[2]]

## Single-cell data importation
data_RNA <- readRDS("R.Data/data_RNA_19.rds") 
data_miRNA <- readRDS("R.Data/data_miRNA_19.rds") 

## Sort miRNAs by mean expression 
mean_miRNAs <- apply(data_miRNA,1,mean)
list_miRNA <- names(sort(mean_miRNAs[mean_miRNAs > -13],decreasing = TRUE))
mean_mRNAs <- apply(data_RNA,1,function(x) mean(x, na.rm =T))
hist(mean_mRNAs)


###########################################################################
## Apply GSEA for miR-92a-3p

miR <- 'hsa-miR-92a-3p'

## save plot in a pdf file
file_output <- 'R.results/Supp_5_bcd_miR-92a-3p_plot_Conservation.pdf'
pdf(file_output, width = 14.07, height = 8)


## GSEA with conserved targets only
apply_GSEA(miR, conservation = 'conserved', thr_exp = 4,
           threshold = 0, display = TRUE)

## GSEA with non conserved targets only
apply_GSEA(miR, conservation = 'non conserved', thr_exp = 4,
           threshold = 0, display = TRUE)

## GSEA with both conserved and non conserved targets
apply_GSEA(miR, conservation = 'both', thr_exp = 4,
           threshold = 0, display = TRUE)

dev.off()


###########################################################################
## Boxplot correlation miR and mRNAs
miRNA <- 'hsa-miR-92a-3p'
thr_exp <- 4

exp_miR <- data_miRNA[which(rownames(data_miRNA) == miRNA ),]
genes_selected <- names(mean_mRNAs[which(mean_mRNAs > thr_exp)])


targets_cons <- targets_selection(miRNA, conservation = 'conserved', 
                                  selection = 'all')
targets_noncons <- targets_selection(miRNA, conservation = 'non conserved', 
                                  selection = 'all')
targets_both <- targets_selection(miRNA, conservation = 'both', 
                             selection = 'all')

non_targets <- genes_selected[-which(genes_selected %in% targets_both)]

## select targets using threshold on expression 
tar_cons <- targets_cons[which(targets_cons %in% genes_selected)]
tar_noncons <- targets_noncons[which(targets_noncons %in% genes_selected)]
tar_both <- targets_both[which(targets_both %in% genes_selected)]

## compute correlations
cor_non_targets <- t(cor(t(exp_miR), t(data_RNA[which(rownames(data_RNA) %in% non_targets),]), 
                         method ='pearson', use = 'pairwise.complete.obs'))
cor_targets_cons <- t(cor(t(exp_miR), t(data_RNA[which(rownames(data_RNA) %in% tar_cons),]), 
                          method ='pearson', use = 'pairwise.complete.obs'))
cor_targets_noncons <- t(cor(t(exp_miR), t(data_RNA[which(rownames(data_RNA) %in% tar_noncons),]), 
                          method ='pearson', use = 'pairwise.complete.obs'))
cor_targets_both <- t(cor(t(exp_miR), t(data_RNA[which(rownames(data_RNA) %in% tar_both),]), 
                     method ='pearson', use = 'pairwise.complete.obs'))


l_non_targets <- length(cor_non_targets)
l_cons <- length(cor_targets_cons)
l_noncons <- length(cor_targets_noncons)
l_both <- length(cor_targets_both)

vec_length <- c(l_non_targets, l_cons, l_noncons, l_both)

## save plot in a pdf file
file_output2 <- 'R.results/Supp_5a_miR-92a-3p_box_plot_correlations_conservation.pdf'
pdf(file_output2, width = 7, height = 7)

colors <- c('#E0E0E0','#CCE5FF','#FFE5CC','#CCCCFF')
boxplot(list(cor_non_targets, cor_targets_cons, cor_targets_noncons, cor_targets_both), frame =FALSE,
        main =paste('Correlations between miR-92a-3p and mRNAs \nwith mean exprssion >', thr_exp),
        xaxt ='n', ylab ='Correlations Pearson', col =colors, ylim = c(-1,1.2))
labels =c('Non targets','conserved targets', 
          'non conserved targets',' both conserved and non \n conserved targets')
abline(h =median(cor_non_targets), col ='red')
text(x = seq(1:4),
     y = par("usr")[3] - 0.1,
     labels = labels,
     xpd = NA,
     srt = 20,
     adj = 0.7)

text(x = seq(1:4),
     y = -0.1,
     labels = paste('n =', vec_length),
     xpd = NA,
     adj = 0.5)

w1 <- wilcox.test(cor_non_targets, cor_targets_cons, 
                  eps = 1e-30, alternative = "two.sided")
w2 <- wilcox.test(cor_non_targets, cor_targets_noncons,  
                  eps = 1e-30, alternative = "two.sided")
w3 <- wilcox.test(cor_non_targets, cor_targets_both,  
                  eps = 1e-30, alternative = "two.sided")

legend (x= 'topleft', title = 'Wilcoxon test', bty ='n',
        legend = paste('p =',c(s3(w1$p.value), s3(w2$p.value), s3(w3$p.value))), 
        fill = colors[2:4])


dev.off()

