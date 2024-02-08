###########################################################################
## Script for Matter Arising 
## Cumulative distribution plot
## Kolmogorov-Smirnov analysis with both type of H0 
## Distribution of H0 for different threshold
## Analysis for miR-92a-3p, miR-26a-3p, miR-125a-3p and let-7i-5p
## Parameters: Both, Expression > 4 RPKM and TCS <= 0
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
hist(mean_miRNAs)


###########################################################################
###########################################################################
## to contradict with figure 4d of the article

pdf(file ='R.results/Supp_1abef_KS_test.pdf', height = 6.5, width =7.5)

for (miRNA in list_miRNA[c(1,7,9,13)]){
  
  name_legend <- str_remove(miRNA, 'hsa-')
  
  high_genes <- names(mean_mRNAs[which(mean_mRNAs > 4)])
  moderate_genes <- names(mean_mRNAs[which(mean_mRNAs <= 4 & mean_mRNAs >0)])
  weak_genes <- names(mean_mRNAs[which(mean_mRNAs <= 0)])
  
  length(weak_genes)+ length(high_genes)+length(moderate_genes)
  
  ## targets
  exp1 <- data_miRNA[which(rownames(data_miRNA) == miRNA),]
  
  
  targets <- targets_selection (miRNA, conservation = 'conserved',
                                selection = 'TCS', threshold = -0.1)
  
  
  high_targets <- targets [which(targets %in% high_genes)]
  moderate_targets <- targets [which(targets %in% moderate_genes)]
  weak_targets <- targets [which(targets %in% weak_genes)]
  
  
  ## H0 
  no_targets <- names(mean_mRNAs)[-which(names(mean_mRNAs) %in% targets)]
  high_no_targets <- no_targets[which(no_targets %in% high_genes)]
  
  
  ## compute correlation
  correlations <- as.matrix(cor(t(data_RNA), t(exp1), 
                                method ='pearson', use = 'pairwise.complete.obs'))
  
  corr_high_targets <- correlations[which(rownames(correlations) %in% high_targets),]
  corr_moderate_targets <- correlations[which(rownames(correlations) %in% moderate_targets),]
  corr_weak_targets <- correlations[which(rownames(correlations) %in% weak_targets),]
  corr_H0_all <- correlations[which(rownames(correlations) %in% no_targets),]
  corr_H0_high <- correlations[which(rownames(correlations) %in% high_no_targets),]
  
  ## Cumulative distribution
  CDF_high_tar <- build_ecdf(corr_high_targets)
  CDF_moderate_tar <- build_ecdf(corr_moderate_targets)
  CDF_weak_tar <- build_ecdf(corr_weak_targets)
  CDF_H0_all <- build_ecdf (corr_H0_all)
  CDF_H0_high <- build_ecdf (corr_H0_high)
  
  ## KS test
  ks_all_H0 <- ks.test(corr_H0_all, corr_high_targets)
  ks_high_H0 <- ks.test(corr_H0_high, corr_high_targets)
  
  ## Plot
  plot(CDF_high_tar[[1]], CDF_high_tar[[2]] , type = 'l', col = 'red', lwd = 2, 
       main = name_legend,
       xlab = paste('Correlation with', name_legend), 
       ylab ='Cumulative distribution',
       xlim = c(-0.8, 0.8), xaxt = 'n')
  axis (1, at = c(-0.8,-0.5,0,0.5,0.8))
  lines(CDF_moderate_tar[[1]], CDF_moderate_tar[[2]], type ='l',col ='green', lwd = 2)
  lines(CDF_H0_all[[1]], CDF_H0_all[[2]], type ='l', col = 'darkblue', lwd = 2, lty =5)
  lines(CDF_H0_high[[1]], CDF_H0_high[[2]], col = 'orange', lwd = 2, lty =5)
  
  legend (x ='bottomright', col = c('red','green','orange','darkblue'),
          bty ='n',  lty = c (1,1,5,5),
          legend = c(paste0('Targets exp > 4 (n = ',length(corr_high_targets),')'),
                     paste0('Targets 0 < exp <= 4 (n = ',length(corr_moderate_targets),')'),
                     paste0('Non targets exp > 4 (n = ',length(corr_H0_high),')'),
                     paste0('Non targets exp > -7 (n = ',length(corr_H0_all),')')))
  
  legend (x ='topleft', fill = c('darkblue', 'orange'), bty ='n',
          legend = c(paste('p =',signif(ks_all_H0$p.value,3)),
                     paste('p =', signif(ks_high_H0$p.value,3))))


}


dev.off()


###############################################################################
## Distribution of correlation with all genes ----
## Identify targets

## Boxplot of distribution
pdf('R.results/Supp_1_bdgh_box_plot_distribution_shifted_H0.pdf', height =6,width =9)
for (miRNA in list_miRNA[c(1,7,9,13)]){
  
  exp1 <- data_miRNA[which(rownames(data_miRNA) == miRNA),]
  
  name_legend <- str_remove(miRNA,'hsa-')
  
  ## predict targets 
  target_miR <- targets_selection(miRNA, conservation = 'both',
                                  selection = 'all')
  
  ## select H0 only
  exp_H0 <- data_RNA[-which(rownames(data_RNA) %in% target_miR),]
  
  print(dim(exp_H0))
  ## compute correlation for H0
  correlations <- as.matrix(cor(t(exp_H0), t(exp1), 
                                method ='pearson', use = 'pairwise.complete.obs'))
  
  ## Histogram
  list_cor <- list()
  vec_thr_exp <- c(-7,-4,0,4,6) 
  
  for (i in 1:length(vec_thr_exp)){
    
    selected_genes <- names(mean_mRNAs[which(mean_mRNAs > vec_thr_exp[i])])
    selected_cor <- correlations[which(rownames(correlations) %in% selected_genes)]  
    list_cor[[length(list_cor) +1]] <- c(selected_cor)
  }
  
  list_cor
  
  boxplot(list_cor, main =paste('Correlations between',miRNA, 'and all non targets genes'),
          xaxt ='n', xlab ='Threshold applied on expression',ylab ='Pearson correlations', 
          ylim = c(-0.8,0.8))
  text(x = seq(1:length(vec_thr_exp)),
       y = par("usr")[3] - 0.1,
       labels = vec_thr_exp,
       xpd = NA,
       srt = 0,
       adj = 0.7)
  abline(h =median(list_cor[[1]]), col ='red')
  
  
}

dev.off()
