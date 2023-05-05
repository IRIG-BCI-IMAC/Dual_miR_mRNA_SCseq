
###########################################################################
###########################################################################
## Script for Matter Arising 
## Cumulative distribution plot
## Kolmogorov-Smirnov analysis with both type of H0 
## Distribution of H0 for different threshold
## Analysis for miR-92a-3p, miR-26a-3p, miR-125a-3p and let-7i-5p
## Parameters: Both, Expression > 4 RPKM and TCS <= 0
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


hist(mean_miRNAs)



mean_mRNA <- apply(data_RNA_19, 1, function(x) mean(x, na.rm =TRUE))
hist(mean_mRNA)
min(mean_mRNA)


exp_miRNA <- data_miRNA_19
exp_mRNA <- data_RNA_19


###########################################################################
###########################################################################
## to contradict with figure 4d of the article

pdf(file ='R.results/Supp_1a_KS_test.pdf', height = 6.5, width =7.5)

for (miRNA in list_miRNA[c(1,7,9,13)]){
  
  name_legend <- str_remove(miRNA, 'hsa-')
  
  high_genes <- names(mean_mRNA[which(mean_mRNA > 4)])
  moderate_genes <- names(mean_mRNA[which(mean_mRNA <= 4 & mean_mRNA >0)])
  weak_genes <- names(mean_mRNA[which(mean_mRNA <= 0)])
  
  length(weak_genes)+ length(high_genes)+length(moderate_genes)
  
  ## targets
  exp1 <- exp_miRNA[which(rownames(exp_miRNA) == miRNA),]
  
  
  targets <- targets_selection (miRNA, conservation = 'conserved',
                                selection = 'TCS', threshold = -0.1)
  
  
  high_targets <- targets [which(targets %in% high_genes)]
  moderate_targets <- targets [which(targets %in% moderate_genes)]
  weak_targets <- targets [which(targets %in% weak_genes)]
  
  
  ## H0 
  weak_no_targets <- names(mean_mRNA)[-which(names(mean_mRNA) %in% targets)]
  moderate_no_targets <- weak_no_targets[which(weak_no_targets %in% moderate_genes)]
  high_no_targets <- weak_no_targets[which(weak_no_targets %in% high_genes)]
  
  
  ## compute correlation
  correlations <- as.matrix(cor(t(exp_mRNA), t(exp1), 
                                method ='pearson', use = 'pairwise.complete.obs'))
  
  corr_high_targets <- correlations[which(rownames(correlations) %in% high_targets),]
  corr_moderate_targets <- correlations[which(rownames(correlations) %in% moderate_targets),]
  corr_weak_targets <- correlations[which(rownames(correlations) %in% weak_targets),]
  corr_H0_weak <- correlations[which(rownames(correlations) %in% weak_no_targets),]
  corr_H0_moderate <- correlations[which(rownames(correlations) %in% moderate_no_targets),]
  corr_H0_high <- correlations[which(rownames(correlations) %in% high_no_targets),]
  
  ## Cumulative distribution
  CDF_high_tar <- build_ecdf(corr_high_targets)
  CDF_moderate_tar <- build_ecdf(corr_moderate_targets)
  CDF_weak_tar <- build_ecdf(corr_weak_targets)
  CDF_H0_weak <- build_ecdf (corr_H0_weak)
  CDF_H0_moderate <- build_ecdf(corr_H0_moderate)
  CDF_H0_high <- build_ecdf (corr_H0_high)
  
  ## KS test
  
  ks_H0_weak <- ks.test(corr_H0_weak, corr_high_targets)
  ks_high_H0 <- ks.test(corr_H0_high, corr_high_targets)
  
  ## Plot
  plot(CDF_high_tar[[1]], CDF_high_tar[[2]] , type ='l', col ='red', 
       main = name_legend,
       xlab =paste('Correlation with', name_legend), 
       ylab ='Cummulatide distribution',
       xlim = c(-0.8, 0.8), xaxt ='n')
  axis (1, at =c(-0.8,-0.5,0,0.5,0.8))
  lines(CDF_moderate_tar[[1]], CDF_moderate_tar[[2]], type ='l',col ='green', lwd = 2)
  lines(CDF_weak_tar[[1]], CDF_weak_tar[[2]],type ='l', col ='cyan', lwd = 2)
  lines(CDF_H0_weak[[1]], CDF_H0_weak[[2]], type ='l', col = 'darkblue', lwd = 2, lty =5)
  lines(CDF_H0_moderate[[1]], CDF_H0_moderate[[2]], type ='l', col = 'darkgreen', lwd = 2, lty =5)
  lines(CDF_H0_high[[1]], CDF_H0_high[[2]], col = 'red4', lwd = 2, lty =5)
  
  legend (x ='bottomright', col = c('red','green','cyan','red4','darkgreen','darkblue'),
          bty ='n',  lty = c (1,1,1,5,5,5),
          legend = c(paste0('Targets exp > 4 (n =',length(corr_high_targets),')'),
                     paste0('Targets 0 < exp <= 4 (n =',length(corr_moderate_targets),')'),
                     paste0('Targets exp <= 0 (n =',length(corr_weak_targets),')'),
                     paste0('Non targets exp > 4 (n =',length(corr_H0_high),')'),
                     paste0('Non targets exp > 0 (n =',length(corr_H0_moderate),')'),
                     paste0('Non targets exp > -7 (n =',length(corr_H0_weak),')')))
  
  legend (x ='topleft', fill = c('darkblue', 'red4'), bty ='n',
          legend = c(paste('p =',signif(ks_H0_weak$p.value,3)),
                     paste('p =', signif(ks_high_H0$p.value,3))))


}


dev.off()


###############################################################################
## Distribution of correlation with all genes ----
## Identify targets


pdf(file ='R.results/Supp_1c_distribution_correlation_H0_by_thr.pdf', height = 6.5, width =7.5)

for (miRNA in list_miRNA[c(1,7,9,13)]){
  
  exp1 <- exp_miRNA[which(rownames(exp_miRNA) == miRNA),]
  
  name_legend <- str_remove(miRNA,'hsa-')
  
  ## predict targets 
  target_miR <- targets_selection(miRNA, conservation = 'both',
                                  selection = 'all')
  
  ## select H0 only
  exp_H0 <- exp_mRNA[-which(rownames(exp_mRNA) %in% target_miR),]
  
  print(dim(exp_H0))
  ## compute correlation for H0
  correlations <- as.matrix(cor(t(exp_H0), t(exp1), 
                                method ='pearson', use = 'pairwise.complete.obs'))
  
  ## Histogram
  hist(correlations, prob = TRUE, ylim = c(0,2), col ='white', border ='darkblue',
       main = paste('Histogram of correlations for', name_legend ,'with non targets'),
       xlab = paste0('Correlation with ',name_legend),
       xlim = c(-0.8, 0.8), xaxt ='n')
  axis(1, at = c(-0.8,-0.5, 0, 0.5, 0.8))
  
  
  vec_thr_exp <- c(-7,0,4,8)  
  vec_colors <- c('darkblue','darkgreen','red4','black')
  vec_type <- c("solid","dotdash","longdash","solid")
  
  vec_lwd <- c(0.8,2.5,2.5,1.8)
  
  for (i in 1:length(vec_thr_exp)){
    
    selected_genes <- names(mean_mRNA[which(mean_mRNA > vec_thr_exp[i])])
    selected_cor <- correlations[which(rownames(correlations) %in% selected_genes)]  
    lines(density(selected_cor), col = vec_colors[i], lwd = vec_lwd[i], lty = vec_type[i])
  }
  
  legend (x ='topleft', legend = paste ('>',vec_thr_exp), title ='Thresold on expression',
          col = vec_colors, bty ='n', lty =vec_type, lwd =1.8)
  
}

dev.off()
