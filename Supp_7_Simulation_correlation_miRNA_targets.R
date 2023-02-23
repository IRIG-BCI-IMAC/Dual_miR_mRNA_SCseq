###########################################################################
## Script for Matter Arising
## Simulation of repression of a miRNA on its targets
## targets with different expression
## miRNA with different efficacy of repression
###########################################################################

set.seed (0)

source ('Functions.R')

nb_cells <- 19

tested_efficacy <- c(1e-1, 1e-2, 1e-3)*2 # miRNA efficacy

l_exp <- c(seq(-7, 8, by = 1)) # target expression log2 scale 
tae <- 2^(l_exp) * 1e3 # target average expression



miRNA_exp <- rnorm(19, mean = 1, sd = 0.1) 
boxplot(miRNA_exp)
mean(miRNA_exp)

results_correlation <- as.data.frame(matrix(ncol = length(l_exp),
                             nrow = 3))# matrix for cor median results
results_p_values <- as.data.frame(matrix(ncol = length(l_exp),
                             nrow = 3))# matrix for p_values median results

results_mean <- as.data.frame(matrix(ncol = length(l_exp),
                             nrow = 3))# matrix for p_values median results

for (i in 1:length(tested_efficacy)){
efficacy <- tested_efficacy[i] # miRNA efficacy to repress mRNA targets

repression <- miRNA_exp * efficacy
repression # repression is proportional to miRNA expression




for(j in 1:length(tae)){
  print(paste('mRNA mean expression =', tae[j]))
  ## create empty vectors
  vec_corr <- c()
  vec_pval <- c()
  vec_mean <- c() 
  
  for (x in 1:1000){
    ## loop for 200 replicates
  mRNAf <- rep(tae[j],19)
  mRNAf
  
  mRNAi <- (mRNAf / efficacy) / ((1-efficacy) / efficacy )
  
  repressed_mRNA <- mRNAi - mRNAi * repression
  
  ## apply poisson distribution on mRNA expression after repression
  vec_counts<- c()
  for ( value in repressed_mRNA){
    counts <- rpois(1, value)
    vec_counts <- c(vec_counts, counts)
  }

  ## compute correlation and p-value
  correlation <- signif(cor (miRNA_exp, vec_counts),2)
  pval <- signif(cor.test(miRNA_exp, vec_counts)$p.value, 2)
  mean_mRNA <- mean(vec_counts)
  
  ## store results for replicates
  vec_corr <- c(vec_corr, correlation)
  vec_pval <- c(vec_pval, mean_mRNA)
  vec_mean <- c(vec_mean, mean_mRNA)
  
  ## print 2 examples of correlations
  if( efficacy %in% c(0.02,0.2) && j == 12 && x == 499){

    file_output1 <- paste('R.results/Supp_7',
                          ifelse(efficacy == 0.02,'b','c'),
                          '_correlation_miRNA_targets.pdf', sep = '')
    pdf(file_output1, width = 6.5, height = 6.5)
    Reg <- lm(vec_counts ~ miRNA_exp)
    plot(miRNA_exp, vec_counts, 
         main = paste('miRNA efficacy =', efficacy,
                      '\n corr =', correlation,'- pval =',pval),
         xlab ='miRNA expression', 
         ylab =paste('mRNA expression'))
    legend (x ='bottomleft', legend = c(paste('corr =', correlation),
                                  paste('p-val =', pval)), bty ='n')
    abline(Reg, col='cyan4')
    
    
    dev.off()
    
    
  }
  
  }
  
  results_correlation[i,j] <- median (vec_corr) 
  results_p_values[i,j] <- median (vec_pval)
  results_mean[i,j] <- median (vec_mean)
}
}

results_correlation
results_mean


###########################################################################
## plot results correlation by mean mRNA expression
file_output2 <- 
  paste('R.results/Supp_7a_Simu_correlation_miRNA_targets.pdf', sep = '')
pdf(file_output2, width = 8, height = 6.5)

vec_col <- c('black', 'red','blue')
plot(1,1, xlim = c(-7,8), ylim =c(-1,0.05), type ='l',
     xlab ='mean (log2 (mRNA expression))',
     ylab = 'Pearson coefficient correlation',
     xaxt ='n')
axis(1, at = -7:8)
grid()
line_col1 <- rgb(0, 0.6, 0.4,0.3)
abline(v = c(0,4), col = line_col1, lty = 2, lwd = 3)

mean_RPKM <- log2(results_mean/1e3) 

for (i in 1:dim(results_correlation)[1]){
  lines(mean_RPKM[i,], results_correlation[i,], type = 'l', 
        col = vec_col[i], lwd = 2)
}

legend( x ='bottomleft', legend = tested_efficacy, fill = vec_col[1:3],
        title ='miRNA repression efficacy', bg ='white')


dev.off()

