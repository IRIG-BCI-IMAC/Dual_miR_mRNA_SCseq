###########################################################################
## Script for Matter Arising 
## GSEA analysis and Kolmogorov-Smirnov test
## Apply on simulation for different size of set of targets
## Results are the median on 500 replicates 
###########################################################################


## Importation
source("Functions.R")# functions importation


###########################################################################
## Genes and targets simulation -------------------------------------------
###########################################################################

## Simulation on artificial distribution
nb_genes <- 20000
nb_targets <- 2000


## Plot simulation distribution
file_output1 <- 'R.results/Supp_2a_Simulation_distribution.pdf'
pdf (file_output1)

x_genes <- seq(-1, 1, length = nb_genes )
y_genes <- dnorm(x_genes, mean = 0, sd = 0.2 )
plot(x_genes,y_genes, type = "l", lwd = 4, 
     xlab = 'Pearson correlation coefficient', ylab = 'Density',
     col = 'grey40')

x_tar <- seq(-1, 1, length = nb_targets )
y_tar <- dnorm(x_tar, mean = -0.03, sd = 0.2 )
lines(x_tar,y_tar, type = "l", lwd = 4, xlab = "", ylab = "", 
      col = "#FB9A99")
legend (x = 'topleft', legend = c(paste('Genes N =',nb_genes ),
                                  paste('Targets N =', nb_targets)), 
        fill = c('grey40', "#FB9A99"))

dev.off()




## Cumulative Distribution Function
simu_genes <- rnorm (nb_genes, mean = 0, sd = 0.2)
simu_targets <- rnorm (nb_targets, mean = -0.03, sd = 0.2)
CDF_genes <- ecdf(simu_genes)
CDF_targets <- ecdf(simu_targets)

plot (CDF_genes, col='grey', main= 'ECDF(x)',
      ylab='Cumulative distribution ', lwd=3)
plot(CDF_targets, col='orange',add=T, lwd=3)
legend(x = 'left', col =c('grey','orange'), cex=0.9, 
       lwd ='3', legend =c('Genes N = 20 000','Targets N = 2 000'))



## Number of targets used for the test
vec_nb_target <- c(100,250,500,750,1000,1500,2000)

###########################################################################
## KS test study ----------------------------------------------------------
###########################################################################

res_KS <- as.matrix(t(vec_nb_target)) 

for (j in 1:500){
  vec_KS_p <- c()
  for (nb in vec_nb_target){
    used_target <- rnorm (nb, mean = -0.03, sd = 0.2)
    resy_KS <- ks.test(used_target, simu_genes, eps = 1e-30, 
                       alternative = 'two.sided')
    p <- resy_KS[[2]]
    vec_KS_p <- c(vec_KS_p,p)
  }
  
  res_KS <- rbind(res_KS,vec_KS_p)
  
}

## use first line as column names
matrix_KS <- renamecols (res_KS)
matrix_KS

median_KS <- apply(matrix_KS,2,median)



plot(vec_nb_target,log10(median_KS), 
     pch=19, type="b", cex.lab = 1.2,
     col = 'seagreen', xaxt = 'n',
     xlab='Nombre de cibles utilis?es', 
     ylab = 'log10(p-valeur)',
     main = paste('log10(p-valeur) obtenue avec le test KS\n',
                  'en fonction du nombre de cibles utilisees'))
abline (h=log10(0.05),col='lightblue')
axis(1, at=vec_nb_target, labels=vec_nb_target)




## Store results
## Create file

file_table <- 'R.results/P-values_simulation_table_test.xlsx'
name_wb <- 'Simulation results'
wb <- createWorkbook(name_wb)
sheet <- 'KS results'
addWorksheet(wb, sheet)

writeData(wb, sheet, matrix_KS)
saveWorkbook(wb, file_table, overwrite = TRUE)

###########################################################################
## GSEA study -------------------------------------------------------------
###########################################################################
res_GSEA <- matrix(vec_nb_target) 

for (j in 1:500){
  print(paste("################",j))
  vec_GSEA_p <- c()
  for (nb in vec_nb_target){
    print(nb)
    used_target <- rnorm (nb, mean = -0.03, sd = 0.2)
    
    
    
    vec <- c()
    for(k in 1:nb_genes){vec <- c(vec, paste('genex',k, sep=''))}
    names(simu_genes) <- vec
    
    
    #######################################################################
    ## For y 
    vec <- c()
    for(i in 1:nb){vec <- c(vec, paste('geney',i, sep=''))}
    names(used_target) <- vec
    
    
    ## Build term2gene
    vec_legend <- rep('test_y', length(used_target))
    
    df_TERM2GENE <- as.data.frame(cbind(vec_legend, names(used_target)))
    
    geneList <- sort(c(used_target,simu_genes), decreasing = TRUE)
    
    resy_GSEA <- GSEA(geneList, TERM2GENE = df_TERM2GENE, verbose = TRUE,
                      pvalueCutoff = 1,maxGSSize = 7000, minGSSize = 9,
                      by = 'fgsea', nPermSimple = 100000, eps = 1e-30 )
    
    resy_GSEA@result[["pvalue"]]
    resy_GSEA@result[["enrichmentScore"]]
    
    
    
    
    
    vec_GSEA_p <- c(vec_GSEA_p, resy_GSEA@result[["pvalue"]])
  }
  
  res_GSEA <- cbind(res_GSEA,vec_GSEA_p)
  
}


## use first line as column names
matrix_GSEA <- renamecols (t(res_GSEA))
matrix_GSEA

median_GSEA <- apply(matrix_GSEA,1,median)


## store results
sheet <- 'GSEA_ results'
wb <- loadWorkbook(file_table)
addWorksheet(wb, sheet)
writeData(wb,sheet, matrix_GSEA)
saveWorkbook(wb, file_table, overwrite = TRUE)





###########################################################################
## plot -------------------------------------------------------------------
###########################################################################

## Number of targets used for the test
vec_nb_target <- c(100,250,500,750,1000,1500,2000)


file_stored_table <-  'R.results/P-values_simulation_table.xlsx'
sheet_names <- excel_sheets(path = file_stored_table)

## KS results
table_KS <- 
  as.data.frame(
    read_excel(file_stored_table, 1))

median_KS <- apply(table_KS,2,median)

## GSEA results
table_GSEA <- 
  as.data.frame(
    read_excel(file_stored_table, 2))

median_GSEA <- apply(table_GSEA,2,median)


file_output2 <- 'R.results/Supp_2b_P-values_simulation_plot.pdf'
pdf (file_output2)

## Build plot ----
plot(vec_nb_target,log10(median_GSEA), 
     xaxt='n', cex=1.2, lty = 2,
     pch=15, type="b", cex.lab = 1.2,
     col = 'coral', ylim = c(-10,-0.8), 
     xlab='Number of used targets', 
     ylab = 'log10(p-value)',
     main = paste('log10(p-value) \nby number of used targets'))
points(vec_nb_target,log10(median_KS), 
       pch=19, type='b', col='seagreen', cex=1.2, lty = 1)
axis(1, at=vec_nb_target, labels=vec_nb_target)
legend (x='topright', legend =c('KS test','GSEA'),
        col = c('seagreen','coral'), pch =c(19,15), cex=1.2, 
        lty=c(1,2))

dev.off()