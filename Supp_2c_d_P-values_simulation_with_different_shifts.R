###########################################################################
## Script for Matter Arising 
## GSEA analysis and Kolmogorov-Smirnov test
## Apply on simulation for different size of set of targets
## Results are the median on 100 replicates 
## Different shift are used for different groups of targets
###########################################################################

## Importation
source("Functions.R")# functions importation

###########################################################################
## Genes and targets simulation -------------------------------------------
###########################################################################

## Simulation of artificial distribution
nb_genes <- 20000
nb_targets <- 2000

## Define shift
shift1 <- -0.1
shift2 <- -0.06
shift3 <- -0.03
shift4 <- -0.01


## Plot distribution
file_output1 <- 'R.results/Supp_2c_Simulation_distribution_shift.pdf'
pdf(file_output1)
plot(0,type = 'n', lwd = 3, xlim = c(-1,1), ylim = c(0,2),
     xlab = 'Pearson correlation coefficient', ylab = 'Density',
     col = 'grey')

vec_target <- c(100, 400, 500, 1000)
vec_shift <- c(shift1, shift2, shift3, shift4)
vec_colors <- rev(brewer.pal(n = 9, name = 'Paired'))[c(1,3,5,9)]

for (i in 1:length(vec_target)){
  x <- seq(-1, 1, length = vec_target[i] )
  y <- dnorm(x, mean = vec_shift[i], sd = 0.2)
  
  lines(x,y, type = "l", lwd = 4, xlab = "", ylab = "", 
        col = vec_colors[i])
  
}
legend (x = 'topleft', legend  = c(paste(shift1,'(top 100)'),
                                   paste(shift2,'(next 400)'),
                                   paste(shift3,'(next 500)'),
                                   paste(shift4,'(next 1000)'),
                                   'H0'),
        title ='Shift', fill = c(vec_colors, 'grey40'))


x_genes <- seq(-1, 1, length = nb_genes )
y_genes <- dnorm(x_genes, mean = 0, sd = 0.2 )
lines(x_genes,y_genes, type = "l", lwd = 4, 
      xlab = 'Pearson correlation coefficient', ylab = 'Density',
      col = 'grey40')

dev.off()





###########################################################################
## Function to study the sets of simulated targets
compute_ks <- function (list_target){
  
  l <- length(list_target)
  the_mean <- mean(list_target)
  h0 <- rnorm (nb_genes, mean = 0, sd = 0.2)
  res.ks <- ks.test(list_target, h0, eps = 1e-30, 
                    alternative = 'two.sided')
  p <- res.ks$p.value
  if (p == 0){
    p <- 2.2e-16
  }
  
  
  return (c(l, the_mean, p ))
}

###########################################################################
## Function to compute targets
select_random <- function (nb_target){
  
  ## Target definition
  if (nb_target <= 100){
    top100 <- rnorm(100, mean = shift1, sd = 0.2)
    set_target <- sample(top100,nb_target)
  }
  if (nb_target <= 500 && nb_target > 100){
    top100 <- rnorm(100, mean = shift1, sd = 0.2)
    top400 <- rnorm(400, mean = shift2, sd = 0.2)
    set_target <- c(top100, sample(top400,nb_target-100))
  }
  if (nb_target <= 1000 && nb_target > 500){
    top100 <- rnorm(100, mean = shift1, sd = 0.2)
    top400 <- rnorm(400, mean = shift2, sd = 0.2)
    top500 <- rnorm(500, mean = shift3, sd = 0.2)
    set_target <- c(top100, top400, sample(top500,nb_target-500))
  }
  if (nb_target <= 2000 && nb_target > 1000){
    top100 <- rnorm(100, mean = shift1, sd = 0.2)
    top400 <- rnorm(400, mean = shift2, sd = 0.2)
    top500 <- rnorm(500, mean = shift3, sd = 0.2)
    top1000 <- rnorm(1000, mean = shift4, sd = 0.2)
    set_target <- c(top100, top400, top500, sample(top1000,nb_target-1000))
  }
  
  return (set_target)
}

###########################################################################
## KS test study -------------------------------------------------------------
###########################################################################

## Compute p-values
vec_nb_target <- c(100,250,500,750,1000,1250,1500,1750,2000)

mat_mean_KS <- as.matrix(t(vec_nb_target))
mat_pval_KS <- as.matrix(t(vec_nb_target))

## Beginning of the loop for replicate
for (i in 1:100){
  
  vec_length <- c()
  vec_mean <- c()
  vec_pval <- c()
  
  ## Beginning of the loop for different number of targets
  for (nb in vec_nb_target){
    list_target <- select_random(nb) 
    res <- compute_ks (list_target)
    vec_length <- c(vec_length, res[1])
    vec_mean <- c(vec_mean, res[2])
    vec_pval <- c(vec_pval, res[3])
    
  }
  vec_length
  
  mat_mean_KS <- rbind(mat_mean_KS, vec_mean)
  mat_pval_KS <- rbind(mat_pval_KS, vec_pval)
}



mat_mean_KS
mat_pval_KS

## use first line as column names
matrix_pval_KS <- renamecols (mat_pval_KS)


## Compute median
median_pval_KS <- apply (mat_pval_KS[-1,], 2, median)
mean_mean_KS <- signif(apply (mat_mean_KS[-1,], 2, mean),2)
names(mean_mean_KS) <- vec_nb_target
  
## save mean shift

file_table <- 'R.results/P-values_simulation_shift_table.xlsx'
name_wb <- 'Simulation results with shift'
wb <- createWorkbook(name_wb)
sheet <- 'KS mean'
addWorksheet(wb, sheetName = sheet )

writeData(wb, sheet, t(mean_mean_KS))
saveWorkbook(wb, file_table, overwrite = TRUE)



## Plot p-values by number of targets used
vec_colors <- rev(brewer.pal(n = 9, name = 'Paired'))

plot (vec_nb_target, log10(median_pval_KS), type = 'b', pch = 19,
      col = vec_colors, xaxt = 'n',
      ylim = c(-17,0),
      xlab= 'Number of targets used', 
      ylab = 'log10(p-value)', 
      main = "KS test Simulation using different shift\n 100 replicates")
axis(1, at=vec_nb_target, labels=vec_nb_target)
abline (h = log10(0.05), col = 'lightblue')
legend (x = 'bottomleft', legend = mean_mean_KS, title = 'Mean shift',
        fill = vec_colors)
legend (x = 'bottomright', legend  = c(paste(shift1,'(top 100)'),
                                   paste(shift2,'(next 400)'),
                                   paste(shift3,'(next 500)'),
                                   paste(shift4,'(next 1000)')),
        title ='Shift', bty = 'n')


## Store KS results
sheet <- 'KS results'
wb <- loadWorkbook(file_table)
addWorksheet(wb, sheet)

writeData(wb, sheet, matrix_pval_KS)
saveWorkbook(wb, file_table, overwrite = TRUE)





###########################################################################
## GSEA study -------------------------------------------------------------
###########################################################################

## Function to study the sets of simulated targets
compute_GSEA <- function (list_target){
  
  l <- length(list_target)
  print(l)
  the_mean <- mean(list_target)
  h0 <- rnorm (nb_genes, mean = 0, sd = 0.2)
  
  vec <- c()
  for(k in 1:nb_genes){vec <- c(vec, paste('genex',k, sep=''))}
  names(h0) <- vec
  
  
  ##########################################################################
  ## For y 
  vec <- c()
  for(i in 1:length(list_target)){vec <- c(vec, paste('geney',i, sep=''))}
  names(list_target) <- vec
  
  ## Build term2gene
  vec_legend <- rep('test_y', length(list_target))
  
  df_TERM2GENE <- as.data.frame(cbind(vec_legend, names(list_target)))
  
  geneList <- sort(c(list_target,h0), decreasing = TRUE)
  
  resy_GSEA <- GSEA(geneList, TERM2GENE = df_TERM2GENE, verbose = TRUE,
                    pvalueCutoff = 1,maxGSSize = 7000, minGSSize = 9,
                    by = 'fgsea', nPermSimple = 100000, eps = 1e-30 )
  
  resy_GSEA@result[["pvalue"]]
  resy_GSEA@result[["enrichmentScore"]]
  

  p <- resy_GSEA@result[["pvalue"]]
  
  return (c(l, the_mean, p ))
}


## Compute p-value for GSEA
vec_nb_target <- c(100,250,500,750,1000,1250,1500,1750,2000)

mat_mean_GSEA <- as.matrix(t(vec_nb_target))
mat_pval_GSEA <- as.matrix(t(vec_nb_target))

## Beginning of the loop for replicate
for (i in 1:100){
  print(paste('##############', i))
  vec_length <- c()
  vec_mean <- c()
  vec_pval <- c()
  
  ## Beginning of the loop for different nb of targets
  for (nb in vec_nb_target){
    list_target <- select_random(nb) 
    res <- compute_GSEA (list_target)
    vec_length <- c(vec_length, res[1])
    vec_mean <- c(vec_mean, res[2])
    vec_pval <- c(vec_pval, res[3])
    
  }
  vec_length
  
  mat_mean_GSEA <- rbind(mat_mean_GSEA, vec_mean)
  mat_pval_GSEA <- rbind(mat_pval_GSEA, vec_pval)
}



mat_mean_GSEA
mat_pval_GSEA

## use first line as column names
matrix_pval_GSEA <- renamecols (mat_pval_GSEA)

## Compute median
median_pval_GSEA <- apply (mat_pval_GSEA[-1,], 2, median)
mean_mean_GSEA <- signif(apply (mat_mean_GSEA[-1,], 2, mean),2)
names(mean_mean_GSEA) <- vec_nb_target

sheet <- 'GSEA mean'
wb <- loadWorkbook(file_table)
addWorksheet(wb, sheet)

writeData(wb, sheet, t(mean_mean_GSEA))
saveWorkbook(wb, file_table, overwrite = TRUE)



vec_colors <- rev(brewer.pal(n = 9, name = 'Paired'))

## Plot p-values by number of targets used
plot (vec_nb_target, log10(median_pval_GSEA), type = 'b', pch = 19,
      col = vec_colors, xaxt = 'n', lty = 2,
      ylim = c(-17,0),
      xlab= 'Number of targets used', 
      ylab = 'log10(p-value)', 
      main = "GSEA Simulation using different shift\n 100 replicates")
axis(1, at=vec_nb_target, labels=vec_nb_target)
abline (h = log10(0.05), col = 'lightblue')
legend (x = 'bottomleft', legend = mean_mean_GSEA, title = 'Mean shift',
        fill = vec_colors)
legend(x='bottomright', legend  = c(paste('N = 100, shift =',shift1),
                                    paste('N = 400, shift =',shift2),
                                    paste('N = 500, shift =',shift3),
                                    paste('N = 1000, shift =',shift4)))


## Store results
sheet <- 'GSEA results'
wb <- loadWorkbook(file_table)
addWorksheet(wb, sheet)

writeData(wb, sheet, matrix_pval_GSEA)
saveWorkbook(wb, file_table, overwrite = TRUE)




###########################################################################
## Plot results -----------------------------------------------------------
###########################################################################

## Number of targets used for the test
vec_nb_target <- c(100,250,500,750,1000,1250,1500,1750,2000)

## Define shift
shift1 <- -0.1
shift2 <- -0.06
shift3 <- -0.03
shift4 <- -0.01

file_stored_table <-  'R.results/P-values_simulation_shift_table.xlsx'
sheet_names <- excel_sheets(path = file_stored_table)

## KS results
table_KS <- 
  as.data.frame(
    read_excel(file_stored_table, 2))

median_KS <- apply(table_KS,2,median)

mean_mean_KS <-  as.data.frame(read_excel(file_stored_table, 1))

## GSEA results
table_GSEA <- 
  as.data.frame(
    read_excel(file_stored_table, 4))

median_GSEA <- apply(table_GSEA,2,median)

mean_mean_GSEA <-  as.data.frame(read_excel(file_stored_table, 3))


mean_mean <- signif(apply(rbind(mean_mean_GSEA, mean_mean_KS),2,mean),2)

file_output2 <- 'R.results/Supp_2d_P-values_simulation_shift_plot.pdf'
pdf (file_output2)

## Build plot ----
vec_colors <- rev(brewer.pal(n = 9, name = 'Paired'))
plot(vec_nb_target,log10(median_GSEA), 
     xaxt = 'n', cex = 1.2, lty = 2,
     pch = 15, type="b", cex.lab = 1.2,
     col = vec_colors, ylim = c(-16,-0.8), 
     xlab='Number of used targets', 
     ylab = 'log10(p-value)',
     main = paste('log10(p-value) \nby number of used targets'))
points(vec_nb_target,log10(median_KS), 
       pch = 19, type = 'b', col = vec_colors, cex = 1.2)
axis(1, at=vec_nb_target, labels = vec_nb_target)
legend (x ='topright', legend = c( 'KS test', 'GSEA'),
        col = 'grey', pch = c(19,15), cex = 1.2, 
        lty = c(1,2))
legend (x = 'topleft', legend  = c(paste(shift1,'(top 100)'),
                                   paste(shift2,'(next 400)'),
                                   paste(shift3,'(next 500)'),
                                   paste(shift4,'(next 1000)')),
                       title ='Shift', bty = 'n')
legend (x = 'bottomright', legend = mean_mean, title = 'Mean shift',
        fill = vec_colors)


dev.off()
