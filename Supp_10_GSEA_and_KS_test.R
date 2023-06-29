###########################################################################
## Script for Matter Arising 
## Comparison of GSEA and Kolmogorov-Smirnov test results
## Parameters: Conserved targets or Both  
## Expression > 4 RPKM and TCS < -0.1 or 0
###########################################################################

## Importation 
source("Functions.R")# functions importation

## Data importation -------------------------------------------------------
## TargetScan v7.1 data importation 
if(!exists('TS'))
  TS <- TargetScan_importation()
TargetScan <- TS[[1]] ; families <- TS[[2]]

## Single-cell data importation
data_imported <- import_SCdata()
data_RNA_19 <- data_imported[[1]] ; data_miRNA_19 <- data_imported[[2]] 

## Sort miRNAs by mean expression 
mean_miRNAs <- apply(data_miRNA_19,1,mean)
list_miRNA <- names(sort(mean_miRNAs[mean_miRNAs > -13],decreasing = TRUE))


###########################################################################
## Compute GSEA and KS test table -----------------------------------------
###########################################################################

###########################################################################
## Build Table for KS vs GSEA 
## conserved, 4, -0.1
## Or both, 4, 0

choice <- 'final'  #  'final' or 'article'

## Parameters definition

if (choice == 'final'){
  conserv <- 'both'
  exp <- 4
  efficacy <- 0
}
if (choice == 'article'){
  conserv <- 'conserved'
  exp <- 4
  efficacy <- -0.1
  
}




###########################################################################
## Compute GSEA table -----------------------------------------------------
## Results preparation 
colnames_mat <- c('miRNA','mean','targets', 'H0', 'pvalue', 'ES')
matrix_result <- t(as.matrix(colnames_mat))

print ("######################")
print (paste ("PARAMETERS : ", conserv, exp, efficacy))


## Beginging of the loop
for (miRNA in list_miRNA){
  
  mean_interest <- mean_miRNAs[which(names(mean_miRNAs) == miRNA)]
  
  res <- apply_GSEA(miRNA, conservation = conserv, thr_exp = exp,
                    selection='TCS', threshold = efficacy, H0 ='selected')
  
  vec_result <- c(miRNA, mean_interest[[1]], res)
  
  matrix_result <- rbind(matrix_result, vec_result)
}

## use first line as column names
results <- renamecols (matrix_result)
results


##########################################################################
## Store results 
file_table <- paste('R.results/GSEA_and_KS_test_',
                    choice,'_table.xlsx', sep ='')
name_wb <- " GSEA results"
sheet <- paste("GSEA", conserv, exp, efficacy)
wb <- createWorkbook(name_wb)
addWorksheet(wb, sheet)
writeData(wb, sheet, results)
saveWorkbook(wb, file_table, overwrite = TRUE)




###########################################################################
## Compute KS test table --------------------------------------------------

## Results preparation 
colnames_mat <- c('miRNA','mean','targets', 'H0', 'pvalue', 'ES')
matrix_result <- t(as.matrix(colnames_mat))

print ("######################")
print (paste ("PARAMETERS : ", conserv, exp, efficacy))


## Beginning of the loop
for (miRNA in list_miRNA){
  
  mean_interest <- mean_miRNAs[which(names(mean_miRNAs) == miRNA)]
  
  res <- KS_test(miRNA, conservation = conserv, thr_exp = exp,
                 selection='TCS', threshold = efficacy)
  
  vec_result <- c(miRNA,mean_interest[[1]],res)
  
  matrix_result <- rbind(matrix_result, vec_result)
}

## use first line as column names
results <- renamecols (matrix_result)
results

##########################################################################
## Store results 
sheet <- paste("KS", conserv, exp, efficacy)
wb <- loadWorkbook(file_table)
addWorksheet(wb, sheet)
writeData(wb, sheet, results)
saveWorkbook(wb, file_table, overwrite = TRUE)



###########################################################################
## Plots ------------------------------------------------------------------
###########################################################################

###########################################################################
## Load GSEA and KS test table
file_stored_table <- paste('R.results/GSEA_and_KS_test_',choice,'_table.xlsx', sep ='')
sheet_names <- excel_sheets(path = file_stored_table)


table <- 
  as.data.frame(
    read_excel(file_stored_table, 1))

miRNA_names <- str_remove(table[,1], 'hsa-')
top10 <- miRNA_names[1:10] 


###########################################################################
## Build boxplot of number of targets for each conditions ----
nb_sheet <- length(sheet_names)

list_target <- list()
vec_mean <- c()

for (x in 1:nb_sheet){
  table <- 
    as.data.frame(
      read_excel(file_stored_table, x))
  
  mean_targets <- mean(as.numeric(table[,3]))
  
  vec_mean <- c(vec_mean, mean_targets)
  list_target <- c(list_target, list(as.numeric(table[,3])) )
  
}

par(mfrow=c(1,1))

boxplot (list_target, xaxt ="n",
         main = 'Boxplot of number of used targets by different conditions')
text(x = 1:nb_sheet,
     y = par("usr")[3] - 0.3,
     labels = sheet_names,
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     cex = 1,
     adj = 0.8)



###########################################################################
## Plots: GSEA and KS test ------------------------------------------------
###########################################################################

BH = TRUE

line_col1 <- rgb(0, 0.6, 0.4,0.3)
line_col2 <- rgb(1, 0.8, 0, 0.5) 

###########################################################################
## Plot GSEA ----
table_GSEA <- 
  as.data.frame(
    read_excel(file_stored_table, 1))

p_GSEA <- as.numeric(table_GSEA[,5])
ES_GSEA <- as.numeric(table_GSEA[,6])

## process p-value and ES
vec_colors_GSEA <- p_value_process(p_GSEA, ES_GSEA, BH = BH)[[2]]
log10_p_GSEA <- p_value_process(p_GSEA, ES_GSEA, BH = BH)[[1]]

## Plot
par(mfrow=c(1,1))
plot (ES_GSEA,log10_p_GSEA, 
      pch = c(rep(19,10), rep(4,length(log10_p_GSEA)-10)), 
      col = vec_colors_GSEA, 
      main = paste('log10(p-value) by ES for', sheet_names[1]))
signif <- which(vec_colors_GSEA != 'grey')
addTextLabels(ES_GSEA[signif],log10_p_GSEA[signif], 
              label = miRNA_names[signif], 
              col.label = ifelse (miRNA_names[signif] %in% top10,
                                  'black', 'grey'))

###########################################################################
## Plot KS ----
table_KS <- 
  as.data.frame(
    read_excel(file_stored_table,2))

p_KS <- as.numeric(table_KS[,5])
D_KS <- as.numeric(table_KS[,6])

## process p-value and ES
vec_colors_KS <- p_value_process(p_KS, D_KS, BH = BH)[[2]]
log10_p_KS <- p_value_process(p_KS, D_KS, BH = BH)[[1]]


## Plot
par(mfrow=c(1,1))
plot (D_KS,log10_p_KS, 
      pch = c(rep(19,10), rep(4,length(log10_p_KS)-10)), 
      col = vec_colors_KS,
      main = paste('log10(p-value) by D for', sheet_names[2]))
signif <- which(vec_colors_KS != 'grey')
addTextLabels(D_KS[signif],log10_p_KS[signif], 
              label = miRNA_names[signif], 
              col.label = ifelse (miRNA_names[signif] %in% top10,
                                  'black', 'grey'))


###########################################################################
## Plot GSEA VS KS test ---------------------------------------------------
###########################################################################
if (choice == 'final'){
  file_output <- paste('R.results/Supp_10_GSEA_vs_KS_',choice,'_plot.pdf', sep ='')
}else {
  file_output <- paste('R.results/GSEA_vs_KS_',choice,'_plot.pdf', sep ='')
}
pdf(file_output)

sign_p_GSEA <- p_value_process(p_GSEA,ES_GSEA, 
                               BH = BH, double_sign = TRUE)[[1]]
sign_p_KS <- -p_value_process(p_KS,D_KS, 
                             BH = BH, double_sign = TRUE)[[1]]


## axis title
if (BH == TRUE){
  ylab.name <- expression (paste ('+/- log10 (',p[BH],') - GSEA')) 
  xlab.name <- expression (paste (' +/-log10 (',p[BH],') - KS'))
} else {
  ylab.name <- paste ('+/-log10 (p-value) - GSEA')
  xlab.name <- paste ('+/-log10 (p-value) - KS')
}

plot (sign_p_KS,sign_p_GSEA, col = vec_colors_GSEA,
      pch = c(rep(19,10), rep(4,length(sign_p_GSEA)-10)),
      xlim = c(-17,13), ylim = c(-17,13),
      lwd = 2, cex = 1.2,
      ylab= ylab.name, xlab= xlab.name,
      main = paste('Results KS vs GSEA\n',
                   conserv,', Expression >',exp,'and TCS <', efficacy))

signif <-  unique(c(which(vec_colors_GSEA != 'grey'), 
                    which(vec_colors_KS != 'grey')))
signiftop10 <- signif[which(signif < 11)]


grid()
thr_p <- log10(0.05)
colors <- c('grey',line_col1,line_col1)
abline(h = c(0, thr_p, -thr_p), col = colors, lwd = 2)
abline(v = c(0, thr_p, -thr_p), col = colors, lwd = 2)
abline(a = 0, b = 1, col = line_col2,lwd = 1)     
addTextLabels(sign_p_KS[signiftop10], sign_p_GSEA[signiftop10], 
              label = miRNA_names[signiftop10], col.label = 'black')
identify(sign_p_KS, sign_p_GSEA, label = miRNA_names )

## miR signif with GSEA
print(length(vec_colors_GSEA[which(vec_colors_GSEA %in% 
                                     c('#33ccff','#000099') )])) # blue
print(length(vec_colors_GSEA[which(vec_colors_GSEA %in% 
                                     c('#DA261F','#69100D') )])) # red


## miR signif with KS
print(length(vec_colors_KS[which(vec_colors_KS %in% 
                                   c('#33ccff','#000099') )])) # blue
print(length(vec_colors_KS[which(vec_colors_KS %in% 
                                   c('#DA261F','#69100D') )])) # red


dev.off()
