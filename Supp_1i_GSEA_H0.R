###########################################################################
## Script for Matter Arising 
## GSEA analysis to study H0 influence
## Use either all genes for H0 or genes with expression above the threshold
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
data_imported <- import_SCdata()
data_RNA_19 <- data_imported[[1]] ; data_miRNA_19 <- data_imported[[2]] 

## Sort miRNAs by mean expression 
mean_miRNAs <- apply(data_miRNA_19,1,mean)
all_miRNAs <- names(sort(mean_miRNAs[mean_miRNAs > -13],decreasing = TRUE))
list_miRNA <- all_miRNAs [-which(all_miRNAs %in% 
                                   c("hsa-miR-183-5p","hsa-miR-140-3p"))]


###########################################################################
## Compute GSEA H0 table ------------------------------------------
###########################################################################

###########################################################################
## Build Table for GSEA analysis using all genes for H0
## With parameters to select targets: both, 4, 0

## Parameters definition
conserv <- 'both'
exp <- 4
efficacy <- 0


###########################################################################
## Results preparation ----
colnames_mat <- c('miRNA','mean','targets', 'H0', 'pvalue', 'ES')
matrix_result <- t(as.matrix(colnames_mat))

print ("######################")
print (paste ("PARAMETERS : H0 = all", conserv, exp, efficacy))


## Beginning of the loop
for (miRNA in list_miRNA){
  print(miRNA)
  mean_interest <- mean_miRNAs[which(names(mean_miRNAs) == miRNA)]
  
  res <- apply_GSEA(miRNA, conservation = conserv, thr_exp = exp,
                    selection='all', threshold = efficacy, H0 ='all')
  
  vec_result <- c(miRNA, mean_interest[[1]], res)
  
  matrix_result <- rbind(matrix_result, vec_result)
}

## use first line as column names
results <- renamecols (matrix_result)
results

## Store results 
file_table <- 'R.results/GSEA_H0_table.xlsx'

name_wb <- "H0 results"
sheet <- paste("GSEA H0 = all",conserv,exp,efficacy)
wb <- createWorkbook(name_wb)
addWorksheet(wb, sheet)
writeData(wb, sheet, results)
saveWorkbook(wb, file_table, overwrite = TRUE)


###########################################################################
###########################################################################
## Apply GSEA with genes above threshold on expression for H0
## Results preparation ----
colnames_mat <- c('miRNA','mean','targets', 'H0', 'pvalue', 'ES')
matrix_result <- t(as.matrix(colnames_mat))

print ("######################")
print (paste ("PARAMETERS H0 =  selected : ", conserv, exp, efficacy))


## Beginging of the loop
for (miRNA in list_miRNA){
  print(miRNA)
  mean_interest <- mean_miRNAs[which(names(mean_miRNAs) == miRNA)]
  
  res <- apply_GSEA(miRNA, conservation = conserv, thr_exp = exp,
                    selection='all', threshold = efficacy, 
                    H0 ='selected')
  
  vec_result <- c(miRNA, mean_interest[[1]], res)
  
  matrix_result <- rbind(matrix_result, vec_result)
}

## use first line as column names
results <- renamecols (matrix_result)
results

## Write in xlsx file
sheet <- paste("GSEA H0 = selected ",conserv,exp,efficacy)
wb <- loadWorkbook(file_table)
addWorksheet(wb, sheet)
writeData(wb, sheet, results)
saveWorkbook(wb, file_table, overwrite = TRUE)



###########################################################################
## Plots ------------------------------------------------------------------
###########################################################################

###########################################################################
## Load GSEA H0 table
file_stored_table <- 'R.results/GSEA_H0_table.xlsx'
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
## Plots: H0 = all vs selected --------------------------------------------
###########################################################################

BH = TRUE

line_col1 <- rgb(0, 0.6, 0.4,0.3)
line_col2 <- rgb(1, 0.8, 0, 0.5) 

## Plot H0 = all ----
table_H0all <- 
  as.data.frame(
    read_excel(file_stored_table, 1))

p_H0all <- as.numeric(table_H0all[,5])
ES_H0all <- as.numeric(table_H0all[,6])

## process p-value and ES
vec_colors_H0all <- p_value_process(p_H0all, ES_H0all, BH = BH)[[2]]
log10_p_H0all <- p_value_process(p_H0all, ES_H0all, BH = BH)[[1]]

## Plot
par(mfrow=c(1,1))
plot (ES_H0all,log10_p_H0all, 
      pch = c(rep(19,10), rep(4,length(log10_p_H0all)-10)), 
      col = vec_colors_H0all, 
      main = paste('log10(p-value) by ES for', sheet_names[1]))
grid ()
abline(h = log10(0.05), col=line_col1)
abline(v = 0, col=line_col1)
signif <- which(vec_colors_H0all != 'grey')
addTextLabels(ES_H0all[signif],log10_p_H0all[signif], 
              label = miRNA_names[signif], 
              col.label = ifelse (miRNA_names[signif] %in% top10, 
                                  'black', 'grey'))


###########################################################################
## Plot H0 = selected ----
table_H0selected <- 
  as.data.frame(
    read_excel(file_stored_table,2))

p_H0selected <- as.numeric(table_H0selected[,5])
ES_H0selected <- as.numeric(table_H0selected[,6])


## process p-value and ES
vec_colors_H0selected <- p_value_process(p_H0selected, ES_H0selected, 
                                       BH = BH)[[2]]
log10_p_H0selected <- p_value_process(p_H0selected, ES_H0selected, 
                                    BH = BH)[[1]]



## Plot
par(mfrow=c(1,1))
plot (ES_H0selected,log10_p_H0selected, 
      pch = c(rep(19,10), rep(4,length(log10_p_H0selected)-10)), 
      col = vec_colors_H0selected,
      main = paste('log10(p-value) by ES for', sheet_names[2]))
grid ()
abline(h = log10(0.05), col=line_col1)
abline(v = 0, col=line_col1)
signif <- which(vec_colors_H0selected != 'grey')
addTextLabels(ES_H0selected[signif],log10_p_H0selected[signif], 
              label = miRNA_names[signif], 
              col.label = ifelse (miRNA_names[signif] %in% top10, 
                                  'black', 'grey'))



###########################################################################
## Plot H0 all VS H0selected -----------------------------------------------
###########################################################################
file_output <- paste('R.results/Supp_1i_GSEA_H0.pdf', sep ='')
pdf(file_output)


sign_p_H0all <- p_value_process(p_H0all,ES_H0all, 
                                  BH = BH, double_sign = TRUE)[[1]]
sign_p_H0selected <- p_value_process(p_H0selected,ES_H0selected, 
                                   BH = BH, double_sign = TRUE)[[1]]


## axis title
if (BH == TRUE){
  ylab.name <- expression (paste ('+/- log10 (',p[BH],') - H0 all')) 
  xlab.name <- expression (paste (' +/-log10 (',p[BH],') - H0 selected'))
} else {
  ylab.name <- paste ('+/-log10 (p-value) - H0 all')
  xlab.name <- paste ('+/-log10 (p-value) - H0 selected')
}



plot (0,0, type ='n',
      xlim = c(-18,13), ylim = c(-18,13),
      ylab= ylab.name, xlab= xlab.name,
      main = 
        paste('Results H0 selected vs H0 all\n',
              'Both, Expression > 4 and TCS < 0'))

signif <-  unique(c(which(vec_colors_H0all != 'grey'), 
                    which(vec_colors_H0selected != 'grey')))
signiftop10 <- signif[which(signif < 11)]

grid()
thr_p <- log10(0.05)
colors <- c('grey',line_col1,line_col1)
abline(h = c(0, thr_p, -thr_p), col = colors, lwd = 2)
abline(v = c(0, thr_p, -thr_p), col = colors, lwd = 2)
abline(a = 0, b = 1, col = line_col2,lwd = 1)     

points (sign_p_H0selected,sign_p_H0all, col = vec_colors_H0selected,
        pch = c(rep(19,10), rep(4,length(sign_p_H0all)-10)), lwd = 2, cex = 1.2,)
addTextLabels(sign_p_H0selected[signiftop10], sign_p_H0all[signiftop10], 
              label = miRNA_names[signiftop10], col.label = 'black')



## miR signif with H0all
print(length(vec_colors_H0all[which(vec_colors_H0all %in% 
                                        c('#33ccff','#000099') )])) # blue
print(length(vec_colors_H0all[which(vec_colors_H0all %in% 
                                        c('#DA261F','#69100D') )])) # red


## miR signif with H0selected
print(length(vec_colors_H0selected[which(vec_colors_H0selected %in% 
                                         c('#33ccff','#000099') )])) # blue
print(length(vec_colors_H0selected[which(vec_colors_H0selected %in% 
                                         c('#DA261F','#69100D') )])) # red


dev.off()
