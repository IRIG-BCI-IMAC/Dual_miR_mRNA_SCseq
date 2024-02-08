###########################################################################
## Script for Matter Arising 
## GSEA analysis to study parameters of method of correlation 
## Use either Pearson or Spearman correlations 
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


###########################################################################
## Compute GSEA Correlation table ------------------------------------------
###########################################################################

###########################################################################
## Build Table for GSEA analysis using Pearson or Spearman correlations 
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
print (paste ("PARAMETERS : Pearson", conserv, exp, efficacy))


## Beginning of the loop
for (miRNA in list_miRNA){
  print(miRNA)
  mean_interest <- mean_miRNAs[which(names(mean_miRNAs) == miRNA)]
  
  res <- apply_GSEA(miRNA, conservation = conserv, thr_exp = exp,
                    selection='all', threshold = efficacy, H0 ='selected')
  
  vec_result <- c(miRNA, mean_interest[[1]], res)
  
  matrix_result <- rbind(matrix_result, vec_result)
}

## use first line as column names
results <- renamecols (matrix_result)
results

## Store results 
file_table <- 'R.results/GSEA_Correlation_table.xlsx'

name_wb <- "Correlation results"
sheet <- paste("GSEA Pearson",conserv,exp,efficacy)
wb <- createWorkbook(name_wb)
addWorksheet(wb, sheet)
writeData(wb, sheet, results)
saveWorkbook(wb, file_table, overwrite = TRUE)


###########################################################################
###########################################################################
## Apply GSEA with Spearman correlations
## Results preparation ----
colnames_mat <- c('miRNA','mean','targets', 'H0', 'pvalue', 'ES')
matrix_result <- t(as.matrix(colnames_mat))

print ("######################")
print (paste ("PARAMETERS Spearman: ", conserv, exp, efficacy))


## Beginging of the loop
for (miRNA in list_miRNA){
  print(miRNA)
  mean_interest <- mean_miRNAs[which(names(mean_miRNAs) == miRNA)]
  
  res <- apply_GSEA(miRNA, conservation = conserv, thr_exp = exp,
                    selection='all', threshold = efficacy, 
                    spearman = TRUE, H0 ='selected')
  
  vec_result <- c(miRNA, mean_interest[[1]], res)
  
  matrix_result <- rbind(matrix_result, vec_result)
}

## use first line as column names
results <- renamecols (matrix_result)
results

## Write in xlsx file
sheet <- paste("GSEA Spearman",conserv,exp,efficacy)
wb <- loadWorkbook(file_table)
addWorksheet(wb, sheet)
writeData(wb, sheet, results)
saveWorkbook(wb, file_table, overwrite = TRUE)



###########################################################################
## Plots ------------------------------------------------------------------
###########################################################################

###########################################################################
## Load GSEA Correlation table
file_stored_table <- 'R.results/GSEA_Correlation_table.xlsx'
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
## Plots: Pearson and Spearman --------------------------------------------
###########################################################################

BH = TRUE

line_col1 <- rgb(0, 0.6, 0.4,0.3)
line_col2 <- rgb(1, 0.8, 0, 0.5) 

## Plot Pearson ----
table_pearson <- 
  as.data.frame(
    read_excel(file_stored_table, 1))

p_pearson <- as.numeric(table_pearson[,5])
ES_pearson <- as.numeric(table_pearson[,6])

## process p-value and ES
vec_colors_pearson <- p_value_process(p_pearson, ES_pearson, BH = BH)[[2]]
log10_p_pearson <- p_value_process(p_pearson, ES_pearson, BH = BH)[[1]]

## Plot
par(mfrow=c(1,1))
plot (ES_pearson,log10_p_pearson, 
      pch = c(rep(19,10), rep(4,length(log10_p_pearson)-10)), 
      col = vec_colors_pearson, 
      main = paste('log10(p-value) by ES for', sheet_names[1]))
grid ()
abline(h = log10(0.05), col=line_col1)
abline(v = 0, col=line_col1)
signif <- which(vec_colors_pearson != 'grey')
addTextLabels(ES_pearson[signif],log10_p_pearson[signif], 
              label = miRNA_names[signif], 
              col.label = ifelse (miRNA_names[signif] %in% top10, 
                                  'black', 'grey'))


###########################################################################
## Plot Spearman ----
table_spearman <- 
  as.data.frame(
    read_excel(file_stored_table,2))

p_spearman <- as.numeric(table_spearman[,5])
ES_spearman <- as.numeric(table_spearman[,6])


## process p-value and ES
vec_colors_spearman <- p_value_process(p_spearman, ES_spearman, 
                                       BH = BH)[[2]]
log10_p_spearman <- p_value_process(p_spearman, ES_spearman, 
                                    BH = BH)[[1]]



## Plot
par(mfrow=c(1,1))
plot (ES_spearman,log10_p_spearman, 
      pch = c(rep(19,10), rep(4,length(log10_p_spearman)-10)), 
      col = vec_colors_spearman,
      main = paste('log10(p-value) by ES for', sheet_names[2]))
grid ()
abline(h = log10(0.05), col=line_col1)
abline(v = 0, col=line_col1)
signif <- which(vec_colors_spearman != 'grey')
addTextLabels(ES_spearman[signif],log10_p_spearman[signif], 
              label = miRNA_names[signif], 
              col.label = ifelse (miRNA_names[signif] %in% top10, 
                                  'black', 'grey'))



###########################################################################
## Plot Pearson VS Spearman -----------------------------------------------
###########################################################################
file_output <- paste('R.results/Supp_17_GSEA_Correlation_plot.pdf', sep ='')
pdf(file_output)


sign_p_pearson <- p_value_process(p_pearson,ES_pearson, 
                                  BH = BH, double_sign = TRUE)[[1]]
sign_p_spearman <- p_value_process(p_spearman,ES_spearman, 
                                   BH = BH, double_sign = TRUE)[[1]]


## axis title
if (BH == TRUE){
  ylab.name <- expression (paste ('+/- log10 (',p[BH],') - pearson')) 
  xlab.name <- expression (paste (' +/-log10 (',p[BH],') - spearman'))
} else {
  ylab.name <- paste ('+/-log10 (p-value) - pearson')
  xlab.name <- paste ('+/-log10 (p-value) - spearman')
}


mini <- min(sign_p_spearman,sign_p_pearson)
maxi <- max(sign_p_spearman,sign_p_pearson)


plot (sign_p_spearman,sign_p_pearson, col = vec_colors_pearson,
      pch = c(rep(19,10), rep(4,length(sign_p_pearson)-10)),
      xlim = c(-17,20), ylim = c(-17,20),
      lwd = 2, cex = 1.2,
      ylab= ylab.name, xlab= xlab.name,
      main = 
        paste('Results spearman vs pearson\n',
              'Both, Expression > 4 and TCS < 0'))

signif <-  unique(c(which(vec_colors_pearson != 'grey'), 
                    which(vec_colors_spearman != 'grey')))
signiftop10 <- signif[which(signif < 11)]


grid()
thr_p <- log10(0.05)
colors <- c('grey',line_col1,line_col1)
abline(h = c(0, thr_p, -thr_p), col = colors, lwd = 2)
abline(v = c(0, thr_p, -thr_p), col = colors, lwd = 2)
abline(a = 0, b = 1, col = line_col2,lwd = 1)     
addTextLabels(sign_p_spearman[signiftop10], sign_p_pearson[signiftop10], 
              label = miRNA_names[signiftop10], col.label = 'black')

## miR signif with pearson
print(length(vec_colors_pearson[which(vec_colors_pearson %in% 
                                        c('#33ccff','#000099') )])) # blue
print(length(vec_colors_pearson[which(vec_colors_pearson %in% 
                                        c('#DA261F','#69100D') )])) # red


## miR signif with spearman
print(length(vec_colors_spearman[which(vec_colors_spearman %in% 
                                         c('#33ccff','#000099') )])) # blue
print(length(vec_colors_spearman[which(vec_colors_spearman %in% 
                                         c('#DA261F','#69100D') )])) # red


dev.off()
