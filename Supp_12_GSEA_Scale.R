###########################################################################
## Script for Matter Arising
## GSEA analysis to study the parameter of scale to copute correlations
## Use either linear scale or log2 scale for miRNA and mRNA expression
## Parameters: Both, Expression > 4 RPKM and TCS < 0.1
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
## Compute GSEA Scale table -----------------------------------------------
###########################################################################

###########################################################################
## Build  Table for Both conserved and non-conserved targets
## with Expression > 4 RPKM ans TCS < 0

## Results preparation for log2 scale -------------------------------------
colnames_mat <- c('miRNA','mean','targets', 'H0', 'pvalue', 'ES')
matrix_result_log2 <- t(as.matrix(colnames_mat))
  
print ("######################")
print (paste ("PARAMETERS : both 4 0"))
  
  
## Beginging of the loop
for (miRNA in list_miRNA){
  print(miRNA)
  mean_interest <- mean_miRNAs[which(names(mean_miRNAs) == miRNA)]
    
  res <- apply_GSEA(miRNA, conservation = 'both', thr_exp = 4,
                      selection='all', threshold = 0, scale ='log2')
    
  vec_result <- c(miRNA, mean_interest[[1]], res)
    
  matrix_result_log2 <- rbind(matrix_result_log2, vec_result)
}
  
## use first line as column names
results_log2 <- renamecols (matrix_result_log2)
results_log2

##########################################################################
## Store results 
file_table <- file_table <- 'R.results/GSEA_Scale_table.xlsx'
name_wb <- "Scale results"
sheet <- paste("GSEA both, 4, 0, log2")
wb <- createWorkbook(name_wb)
addWorksheet(wb, sheet)
writeData(wb, sheet, results_log2)
saveWorkbook(wb, file_table, overwrite = TRUE)




##########################################################################
##########################################################################
## Results preparation for linear scale ----------------------------------
colnames_mat <- c('miRNA','mean','targets', 'H0', 'pvalue', 'ES')
matrix_result_linear <- t(as.matrix(colnames_mat))
  
print ("######################")
print (paste ("PARAMETERS : both 4 0"))
  
  
## Beginging of the loop
for (miRNA in list_miRNA){
  print(miRNA)
  mean_interest <- mean_miRNAs[which(names(mean_miRNAs) == miRNA)]
    
  res <- apply_GSEA(miRNA, conservation = 'both', thr_exp = 4,
                    selection='all', threshold = 0, scale = 'linear')
    
  vec_result <- c(miRNA, mean_interest[[1]], res)
    
  matrix_result_linear <- rbind(matrix_result_linear, vec_result)
}

## use first line as column names
results_linear <- renamecols (matrix_result_linear)
results_linear

## Write in xlsx file
sheet <- paste("GSEA both, 4, 0, linear")
wb <- loadWorkbook(file_table)
addWorksheet(wb, sheet)
writeData(wb, sheet, results_linear)
saveWorkbook(wb, file_table, overwrite = TRUE)



###########################################################################
## Plots ------------------------------------------------------------------
###########################################################################

###########################################################################
## Load GSEA Scale table
file_stored_table <- "R.results/GSEA_Scale_table.xlsx"
sheet_names <- excel_sheets(path = file_stored_table)


table <- 
  as.data.frame(
    read_excel(file_stored_table, 2))

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
         main = 
           'Boxplot of number of used targets by different conditions')
text(x = 1:nb_sheet,
     y = par("usr")[3] - 0.3,
     labels = sheet_names,
     xpd = NA,
     srt = 35,
     cex = 1,
     adj = 0.8)


###########################################################################
## Plots: Log2 and Linear -------------------------------------------------
###########################################################################

BH = TRUE

line_col1 <- rgb(0, 0.6, 0.4,0.3)
line_col2 <- rgb(1, 0.8, 0, 0.5) 

###########################################################################
## Plot log2 ----
table_log2 <- 
  as.data.frame(
    read_excel(file_stored_table, 1))

p_log2 <- as.numeric(table_log2[,5])
ES_log2 <- as.numeric(table_log2[,6])


## process p-value and ES
vec_colors_log2 <- p_value_process(p_log2, ES_log2, BH = BH)[[2]]
log10_p_log2 <- p_value_process(p_log2, ES_log2, BH = BH)[[1]]

## Plot
par(mfrow=c(1,1))
plot (ES_log2,log10_p_log2, 
      pch = c(rep(19,10), rep(4,length(log10_p_log2)-10)), 
      col = vec_colors_log2, 
      main = paste('log10(p-value) by ES for', sheet_names[1]),)
grid ()
abline(h = log10(0.05), col=line_col1)
abline(v = 0, col=line_col1)
signif <- which(vec_colors_log2 != 'grey')
addTextLabels(ES_log2[signif],log10_p_log2[signif], 
              label = miRNA_names[signif], 
              col.label = ifelse (miRNA_names[signif] %in% top10,
                                  'black', 'grey'))

###########################################################################
## Plot linear ----
table_linear <- 
  as.data.frame(
    read_excel(file_stored_table,2))

p_linear <- as.numeric(table_linear[,5])
p_linear[which(p_linear == 0)] <- 2.2e-16
ES_linear <- as.numeric(table_linear[,6])

## process p-value and ES
vec_colors_linear <- p_value_process(p_linear, ES_linear, BH = BH)[[2]]
log10_p_linear <- p_value_process(p_linear, ES_linear, BH = BH)[[1]]



## Plot
par(mfrow=c(1,1))
plot (ES_linear,log10_p_linear, 
      pch = c(rep(19,10), rep(4,length(log10_p_linear)-10)), 
      col = vec_colors_linear,
      main = paste('log10(p-value) by ES for', sheet_names[2]))
grid ()
abline(h = log10(0.05), col=line_col1)
abline(v = 0, col=line_col1)
signif <- which(vec_colors_linear != 'grey')
addTextLabels(ES_linear[signif],log10_p_linear[signif], 
              label = miRNA_names[signif], 
              col.label = ifelse (miRNA_names[signif] %in% top10,
                                  'black', 'grey'))


###########################################################################
## Plot log2 VS linear ----------------------------------------------------
###########################################################################

file_output <- 'R.results/Supp_16_GSEA_Scale_plot.pdf'
pdf(file_output)

sign_p_log2 <- p_value_process(p_log2,ES_log2, 
                               BH = BH, double_sign = TRUE)[[1]]
sign_p_linear <- p_value_process(p_linear,ES_linear, 
                                 BH = BH, double_sign = TRUE)[[1]]


## axis title
if (BH == TRUE){
  ylab.name <- expression (paste ('+/- log10 (',p[BH],') - log2')) 
  xlab.name <- expression (paste (' +/-log10 (',p[BH],') - linear'))
} else {
  ylab.name <- paste ('+/-log10 (p-value) - log2')
  xlab.name <- paste ('+/-log10 (p-value) - linear')
}

plot (sign_p_linear,sign_p_log2, col = vec_colors_log2,
      pch = c(rep(19,10), rep(4,length(sign_p_log2)-10)),
      xlim = c(-17,13), ylim = c(-17,13),
      lwd = 2, cex = 1.2,
      ylab= ylab.name, xlab= xlab.name,
      main = paste('Results linear vs log2
      Both, expression > 4 and TCS <= 0'))

signif <-  unique(c(which(vec_colors_log2 != 'grey'), 
                    which(vec_colors_linear != 'grey')))
signiftop10 <- signif[which(signif < 11)]


grid()
thr_p <- log10(0.05)
colors <- c('grey',line_col1,line_col1)
abline(h = c(0, thr_p, -thr_p), col = colors, lwd = 2)
abline(v = c(0, thr_p, -thr_p), col = colors, lwd = 2)
abline(a = 0, b = 1, col = line_col2,lwd = 1)     
addTextLabels(sign_p_linear[signiftop10], sign_p_log2[signiftop10], 
              label = miRNA_names[signiftop10], col.label = 'black')

## miR signif with log2
print(length(vec_colors_log2[which(vec_colors_log2 %in% 
                                     c('#33ccff','#000099') )])) # blue
print(length(vec_colors_log2[which(vec_colors_log2 %in% 
                                     c('#DA261F','#69100D') )])) # red


## miR signif with linear
print(length(vec_colors_linear[which(vec_colors_linear %in% 
                                       c('#33ccff','#000099') )])) # blue
print(length(vec_colors_linear[which(vec_colors_linear %in% 
                                       c('#DA261F','#69100D') )])) # red
#identify(sign_p_linear,sign_p_log2, 
         #label = miRNA_names)


dev.off()
