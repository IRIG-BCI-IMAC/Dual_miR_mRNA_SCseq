###########################################################################
## Script for Matter Arising
## GSEA analysis to study the threshold to selected genes
## Parameters : Conservation = both
## Expression = -7, -4, 0, 4, 6 or 8  and TCS  <= 0
## H0 = non targets with expression > threshold 
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
## Compute GSEA Expression table ------------------------------------------
###########################################################################

###########################################################################
## BuildTable for GSEA with different expression of threshold 

## Parameters definition
conserv <- 'both' 
param_exp <- c(8,6,4,0,-4,-7)
efficacy <- -0.1

## Create file
file_table <- 'R.results/GSEA_Expresion_table.xlsx'
name_wb <- 'Expression results'
wb <- createWorkbook(name_wb)
saveWorkbook(wb, file_table, overwrite = FALSE)


for (exp in param_exp){
  ## Results preparation 
  colnames_mat <- c('miRNA','mean','targets', 'H0', 'pvalue', 'ES')
  matrix_result <- t(as.matrix(colnames_mat))
  
  print ("######################")
  print (paste ("PARAMETERS : ", conserv, exp, efficacy))
  
  
  ## Beginning of the loop
  for (miRNA in list_miRNA){
    
    print(miRNA)
    mean_interest <- mean_miRNAs[which(names(mean_miRNAs) == miRNA)]
    
    res <- apply_GSEA(miRNA, conservation = conserv, thr_exp = exp,
                      selection='TCS', threshold = efficacy, H0='selected' )
    
    vec_result <- c(miRNA, mean_interest[[1]], res)
    
    matrix_result <- rbind(matrix_result, vec_result)
  }
  
  ## use first line as column names
  results <- renamecols (matrix_result)
  results
  
  
  ## Write in xlsx file
  sheet <- paste("GSEA", conserv ,exp, efficacy)
  wb <- loadWorkbook(file_table)
  addWorksheet(wb, sheet)
  writeData(wb, sheet, results)
  saveWorkbook(wb, file_table, overwrite = TRUE)
  
}

## delete the first sheet
wb = loadWorkbook(file_table)
removeWorksheet(wb, 1)
saveWorkbook(wb, file_table, overwrite = TRUE)




###########################################################################
## Plots ------------------------------------------------------------------
###########################################################################

###########################################################################
## Load GSEA Expression table
file_stored_table <- 'R.results/GSEA_Expresion_table.xlsx'
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

par(mfrow = c(1,1))

boxplot (list_target, xaxt ="n",
         main = 
           'Boxplot of number of used targets by different conditions')
text(x = 1:nb_sheet,
     y = par("usr")[3] - 0.3,
     labels = sheet_names,
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     cex = 1,
     adj = 0.8)



###########################################################################
## Plots: Expression thresholds -------------------------------------------
###########################################################################
file_output <-"R.results/Supp_6_GSEA_Expression_plot.pdf"
pdf(file_output, width = 6.5, height = 5)


BH = TRUE 


## Fit results for TCS threshold == 0
table <- 
  as.data.frame(
    read_excel(file_stored_table, 3))

p <- as.numeric(table[,5])
ES <- as.numeric(table[,6])

## process p-value and ES
vec_colors <- p_value_process(p, ES, BH = BH)[[2]]
log10_p <- p_value_process(p, ES, BH = BH)[[1]]


## model for miR with a positive enrichment
pos_miR <- which(vec_colors  %in% c('#DA261F','#69100D') )
pos_model <- cbind(ES[pos_miR],log10_p[pos_miR])
r_pos <- princomp(pos_model)
b_pos <- r_pos$loadings[2,1] / r_pos$loadings[1,1]
a_pos <- r_pos$center[2] - b_pos * r_pos$center[1]



## model for miR with a negative enrichment
neg_miR <- which(vec_colors  %in% c('#33ccff','#000099') )
neg_model <- cbind(ES[neg_miR],log10_p[neg_miR])
r_neg <- princomp(neg_model)
b_neg <- r_neg$loadings[2,1] / r_neg$loadings[1,1]
a_neg <- r_neg$center[2] - b_neg * r_neg$center[1]





line_col1 <- rgb(0, 0.6, 0.4,0.3)
line_col2 <- rgb(1, 0.8, 0, 0.5) 
ind_cex <- 0.7

#par(mfrow=c(3,2))
for (x in 1:nb_sheet){
  table <- 
    as.data.frame(
      read_excel(file_stored_table, x))
  
  p <- as.numeric(table[,5])
  ES <- as.numeric(table[,6])
  
  ## process p-value and ES
  vec_colors <- p_value_process(p, ES, BH = BH)[[2]]
  log10_p <- p_value_process(p, ES, BH = BH)[[1]]
  
  
  signif <- which(vec_colors != 'grey')
  signiftop10 <- which (which(miRNA_names %in% top10) %in% signif)
  high_signif <- which(vec_colors %in% c('#000099','#69100D','#33ccff') )
  
  to_show <- unique(c(signiftop10))
  ## axis title
  if (BH == TRUE){
    ylab.name <- expression (paste ('log10 (',p[BH],')')) 
  } else {
    ylab.name <- paste ('log10 (p-value)')
  }
  
  mini <- min(log10_p, na.rm =TRUE)
  maxi <- max(log10_p, na.rm =TRUE)
  
  print(c(mini, maxi))
  
  
  plot (ES,log10_p, 
        pch = c(rep(19,10), rep(4,length(ES)-10)),
        col = vec_colors,
        ylim = c(-13,0), xlim = c(-0.6,0.6),
        main = paste(sheet_names[x],
                     "\n mean number of targets =", round(vec_mean[x])) ,
        ylab = ylab.name)
  addTextLabels(ES[to_show],log10_p[to_show], 
                label = miRNA_names[to_show], col.label = "black")
  grid()
  lines (x = c(par('usr')[1],par('usr')[2]), 
         y = c(log10(0.05),log10(0.05)), col = 'lightblue', lty =2)
  lines(x = c(0,0), y = c(par('usr')[3],par('usr')[4]), 
        col ='lightblue', lty =2)
  abline(a_neg,b_neg, col = line_col1, cex = ind_cex)
  abline(a_pos,b_pos, col = line_col1, cex = ind_cex)
  
  
  print(length(vec_colors[which(vec_colors != 'grey')]))
  print(length(vec_colors
               [which(vec_colors %in% c('#33ccff','#000099') )])) # blue
  print(length(vec_colors
               [which(vec_colors %in% c('#DA261F','#69100D') )])) # red
  
}

dev.off()

