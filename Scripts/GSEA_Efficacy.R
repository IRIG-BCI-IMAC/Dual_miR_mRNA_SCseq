###########################################################################
## Script for Matter Arising
## GSEA analysis to study the parameters of efficacy (TCS) 
## Parameters: Conservation = both and expression > 4 RPKM 
## and TCS < -0.5, -0..4, -0.3, -0.2, -0.1, 0
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
## Compute GSEA Efficacy table --------------------------------------------
###########################################################################

###########################################################################
## Build  Table for Efficacy between -0.5 and 0
## Conservation = both
## Expression = 4 RPKM
## Parameters definition

conserv <-'both' 
exp <- 4
param_efficacy <- seq (-0.5, 0, 0.1)

## Create file
name_wb <- 'Efficacy results'
file_table <- "R.results/GSEA_Efficacy_table.xlsx"
wb <- createWorkbook(name_wb)
saveWorkbook(wb, file_table, overwrite = FALSE)


for (efficacy in param_efficacy){
  
  ## Results preparation ----
  colnames_mat <- c('miRNA','mean','targets', 'H0', 'pvalue', 'ES')
  matrix_result <- t(as.matrix(colnames_mat))
  
  print ("######################")
  print (paste ("PARAMETERS : ", conserv, exp, efficacy))
  
  
  ## Beginging of the loop
  for (miRNA in list_miRNA){
    
    print(miRNA)
    
    mean_interest <- mean_miRNAs[which(names(mean_miRNAs) == miRNA)]
    
    res <- apply_GSEA(miRNA, conservation = conserv, thr_exp = exp,
                      selection='TCS', threshold = efficacy)
    
    vec_result <- c(miRNA, mean_interest[[1]], res)
    
    matrix_result <- rbind(matrix_result, vec_result)
  }
  
  ## use first line as column names
  results <- renamecols (matrix_result)
  results
  
  ## Write in xlsx file
  sheet <- paste("GSEA", conserv, exp, efficacy)
  wb <- loadWorkbook(file_table)
  addWorksheet(wb, sheet)
  writeData(wb, sheet, results)
  saveWorkbook(wb, file_table, overwrite = TRUE)
  
}

wb = loadWorkbook(file_table)
removeWorksheet(wb, 1)
saveWorkbook(wb, file_table, overwrite = TRUE)


###########################################################################
## Plots ------------------------------------------------------------------
###########################################################################

###########################################################################
## Load GSEA Efficacy table
file_stored_table <- "R.results/GSEA_Efficacy_table.xlsx"
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
     ## Rotate the labels by 35 degrees.
     srt = 35,
     cex = 1,
     adj = 0.8)



###########################################################################
## Plots: Efficacy thresholds ----------------------------------------------
###########################################################################
file_output <- "R.results/Supp_6_GSEA_Efficacy_plot.pdf"
pdf(file_output, width = 6.5, height = 5)

BH = TRUE 

## Fit results for TCS threshold == 0
table <- 
  as.data.frame(
    read_excel(file_stored_table, 6))


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



## define colors
line_col1 <- rgb(0, 0.6, 0.4,0.3)
line_col2 <- rgb(1, 0.8, 0, 0.5) 
ind_cex <- 0.7


#par(mfrow=c(3,2))
for (x in 1:nb_sheet){
  table <- 
    as.data.frame(
      read_excel(file_stored_table, x))
  
  
  p <- as.numeric(table[,5])
  ES<- as.numeric(table[,6])
 
  
  ## process p-value and ES
  vec_colors <- p_value_process(p, ES, BH = BH)[[2]]
  log10_p <- p_value_process(p, ES, BH = BH)[[1]]
  
  
  signif <- which(vec_colors != 'grey')
  signiftop10 <- which (which(miRNA_names %in% top10) %in% signif)
  
  ## axis title
  if (BH == TRUE){
    ylab.name <- expression (paste ('log10 (',p[BH],')')) 
  } else {
    ylab.name <- paste ('log10 (p-value)')
  }
  
  
  
  plot (ES,log10_p, 
        pch = c(rep(19,10), rep (4,length(ES)-10)), 
        col = vec_colors, cex.axis = ind_cex, cex.main = ind_cex,
        ylim = c(-27,1), xlim = c(-0.7,0.7),
        main = paste(sheet_names[x],
                     "\n mean nb of targets =", round(vec_mean[x])) ,
        ylab = ylab.name)
  addTextLabels(ES[signiftop10 ],log10_p[signiftop10 ], 
                label = miRNA_names[signiftop10 ], 
                col.label = "black", cex.label = ind_cex)
  grid()
  lines (x = c(par('usr')[1],par('usr')[2]), 
         y = c(log10(0.05),log10(0.05)), 
         col = line_col1, lty =2, cex = ind_cex)
  lines(x = c(0,0), y = c(par('usr')[3],par('usr')[4]), lty = 2, 
        col = line_col1, cex = ind_cex)
  abline(a_pos,b_pos ,col = line_col1, cex = ind_cex)
  abline(a_neg,b_neg ,col = line_col1, cex = ind_cex)
  
  
  print(length(vec_colors[which(vec_colors != 'grey')]))
  print(length(vec_colors
               [which(vec_colors %in% c('#33ccff','#000099') )])) # blue
  print(length(vec_colors
               [which(vec_colors %in% c('#DA261F','#69100D') )])) # red
  
}

dev.off()





###########################################################################
###########################################################################
## Visualize using sd by mean graph ---------------------------------------
###########################################################################

## Reduce dataset to interesting miRNAs
data_miRNA <- data_miRNA_19[which(rownames(data_miRNA_19) %in% 
                                    list_miRNA),]
df_miR <- data_miRNA[order(-apply(data_miRNA, 1, mean)),]

## Compute mean and sd 
mean_miR <- apply(df_miR, 1, mean)
sd_miR <- apply(df_miR, 1, sd)


line_col1 <- rgb(0, 0.6, 0.4,0.4)
line_col2 <- rgb(1, 0.8, 0, 0.5)


file_output2 <-"R.results/GSEA_Efficacy_sd_by_mean.pdf"
pdf(file_output2)

BH = TRUE 

for (x in 1:nb_sheet){
  table <- 
    as.data.frame(
      read_excel(file_stored_table, x))
  
  p <- as.numeric(table[,5])
  ES <- as.numeric(table[,6])
  
  ## process p-value and ES
  vec_colors <- p_value_process(p, ES, BH = BH)[[2]]
  log10_p <- p_value_process(p, ES, BH = BH)[[1]]
  
  ## identify significant miRNAs
  signif <- which(vec_colors != 'grey')
  signiftop10 <- which (which(miRNA_names %in% top10) %in% signif)
  
  
  ## axis title
  if (BH == TRUE){
    ylab.name <- expression (paste ('log10 (',p[BH],')')) 
  } else {
    ylab.name <- paste ('log10 (p-value)')
  }
  
  plot (mean_miR, sd_miR, 
        main = paste('Sd by mean 
      colored with results from',sheet_names[x]),
        xlim = c(-13,0), 
        pch = c(rep(19,10), rep(4,length(ES)-10)),
        col = vec_colors,
        #main = 'miRNA expression variability across 19 single cells',
        xlab = 'Mean log2 (miRNA expression)', 
        ylab = 'SD log2 (miRNA expression)')
  panel.smooth(mean_miR, sd_miR, col = NULL, col.smooth = line_col2)
  
  interest <- signif [which (signif < 30)]
  addTextLabels(mean_miR[interest],sd_miR[interest], 
                label = miRNA_names[interest], 
                col.label = vec_colors[interest])
  
  
}

dev.off()


###########################################################################
## Plot figure 1b alone  --------------------------------------------------
file_output3 <- "R.results/Figure_1b_GSEA_Efficacy_sd_by_mean.pdf"
pdf(file_output3)

BH = TRUE 

x <- 6
table <- 
  as.data.frame(
    read_excel(file_stored_table, x))
  
p <- as.numeric(table[,5])
ES <- as.numeric(table[,6])
  
## process p-value and ES
vec_colors <- p_value_process(p, ES, BH = BH)[[2]]
log10_p <- p_value_process(p, ES, BH = BH)[[1]]
  
## identify significant miRNAs
signif <- which(vec_colors != 'grey')
signiftop10 <- which (which(miRNA_names %in% top10) %in% signif)
  
  
## axis title
if (BH == TRUE){
  ylab.name <- expression (paste ('log10 (',p[BH],')')) 
} else {
  ylab.name <- paste ('log10 (p-value)')
}
  
plot (mean_miR, sd_miR, 
      main = paste('Sd by mean 
      colored with results from',sheet_names[x]),
      xlim = c(-13,0), 
      pch = c(rep(19,10), rep(4,length(ES)-10)),
      col = vec_colors,
      #main = 'miRNA expression variability across 19 single cells',
      xlab = 'Mean log2 (miRNA expression)', 
      ylab = 'SD log2 (miRNA expression)')
panel.smooth(mean_miR, sd_miR, col = NULL, col.smooth = line_col2)
  
interest <- signif [which (signif < 30)]
addTextLabels(mean_miR[interest],sd_miR[interest], 
              label = miRNA_names[interest], 
              col.label = vec_colors[interest])

dev.off()


