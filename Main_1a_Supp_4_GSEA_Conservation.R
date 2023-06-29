###########################################################################
## Script for Matter Arising
## GSEA analysis to study the parameter of target conservation
## Parameters : Conserved or conserved + non-conserved targets
## Expression > 4 RPKM and TCS < 0.1
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
list_miRNA <- names(sort(mean_miRNAs[mean_miRNAs > -13],decreasing = TRUE))

###########################################################################
## Compute GSEA Conservation table ----------------------------------------
###########################################################################

###########################################################################
## Build  Table for Both or conserved, 
## with Expression > 4 RPKM ans TCS < -0.1

## Parameters definition
conservation <- c('both', 'conserved')
exp <- 4
efficacy <- -0.1

## Create file
file_table <- 'R.results/GSEA_Conservation_table.xlsx'
name_wb <- 'Conservation results'
wb <- createWorkbook(name_wb)
saveWorkbook(wb, file_table, overwrite = FALSE)

for (conserv in conservation){
  
  
  ## Results preparation 
  colnames_mat <- c('miRNA','mean','targets', 'H0', 'pvalue', 'ES')
  matrix_result <- t(as.matrix(colnames_mat))
  
  print ("######################")
  print (paste ("PARAMETERS : ", conserv, exp, efficacy))
  
  
  ## Beginging of the loop
  for (miRNA in list_miRNA){
    print(miRNA)
    mean_interest <- mean_miRNAs[which(names(mean_miRNAs) == miRNA)]
    
    res <- apply_GSEA(miRNA, conservation = conserv, thr_exp = exp,
                      selection = 'TCS', threshold = efficacy,
                      H0 = 'selected')
    
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

## remove the empty sheet
wb <- loadWorkbook(file_table)
removeWorksheet(wb, 1)
saveWorkbook(wb, file_table, overwrite = TRUE)


###########################################################################
## Plots ------------------------------------------------------------------
###########################################################################

###########################################################################
## Load GSEA Conservation table
file_stored_table <- "R.results/GSEA_Conservation_table.xlsx"
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
## Plots: Conservation Both or Conserved ----------------------------------
###########################################################################
file_output <- "R.results/Supp_4_GSEA_Conservation_plot.pdf"
pdf(file_output, width = 6.5, height = 5)


BH = TRUE 


## Fit results for condition with both conserved and non conserved targets 
table <- 
  as.data.frame(
    read_excel(file_stored_table, 1))


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


for (x in 1:nb_sheet){
  table <- 
    as.data.frame(
      read_excel(file_stored_table, x))
  
  p_values <- as.numeric(table[,5])
  ES_res <- as.numeric(table[,6])
  
  
  ## process p-value and ES
  vec_colors <- p_value_process(p_values, ES_res, BH = BH)[[2]]
  log10_p <- p_value_process(p_values, ES_res, BH = BH)[[1]]
  
  
  
  signif <- which(vec_colors != 'grey')
  signiftop10 <- which (which(miRNA_names %in% top10) %in% signif)
  
  
  
  ## axis title
  if (BH == TRUE){
    ylab.name <- expression (paste ('log10 (',p[BH],')')) 
  } else {
    ylab.name <- paste ('log10 (p-value)')
  }
  
  
  plot (ES_res,log10_p, 
        pch = c(rep(19,10), rep(4,length(ES)-10)),
        col = vec_colors,
        ylim = c(-17,0), xlim = c(-0.6,0.6),
        main = paste(sheet_names[x],
                     "\n mean number of targets =", round(vec_mean[x])) ,
        ylab = ylab.name, xlab ='ES')
  addTextLabels(ES_res[signiftop10],log10_p[signiftop10], 
                label = miRNA_names[signiftop10], col.label = "black")
  grid()
  lines (x = c(par('usr')[1],par('usr')[2]), 
         y = c(log10(0.05),log10(0.05)), col = 'lightblue', lty =2)
  lines(x = c(0,0), 
        y = c(par('usr')[3],par('usr')[4]), col ='lightblue', lty =2)
  abline(a_neg,b_neg, col = line_col1, cex = ind_cex)
  abline(a_pos,b_pos, col = line_col1, cex = ind_cex)
  
  
  print(length(vec_colors[which(vec_colors != 'grey')]))
  print(length(vec_colors
               [which(vec_colors %in% c('#33ccff','#000099') )])) # blue
  print(length(vec_colors
               [which(vec_colors %in% c('#DA261F','#69100D') )])) # red
  
}

dev.off()



###########################################################################
## Plots: Both and Conserved ----------------------------------------------
###########################################################################

## Choose to apply BH correction 
BH = TRUE 

line_col1 <- rgb(0, 0.6, 0.4,0.4)
line_col2 <- rgb(1, 0.8, 0, 0.5)

###########################################################################
## Plot Both ----
table_both <- 
  as.data.frame(
    read_excel(file_stored_table, 1))

p_both <- as.numeric(table_both[,5])
ES_both <- as.numeric(table_both[,6])


vec_colors_both <- p_value_process(p_both, ES_both, BH = BH)[[2]]
log10_p_both <- p_value_process(p_both, ES_both, BH = BH)[[1]]

## Plot
par(mfrow=c(1,1))
plot (ES_both,log10_p_both, 
      pch = c(rep(19,10), rep(4,length(log10_p_both)-10)), 
      col = vec_colors_both, 
      xlab = 'ES',
      ylab = paste('log10(p-value)',ifelse(BH == TRUE, '- BH','')),
      main = paste('log10(p-value) by ES for', sheet_names[1]))
grid()
thr_p <- log10(0.05)
abline(h = thr_p, col = line_col1, lwd = 2)
abline(v = 0, col = line_col1, lwd = 2)
signif <- which(vec_colors_both != 'grey')
addTextLabels(ES_both[signif],log10_p_both[signif], 
              label = miRNA_names[signif], 
              col.label = ifelse (miRNA_names[signif] %in% top10, 
                                  'black', 'grey'))





###########################################################################
## Plot Conserved ----
table_cons <- 
  as.data.frame(
    read_excel(file_stored_table,2))

p_cons <- as.numeric(table_cons[,5])
ES_cons <- as.numeric(table_cons[,6])

## process p-value and ES
vec_colors_cons <- p_value_process(p_cons, ES_cons, BH = BH)[[2]]
log10_p_cons <- p_value_process(p_cons, ES_cons, BH = BH)[[1]]

## Plot
par(mfrow=c(1,1))
plot (ES_cons,log10_p_cons, 
      pch = c(rep(19,10), rep(4,length(log10_p_cons)-10)), 
      col = vec_colors_cons, cex = 1.8,
      xlab = 'ES',
      ylab = paste('log10(p-value)',ifelse(BH == TRUE, '- BH','')),
      main = paste('log10(p-value) by ES for', sheet_names[2]))
grid()
thr_p <- log10(0.05)
abline(h = thr_p, col = line_col1, lwd = 2)
abline(v = 0, col = line_col1, lwd = 2)
signif <- which(vec_colors_cons != 'grey')
addTextLabels(ES_cons[signif],log10_p_cons[signif], 
              label = miRNA_names[signif], 
              col.label = ifelse (miRNA_names[signif] %in% top10,
                                  'black', 'grey'))




###########################################################################
## Plot Both VS Conserved targets -----------------------------------------
###########################################################################

file_output2 <- 'R.results/Main_1_GSEA_Conserved_vs_Both.pdf'
pdf(file_output2)
## Plot Both vs Conserved ----
sign_p_both <- p_value_process(p_both,ES_both, BH = BH, 
                               double_sign = TRUE)[[1]]
sign_p_cons <- p_value_process(p_cons,ES_cons, BH = BH, 
                               double_sign = TRUE)[[1]]

colors <- c('grey',line_col1,line_col1)
## axis title
if (BH == TRUE){
  ylab.name <- expression (paste ('+/- log10 (',p[BH],
                                  ') - Both conserved and non conserved targets')) 
  xlab.name <- expression (paste ('+/- log10 (',p[BH],
                                  ') - Conserved targets'))
} else {
  ylab.name <- paste ('+/- log10 (p-value)',
                      '- Both conserved and non conserved targets')
  xlab.name <- paste ('+/- log10 (p-value) - Conserved targets')
}


plot(1,1, type ='n',
     xlim = c(-13,13), ylim = c(-13,13),
     ylab = ylab.name, xlab = xlab.name, 
     main = 'Condition both vs conserved\nexpression > 4 and TCS < -0.1')

thr_p <- log10(0.05)
abline(h = c(0, thr_p, -thr_p), col = colors, lwd = 2)
abline(v = c(0, thr_p, -thr_p), col = colors, lwd = 2)
abline(a = 0, b = 1, col = line_col2,lwd = 1)  
grid()
points (sign_p_cons,sign_p_both, col = vec_colors_both,
      pch = c(rep(19,10), rep(4,length(sign_p_both)-10)),
      lwd = 2, cex = 1.2)


signif_cons <-  which(vec_colors_cons != 'grey')
signif_both <- which(vec_colors_both != "grey") 
#signif <- unique(c(signif_cons, signif_both[which(signif_both < 11)]))

miR2highlight <- c('miR-92a-3p','miR-25-3p','miR-30d-5p',
                     'miR-182-5p','miR-27b-3p','miR-125a-5p','miR-26a-5p')
signif <- which(miRNA_names %in% miR2highlight)


   
addTextLabels(sign_p_cons[signif],sign_p_both[signif], 
              label = miRNA_names[signif], 
              col.label ='black')
legend(x='topleft', legend = c('top 10 most expressed microRNAs',
                               'Others microRNAs'),
       col = 'black', pch =c(19,4), cex = 0.8)

dev.off()
