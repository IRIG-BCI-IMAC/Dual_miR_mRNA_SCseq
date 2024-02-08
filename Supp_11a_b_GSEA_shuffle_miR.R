###########################################################################
## GSEA analysis for negative control with shuffle microRNAs
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
mean_mRNAs <- apply(data_RNA,1,function(x) mean(x, na.rm =TRUE))
hist(mean_mRNAs)

folder_option <- 'shuffle'    ## 'shuffle' or 'shuffle_only_seed'
###########################################################################
## Build  Table for Both or conserved, 
## with Expression > 4 RPKM ans TCS < -0.1

## Parameters definition
conserv <- 'both'
thr_exp <- c(-7,0,4)
efficacy <- 0

## Create file
file_table <- paste0("R.results/GSEA_",folder_option,"_table.xlsx")

name_wb <- 'Conservation results'
wb <- createWorkbook(name_wb)
saveWorkbook(wb, file_table, overwrite = FALSE)

for (exp in thr_exp){
  
  ## Results preparation 
  colnames_mat <- c('miRNA','mean','targets', 'H0', 'pvalue', 'ES')
  matrix_result <- t(as.matrix(colnames_mat))
  
  print ("######################")
  print (paste ("PARAMETERS : ", conserv, exp, efficacy))
  
  
  ## Beginging of the loop
  for (i in 1:100){
    miRNA <- 'hsa-miR-92a-3p'
    miRNA_for_targets <- paste0('shuffle_',i)
    print(miRNA_for_targets)
    mean_interest <- mean_miRNAs[which(names(mean_miRNAs) == miRNA)]
    
    res <- apply_GSEA(miRNA, conservation = conserv, thr_exp = exp,
                      miRNA_for_targets = miRNA_for_targets,
                      selection = 'all', threshold = efficacy,
                      H0 = 'selected', targets_option = 'perl', 
                      folder = folder_option )
    
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
## Load GSEA Shuffle table
file_stored_table <- paste0("R.results/GSEA_",folder_option,"_table.xlsx")
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
## Plots: Shuffle results -------------------------------------------------
###########################################################################
file_output <- paste0("R.results/Supp_11_",folder_option,"_plot.pdf")
pdf(file_output, width = 6.5, height = 5)

BH = FALSE


## Fit results for condition with both conserved and non conserved targets 
table <- 
  as.data.frame(
    read_excel(file_stored_table, 1))


p <- as.numeric(table[,5])
ES <- as.numeric(table[,6])

## process p-value and ES
vec_colors <- p_value_process(p, ES, BH = BH)[[2]]
log10_p <- p_value_process(p, ES, BH = BH)[[1]]



line_col1 <- rgb(0, 0.6, 0.4,0.3)
line_col2 <- rgb(1, 0.8, 0, 0.5) 
ind_cex <- 0.7


x =1
  
  ## compute miR-92a-3p GSEA
  sheet_names[x]
  conserv <- str_split(sheet_names[x],' ')[[1]][2]
  exp <- str_split(sheet_names[x],' ')[[1]][3]

  res_92 <- apply_GSEA('hsa-miR-92a-3p', conservation = conserv, thr_exp = exp,
                       selection = 'all', threshold = efficacy,
                       H0 = 'selected', targets_option = 'online' )
  p92 <- res_92[[3]]
  ES92 <- res_92[[4]]
  

  
  ## load ES and p-value
  table <- 
    as.data.frame(
      read_excel(file_stored_table, x))
  
  p_values <- c(as.numeric(table[,5]), p92)
  ES_res <- c(as.numeric(table[,6]),ES92)
  
  
  ## process p-value and ES
  vec_colors <- p_value_process(p_values, ES_res, BH = BH)[[2]]
  log10_p <- p_value_process(p_values, ES_res, BH = BH)[[1]]
  
  vec_colors[length(vec_colors)] <- '#000099'
  
  signif <- which(vec_colors != 'grey')

  
  
  
  ## axis title
  if (BH == TRUE){
    ylab.name <- expression (paste ('log10 (',p[BH],')')) 
  } else {
    ylab.name <- paste ('log10 (p-value)')
  }
  
  
  plot (ES_res,log10_p, 
        pch = c(rep(4,length(ES_res)-1),5),
        col = vec_colors,
        ylim = c(-22,0), xlim = c(-0.6,0.6),
        main = paste(sheet_names[x],
                     "\n mean number of targets =", round(vec_mean[x])) ,
        ylab = ylab.name, xlab ='ES')
  addTextLabels(ES_res[length(ES_res)],log10_p[length(ES_res)], 
                label = 'miR-92a-3p', col.label = "black")
  grid()
  lines (x = c(par('usr')[1],par('usr')[2]), 
         y = c(log10(0.05),log10(0.05)), col = 'lightblue', lty =2)
  lines(x = c(0,0), 
        y = c(par('usr')[3],par('usr')[4]), col ='lightblue', lty =2)

  
  
  print(length(vec_colors[which(vec_colors != 'grey')]))
  print(length(vec_colors
               [which(vec_colors %in% c('#33ccff','#000099') )])) # blue
  print(length(vec_colors
               [which(vec_colors %in% c('#DA261F','#69100D') )])) # red
  


dev.off()
