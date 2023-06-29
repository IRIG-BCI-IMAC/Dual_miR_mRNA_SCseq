###########################################################################
## Script for Matter Arising
## Negative Control analysis
## Correlations are compute using hsa-miR-92a-3p expression for every cases
## Set of predicted targets are the predicted targets for each miRNAs
## Parameters: conserved, expression > -7 RPKM and TCS < -0.1
## or both, expression > 4 and TCS <= 0
###########################################################################

## Importation 
source("Functions.R")# functions importation

## Data importation -------------------------------------------------------
## TargetScan v7.1 data importation 
if (!exists('TS'))
  TS <- TargetScan_importation()
TargetScan <- TS[[1]] ; families <- TS[[2]]

## Single-cell data importation
data_imported <- import_SCdata()
data_RNA_19 <- data_imported[[1]] ; data_miRNA_19 <- data_imported[[2]] 

## Sort miRNAs by mean expression 
mean_miRNAs <- apply(data_miRNA_19,1,mean)
list_miRNA <- names(sort(mean_miRNAs[mean_miRNAs > -13],decreasing = TRUE))


###########################################################################
## Compute GSEA Negative Control Table ------------------------------------
###########################################################################
###########################################################################
## Build Table Negative Control
## conserved, 4, -0.1
## Or both, 4, 0

choice <- 'article'  # 'final' or 'article'

## Parameters definition

if (choice == 'final'){
  conserv <- 'both'
  exp <- 4
  efficacy <- 0
}
if (choice == 'article'){
  conserv <- 'conserved'
  exp <- -7
  efficacy <- -0.1
  
}


##########################################################################
## Results preparation ----
colnames_mat <- c('miRNA','mean','targets', 'H0', 'pvalue', 'ES')
matrix_result <- t(as.matrix(colnames_mat))

print ("######################")
print (paste ("PARAMETERS : ", conserv, exp, efficacy))


## Beginning of the loop
for (miRNA in list_miRNA){
  print(miRNA)
  mean_interest <- mean_miRNAs[which(names(mean_miRNAs) == miRNA)]
  
  res <- apply_GSEA('hsa-miR-92a-3p',
                           miRNA_for_targets = miRNA, 
                           conservation = conserv, thr_exp = exp,
                           selection='TCS', threshold = efficacy)
  
  vec_result <- c(miRNA, mean_interest[[1]], res)
  
  matrix_result <- rbind(matrix_result, vec_result)
}

## use first line as column names
results <- renamecols (matrix_result)
results


##########################################################################
## Store results 
file_table <- paste('R.results/GSEA_Negativ_Control_',
                    choice,'_table.xlsx', sep ='')

name_wb <- "Negative Control results"
sheet <- paste("GSEA",conserv,exp,efficacy)
wb <- createWorkbook(name_wb)
addWorksheet(wb, sheet)
writeData(wb, sheet, results)
saveWorkbook(wb, file_table, overwrite = TRUE)






###########################################################################
## Plots ------------------------------------------------------------------
###########################################################################

###########################################################################
## Load GSEA Negative Control table 
file_stored_table <-  paste('R.results/GSEA_Negativ_Control_'
               ,choice,'_table.xlsx', sep ='')
sheet_names <- excel_sheets(path = file_stored_table)


table <- 
  as.data.frame(
    read_excel(file_stored_table, 1))

miRNA_names <- str_remove(table[,1], 'hsa-')
top10 <- miRNA_names[1:10] 


###########################################################################
## Build boxplot of number of targets ----
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
## Plot: Negative Control -------------------------------------------------
###########################################################################
if(choice == 'article'){
  file_output <- paste("R.results/Supp_9b_GSEA_Negative_Control_",
                choice,"_plot.pdf", sep = '')
} else (
  file_output <- paste("R.results/GSEA_Negative_Control_",
                choice,"_plot.pdf", sep = '')
)
pdf (file_output, width = 7.5, height = 6)


BH = FALSE
line_col1 <- rgb(0, 0.6, 0.4,0.3)
line_col2 <- rgb(1, 0.8, 0, 0.5) 

table <- 
  as.data.frame(
    read_excel(file_stored_table, 1))
  
p <- as.numeric(table[,5])
ES <- as.numeric(table[,6])

  
## process p-value and ES
vec_colors <- p_value_process(p, ES, BH = BH)[[2]]
log10_p <- p_value_process(p, ES, BH = BH)[[1]]
  
highlight <- which(miRNA_names %in% c('miR-92a-3p', 'miR-25-3p',
                                 'miR-92b-3p', 'miR-125a-5p'))
  
## axis title
if (BH == TRUE){
  ylab.name <- expression (paste ('log10 (',p[BH],')')) 
} else {
  ylab.name <- paste ('log10 (p-value)')
}


## plot
plot (ES,log10_p, 
      pch = 4, cex = 1,
      col = vec_colors,
      ylim = c(-7,0), xlim = c(-0.5,0.5),
      main = paste("Negative Control",sheet_names[1],
                   "\n mean number of targets =", round(vec_mean[x])) ,
      ylab = ylab.name)
addTextLabels(ES[highlight],log10_p[highlight], label = miRNA_names[highlight], 
              col.label = "black", cex.label = 0.8)

grid()
lines (x = c(par('usr')[1],par('usr')[2]), 
       y = c(log10(0.05),log10(0.05)), col = line_col1, lty =2, lwd = 3)
lines(x = c(0,0), 
      y = c(par('usr')[3],par('usr')[4]), col =line_col1, lty = 2, lwd = 3)
length(which(vec_colors != 'grey'))
length(which(vec_colors == 'grey'))
length(which(!is.na(vec_colors)))

dev.off()

