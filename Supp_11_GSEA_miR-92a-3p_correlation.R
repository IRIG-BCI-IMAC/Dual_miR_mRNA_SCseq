###########################################################################
## Script for Matter Arising 
## GSEA analysis to study GSEA results and correlation with miR-92a-3p
## Parameters: Both, Expression > 4 RPKM and TCS <= 0
## Pearson correlation between miR-92a-3p expression and 
## other miR expression
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

genes <- rownames(data_RNA_19)


###########################################################################
## Compute Correlation with miR-92a-3p
exp_92 <- data_miRNA_19[which(rownames(data_miRNA_19) == 'hsa-miR-92a-3p'),]


used_miRNA_exp <- data_miRNA_19[which(rownames(data_miRNA_19) 
                                      %in% list_miRNA),]

dim(used_miRNA_exp)

corr_with_92 <- t(as.matrix(cor(t(exp_92), t(used_miRNA_exp), 
                                method ='pearson')))

mat_cor_92 <- as.matrix(cbind(str_remove(rownames(corr_with_92),'hsa-'), corr_with_92))
colnames(mat_cor_92) <- c('miRNA_names','corr')


###########################################################################
## Load GSEA Expression table
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
## Plots: Collect GSEA results for both 4 0 -------------------------------
###########################################################################
BH = TRUE 

ind_sheet <- which(sheet_names == 'GSEA both 4 0')
## Fit results for TCS threshold == 0
table <- 
  as.data.frame(
    read_excel(file_stored_table, ind_sheet))

p <- as.numeric(table[,5])
ES <- as.numeric(table[,6])


###########################################################################
## Plot GSEA results by correlation with miR-92

file_output <-"R.results/Supp_11_GSEA_Expression_by_cor_92.pdf"
pdf(file_output, width = 7.5, height = 7.5)
BH = TRUE 

line_col1 <- rgb(0, 0.6, 0.4,0.3)
line_col2 <- rgb(1, 0.8, 0, 0.5) 
ind_cex <- 0.7
 
  
## process p-value and ES
vec_colors <- p_value_process(p, ES, BH = BH)[[2]]
mat_colors <- cbind (miRNA_names, vec_colors)
log10_p <- p_value_process(p, ES, BH = BH, double_sign = TRUE)[[1]]
mat_log10_p <- cbind (miRNA_names, log10_p)
  
mat_plot <- merge(mat_log10_p, mat_cor_92, by = "miRNA_names")
mat_plot <- merge(mat_plot, mat_colors, by = "miRNA_names")
mat_plot$corr <- as.numeric(mat_plot$corr)
mat_plot$log10_p <- as.numeric(mat_plot$log10_p)

  
signif <- which(mat_plot$vec_colors != 'grey')
top10_index <- which(mat_plot$miRNA_names %in% top10)
signiftop10 <- signif[which(signif %in% top10_index)]

miR2highlight <- c('miR-92a-3p','miR-92b-3p','miR-25-3p',
                     'miR-30d-5p', 'miR-182-5p',
                     'miR-125a-5p','miR-26a-5p')  

to_show <- which(mat_plot$miRNA_names %in% miR2highlight)

## axis title
if (BH == TRUE){
  ylab.name <- expression (paste ('log10 (',p[BH],') of GSEA results')) 
} else {
  ylab.name <- paste ('log10 (p-value) of GSEA results')
}
  
vec_pch <- rep(4,dim(mat_plot)[1])
vec_pch[top10_index] <- 19
  
plot (mat_plot$corr,mat_plot$log10_p, 
      pch =vec_pch,
      col = mat_plot$vec_colors , cex =1.8,
      ylim = c(-17,13), xlim = c(-1,1),
      main = paste(sheet_names[x],
                   "\n mean number of targets =", round(vec_mean[x])) ,
      ylab = ylab.name, 
      xlab ='Correlation with miR-92a-3p')
  

  
grid()
abline(h = c(log10(0.05),-log10(0.05)), col = 'seagreen', lty =2)
abline(v = 0, col = 'seagreen', lty =2, lwd =3)


## model for results
model_fit <- cbind(mat_plot$corr,mat_plot$log10_p)
r_model <- princomp(model_fit)
b_model <- r_model$loadings[2,1] / r_model$loadings[1,1]
a_model <- r_model$center[2] - b_model * r_model$center[1]
abline(a_model,b_model ,col = 'orange', cex = 1.5)

addTextLabels(as.numeric(mat_plot$corr[to_show]),
              as.numeric(mat_plot$log10_p[to_show]), 
              label = mat_plot$miRNA_names[to_show], col.label = "black")
  
print(length(vec_colors[which(vec_colors != 'grey')]))
print(length(vec_colors
             [which(vec_colors %in% c('#33ccff','#000099') )])) # blue
print(length(vec_colors
             [which(vec_colors %in% c('#DA261F','#69100D') )])) # red
  


dev.off()

