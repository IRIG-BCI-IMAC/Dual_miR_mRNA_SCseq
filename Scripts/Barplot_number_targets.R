###########################################################################
## Script for Matter Arising 
## Barplot of number of targets predicted by TargetScan
## Both conserved and non conserved targets
## Using all targets or only targets with Expression(log2 scale) > 4 RPKM 
###########################################################################

## Importation
source("Functions.R")# functions importation


## Data importation -------------------------------------------------------
## TargetScan v7.1 data importation 
if (!exists("TS"))
  TS <- TargetScan_importation() 
TargetScan <- TS[[1]] ; families <-TS[[2]]

## Single-cell data importation
data_imported <- import_SCdata()
data_RNA_19 <- data_imported[[1]] ; data_miRNA_19 <- data_imported[[2]] 

genes <- rownames(data_RNA_19)
mean_mRNAs <- apply(data_RNA_19, 1, 
                    function(x) mean(x, na.rm = TRUE)) 
highly_genes <- rownames(data_RNA_19[which(mean_mRNAs > 4),])
length(highly_genes)

## Sort miRNAs by mean expression 
mean_miRNAs <- apply(data_miRNA_19,1,mean)
all_miRNAs <- names(sort(mean_miRNAs[mean_miRNAs > -13],decreasing = TRUE))
list_miRNA <- all_miRNAs [-which(all_miRNAs %in% 
                                   c("hsa-miR-183-5p","hsa-miR-140-3p"))]

###########################################################################
## Compute number of predicted targets ------------------------------------
###########################################################################


## Select only highly expressed genes or all genes 
highly <- FALSE  #'FALSE' for all or 'TRUE' for highly

if (highly == TRUE){
  mean_mRNAs <- apply(data_RNA_19, 1, 
                      function(x) mean(x, na.rm = TRUE)) 
  genes <- rownames(data_RNA_19[which(mean_mRNAs > 4),])
  length(genes)
  
}
  
## Build result matrix
colnames_mat <- c('miRNA', 'nb Conserved targets', 'nb Non conserved targets') 

matrix_result <- t(as.matrix(colnames_mat))

for (miRNA in list_miRNA){
  print(miRNA)
  
  Target_Conserved <- unique(targets_selection 
                             (miRNA = miRNA, 
                              conservation ='conserved',
                              selection = 'all'))
  used_conserved <- Target_Conserved[which(Target_Conserved %in% genes)]
  Target_Non_conserved <- unique(targets_selection 
                                 (miRNA = miRNA, 
                                  conservation ='non conserved',
                                  selection='all'))
  used_non_conserved <- Target_Non_conserved [which(
    Target_Non_conserved %in% genes)]
  
  vec_nb_targets <- c(miRNA, length(used_conserved), 
                      length(used_non_conserved))
  matrix_result <- rbind(matrix_result, vec_nb_targets)
}

df <- as.data.frame(matrix_result[-1,-1])
colnames(df) <- matrix_result[1,-1]
rownames(df) <- str_remove(matrix_result[-1,1],'hsa-')
df[,1] <- as.numeric(df[,1])
df[,2] <- as.numeric(df[,2])

## look at results
len <- length(df[,1])
max_C <- max(df[,1])
mean_C <- round(mean(df[,1]))
median_C <- median(df[,1])

max_NC <- max(df[,2])
mean_NC <- round(mean(df[,2]))
median_NC <- median(df[,2])

## add a column for Both
df[,3] <- df[,2]+df[,1]
colnames(df)[3] <- 'Both'

## Compute percentage of non conserved targets
percentage <- df[2]/df[3] * 100
min (percentage)
max(percentage)
mean(unlist(percentage))
median(unlist(percentage))

## Improve miRNA names
miRNA_name <- str_remove(list_miRNA,'hsa-')
miRNA_name <- str_remove (miRNA_name, "\\.1")


###########################################################################
## Order miRNAs according to results --------------------------------------
###########################################################################

choice <- 'conserved'   # order by 'conserved' or  'both'

  if (choice == 'conserved'){
    ## To order by nb of conserved targets 
    data_res <- t(order_dataset(df))
    
    option <- "number of conserved targets"
  }
  if (choice == 'both'){
    ## To order by nb of total targets 
    data_res <- t(df[order(df[,3],decreasing=FALSE),])
    option <- "number of targets"
  }


###########################################################################
## Plot -------------------------------------------------------------------
###########################################################################

## xaxis labels
highlight <- c(str_remove(all_miRNAs[1:10],"hsa-"))
names_miRNAs <- c()
for (name in colnames(data_res)){
  new_name <- ifelse(name %in% highlight, name, '')
  names_miRNAs <- c(names_miRNAs, new_name)
}
names_miRNAs
names_miRNAs[59] <- names_miRNAs[57]  
names_miRNAs[57] <- "" 

## Save plot
if (highly == FALSE){
  file_output <- 
    paste("R.results/Supp_3a_Barplot_number_targets.pdf", sep ='')
} else{
  file_output <- 
    paste("R.results/Supp_3b_Barplot_number_targets_highly.pdf", sep ='')
}

pdf (file_output, width = 11, height = 8)

## plot ----
par(mar = c(7,5,5,1))

plt <- barplot(data_res[1:2,], 
               col = c("dodgerblue4", "coral2") , 
               border = "white", space = 0.04, 
               ylim = c(0,100 + max(data_res[3,])),
               font.axis=4, xaxt = 'n',
               xlab = "", ylab = "Number of targets",
               main = 
                 paste('Number of targets predicted by TargetScan v7.1\n',
                       'for miRNAs of the dataset'),
               las=2)


text(plt, par("usr")[1], labels = names_miRNAs, srt = 50,
     adj = 1.12, xpd = TRUE, cex = 0.9, col = 'black' ) 

text(x= ifelse(highly == FALSE,4000,800), 
     paste("n =", length(data_res[1,]))) 

legend(x='topleft', fill = c("coral2", "dodgerblue4"),
       legend=c('Non conserved targets', 'Conserved targets'))

dev.off()







