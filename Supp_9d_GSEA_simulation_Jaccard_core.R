###########################################################################
## Script for Matter Arising 
## Simulation for negative control with miR-92a-3p
## Generated sets of 570 targets 
## With different Jaccard index with miR-92a-3p core targets
## Apply GSEA with miR-92a-3p expression for 100 replicates
###########################################################################

set.seed(0)

## Importation
source("Functions.R")# functions importation


## Data importation -------------------------------------------------------
## TargetScan v7.1 data importation 
if(!exists('TS'))
  TS <- TargetScan_importation()
TargetScan <- TS[[1]] ; families <- TS[[2]]


## Single-cell data importation
data_imported <- import_SCdata()
data_RNA_19 <- data_imported[[1]] ; data_miRNA_19 <- data_imported[[2]] 

genes <- rownames(data_RNA_19)


## level of expression of mRNAs
mean_mRNAs <- apply(data_RNA_19, 1, 
                    function(x) mean(x, na.rm = TRUE)) 

## Sort miRNAs by mean expression 
mean_miRNAs <- apply(data_miRNA_19,1,mean)
list_miRNA <- names(sort(mean_miRNAs[mean_miRNAs > -13], decreasing = TRUE))


###########################################################################
## Conditions choice  -----------------------------------------------------
###########################################################################

choice <- 'article'  #  article' or 'article4' or 'final'

## Parameters definition
if (choice == 'article'){
  conservation <- 'conserved'
  thr_exp <- -7
  selection <- 'TCS'
  threshold <- -0.1
  
}
if (choice == 'article4'){
  conservation <- 'conserved'
  thr_exp <- 4
  selection <- 'TCS'
  threshold <- -0.1
  
}
if (choice == 'final'){
  conservation <- 'both'
  thr_exp <- 4
  selection <- 'all'
  threshold <- 0
}

## select genes according to the threshold used on targets expression
select_genes <- names(mean_mRNAs)[which(mean_mRNAs > -7)]

## Select targets of miR-92a-3p ----
targets92TS <- targets_selection('hsa-miR-92a-3p', 
                                 conservation = 'both',
                                 selection = 'all', threshold = -0.1)
targets_92 <- targets92TS [which(targets92TS %in% select_genes)] 
length(targets_92)

## Select core targets of miR-92a-3p
res_GSEA <- apply_GSEA(miRNA ='hsa-miR-92a-3p', conservation = 'both',
                       thr_exp = -7, selection='all', threshold = 0)

core92 <- res_GSEA[[5]][[1]]
length(core92)


###########################################################################
## Compute targets and Jaccard index --------------------------------------
###########################################################################
colnames_mat <- c('miRNA', 'Nb Targets','common' ,'Jaccard index',
                  'common core' ,'Jaccard core')

matrix_res <- t(as.matrix(colnames_mat))
## browse the list of miRNAs
for (miRNA in list_miRNA){
  name_legend <- str_remove (miRNA, 'hsa-')
  
  ## Number of targets
  targetsTS <- targets_selection(miRNA, 
                                 conservation = conservation,
                                 selection =selection, 
                                 threshold = threshold)
  targets <- targetsTS [which(targetsTS %in% select_genes)] 
  
  ## Number of common targets
  common_targets <- targets[which(targets %in% targets_92)]
  print(length(common_targets)/ length(targets) * 100)
  
  ## Compute Jaccard index
  jaccard_index <- jaccard(targets, targets_92)
  print(jaccard_index)
  
  ## Number of common targets with core 
  common_core <- targets[which(targets %in% core92)]
  print(length(common_core)/ length(targets) * 100)
  
  ## Compute Jaccard index
  jaccard_core <- jaccard(targets, core92)
  print(jaccard_core)
  
  ## Store results
  res <- c(name_legend, length(targets), length(common_targets),
           jaccard_index, length(common_core), jaccard_core)  
  matrix_res <- rbind (matrix_res,res)
}

mat <- renamecols(matrix_res)
mat

## Delete miR without targets
indnon0 <- which (mat[,2] != 0)
mat_reduce <- mat[indnon0,]
dim(mat_reduce)[1]





##########################################################################
##########################################################################
## Analysis preparation

print ("######################")

nb_targets <- as.numeric(mat_reduce[,2])
high_nb_targets <- nb_targets[which(nb_targets > 100)]
median_nb_targets <- round(median(high_nb_targets))
median_nb_targets


## Jaccard index to test 
simu_jaccard <- seq(0,0.07,0.01)
simu_jaccard

## Number of core tragets needed in simulations 
number_core_targets <- round(simu_jaccard * 
                               (median_nb_targets + length(core92))) 
number_core_targets

## other parameters
exp <- -7
miRNA  <- 'hsa-miR-92a-3p'
thr_exp <- -7
display <- FALSE
spearman <- FALSE 


## select genes according to their level of expression
selected_genes <- names(mean_mRNAs)[which(mean_mRNAs > thr_exp)]

##  miRNA information ----
ind_miRNA <- which(rownames(data_miRNA_19) == miRNA)[1]
interest_miRNA <- data_miRNA_19[ind_miRNA,]
SD_interet <- sd(interest_miRNA)


vec_ES <- c()
vec_p <- c()
vec_len_H1 <- c()
vec_len_H0 <- c()
vec_common <- c()
vec_core <- c()

## genes which are not miR-92a-3p targets
no_92_tar <- selected_genes [-which(selected_genes %in% targets_92)]

## Results preparation ----
colnames_mat <- c('len_core','Jaccard','targets', 'H0', 'pvalue', 'ES')
matrix_result <- t(as.matrix(colnames_mat))

###########################################################################
## Beginning of the loop
for (len_core in number_core_targets){

    print("#######################################")
    print(paste('Number of core targets :',len_core))
    
    k_ES <- c()
    k_p <- c()
    k_len_H1 <- c()
    k_len_H0 <- c()
    k_common <- c()
    k_core <- c()
    
    for (k in 1:100){
      
      len_no_tar <- median_nb_targets - len_core
      
      used_no_92_tar <- sample(no_92_tar, len_no_tar) 
      used_core92 <- sample(core92, len_core) 
      
      used_targets <- c(used_core92, used_no_92_tar)
      len_H1 <- length(used_targets)
      print(paste('targets =',len_H1))
      len_core <- length(which(used_targets %in% core92))
      print(paste('Common with core miR-92a-3p =', 
                  length(which(used_targets %in% core92))))
      
      ## build the TERM2GENE data frame
      name_legend <- str_replace(miRNA, 'hsa-', '')
      l <- length(used_targets)
      vec_legend <- rep(name_legend, l)
      
      df_TERM2GENE <- as.data.frame(cbind(vec_legend, used_targets))
      
      
      
      ## Delete unselected targets from the dataset -----------------------
      
      unselected_targets <- 
        targets_92[-which(targets_92 %in% used_targets)]
      
      to_be_del <- which(rownames(data_RNA_19) %in% unselected_targets)
      
      if (length(to_be_del) > 0){
        data_RNA_19_reduce <- data_RNA_19[-to_be_del,]
      }
      if (length(to_be_del) == 0 ){
        data_RNA_19_reduce <- data_RNA_19
      }
      
      print(dim(data_RNA_19_reduce))
      
      
      ## Correlation between the miRNA and mRNAs --------------------------
      if(spearman == TRUE){
        correlation <- cor(t(interest_miRNA), t(data_RNA_19_reduce), 
                        method = 'spearman', use = "pairwise.complete.obs")
      } else {
        correlation <- cor(t(interest_miRNA), t(data_RNA_19_reduce), 
                       method = 'pearson', use = "pairwise.complete.obs")
      }
      
      corr_miRNA <- t(correlation)
      len_H0 <- dim(data_RNA_19_reduce)[1]
      ## Build the geneList -----------------------------------------------
      ## Ranked mRNAs due to correlation coefficient 
      mat_corr_miRNA <- as.data.frame(corr_miRNA)
      colnames(mat_corr_miRNA) <-'corr_miRNA'
      gene_df <- mat_corr_miRNA %>% mutate(rank = rank(corr_miRNA,  
                            ties.method = "random"))%>%arrange(desc(rank))
      geneList <- gene_df$corr_miRNA
      names(geneList) <- rownames(gene_df)
      
      print( paste('H0 =', length(geneList)))
      
      if (dim(df_TERM2GENE)[1] > 3){
        
        ## Apply GSEA -----------------------------------------------------
        res_GSEA <- GSEA(geneList, TERM2GENE = df_TERM2GENE, 
                         verbose = TRUE, pvalueCutoff = 1,
                         maxGSSize = 7000, minGSSize = 3,
                         by = 'fgsea', nPermSimple = 100000, eps = 1e-30 )
        
        
        ## Print results
        p <- res_GSEA@result[["pvalue"]]
        ES <- res_GSEA@result[["enrichmentScore"]]
        core <- str_split (res_GSEA@result[['core_enrichment']], '/')
        
        print(p)
        print (ES)
        
      }else {
        p <- "NA"
        ES <- "NA"
        core <- "NA"
      }
      
      
      jaccard <- signif(len_core / (median_nb_targets + length(core92)),1)
      
      vec_result <- c(len_core, jaccard, len_H1, len_H0, p, ES)
      
      matrix_result <- rbind(matrix_result, vec_result) 
      
    }
    
}

df <- renamecols (matrix_result)



##########################################################################
## Store results 
file_table <- paste('R.results/GSEA_simulation_core_100_seed_',
                    choice,'_table.xlsx', sep ='')

name_wb <- "Negative Control results"
sheet <- paste("GSEA simulation")
wb <- createWorkbook(name_wb)
addWorksheet(wb, sheet)
writeData(wb, sheet, df)
saveWorkbook(wb, file_table, overwrite = TRUE)







###########################################################################
## Plots ------------------------------------------------------------------
###########################################################################

###########################################################################
## Load GSEA Negative Control table 
file_stored_table <-  paste('R.results/GSEA_simulation_core_100_seed_',
                            choice,'_table.xlsx', sep ='')
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
## Plot: Boxplot of p-value by jaccard index ------------------------------
###########################################################################

BH = FALSE

## Fit results for condition with miR-92a-3p targets
table <- 
  as.data.frame(
    read_excel(file_stored_table, 1))


p <- as.numeric(table[,5])
ES <- as.numeric(table[,6])

## process p-value and ES
log10_p <- p_value_process(p, ES, BH = BH)[[1]]

jaccard <- as.numeric(table[,2])


data_table <- data.frame(jaccard,log10_p)



## Boxplot
file_output <- 'R.results/Supp_9d_Boxplot_p-value_GSEA_simu.pdf'
pdf(file_output, width = 7.5, height = 6)
boxplot(data_table$log10_p ~ data_table$jaccard ,
        xlab = 'Jaccard index',
        ylab = 'log10 (p-value)', 
        col = 'lightblue')

dev.off()



