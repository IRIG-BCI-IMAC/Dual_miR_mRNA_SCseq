###########################################################################
## Function file 
## Contain functions called in other scripts
###########################################################################

## Packages installation
if (!'ggplot2' %in% installed.packages()[, "Package"]){
  install.packages("ggplot2", dependencies = TRUE)
}
##
if (!"tidyverse" %in% installed.packages()[, "Package"]){
  install.packages("tidyverse", dependencies = TRUE)
}
##
if (!"dplyr" %in% installed.packages()[, "Package"]){
  install.packages("dplyr", dependencies = TRUE)
}
##
if (!"devtools" %in% installed.packages()[, "Package"]){
  install.packages("devtools", dependencies = TRUE)
}
##
if (!"basicPlotteR" %in% installed.packages()[, "Package"]){
  devtools::install_github("JosephCrispell/basicPlotteR")
}
##
if (!"clusterProfiler" %in% installed.packages()[, "Package"]){
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("clusterProfiler")
}
##
if (!"enrichplot" %in% installed.packages()[, "Package"]){
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("enrichplot")
}
##
if (!"gridExtra" %in% installed.packages()[, "Package"]){
  install.packages("gridExtra")
}
##
if (!"RColorBrewer" %in% installed.packages()[, "Package"]){
  install.packages("RColorBrewer")
}
##
if (!"openxlsx" %in% installed.packages()[, "Package"]){
  install.packages("openxlsx", dependencies = TRUE)
}



## Useful importation
library("stringr") # to manipulate string
library("basicPlotteR") # to add labels without overlapping
library("readxl") # read table file
library("dplyr") # to use the arrange function
library("clusterProfiler")# to apply GSEA
library("enrichplot")# to visualize GSEA plot 
library("ggplot2") # to plot results
library("gridExtra")# to apply grid function on graph
library("openxlsx") # to write on xlsx file
library("RColorBrewer") # to use palette 
source("Function_GSEA_plot.R")# to display GSEA plot with annotation

###########################################################################
###########################################################################
## Function to import and prepare data ------------------------------------
## data from Wang et al., 2019 article
###########################################################################
data2rds <- function (){
  
  ## mRNA data
  data_RNA <- 
    read.table("R.Data/GSE114071_NW_K562_RNA_processed.gct", 
               sep = '\t', skip = 1, header = TRUE)
  saveRDS (data_RNA, "R.Data/data_mRNA_19.rds")
  
  ## miRNA data
  data_miRNA <- 
    read.table("R.Data/GSE114071_NW_K562_miRNA_log2.gct", 
               sep = '\t', header = TRUE )
  saveRDS (data_miRNA, "R.Data/data_miRNA_19.rds")
  
}

###########################################################################
###########################################################################
## Function to import and prepare data ------------------------------------
## from Wang et al., 2019 article
###########################################################################
import_SCdata <- function(){
  
  ## Single-cell data importation 
  data_RNA <- readRDS("R.Data/data_mRNA_19.rds")
  data_miRNA <- readRDS("R.Data/data_miRNA_19.rds")
  
  ## Preparation for miRNAs data 
  library(stringr)
  names_miRNAs_SC <- data_miRNA$Name
  data_miRNA$Name <- str_replace(data_miRNA$Name, 'hsa.', '')
  data_miRNA$Name <- str_replace_all(data_miRNA$Name,
                                     '[:punct:]', '-')
  prettier_names_miRNAs_SC <- data_miRNA$Name
  data_miRNA_19 <- data_miRNA[,c(-1,-2,-22,-23)]
  
  
  ## Preparation for mRNA data
  row_names <- data_RNA$NAME
  row_names
  .rowNamesDF(data_RNA, make.names = TRUE) <- row_names
  data_RNA$NAME <- NULL
  data_RNA$Description <- NULL
  data_RNA_19 <- data_RNA[,c(-6,-21,-22)]
  
  ## Rename columns to be the same in the two data set
  namescol <- colnames(data_miRNA_19)
  colnames(data_RNA_19) <- namescol
  
  
  ## Delete mRNAs with more than 5 NAs 
  count_NA <- rowSums(is.na(data_RNA_19))
  to_be_del <- which(count_NA > 5)
  data_RNA_19_reduce_1 <- data_RNA_19[-to_be_del,]
  
  ## Apply a minimal expression level at 1e-2
  data_RNA_19_reduce_1[data_RNA_19_reduce_1 < 0.01 ] <- 0.01
  
  ## Delete mRNAs with all values == 1e-2
  no_exp <- c()
  for (x in 1:dim(data_RNA_19_reduce_1)[1]){
    line <- data_RNA_19_reduce_1[x,]
    l <- length(which(line == 0.01))
    NAS <- length(which(is.na(line)))
    if (l+NAS==19){
      no_exp <- c(no_exp, x)
    }
  }
  length(no_exp)
  data_RNA_19_reduce <- data_RNA_19_reduce_1[-no_exp,]
  dim(data_RNA_19_reduce)
  
  
  ## keep only one loci for each miRNA
  data_miRNA_19 <- keep_one_loci(data_miRNA = data_miRNA_19, 
                                 names_miRNA = names_miRNAs_SC)
  
  rownames(data_miRNA_19)[which(rownames(data_miRNA_19) 
                                == 'hsa-miR-101-3p')] <- 'hsa-miR-101-3p.1'
  rownames(data_miRNA_19)[which(rownames(data_miRNA_19) 
                                == 'hsa-miR-126-3p')] <- 'hsa-miR-126-3p.1'
  rownames(data_miRNA_19)[which(rownames(data_miRNA_19) 
                                == 'hsa-miR-183-5p')] <- 'hsa-miR-183-5p.1'
  rownames(data_miRNA_19)[which(rownames(data_miRNA_19) 
                                == 'hsa-miR-140-3p')] <- 'hsa-miR-140-3p.2'
  
  data_RNA_19_reduce <- log2(data_RNA_19_reduce)
  
  return (list(data_RNA_19_reduce, data_miRNA_19))
}





###########################################################################
###########################################################################
## Function to select only one loci for each miRNA ------------------------
## The choosen loci is the one with the higher median 
###########################################################################
keep_one_loci <- function (data_miRNA = data_miRNA_19, 
                           names_miRNA = names_miRNAs_SC){
  
  ## build results dataframe
  new_data <- data.frame(matrix(ncol = 19, nrow = 0)) 
  colnames(new_data) <- colnames(data_miRNA)
  
  ## for all miRNAs
  for (miRNA in unique (names_miRNA)){
    ind <- which(names_miRNA == miRNA)
    ## if there are different loci
    if (length(ind) > 1){
      lines <- data_miRNA[ind,]
      mean_miR <- apply (lines, 1, mean)
      indice <- which(mean_miR == max(mean_miR))[1]
      max_line <- data_miRNA[ind[indice],]
    }else ( # if there is only one loci
      max_line <- data_miRNA[ind,]
    )
    new_data <- rbind (new_data, max_line)
    
  }
  rownames(new_data) <- unique (names_miRNA)
  
  return (new_data)
}

###########################################################################
###########################################################################
## Select human information only from TargetScan files --------------------
## the txt files are the one download from TargetScan.org
###########################################################################
targetscanHuman <- function (){
  for (version in c("7.1", "8.0")){
    version <- paste ("_",version, sep = '')
    
    ## import families information
    families <- read.table (paste("R.data/miR_Family_Info",version,".txt", sep =''), 
                            sep = '\t',header = TRUE)
    Human_families <- families [which(families$Species.ID == 9606),] 
    
    
    ## Non conserved predictions
    all <-  read.table (paste("R.data/Summary_Counts.all_predictions",version,".txt", sep =''), 
                        sep = '\t',header = TRUE)
    Human_all <- all [which(all$Species.ID == 9606),] 
    
    
    ## Conserved predictions
    default <-  read.table (paste("R.data/Summary_Counts.default_predictions",version,".txt", sep =''), 
                            sep = '\t',header = TRUE)
    Human_default <- default [which(default$Species.ID == 9606),] 
    
    
    ## merge Conserved and Non conserved prediction
    both <- rbind (Human_default, Human_all)
    
    write.table (Human_families, paste("R.data/Human_miR_Family_Info",version,".txt", sep =''))
    
    write.table (both, paste("R.data/Human_Summary_Counts.both_predictions",version,".txt", sep =''))
    
    
  }
}

###########################################################################
###########################################################################
## Transform TargetScan txt with human information files into rds file ----
###########################################################################
targetscan2rds <- function (){
  for (version in c("7.1", "8.0")){
    version <- paste ("_",version, sep = '')
    

    name_prediction <- paste ("R.data/Human_Summary_Counts.both_predictions",
                              version,".txt", sep ='')
    data_TargetScan <-
      read.table(name_prediction,
                 sep = " ", header = TRUE)
    
    saveRDS(data_TargetScan,
            file = paste("R.data/target_scan",version,".rds", sep =''))
    
  }
}

###########################################################################
###########################################################################
## Function to import TargetScan data -------------------------------------
## The option conservation is used to select conserved targets only 
## or conserved and non conserved together
## conservation = 'all' for both conserved and non conserved targets
## conservation = 'conserved' for conserved targets  only 
###########################################################################

TargetScan_importation <- function( version = "7.1"){
  
  version <- paste ("_",version, sep = '') 
  
  ## TargetScan v8 importation of miRNAs family 
  name_family <- paste ("R.data/Human_miR_family_Info",version,".txt",
                        sep = '')
  families <- read.table(name_family, sep = " ", header = TRUE)
  
  ## prepare families data
  families[,6] <- as.numeric(families[,6])
  
  message('miRNAs families have been imported ....')
  
  ## TargetScan v8 data importation for conserved and non conserved 
  name_prediction <- paste ("R.data/target_scan",
                            version,".rds", sep ='')
  data_TargetScan <- 
    readRDS(name_prediction)
  
  ## prepare TargetScan data
  data_TargetScan[,15] <- as.numeric(data_TargetScan[,15])
  data_TargetScan[,16] <- as.numeric(data_TargetScan[,16])
  
  return (list(data_TargetScan, families))
  
}



###########################################################################
###########################################################################
## Function to select targets ---------------------------------------------
## conservation = 'both' to select conserved and non conserved targets
## conservation = 'conserved" to select conserved targets only 
## selection = 'all' to select all targets
## selection = 'CWCS' to apply a threshold on 
## Cummulative Weighted Contexte ++ Score 
## selection = 'TCS' to apply a threshold on Total Context ++ Score 
## Threshold is a numeric variable to apply if selection = 'CWCS' or 'TCS'
###########################################################################
targets_selection <- function (miRNA='hsa-miR-92a-3p', conservation='both',
                               selection='all', threshold=0){
  # if (miRNA == 'hsa-miR-140-3p'){
  #   miRNA <- 'hsa-miR-140-3p.2'
  # }
  # if (miRNA == 'hsa-miR-183-5p'){
  #   miRNA <- 'hsa-miR-183-5p.1'
  # }
  ## Conservation selection
  if (conservation ==  'both' ){
    data_TargetScan <- TargetScan
  }else if (conservation == 'conserved'){
    data_TargetScan <- 
      TargetScan[which(TargetScan$Total.num.conserved.sites > 0),]
  }else if (conservation == 'non conserved'){
    data_TargetScan <- 
      TargetScan[-which(TargetScan$Total.num.conserved.sites > 0),]
  }
  ## find miRNA family 
  family <- families$Seed.m8[which(families$MiRBase.ID == miRNA)]
  
  ## selection
  if (selection == 'all'){
    row_selection <- 
      data_TargetScan[which(data_TargetScan$miRNA.family == family), ]
  }else if (selection == 'CWCS'){
    
    rows <- data_TargetScan[which(data_TargetScan$miRNA.family == family), ]
    
    indNA <- which(is.na(rows$Cumulative.weighted.context...score))
    rows$Cumulative.weighted.context...score[indNA] <- 0
    
    if (threshold == 0) {
      row_selection <- 
        rows[which(rows$Cumulative.weighted.context...score <= threshold), ]
    }else {
      row_selection <- 
        rows[which(rows$Cumulative.weighted.context...score < threshold), ]
    }
    
    
  }else if (selection == 'TCS'){
    
    rows <- data_TargetScan[which(data_TargetScan$miRNA.family == family), ]
    indNA <- which(is.na(rows$Total.context...score))
    rows$Total.context...score[indNA] <- 0
    if (threshold == 0){
      row_selection <- rows[which(rows$Total.context...score <= threshold), ]
    } else {
      row_selection <- rows[which(rows$Total.context...score < threshold), ]
      
    }
    
  }
  
  return (unique(row_selection$Gene.Symbol))
  
}


###########################################################################
###########################################################################
## Function to apply GSEA on a microRNA -----------------------------------
## parameters to choose targets conditions of selection 
###########################################################################
apply_GSEA <- function (miRNA='hsa-miR-92a-3p', miRNA_for_targets=NULL,
                        conservation='both', thr_exp=4,
                        selection='all', threshold=0, display=FALSE, 
                        spearman=FALSE, H0='selected', scale ='log2'){
  
  
  ## level of expression of mRNAs
  if (scale == 'log2'){
  mean_mRNAs <- apply(data_RNA_19, 1, 
                      function(x) mean(x, na.rm = TRUE)) 
  }
  
  if (scale == 'linear'){
    
    ## Single-cell data importation reload in linear
    data_imported <- import_SCdata()
    data_RNA_19 <- 2^data_imported[[1]] ; data_miRNA_19 <- 2^data_imported[[2]] 
    
    ## level of expression of mRNAs
    mean_mRNAs <- apply(log2(data_RNA_19), 1, 
                        function(x) mean(x, na.rm = TRUE)) 
    
  }
  
  ## selected genes which RPKM > threshold of expression
  selected_genes <- names(mean_mRNAs)[which(mean_mRNAs > thr_exp)]
  
  ##  miRNA information ----
  ind_miRNA <- which(rownames(data_miRNA_19) == miRNA)[1]
  interest_miRNA <- data_miRNA_19[ind_miRNA,]
  SD_interet <- sd(interest_miRNA)
  

  ## Targets Selection ----------------------------------------------------
  if (is.null(miRNA_for_targets)){
    Target_mir <- targets_selection(miRNA, 
                                    conservation = conservation,
                                    selection = selection, 
                                    threshold = threshold) 
  } else {
    Target_mir <- targets_selection(miRNA_for_targets, 
                                    conservation = conservation,
                                    selection = selection, 
                                    threshold = threshold)  
  }
  
  used_targets <- Target_mir[which(Target_mir %in% selected_genes)]
  len_H1 <- length(used_targets)
  print(paste('targets =',len_H1))
  
  ## build the TERM2GENE data frame
  name_legend <- str_replace(miRNA, 'hsa-', '')
  l <- length(used_targets)
  vec_legend <- rep(name_legend, l)
  
  df_TERM2GENE <- as.data.frame(cbind(vec_legend, used_targets))
  
  ## Compute all targets for the miR with no conditions ----
  if (is.null(miRNA_for_targets)){
    all_targets <- targets_selection(miRNA, conservation = 'both',
                                     selection = 'all')
  } else {
    all_targets <- targets_selection(miRNA_for_targets, 
                                     conservation = 'both',
                                     selection = 'all')  
  }
  
  ## Prepare to delete unselected targets ----
  unselected_targets <- all_targets[-which(all_targets %in% used_targets)]
  
  if(H0 == 'all'){
    to_be_del <- which(rownames(data_RNA_19) %in% unselected_targets)
  }
  if (H0 == 'selected'){
    unselected_genes <- names(mean_mRNAs)[-which(names(mean_mRNAs) %in% selected_genes)]
    print('###')
    print(length(selected_genes))
    print(length(unselected_genes))
    to_be_del_gene <- which(rownames(data_RNA_19) %in% unselected_genes)
    to_be_del_targets <- which(rownames(data_RNA_19) %in% unselected_targets)
    to_be_del <- unique(c(to_be_del_gene,to_be_del_targets))
  }
  
  ## Delete unused genes (unslected targets and unselcted genes)
  if (length(to_be_del) > 0){
    data_RNA_19_reduce <- data_RNA_19[-to_be_del,]
  }
  if (length(to_be_del) == 0 ){
    data_RNA_19_reduce <- data_RNA_19
  }
  
  
  print(paste('genes used:',dim(data_RNA_19_reduce)[1]))
  
  
  ## Correlation between the miRNA and mRNAs ------------------------------
  if(spearman == TRUE){
    correlation <- cor(t(interest_miRNA), t(data_RNA_19_reduce), 
                       method = 'spearman', use = "pairwise.complete.obs")
  } else {
    correlation <- cor(t(interest_miRNA), t(data_RNA_19_reduce), 
                       method = 'pearson', use = "pairwise.complete.obs")
  }
  
  corr_miRNA <- t(correlation)
  len_H0 <- dim(data_RNA_19_reduce)[1] - len_H1
  ## Build the geneList -----------------------------------------------------
  ## Ranked mRNAs due to correlation coefficient 
  mat_corr_miRNA <- as.data.frame(corr_miRNA)
  colnames(mat_corr_miRNA) <-'corr_miRNA'
  gene_df <- mat_corr_miRNA %>% mutate(rank = rank(corr_miRNA,  
                                ties.method = "random"))%>%arrange(desc(rank))
  geneList <- gene_df$corr_miRNA
  names(geneList) <- rownames(gene_df)
  
  print( paste('H0 =', len_H0))
  
  if (dim(df_TERM2GENE)[1] > 3){
    
    
    ## Apply GSEA -------------------------------------------------------------
    res_GSEA <- GSEA(geneList, TERM2GENE = df_TERM2GENE, verbose = TRUE,
                     pvalueCutoff = 1,maxGSSize = 7000, minGSSize = 3,
                     by = 'fgsea', nPermSimple = 100000, eps = 1e-30 )
    
    
    ## Print results
    p <- res_GSEA@result[["pvalue"]]
    ES <- res_GSEA@result[["enrichmentScore"]]
    core <- str_split (res_GSEA@result[['core_enrichment']], '/')
    
    ## GSEA plot
    if (length(res_GSEA@result[["enrichmentScore"]]) != 0 &
        display == TRUE){
      
      print('display')
      title <- paste('GSEA plot for ', name_legend,
                     ifelse (is.null(miRNA_for_targets),'', 
                             paste(' with ',miRNA_for_targets, "targets")),
                     conservation, thr_exp, threshold,
                     "\nEnrichment Score = ",
                     signif(res_GSEA@result[["enrichmentScore"]], 
                            digits = 3),
                     "\nNumber of targets = ", len_H1)
      
      # , 
      #       '\n',nb_used_targets, ' targets from',
      #       ifelse(conservation == 'both',
      #        ' conserved and non conserved',' conserved'),
      #       ifelse(selection == 'all',
      #              ' all targets', paste (' and ',selection,' < ',
      #                             threshold)),
      #       sep = '') 
      
      
      plot <- mygseaplot2(res_GSEA,geneSetID = name_legend,
                          title = title
                          , subplots = c(1,2,3))
      # ,pvalue_table = FALSE 
      return (plot)
      
    }
  }else {
    p <- "NA"
    ES <- "NA"
    core <- "NA"
  }
  
  return (list(len_H1, len_H0, p,ES, core))
}

  
  



###########################################################################
## Function to apply KS test on a microRNA --------------------------------
## parameters to choose targets conditions of selection 
###########################################################################
KS_test <- function (miRNA = "hsa-miR-92a-3p", 
                     conservation = 'both',
                     thr_exp = 4,
                     selection='all', threshold=-0){
  
  ## miRNAs informations ----
  print(paste('CDF plot for' ,miRNA))
  name_legend <- str_replace(miRNA,'hsa.','')
  
  ind_miRNA <- which(rownames(data_miRNA_19) == miRNA)
  
  interest_miRNA <- data_miRNA_19[ind_miRNA,]
  
  median_interet <- signif(median(as.numeric(interest_miRNA)),digits = 2)
  median_interet
  SD_interet <- sd(interest_miRNA)
  
  
  ## Select genes which RPKM > threshold of expression
  mean_mRNAs <- apply(data_RNA_19, 1, 
                      function(x) mean(x, na.rm = TRUE)) 
  
  selected_genes <- names(mean_mRNAs)[which(mean_mRNAs > thr_exp)]
  

  highly_genes <- names(mean_mRNAs)[which(mean_mRNAs > thr_exp)]
  
  
  ## Targets Selection ----
  Target_mir <- targets_selection(miRNA, 
                                 conservation = conservation,
                                 selection = selection, 
                                 threshold = threshold)   
  
  
  nb_used_targets <- length(which(Target_mir %in% selected_genes))
  used_targets <- Target_mir[which (Target_mir%in% selected_genes)]
  
  ## Targets selection without conditions ----
  target_mir_all <- targets_selection(miRNA, conservation = 'both',
                                      selection = 'all', threshold = 0)   
  
  non_targets <- selected_genes[-which(selected_genes 
                                              %in% target_mir_all)]
  
  
  ## Correlation between the miRNA and mRNAs ---- 
  
  correlation <- cor(t(interest_miRNA), t(data_RNA_19), 
                     method = 'pearson',use = "pairwise.complete.obs")
  corr_miRNA <- as.vector(correlation)
  names(corr_miRNA) <- colnames(correlation)
  summary(is.na(corr_miRNA)) # Look for NAs 
  
  
  
  ##  H1 Targets mean log2(RPKM) > 4         red
  corr_targets_selected <- sort(corr_miRNA [which(names(corr_miRNA) %in% used_targets)])
  len_H1 <- length(corr_targets_selected)
  
  ## H0 Non targets                         dark blue
  corr_non_targets <- sort(corr_miRNA [which(names(corr_miRNA) 
                                            %in% non_targets)])
  len_H0 <- length(corr_non_targets)
  
  if (len_H1 > 0){
    
    ## Cumulative distribution ----
    
    CDFH1 <- ecdf(corr_targets_selected)
    CDFH0 <- ecdf(corr_non_targets)
    
    
    ## Kolmogorov-Smirnov test ----
    KS1 <- plot_KStest(corr_targets_selected, corr_non_targets)
    p1 <- KS1[1] ; res1 <- KS1[2]
    
    ## the sign of D will be the relative position of vec2 compare to vec1
    KS <- KStest(corr_targets_selected, corr_non_targets, option ="D")
    D <- KS[2]
    p <- KS[1]
    if (p == 0){
      p = 2.2e-16
    }
    
    ## Plot ----
    title <-paste('plot for', name_legend, 
                  'predicted targets\nmedian log2(miRNA) = ',
                  median_interet) 
    
    
    plot(CDFH0, col = 'orange', cex = 0.3,xaxt = 'n', xlim = c(-0.95,0.8),
         main = title,
         ylab ='Cumulative distribution', 
         xlab = paste('Correlation with', name_legend) )
    plot(CDFH1, col = 'red', add = TRUE, cex = 0.3)
    abline(v = 0, col = 'black', lty = 2)
    axis(side = 1, at = seq(-1,1,0.2))
    
    
    
    legend(x = 'topleft', bg = "white", 
           legend = c(paste('Targets highly (n =', len_H1,')'),
                      paste('Non targets (n =', len_H0,')')),
           fill = c('red', 'orange'))
    
    legend(x = 'right',
           title = 'Compare to Non targets',
           legend = c(paste(p1, ' ', res1)),
           fill = c('red'))
    
    return (c(len_H1, len_H0, p, D))
    
  }
  else {return (c(len_H1, len_H0, 'NA','NA'))} 
}

###########################################################################
###########################################################################
## Function for Kolmogorov-Smirnov test to build a plot -------------------
## Function to return an annotation using the p-value 
## to represent the significance 
###########################################################################
plot_KStest<- function (vector1, vector2){
  test <- ks.test(vector1, vector2, exact = TRUE, tol = 1e-30, 
                  simulate.p.value=TRUE, B=3000)
  p <- test[['p.value']]
  print(p)
  if (p < 0.0005){
    res <- '***'
  }else if(p < 0.005){
    res <- '**'
    
  }else if (p < 0.05){
    res <- '*'  
  }else {
    res <- 'ns' # non significant
  }
  p <- signif(test[['p.value']], digits = 2)
  
  return(c(p,res))
  
}

###########################################################################
###########################################################################
## Function for Kolmogorov-Smirnov test -----------------------------------
## function to return the p-value with a sign 
## depending on the side of the enrichment 
## The sign will be the relative side of vector 2 compare to vector 1 
###########################################################################


KStest <- function (vec1, vec2, option ='D'){
  library('dgof')
  KS <- ks.test(vec1, vec2, eps = 1e-30, alternative = 'two.sided')
  D <- KS[["statistic"]][["D"]]
  p <- KS[['p.value']]
  p <- signif(p, digits = 3)
  
  
  ## Difference between Cumulative Distribution Function 
  diff <- curve(decdf(x,vec1,vec2),
                from = min(vec1,vec2), 
                to = max(vec1,vec2))
  
  res <- signif(sum(diff[["y"]]),3)
  
  if (option == 'D'){
    res <- signif(ifelse (max(diff$y) > abs(min(diff$y)), -D ,D),3)
  }
  
  indmax <- diff$x[which(diff$y == max(diff$y))]
  
  
  result <- c(p,res,indmax)
  return(result)
}

###########################################################################
###########################################################################
## Function for difference between CDF ------------------------------------
###########################################################################
decdf <- function(x, vec1, vec2)  ecdf(vec2)(x) - ecdf(vec1)(x)


###########################################################################
###########################################################################
## Function to build ECDF for graph ---------------------------------------
###########################################################################
build_ecdf <- function (vec_correlation){
  
  step <- seq(min(vec_correlation),max(vec_correlation),0.01)
  
  
  empy <- c() ## to store empirical distributuion value
  
  for (value in step){
    emp <- length(which(vec_correlation < value))/length(vec_correlation)
    empy<- c(empy,emp)
  }
  
  step <- c(-1, step, 1)
  empy <- c(0, empy, 1)
  
  return (list(step, empy))
  
} 


###########################################################################
###########################################################################
## Function to process p-values -------------------------------------------
## return color code vector and log10 p-values
## if BH = TRUE, the Bonferoni Hoechberg correction is applied 
## if double_sign = TRUE, the returned p-values will have the sign of the 
## corresponding enrichment score
###########################################################################
p_value_process <- function (vec_pvalue, vec_res, BH = TRUE, 
                             double_sign = FALSE){
  
  
  ## Apply min and get rid of NA
  vec_pvalue[vec_pvalue == '0'] <- '+2.2e-16'
  vec_pvalue[vec_pvalue == ''] <- NA
  
  vec_res[vec_res == ''] <- NA
  
  
  ## get sign  
  vec_res <- as.numeric(vec_res)
  vec_sign <- c()
  for (x in seq_along(vec_res)){
    sign <- ifelse(vec_res[x]>0,'+','-')
    vec_sign <- c(vec_sign,sign)
  }
  
  if (BH == TRUE) {
    ## Apply BH correction 
    vec_abs <- abs(as.numeric(vec_pvalue))
    vec_BH <- p.adjust(vec_abs, method = 'hochberg')
    
    vec_sign_BH <- c()
    for (x in seq_along(vec_BH)){
      p_BH <- paste(vec_sign[x], vec_BH[x], sep ='')
      vec_sign_BH <- c(vec_sign_BH, p_BH)
    }
    
    ## vec colors for plot
    vec_colors <- pvalue2colors (vec_sign_BH)
    
    ## log10 transformation
    vec_log10 <- log10(vec_BH)
    
  } else {
    vec_p <- c()
    for (x in seq_along(vec_pvalue)){
      p <- paste(vec_sign[x], vec_pvalue[x], sep ='')
      vec_p <- c(vec_p, p)
    }
    
    ## vec colors for plot
    vec_colors <- pvalue2colors (vec_p)
    
    ## log10 transformation
    vec_log10 <- log10(as.numeric(vec_pvalue))
    
  }
  
  if (double_sign == TRUE){
    log10_p <- c()
    for (x in seq_along(vec_log10)){
      p <- paste(vec_sign[x], abs(vec_log10[x]), sep ='')
      log10_p <- c(log10_p, p)
    }
    log10_p <- as.numeric(log10_p)
    log10_p[is.na(log10_p)] <- 0
    
  } else{
    log10_p <- vec_log10
  }
  
  return (list(log10_p, vec_colors))
}



###########################################################################
###########################################################################
## Function to get a vector of color using p value vector -----------------
## the sign on the p-value reflect the side of enrichment
###########################################################################
pvalue2colors <- function (vec_pvalue){
  
  vec_colors = c()
  
  for (x in seq_along(vec_pvalue)) {
    
    if (class(vec_pvalue) == "character" ){
      sign <- substr(vec_pvalue[x],1,1)
      if ( sign == '+'){
        palette <- c('#DA261F','#69100D')
      }else {
        palette <- c('#33ccff','#000099')
      }
    }
    
    p <- abs(as.numeric(vec_pvalue[x]))
    
    if (is.na(p) == TRUE){
      color <- NA
    }else if (p < 0.0005){
      color <- palette[2]
    } else if (p < 0.05){
      color <- palette[1]
    }else {
      color <- 'grey'
    }
    
    vec_colors <- c(vec_colors, color)
    
  }
  return (vec_colors)
}



###########################################################################
###########################################################################
## Function to rename columns of a matrix using the first row -------------
###########################################################################
renamecols <- function (matrix_res){
  col.names <- matrix_res[1,]
  colnames(matrix_res) <- col.names
  
  return (matrix_res[-1,]) 
  
}


###########################################################################
###########################################################################
## Function to order data set in term of number of targets ----------------
## use for bar plot of conserved and non-conserved targets
###########################################################################
order_dataset <- function (dataset){
  
  output <- dataset[FALSE,]
  cons <- unique(dataset[,1])
  
  ordered_cons <- cons[order(cons)]
  
  for (value in ordered_cons){
    data_occurence <- dataset[which(dataset[,1]== value),]
    
    data_ordered <- data_occurence [order
                                    (data_occurence$`nb Non conserved targets`,
                                      decreasing = FALSE),]
    
    output <- rbind(output, data_ordered)
    
    
  }
  
  return (output)
  
}

###########################################################################
###########################################################################
## Function to compute Jaccard index --------------------------------------
###########################################################################
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
