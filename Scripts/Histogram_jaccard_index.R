###########################################################################
## Script for Matter Arising 
## Study Jaccard index between miRNAs and with core targets of miR-92a-3p
## Different conditions
## Histogram of jaccard index
## Real jaccard index are compare with random expected overlap (H0)
###########################################################################

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
all_miRNAs <- names(sort(mean_miRNAs[mean_miRNAs > -13],decreasing = TRUE))
list_miRNA <- all_miRNAs [-which(all_miRNAs %in% 
                                   c("hsa-miR-183-5p","hsa-miR-140-3p"))]


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

## compute core targets of miR-92a-3p
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
                                 selection = selection, 
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

###########################################################################
## Plots histogram of jaccard index ---------------------------------------
###########################################################################
jaccard_choice <- 'core targets'   # or 'targets"

family <- which(mat_reduce[,1] %in% c('miR-92a-3p', 'miR-92b-3p',
                                      'miR-25-3p'))

if (jaccard_choice == 'core targets'){
  vec_jaccard <- as.numeric(mat_reduce[-family,6])
}
if (jaccard_choice == 'targets'){
  vec_jaccard <- as.numeric(mat_reduce[-family,4])
}
  
hist(vec_jaccard, freq = FALSE, breaks = 20,
     col = "lightcoral", 
     xlab= paste('Jaccard index between miR targets and miR-92a-3p',
                 jaccard_choice),
     main = paste('Jaccard index\nfor conditions:',
                  conservation, thr_exp, threshold))
abline(v = quantile(vec_jaccard)[2:4], col=c('cyan3',"blue",'cyan3'),
       lwd=3, lty=2)
legend (x='topright', legend = paste('N =', length(vec_jaccard)), 
        bty = 'n')



## histogram with density
hist(vec_jaccard, freq = FALSE, breaks = 20, 
     col = "lightcoral", 
     xlab= paste('Jaccard index between miR targets and miR-92a-3p'
                 ,jaccard_choice),
     main = paste('Jaccard index \nfor conditions:',
                  conservation, thr_exp, threshold))
legend (x='topright', legend = paste('N =', length(vec_jaccard)), 
        bty = 'n')
lines(density (vec_jaccard), col = 'orange', lwd =3)


###########################################################################
## Compute Null Hypothesis ------------------------------------------------
###########################################################################

size_ensemble <- length(select_genes)
nb_miRNA <- as.numeric(mat_reduce[-family,2])

if (jaccard_choice ==  'core targets'){
  
  nb_ref <- length(core92)
  
}
if (jaccard_choice ==  'targets'){
  
  nb_ref <- length(targets_92)
  
}




## compute random histogram of percentage of common targets

vec_random_jaccard <- c()

for (nb in nb_miRNA){
  
  vec_miRx <- c()
  for (j in 1:100){
    targets_miRx <- sample (select_genes, nb )
    targets_miR92 <- sample (select_genes, nb_ref)
    common <- length(which (targets_miRx %in% targets_miR92))/ nb *100
    
    jaccard_i <- jaccard (targets_miRx, targets_miR92)
    vec_miRx <- c(vec_miRx, jaccard_i)
  }
  vec_random_jaccard <- c(vec_random_jaccard, vec_miRx)
  
}

vec_random_jaccard
hist(vec_random_jaccard)

hist(nb_miRNA, breaks = 15)



###########################################################################
## Plot ------------------------------------------------------------------
###########################################################################
if (choice == 'article'){
  file_output <- paste('R.results/Supp_8a_histogram_jaccard_',
                       jaccard_choice,'_',
                       choice,'.pdf', sep ="")
}else {
  file_output <- paste('R.results/histogram_jaccard_',jaccard_choice,'_',
                       choice,'.pdf', sep ="")
}

pdf(file_output, width = 9, height = 8)

d <- density(vec_random_jaccard)

## plot Both histogram in the same graph
hist(vec_jaccard, freq = FALSE, breaks = 15, 
     xlim = c(0,max(vec_jaccard)+0.05),
     ylim = c(0,max(d$y)+5), density = 30, angle = 20,
     col = "lightcoral", 
     xlab= paste('Jaccard index between miR targets and miR-92a-3p',
                 jaccard_choice),
     main = paste('Jaccard index\nfor conditions:',
                  conservation, thr_exp, threshold))
hist(vec_random_jaccard, add = TRUE, freq = FALSE, col = 'grey', 
     density =30, angle = 170,
     lwd = 3, breaks = ifelse(jaccard_choice == 'core',15,10))
lines(density(vec_random_jaccard), col = 'grey52', lwd = 3)
lines(density (vec_jaccard), col = '#a35437', lwd =3)

res_ks <- ks.test(vec_random_jaccard, vec_jaccard)
legend (x='topright',  bty = 'n', cex = 1,
        legend = paste('N =', length(vec_jaccard),
              '\np-value KS test', ifelse(res_ks$p.value == 0,
                                          '< 2.2e-16', 
                          paste('=', signif(res_ks$p.value,2)))))
legend (x= 'topleft',
        legend = c('Jaccard index','Null hypothesis (H0)'),
        fill = c('lightcoral','grey'), bty = 'n', cex = 1)



dev.off()


max (vec_random_jaccard)

high_jaccard <- 
   mat_reduce[which(as.numeric(mat_reduce[,6]) > max (vec_random_jaccard)),1]



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

## define parameter for the graph
BH = FALSE
line_col1 <- rgb(0, 0.6, 0.4,0.3)
line_col2 <- rgb(1, 0.8, 0, 0.5) 


## get p-value and ES score
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


ind_high_jaccard <- which(miRNA_names %in% high_jaccard)

## plot
plot (ES,log10_p, 
      pch = 4, cex = 1,
      col = vec_colors,
      ylim = c(-7,0), xlim = c(-0.5,0.5),
      main = paste("Negative Control",sheet_names[1],
                   "\n mean number of targets =", round(vec_mean[x])) ,
      ylab = ylab.name)
addTextLabels(ES[highlight],log10_p[highlight], 
              label = miRNA_names[highlight], 
              col.label = "black", cex.label = 0.8)

addTextLabels(ES[ind_high_jaccard],log10_p[ind_high_jaccard], 
              label = miRNA_names[ind_high_jaccard], 
              col.label = "red", cex.label = 0.8)

grid()
lines (x = c(par('usr')[1],par('usr')[2]), 
       y = c(log10(0.05),log10(0.05)), col = line_col1, lty =2, lwd = 3)
lines(x = c(0,0), 
      y = c(par('usr')[3],par('usr')[4]), col =line_col1, lty = 2, lwd = 3)
length(which(vec_colors != 'grey'))
length(which(vec_colors == 'grey'))
length(which(!is.na(vec_colors)))




data_p <- as.data.frame(log10_p)
data_p[ind_high_jaccard,2] <-"high"
data_p[-ind_high_jaccard,2] <-"low"

boxplot(data_p$log10_p ~ data_p$V2,
     col =c('lightblue', 'grey'), xlab ="jaccard index core targets",
     ylab ='log10(p-value) Negative control')
legend ( x = 'bottom', fil = c('lightblue', 'grey'), 
         legend = c('high jaccard index n= 27',
                    'jaccard index in H0 n =30'))

###########################################################################
## Plot log10(p-value) by Jaccard index
## compute correlation

vec_jaccard


## delete results from the miR-92a-3p family 
## and from miR without targets at all
ind_family92 <- which (miRNA_names %in% c('miR-92a-3p',
                                          'miR-92b-3p','miR-25-3p'))
ind_no_targets <- which(list_target[[1]] == 0)
ind_to_del <- c(ind_family92, ind_no_targets )


log10_p_red <- log10_p[-ind_to_del]
log10_p_red[is.na(log10_p_red)] <- 0
vec_colors_red <- vec_colors[-ind_to_del]


###########################################################################
### Define vectors for plot
vec_names <- mat_reduce[-family,1]
log10_p_red
vec_nb_targets <- as.numeric(mat_reduce[-family,2])
vec_percent_common <- 
  as.numeric(mat_reduce[-family,3])/ as.numeric(mat_reduce[-family,2])*100
vec_jaccard <-  as.numeric(mat_reduce[-family,4])
vec_percent_core <- 
  as.numeric(mat_reduce[-family,5])/ as.numeric(mat_reduce[-family,2])*100
vec_jaccard_core <- as.numeric(mat_reduce[-family,6])
vec_colors_red[is.na(vec_colors_red)] <- 'grey'


###
file_output2 <- paste('R.results/ Plot_hypothesis_NC_',choice,'.pdf')
pdf (file_output2, width = 9, height = 8)


list_vec <- list(vec_nb_targets, vec_percent_common, vec_jaccard,
                 vec_percent_core, vec_jaccard_core, log10_p_red)

name_lab <- c('Number of targets',
              'Common targets %','Jaccard index common targets',
              'Core targets % ','Jaccard index core targets', 
              'log10(p-value) Negative Control')

for (i in 1:6){
  for (j in i:6){
    
    if (i != j){
    print(paste(i,'-',j))
    
      ## plot
      plot (list_vec[[i]], list_vec[[j]], pch = 19,
            col = vec_colors_red,
            xlab = name_lab[i],
            ylab = name_lab[j])  
      lines (x = c(par('usr')[1],par('usr')[2]), 
             y = c(log10(0.05),log10(0.05)), col = line_col1, 
             lty =2, lwd = 3)
      
      x_text <- quantile(list_vec[[i]])[2]
      y_text <- quantile(list_vec[[j]])[2]
      
      correlation <- round(cor (list_vec[[i]], list_vec[[j]],
                                    method = 'pearson'),3)
      text (x= x_text, y =y_text, label=paste('R =',correlation))
      
      Reg <- lm( list_vec[[j]] ~ list_vec[[i]])
      abline(Reg, col='cyan4')
      
    }
  }
}


dev.off()




## plot correlation between log10(p-value) and jaccard index core targets

file_output3 <- 'R.results/Supp_8c_Correlation_p-value_jaccard_core.pdf'
pdf(file_output3, width = 7, height = 7)

plot (vec_jaccard_core, log10_p_red, pch = 4,
      col = vec_colors_red,
      xlab = name_lab[5],
      ylab = name_lab[6])  
lines (x = c(par('usr')[1],par('usr')[2]), 
       y = c(log10(0.05),log10(0.05)), col = line_col1, lty =2, lwd = 3)



correlation <- round(cor (vec_jaccard_core, log10_p_red,
                          method = 'pearson'),3)
legend (x = 'bottomleft', legend = paste('R =',correlation), bty ='n')

Reg <- lm(log10_p_red ~  vec_jaccard_core)
abline(Reg, col='cyan4')

ind_name <- which(vec_names == 'miR-125a-5p')
addTextLabels(vec_jaccard_core[ind_name], log10_p_red[ind_name], 
              label = vec_names[ind_name], col.label = 'black')

dev.off()


