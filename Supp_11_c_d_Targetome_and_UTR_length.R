###########################################################################
## Script for Matter Arising
## Compute 3'UTR length for all mRNAs and miRNA targetome
## 3'UTR and mRNA expression in cells grouped by miR-92a-3p expression
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

###########################################################################
## Compute Targetome
l <- dim(data_RNA)[1]
targetome <- matrix(rep(0,l ), nrow =l , ncol =1)
rownames(targetome) <- rownames(data_RNA)
dim(targetome)

## for each miR of the list compute their predicted targets
for (miR in list_miRNA){
  print(miR)
  targets <- targets_selection(miRNA = miR, conservation = 'both',
                               selection ='all')
  for (gene in targets){
    targetome[which(rownames(targetome) == gene),1] <- targetome[which(rownames(targetome) == gene),1] +1
    
  }
}


length(targetome[,1])

data_targetome <- as.matrix(cbind(rownames(targetome), targetome))
colnames(data_targetome) <- c('Gene.Symbol','targetome')



########################################################################
## Import 3' UTR length
human_UTR <- read.table('R.data/TargetScan/UTR_Sequences_human_with_lenght.txt',
                        header = TRUE)

#####################################
## keep only one length of each genes
genes_UTR_red <- matrix(NA, ncol =6, nrow =0)
colnames(genes_UTR_red) <- colnames(human_UTR)

## keep highest length
for(gene in unique(human_UTR$Gene.Symbol)){
  index <- which(human_UTR$Gene.Symbol == gene )
  index_max_length <- index[which(human_UTR$length[index] == max(human_UTR$length[index]))][1]
  genes_UTR_red <- rbind(genes_UTR_red, human_UTR[index_max_length,])
} 

dim(genes_UTR_red)

data_length_UTR <- genes_UTR_red[,c(3,6)]
data_targetome  <- cbind(rownames(targetome), targetome)
colnames(data_targetome) <- c('Gene.Symbol','number')

merge_data <- merge(data_length_UTR, data_targetome, all.x = TRUE, all.y =TRUE)


## Save figure 
file_output <- paste0("R.results/Supp_11_c_targetome_by_UTR_length.pdf")
pdf(file_output, width = 6, height = 6)

plot(merge_data$length, merge_data$number, main ="Number of miR targeting a gene by the 3' UTR length" , 
     xlab =" gene 3'UTR length", ylab ="number of miR targeting a gene")

dev.off()


###########################################################################
## Cluster cells using miR-92a-3p expression
miRNA <- 'hsa-miR-92a-3p'
exp_miRNA <- 2^data_miRNA[which(rownames(data_miRNA) == miRNA),]
hist(as.numeric(exp_miRNA), xlab ='2^ miR-92a-3p expression', breaks =10, 
     main ='histogram of miR-92a-3p expression')


low <- which(exp_miRNA < 0.24)
high <- which(exp_miRNA > 0.24)

length(low)
length(high)


##########################################################################
## Log fold change of means between the two groups
mean_low <- apply(data_RNA[,low],1 ,function(x) mean(x,na.rm =T))
mean_high <- apply(data_RNA[,high],1 ,function(x) mean(x,na.rm =T))

LFC <- as.data.frame(mean_low - mean_high)
max(LFC, na.rm =T)

targets_92 <- targets_selection('hsa-miR-92a-3p')
LFC_targets <- LFC[which(rownames(LFC) %in% targets_92),]
LFC_non_targets <- LFC[-which(rownames(LFC) %in% targets_92),]

data_LFC <- cbind(rownames(LFC), LFC)
colnames(data_LFC) <- c('Gene.Symbol','LFC')


LFC_UTR <- merge(data_LFC, human_UTR, all.x = TRUE, all.y =FALSE)
dim(data_LFC)
dim(human_UTR)
dim(LFC_UTR)


less_in_high <- rownames(LFC)[which(LFC >0)]
less_in_low <- rownames(LFC)[which(LFC <0)]
length(less_in_low)
length(less_in_high)


UTR_less_in_high <- human_UTR$length [human_UTR$Gene.Symbol %in% less_in_high]
UTR_less_in_low <- human_UTR$length [human_UTR$Gene.Symbol %in% less_in_low]


## Save figure 
file_output <- paste0("R.results/Supp_11_d_KS_test_UTR_length_by_groups.pdf")
pdf(file_output, width = 6, height = 6)


## Cummulative distribution with log scale
ks_res <- ks.test(UTR_less_in_low,UTR_less_in_high)
plot(sort(UTR_less_in_low) , 1-ecdf(UTR_less_in_low)(sort(UTR_less_in_low, decreasing =T) ), 
     log="x", col ='darkviolet', xlim = c(50,10000), cex =0.1, 
     main =paste('Cummulative distributon of UTR length\n p val KS =',s3(ks_res$p.value)),
     ylab ='Cummulative distribution', xlab ="3'UTR length (number of bases)")
points(sort(UTR_less_in_high) , 1-ecdf(UTR_less_in_high)(sort(UTR_less_in_high, decreasing = T) ),
     col ='chartreuse4', cex =0.2)
legend (x ='topleft',col =c('darkviolet','chartreuse4'), lwd =2, cex =1,
        legend =c('Log Fold Change < 0',
                  'Log Fold change > 0'))

dev.off()
