###########################################################################
## Script for Matter Arising
## Prepare downloaded dataset from Wang et al. paper
###########################################################################

## Importation ------------------------------------------------------------
source("Functions.R")# functions importation

## Single-cell data importation 
data_RNA <- 
  read.table("R.Data/GSE114071_NW_K562_RNA_processed.gct", 
             sep = '\t', skip = 1, header = TRUE)
data_miRNA <- 
  read.table("R.Data/GSE114071_NW_K562_miRNA_log2.gct", 
             sep = '\t', header = TRUE )

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

dim(data_RNA_19_reduce)
dim(data_miRNA_19)
class(data_RNA_19_reduce)
class(data_miRNA_19)

## save data in RDS files
saveRDS (data_RNA_19_reduce, "R.Data/data_RNA_19.rds")
saveRDS (data_miRNA_19, "R.Data/data_miRNA_19.rds")
