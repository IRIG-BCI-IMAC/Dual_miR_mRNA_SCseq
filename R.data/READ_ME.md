Data used for analysis:
miRNA and mRNA seq data from Wang et al. paper available on GEO under series [GSE114071](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114071) 
Those data are in rds file in this folder.


TargetScan data were downloaded from the TargetScan website https://www.targetscan.org/vert_80/ and https://www.targetscan.org/vert_71/ respectively for v8.0 and v7.1. 

For the v8.0 we used human data from both default predictions and all predictions files: 
https://www.targetscan.org/vert_80/vert_80_data_download/Summary_Counts.default_predictions.txt.zip  
and https://www.targetscan.org/vert_80/vert_80_data_download/Summary_Counts.all_predictions.txt.zip

For the v7.1 we used human data from both default predictions and all predictions files : 
https://www.targetscan.org/vert_71/vert_71_data_download/Summary_Counts.default_predictions.txt.zip  
and https://www.targetscan.org/vert_71/vert_71_data_download/Summary_Counts.all_predictions.txt.zip


TargetScan family information were downloaded for v8.0 and v7.1 respectively from:
https://www.targetscan.org/vert_80/vert_80_data_download/miR_Family_Info.txt.zip  
https://www.targetscan.org/vert_71/vert_71_data_download/miR_Family_Info.txt.zip

For all TargetScan data, only information relating to human has been kept.
This selection was made using the function `targetscan2rds` in the `Functions.R` file. 
