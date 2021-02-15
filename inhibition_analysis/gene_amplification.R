#myc gene copy number dataset
#pre-processed in Python
myc_cn <- read.csv('C:/Users/james/Desktop/Sumoylation_Analysis/project/inhibition_analysis/myc_copy_number.csv', fileEncoding="UTF-8-BOM") 
inhibition <- read.csv("inhibition_data.csv", fileEncoding="UTF-8-BOM") #inhibition data 

#now process inhibition data by eliminating all spaces and special characters
inhibition$cell <- gsub(" ", "", inhibition$cell)
inhibition$cell <- gsub("-", "", inhibition$cell)
inhibition$cell <- gsub("\\.", "", inhibition$cell) #to remove dot, add 2 backslashes
inhibition$cell <- toupper(inhibition$cell)

inhibition

myc_cn_AUC <- merge(x = myc_cn, 
                y = inhibition,
                by = "cell")

distinct(myc_cn_AUC, cell)

#drop duplicate column "type"
myc_cn_AUC <- myc_cn_AUC[, !names(myc_cn_AUC) %in% c("type")]

head(myc_cn_AUC)

mean(myc_cn_AUC$MYC)
min(myc_cn_AUC$MYC)
max(myc_cn_AUC$MYC)

#Set different thresholds for MYC copy number and observe the mean AUC values. 
above_3cn_myc <- myc_cn_AUC[myc_cn_AUC$MYC > 3, ]
#there is only 1 observation

above_2cn_myc <- myc_cn_AUC[myc_cn_AUC$MYC >= 2, ]
below_2cn_myc <- myc_cn_AUC[myc_cn_AUC$MYC < 2, ]
mean(above_2cn_myc$AUC)
mean(below_2cn_myc$AUC)



