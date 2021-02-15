#Pathway Analysis

library(BiocManager)
library(DESeq2)
library(pathview)
library(gage)
library(gageData)
library("AnnotationDbi")
library("org.Hs.eg.db")

columns(org.Hs.eg.db)


#import mutation info
ccl_mut <- read.csv("C:/Users/james/Desktop/Sumoylation_Analysis/data/CCLE_mutations.csv") 
cclID <- read.csv("C:/Users/james/Desktop/Sumoylation_Analysis/data/sample_info.csv") #list of all cell lines, ID and lineage
inhibition <- read.csv("C:/Users/james/Desktop/Sumoylation_Analysis/project/inhibition/inhibition_data.csv", fileEncoding="UTF-8-BOM") #inhibition data 

#dataframe with the cancer cell line names
ccl_df <- data.frame(ID = cclID$DepMap_ID, ccl = cclID$stripped_cell_line_name, lineage = cclID$lineage)


#change the column name of DepMap_ID to ID
colnames(ccl_mut)[which(names(ccl_mut) == "DepMap_ID")] <- "ID"
colnames(ccl_df)[which(names(ccl_df) == "DepMap_ID")] <- "ID"

#subset 3 columns only
mut <- subset(ccl_mut, select = c("ID", "Hugo_Symbol", "Variant_Classification"))

#import the auc data to manipulate
mut_df <- merge(x = ccl_df, 
                y = mut,
                by = "ID")

#eliminate silent mutation
mut_df <- mut_df[mut_df$Variant_Classification != "Silent", ]
mut_df <- mut_df[which(mut_df$lineage == "colorectal" | mut_df$lineage == "pancreas"), ]


#now process inhibition data by eliminating all spaces and special characters
inhibition$cell <- gsub(" ", "", inhibition$cell)
inhibition$cell <- gsub("-", "", inhibition$cell)
inhibition$cell <- gsub("\\.", "", inhibition$cell) #to remove dot, add 2 backslashes

colnames(inhibition)[which(names(inhibition) == "cell")] <- "ccl"

final_df <- merge(x = mut_df,
                  y = inhibition,
                  by = "ccl")

final_df$KEGG <- NA

colo_in <- final_df[final_df$lineage == "colorectal", ]
panc_in <- final_df[final_df$lineage == "pancreas", ]


for(i in 1:nrow(final_df)){
  ourID <- final_df$Hugo_Symbol[i]
  KEGG_ID <- mapIds(org.Hs.eg.db,
                    keys=as.character(ourID), #Our genenames
                    keytype="ALIAS",        # The format of our genenames
                    column="PATH",          # The new format we want to add
                    multiVals="first")
  final_df$KEGG[i] <- KEGG_ID
}


mapIds(org.Hs.eg.db,
       keys="TTN", # Our genenames
       keytype="ALIAS",        # The format of our genenames
       column="PATH",          # The new format we want to add
       multiVals="first")

data(kegg.sets.hs)
head(kegg.sets.hs, 2)

pathview(panc_in$AUC, pathway.id="05410")
