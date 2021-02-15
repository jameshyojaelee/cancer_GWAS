####################################################################################################################################
#classification model (glm) and SVM

library(dplyr)
library(ggplot2)

#gene expression data for just protein coding genes
exp <- read.csv("C:/Users/james/Desktop/Sumoylation_Analysis/data/CCLE_expression.csv", fileEncoding="UTF-8-BOM") 
colnames(exp)[1] <- "ID"

cclID <- read.csv("C:/Users/james/Desktop/Sumoylation_Analysis/data/sample_info.csv") #list of all cell lines, ID and lineage

ccl_df <- data.frame(ID = cclID$DepMap_ID, ccl = cclID$stripped_cell_line_name, lineage = cclID$lineage)

exp_merged <- merge(x = ccl_df, 
                    y = exp,
                    by = "ID")

names(exp_merged) <- sub("\\..*", "", names(exp_merged))

UBA2_exp <- data.frame("ID" = exp_merged[, "ID"], "ccl" = exp_merged[, "ccl"], "lineage"= exp_merged[, "lineage"] ,"UBA2_exp" = exp_merged[, "UBA2"])


gene_effect <- read.csv("C:/Users/james/Desktop/Sumoylation_Analysis/data/Achilles_gene_effect.csv") #CRISPR(Avana) data
names(gene_effect) <- sub("\\..*", "", names(gene_effect))
names(gene_effect)[1] <- "ID"

UBA2_GE <- data.frame("ID" = gene_effect[, "ID"],"UBA2_KO" = gene_effect[, "UBA2"])

df <- merge(x = UBA2_exp,
            y= UBA2_GE,
            by = "ID")

ggplot(data = df, aes(x = UBA2_KO, y= UBA2_exp, color=lineage)) + geom_point()


