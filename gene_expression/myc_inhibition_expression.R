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
names(exp_merged)


#c-MYC subset
MYC_exp <- exp_merged[names(exp_merged)=="MYC"]

MYC_exp <- exp_merged %>% select(ID, ccl, lineage, MYC)

panc_MYC_exp <- MYC_exp[MYC_exp$lineage == "pancreas", ]
colo_MYC_exp <- MYC_exp[MYC_exp$lineage == "colorectal", ]


#inhibition data 
inhibition <- read.csv("inhibition_data.csv", fileEncoding="UTF-8-BOM") 
#now process inhibition data by eliminating all spaces and special characters
inhibition$cell <- gsub(" ", "", inhibition$cell)
inhibition$cell <- gsub("-", "", inhibition$cell)
inhibition$cell <- gsub("\\.", "", inhibition$cell) #to remove dot, add 2 backslashes
inhibition$cell <- toupper(inhibition$cell)
colnames(inhibition)[which(names(inhibition) == "cell")] <- "ccl"
inhibition <- inhibition[,-2]
inhibition <- inhibition[,-2]

panc_MYC_exp <- merge(x= panc_MYC_exp, 
                      y = inhibition,
                      by = "ccl")

colo_MYC_exp <- merge(x= colo_MYC_exp, 
                      y = inhibition,
                      by = "ccl")

both_MYC_exp <- rbind(panc_MYC_exp, colo_MYC_exp)

summary(both_MYC_exp)

both_MYC_exp <- both_MYC_exp[-which(both_MYC_exp$AUC == 0), ]


library(MASS)

mod1 <- lm(MYC ~ AUC, data=both_MYC_exp)
summary(mod1)

mod3 <- lm(MYC ~ AUC + I(AUC^2) + I(AUC^3), data=both_MYC_exp)
summary(mod3)

mod5 <- lm(MYC ~ AUC + I(AUC^2) + I(AUC^3 ) + I(AUC^4) + I(AUC^5), data=both_MYC_exp)
summary(mod5)

ldat <- data.frame(both_MYC_exp,
                   mod1 = predict(mod1),
                   mod3 = predict(mod3),
                   mod5 = predict(mod5))

ggplot(both_MYC_exp, aes(x = AUC, y = MYC)) + 
  geom_point() + 
  geom_line(color = "red", data = ldat, aes(x = AUC, y= mod5)) + 
  labs(title = "AUC vs MYC Expression Level (Polynomial ^5)", x = "inhibition AUC", y= "MYC expression level (log2)") + 
  theme_bw()












