#inhibition analysis using TAK-981 inhibition dataset on pancreatic and colorectal cancer cells


library(dplyr)


#import mutation info
ccl_mut <- read.csv("C:/Users/james/Desktop/Sumoylation_Analysis/data/CCLE_mutations.csv") 
cclID <- read.csv("C:/Users/james/Desktop/Sumoylation_Analysis/data/sample_info.csv") #list of all cell lines, ID and lineage
inhibition <- read.csv("inhibition_data.csv", fileEncoding="UTF-8-BOM") #inhibition data 

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

nrow(distinct(mut_df, ccl, .keep_all = TRUE))
#total of 1694 unique cell lines

#extract colorectal and pancreatic cancer only
#mut_df <- mut_df[which(mut_df$lineage == "colorectal" | mut_df$lineage == "pancreas"), ]
#nrow(distinct(mut_df, ccl, .keep_all = TRUE))
#total of 130 unique colorectal/pancreatic cell lines

#now process inhibition data by eliminating all spaces and special characters
inhibition$cell <- gsub(" ", "", inhibition$cell)
inhibition$cell <- gsub("-", "", inhibition$cell)
inhibition$cell <- gsub("\\.", "", inhibition$cell) #to remove dot, add 2 backslashes

inhibition$cell <- toupper(inhibition$cell)


colnames(inhibition)[which(names(inhibition) == "cell")] <- "ccl"


final_df <- merge(x = mut_df,
                  y = inhibition,
                  by = "ccl")

final_df <- unique(final_df)
distinct(final_df, ccl, .keep_all = TRUE)
nrow(distinct(final_df, ccl, .keep_all = TRUE))
#30 unique cell lines
#DLD1 is the only cell line that doesn't exist in the BROAD institute's dataset
#COLO741 is a colon carcinoma and therefore categorized as "skin" cancer. 

#all the unique cell lines
distinct(final_df, lineage)

colo_in <- final_df[final_df$lineage == "colorectal", ]
nrow(distinct(colo_in, ccl))

#need to include skin cancer because of COLO741 being categorized as colon carcinoma
panc_in <- final_df[final_df$lineage == "pancreas" | final_df$lineage == "skin", ]
nrow(distinct(panc_in, ccl))



#list top 5 common mutations
sort(table(final_df$Hugo_Symbol),decreasing=TRUE)[1:5]

#divide them into two groups: colorectal and pancreatic
sort(table(colo_in$Hugo_Symbol),decreasing=TRUE)[1:10]
sort(table(panc_in$Hugo_Symbol),decreasing=TRUE)[1:10]



##################################################  Pancreatic Cancer ################################################## 


#common mutation with below average AUC value
below_avg_panc_AUC <- panc_in[panc_in$AUC < 400, ]
nrow(distinct(below_avg_panc_AUC, ccl))
#list the 10 most common mutations
panc_below_avg <- sort(table(below_avg_panc_AUC$Hugo_Symbol),decreasing=TRUE)[1:10]
panc_below_avg

#common mutation with above average AUC value
above_avg_panc_AUC <- panc_in[panc_in$AUC >= 400, ]
nrow(distinct(above_avg_panc_AUC, ccl))
#list the 10 most common mutations
panc_above_avg <- sort(table(above_avg_panc_AUC$Hugo_Symbol),decreasing=TRUE)[1:10]
panc_above_avg

# statistical significance of difference in AUC values between cells with that mutation and cells without the mutation

ttest_df <- panc_in
ttest_df$Hugo_Symbol <- as.character(ttest_df$Hugo_Symbol)
ttest_df$Hugo_Symbol[ttest_df$Hugo_Symbol != "TP53"] <- "Lacks Mut"
ttest_df$Hugo_Symbol[ttest_df$Hugo_Symbol == "TP53"] <- "Has Mut"
distinct(ttest_df, Hugo_Symbol)

#save t-test results as a list
t <- list()
t[[1]] <- t.test(AUC ~ Hugo_Symbol, data = ttest_df)

#save p-value and means of the two groups
tsum <- sapply(t, function(x) {
          c(x$estimate[1],
            x$estimate[2],
            p.value = x$p.value)
        })
tsum

#create a function to repeat t-test for other genes
ttest <- function(M, df){
  ttest_df <- df
  ttest_df$Hugo_Symbol <- as.character(ttest_df$Hugo_Symbol)
  ttest_df$Hugo_Symbol[ttest_df$Hugo_Symbol != M] <- "Lacks Mut"
  ttest_df$Hugo_Symbol[ttest_df$Hugo_Symbol == M] <- "Has Mut"

  t <- list()
  t[[1]] <- t.test(AUC ~ Hugo_Symbol, data = ttest_df)
  tsum <- sapply(t, function(x) {
                  c(x$estimate[1],
                    x$estimate[2],
                    p.value = x$p.value)
                })
  return(tsum)
}

ttest("KRAS", panc_in)

#create dataframe with top 10 mutations (with below 400 AUC)
panc_below_avg_df <- as.data.frame(panc_below_avg, stringsAsFactors=FALSE)
colnames(panc_below_avg_df) <- c('mutation', 'freq')

panc_below_avg_df$Mut_AUC_mean <- NA
panc_below_avg_df$others_AUC_mean <- NA
panc_below_avg_df$p.value <- NA

#for loop to fill in the dataframe 
for (i in 1:nrow(panc_below_avg_df)) {
  print(panc_below_avg_df$mutation[i])
  t <- ttest(panc_below_avg_df$mutation[i], panc_in)
  panc_below_avg_df$Mut_AUC_mean[i] <- t[1]
  panc_below_avg_df$others_AUC_mean[i] <- t[2]
  panc_below_avg_df$p.value[i] <- round(t[3], digits=3)
}

write.csv(panc_below_avg_df, "PANC_mutation.csv", row.names=FALSE)



##################################################  Colorectal Cancer ################################################## 

#common mutation with below average AUC value
below_avg_colo_AUC <- colo_in[colo_in$AUC < 350, ]
nrow(distinct(below_avg_colo_AUC, ccl))
#list the 10 most common mutations
colo_below_avg <- sort(table(below_avg_colo_AUC$Hugo_Symbol),decreasing=TRUE)[1:10]
colo_below_avg

#common mutation with above average AUC value
above_avg_colo_AUC <- colo_in[colo_in$AUC >= 350, ]
nrow(distinct(above_avg_colo_AUC, ccl))
#list the 10 most common mutations
colo_above_avg <- sort(table(above_avg_colo_AUC$Hugo_Symbol),decreasing=TRUE)[1:10]
colo_above_avg


#use ttest function that we created during pancreas analysis
ttest("KRAS", colo_in)

#create dataframe with top 10 mutations (with below 400 AUC)
colo_below_avg_df <- as.data.frame(colo_below_avg, stringsAsFactors=FALSE)
colnames(colo_below_avg_df) <- c('mutation', 'freq')

colo_below_avg_df$Mut_AUC_mean <- NA
colo_below_avg_df$others_AUC_mean <- NA
colo_below_avg_df$p.value <- NA

#for loop to fill in the dataframe 
for (i in 1:nrow(colo_below_avg_df)) {
  print(colo_below_avg_df$mutation[i])
  tsum <- ttest(colo_below_avg_df$mutation[i], colo_in)
  colo_below_avg_df$Mut_AUC_mean[i] <- tsum[1]
  colo_below_avg_df$others_AUC_mean[i] <- tsum[2]
  colo_below_avg_df$p.value[i] <- round(tsum[3], digits=3)
}
colo_below_avg_df
write.csv(colo_below_avg_df, "COLO_mutation.csv", row.names=FALSE)



######################################################################################################################


# Additional mutation analysis
# 1. RAD
# 2. XRCC
# 3. BRCA
# 4. ATR


#create subset for RAD genes within pancreatic cells. 
#subset of any gene that starts with RAD
panc_RAD <- droplevels(panc_in[grepl("^RAD", panc_in$Hugo_Symbol), ])
mut_list <- panc_RAD$Hugo_Symbol

panc_RAD_ttest <- panc_in
panc_RAD_ttest$Hugo_Symbol <- as.character(panc_RAD_ttest$Hugo_Symbol)

#change all RAD family genes to RAD
panc_RAD_ttest$Hugo_Symbol[grep(paste0("^", mut_list, collapse="|"), panc_RAD_ttest$Hugo_Symbol)] <- "RAD"
panc_RAD_ttest$Hugo_Symbol[panc_RAD_ttest$Hugo_Symbol != "RAD"] <- "Lacks Mut"
distinct(panc_RAD_ttest, Hugo_Symbol)

t.test(AUC ~ Hugo_Symbol, data=panc_RAD_ttest)



#BRCA in pancreatic
panc_BRCA <- droplevels(panc_in[grepl("^BRCA", panc_in$Hugo_Symbol), ])
mut_list <- panc_BRCA$Hugo_Symbol

panc_BRCA_ttest <- panc_in
panc_BRCA_ttest$Hugo_Symbol <- as.character(panc_BRCA_ttest$Hugo_Symbol)

#change all RAD family genes to RAD
panc_BRCA_ttest$Hugo_Symbol[grep(paste0("^", mut_list, collapse="|"), panc_BRCA_ttest$Hugo_Symbol)] <- "BRCA"
panc_BRCA_ttest$Hugo_Symbol[panc_BRCA_ttest$Hugo_Symbol != "BRCA"] <- "Lacks Mut"
distinct(panc_BRCA_ttest, Hugo_Symbol)

t.test(AUC ~ Hugo_Symbol, data=panc_BRCA_ttest)



#MUC in pancreatic
panc_MUC <- droplevels(panc_in[grepl("^MUC", panc_in$Hugo_Symbol), ])
mut_list <- as.vector(panc_MUC$Hugo_Symbol)

panc_MUC_ttest <- panc_in
panc_MUC_ttest$Hugo_Symbol <- as.character(panc_MUC_ttest$Hugo_Symbol)

#change all RAD family genes to RAD
panc_MUC_ttest$Hugo_Symbol[grep(paste0("^", mut_list, collapse="|"), panc_MUC_ttest$Hugo_Symbol)] <- "MUC"
panc_MUC_ttest$Hugo_Symbol[panc_MUC_ttest$Hugo_Symbol != "MUC"] <- "Lacks Mut"
distinct(panc_MUC_ttest, Hugo_Symbol)

t.test(AUC ~ Hugo_Symbol, data=panc_MUC_ttest)


#MYC in pancreatic
panc_MYC <- droplevels(panc_in[grepl("^MYC", panc_in$Hugo_Symbol), ])
mut_list <- as.vector(panc_MYC$Hugo_Symbol)

panc_MYC_ttest <- panc_in
panc_MYC_ttest$Hugo_Symbol <- as.character(panc_MYC_ttest$Hugo_Symbol)

#change all RAD family genes to RAD
panc_MYC_ttest$Hugo_Symbol[grep(paste0("^", mut_list, collapse="|"), panc_MYC_ttest$Hugo_Symbol)] <- "MYC"
panc_MYC_ttest$Hugo_Symbol[panc_MYC_ttest$Hugo_Symbol != "MYC"] <- "Lacks Mut"
distinct(panc_MYC_ttest, Hugo_Symbol)

t.test(AUC ~ Hugo_Symbol, data=panc_MYC_ttest)


#FANC in pancreatic
panc_FANC <- droplevels(panc_in[grepl("^FANC", panc_in$Hugo_Symbol), ])
mut_list <- as.vector(panc_FANC$Hugo_Symbol)

panc_FANC_ttest <- panc_in
panc_FANC_ttest$Hugo_Symbol <- as.character(panc_FANC_ttest$Hugo_Symbol)

#change all RAD family genes to RAD
panc_FANC_ttest$Hugo_Symbol[grep(paste0("^", mut_list, collapse="|"), panc_FANC_ttest$Hugo_Symbol)] <- "FANC"
panc_FANC_ttest$Hugo_Symbol[panc_FANC_ttest$Hugo_Symbol != "FANC"] <- "Lacks Mut"
distinct(panc_FANC_ttest, Hugo_Symbol)

t.test(AUC ~ Hugo_Symbol, data=panc_FANC_ttest)



#Include other genes in BRCA group and re-do t-test

#BRCA + ATR
panc_BRCA <- droplevels(panc_in[grepl("^BRCA", panc_in$Hugo_Symbol), ])
mut_list <- panc_BRCA$Hugo_Symbol

panc_BRCA_ttest <- panc_in
panc_BRCA_ttest$Hugo_Symbol <- as.character(panc_BRCA_ttest$Hugo_Symbol)

panc_BRCA_ttest$Hugo_Symbol[grep(paste0("^", mut_list, collapse="|"), panc_BRCA_ttest$Hugo_Symbol)] <- "Has_Mut"
panc_BRCA_ttest$Hugo_Symbol[panc_BRCA_ttest$Hugo_Symbol == "ATR"] <- "Has_Mut"
panc_BRCA_ttest$Hugo_Symbol[panc_BRCA_ttest$Hugo_Symbol != "Has_Mut"] <- "Lacks Mut"
distinct(panc_BRCA_ttest, Hugo_Symbol)

t.test(AUC ~ Hugo_Symbol, data=panc_BRCA_ttest)


#BRCA + RAD family
panc_BRCA <- droplevels(panc_in[grepl("^BRCA", panc_in$Hugo_Symbol), ])
mut_list <- panc_BRCA$Hugo_Symbol

panc_BRCA_ttest <- panc_in
panc_BRCA_ttest$Hugo_Symbol <- as.character(panc_BRCA_ttest$Hugo_Symbol)

panc_BRCA_ttest$Hugo_Symbol[grep(paste0("^", mut_list, collapse="|"), panc_BRCA_ttest$Hugo_Symbol)] <- "Has_Mut"
panc_BRCA_ttest$Hugo_Symbol[grepl("^RAD", panc_in$Hugo_Symbol)] <- "Has_Mut"
panc_BRCA_ttest$Hugo_Symbol[panc_BRCA_ttest$Hugo_Symbol != "Has_Mut"] <- "Lacks Mut"
distinct(panc_BRCA_ttest, Hugo_Symbol)

t.test(AUC ~ Hugo_Symbol, data=panc_BRCA_ttest)


#BRCA + FANC family
panc_BRCA <- droplevels(panc_in[grepl("^BRCA", panc_in$Hugo_Symbol), ])
mut_list <- panc_BRCA$Hugo_Symbol

panc_BRCA_ttest <- panc_in
panc_BRCA_ttest$Hugo_Symbol <- as.character(panc_BRCA_ttest$Hugo_Symbol)

panc_BRCA_ttest$Hugo_Symbol[grep(paste0("^", mut_list, collapse="|"), panc_BRCA_ttest$Hugo_Symbol)] <- "Has_Mut"
panc_BRCA_ttest$Hugo_Symbol[grepl("^FANCI", panc_in$Hugo_Symbol)] <- "Has_Mut"
panc_BRCA_ttest$Hugo_Symbol[panc_BRCA_ttest$Hugo_Symbol != "Has_Mut"] <- "Lacks Mut"
distinct(panc_BRCA_ttest, Hugo_Symbol)

t.test(AUC ~ Hugo_Symbol, data=panc_BRCA_ttest)


#BRCA + FANC + ATR + RAD family
panc_all <- droplevels(panc_in[grepl("^BRCA", panc_in$Hugo_Symbol), ])
mut_list <- panc_BRCA$Hugo_Symbol

panc_all_ttest <- panc_in
panc_all_ttest$Hugo_Symbol <- as.character(panc_all_ttest$Hugo_Symbol)

panc_all_ttest$Hugo_Symbol[grep(paste0("^", mut_list, collapse="|"), panc_all_ttest$Hugo_Symbol)] <- "Has_Mut"
panc_all_ttest$Hugo_Symbol[grepl("^FANC", panc_in$Hugo_Symbol)] <- "Has_Mut"
panc_all_ttest$Hugo_Symbol[grepl("^RAD", panc_in$Hugo_Symbol)] <- "Has_Mut"
panc_all_ttest$Hugo_Symbol[grepl("^ATR", panc_in$Hugo_Symbol)] <- "Has_Mut"
panc_all_ttest$Hugo_Symbol[panc_all_ttest$Hugo_Symbol != "Has_Mut"] <- "Lacks Mut"
distinct(panc_all_ttest, Hugo_Symbol)

t.test(AUC ~ Hugo_Symbol, data=panc_all_ttest)


# combination group plot (BRCA+FANC+RAD+ATR)
library(ggplot2)

enrich_title <- paste("AUC of cell lines with mutations in", colorectal$lineage ,"cancer cells \n (",
                      length(which(colorectal$UBA2_gene_effect < colorectal_average)),
                      "cell lines out of ",length(colorectal$UBA2_gene_effect)," with UBA2 score <" ,
                      round(colorectal_average, digits = 3), ")")


combo_boxplot <- ggplot(panc_all_ttest, aes(x=Hugo_Symbol, y=AUC, fill=Hugo_Symbol)) + 
  ggtitle("AUC of cell lines with mutations in \nBRCA+FANC+RAD+ATR \n(p-value = 0.011)") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  geom_boxplot() + 
  theme_bw() +
  theme(plot.title = element_text(size=12))

combo_boxplot





############################################################################################################

mut_list <- c("KRAS", "FANC", "ATR", "RAD", "BRCA")
plot_group <- panc_in[grep(paste0("^", mut_list, collapse="|"), panc_in$Hugo_Symbol), ]
plot_group$Hugo_Symbol <- as.character(plot_group$Hugo_Symbol)

plot_group$Hugo_Symbol[grepl("^FANC", plot_group$Hugo_Symbol)] <- "FANC"
plot_group$Hugo_Symbol[grepl("^ATR", plot_group$Hugo_Symbol)] <- "ATR"
plot_group$Hugo_Symbol[grepl("^RAD", plot_group$Hugo_Symbol)] <- "RAD"
plot_group$Hugo_Symbol[grepl("^BRCA", plot_group$Hugo_Symbol)] <- "BRCA"

plot_group$Hugo_Symbol <- as.factor(plot_group$Hugo_Symbol)

#fix the ordering
plot_group$Hugo_Symbol <- factor(plot_group$Hugo_Symbol, 
                                 levels = c("KRAS", "FANC", "RAD","ATR","BRCA"), ordered = TRUE)


library(ggplot2)
boxplot <- ggplot(plot_group, aes(x=Hugo_Symbol, y=AUC, fill=Hugo_Symbol)) + 
  ggtitle("AUC distributions of cells with significant mutations") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  geom_boxplot() + theme_bw()

boxplot



#p-value dot plot

pvalue_group <- read.csv("C:/Users/james/Desktop/Sumoylation_Analysis/project/inhibition_analysis/mutation_pvalues.csv", fileEncoding="UTF-8-BOM")


pvalue_group$gene <- factor(pvalue_group$gene, 
                                 levels = c("KRAS", "FANC", "RAD","ATR","BRCA", "BRCA+ATR", "BRCA+RAD", "BRCA+FANC", "BRCA+FANC+RAD+ATR"), ordered = TRUE)

pvalue_plot <- ggplot(pvalue_group, aes(x=gene, y=p.value)) + 
  ggtitle("P-values of significant mutations") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  geom_point() + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank())

pvalue_plot







############################################################################################################


#RAD in colorectal
colo_RAD <- colo_in[grepl("^RAD51", colo_in$Hugo_Symbol), ]
mut_list <- as.vector(colo_RAD$Hugo_Symbol)
#remove RADIL since it's not a RAD gene
mut_list <- mut_list[mut_list != "RADIL"]

colo_RAD_ttest <- colo_in
colo_RAD_ttest$Hugo_Symbol <- as.character(colo_RAD_ttest$Hugo_Symbol)

#change all RAD family genes to RAD
colo_RAD_ttest$Hugo_Symbol[grep(paste0("^", mut_list, collapse="|"), colo_RAD_ttest$Hugo_Symbol)] <- "RAD"
colo_RAD_ttest$Hugo_Symbol[colo_RAD_ttest$Hugo_Symbol != "RAD"] <- "Lacks Mut"
distinct(colo_RAD_ttest, Hugo_Symbol)

t.test(AUC ~ Hugo_Symbol, data=colo_RAD_ttest)



#XRCC in colorectal 
colo_XRCC <- colo_in[grepl("^XRCC", colo_in$Hugo_Symbol), ]
mut_list <- as.vector(colo_XRCC$Hugo_Symbol)

colo_XRCC_ttest <- colo_in
colo_XRCC_ttest$Hugo_Symbol <- as.character(colo_XRCC_ttest$Hugo_Symbol)

#change all RAD family genes to RAD
colo_XRCC_ttest$Hugo_Symbol[grep(paste0("^", mut_list, collapse="|"), colo_XRCC_ttest$Hugo_Symbol)] <- "XRCC"
colo_XRCC_ttest$Hugo_Symbol[colo_XRCC_ttest$Hugo_Symbol != "XRCC"] <- "Lacks Mut"
distinct(colo_XRCC_ttest, Hugo_Symbol)

t.test(AUC ~ Hugo_Symbol, data=colo_XRCC_ttest)



#ATR in colorectal 
colo_ATR <- colo_in[grepl("^ATR", colo_in$Hugo_Symbol), ]
mut_list <- as.vector(colo_ATR$Hugo_Symbol)

#delete ATRNL1 since its not a ATR gene
mut_list <- mut_list[mut_list != "ATRNL1"]

#
mut_list <- mut_list[mut_list != "ATRN"]

#
mut_list <- mut_list[mut_list != "ATRX"]

colo_ATR_ttest <- colo_in
colo_ATR_ttest$Hugo_Symbol <- as.character(colo_ATR_ttest$Hugo_Symbol)

#change all RAD family genes to RAD
colo_ATR_ttest$Hugo_Symbol[grep(paste0("^", mut_list, collapse="|"), colo_ATR_ttest$Hugo_Symbol)] <- "ATR"
colo_ATR_ttest$Hugo_Symbol[colo_ATR_ttest$Hugo_Symbol != "ATR"] <- "Lacks Mut"
distinct(colo_ATR_ttest, Hugo_Symbol)

t.test(AUC ~ Hugo_Symbol, data=colo_ATR_ttest)




#BRCA in colorectal 
colo_BRCA <- colo_in[grepl("^BRCA", colo_in$Hugo_Symbol), ]
mut_list <- as.vector(colo_BRCA$Hugo_Symbol)


colo_BRCA_ttest <- colo_in
colo_BRCA_ttest$Hugo_Symbol <- as.character(colo_BRCA_ttest$Hugo_Symbol)

#change all RAD family genes to RAD
colo_BRCA_ttest$Hugo_Symbol[grep(paste0("^", mut_list, collapse="|"), colo_BRCA_ttest$Hugo_Symbol)] <- "BRCA"
colo_BRCA_ttest$Hugo_Symbol[colo_BRCA_ttest$Hugo_Symbol != "BRCA"] <- "Lacks Mut"
distinct(colo_BRCA_ttest, Hugo_Symbol)

t.test(AUC ~ Hugo_Symbol, data=colo_BRCA_ttest)



#MUC in colorectal 
colo_MUC <- colo_in[grepl("^MUC", colo_in$Hugo_Symbol), ]
mut_list <- as.vector(colo_MUC$Hugo_Symbol)


colo_MUC_ttest <- colo_in
colo_MUC_ttest$Hugo_Symbol <- as.character(colo_MUC_ttest$Hugo_Symbol)

#change all RAD family genes to RAD
colo_MUC_ttest$Hugo_Symbol[grep(paste0("^", mut_list, collapse="|"), colo_MUC_ttest$Hugo_Symbol)] <- "MUC"
colo_MUC_ttest$Hugo_Symbol[colo_MUC_ttest$Hugo_Symbol != "MUC"] <- "Lacks Mut"
distinct(colo_MUC_ttest, Hugo_Symbol)

t.test(AUC ~ Hugo_Symbol, data=colo_MUC_ttest)

