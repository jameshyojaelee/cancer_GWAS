#pancreas 

inhibition_scores <- read.csv("pancreatic_inhibition.csv")


AUC <- data.frame(ccl = inhibition_scores$ï..ccl, AUC = inhibition_scores$AUC)

#Merge cell line info and inhibition AUC score
AUC_merged <- merge(x = ccl_df, 
                     y = AUC,
                     by = "ccl")
AUC_merged <- droplevels(AUC_merged)


pancreatic_AUC <- AUC_merged[which(AUC_merged$lineage == "pancreas"),]
pancreatic_AUC <- pancreatic_AUC[order(pancreatic_AUC$AUC), ]
pancreatic_scores <-pancreatic_AUC$AUC
pancreatic_average <- mean(pancreatic_scores)
pancreatic_sd <- sd(pancreatic_scores)
min(pancreatic_scores)
max(pancreatic_scores)
length(pancreatic_AUC$AUC)
length(which(pancreatic_AUC$AUC < pancreatic_average))

enrich_title <- paste("UBA2 gene effect in", pancreatic_AUC$lineage ,"cancer cells \n (",
                      length(which(pancreatic_AUC$AUC < pancreatic_average)),
                      "cell lines out of ",length(pancreatic_AUC$AUC)," with UBA2 score <" ,
                      round(pancreatic_average, digits = 3), ")")

library(ggplot2)
library("scales")
enrich <- ggplot(data = pancreatic_AUC , aes(x = AUC, y = lineage, width= 100)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = AUC)) +
  scale_fill_gradientn(colours = c("red","red","white","white"),
                       guide = "colorbar", limits=c(min(pancreatic_scores),800)) +
  scale_x_continuous(limits=c(min(pancreatic_scores),1000))+
  theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=11, angle = 90, vjust = 3, hjust= 1),
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())
enrich

#################################################################################################
########example with pancreatic cells
AUC_merged <- merge(x = AUC_merged, 
                     y = ccl_mut,
                     by = "ID")

#calculate the number of unique ccl
length(unique(AUC_merged$ccl))
#turns out to be 13 ccls

#subset of pancreatic cell lines with below threshold gene effect score
pancreatic_valid_mut <- AUC_merged[which(AUC_merged$AUC <= pancreatic_average), ]
pancreatic_valid_mut <- droplevels(pancreatic_valid_mut)
#subset of cell lines with above threshold
pancreatic_invalid_mut <- AUC_merged[which(AUC_merged$AUC > pancreatic_average), ]
pancreatic_invalid_mut <- droplevels(pancreatic_invalid_mut)

#eliminate mutations that also occur in invalid
#pancreatic_valid_mut1 <- pancreatic_valid_mut[!duplicated(pancreatic_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
pancreatic_mut_final <- anti_join(pancreatic_valid_mut, pancreatic_invalid_mut, 
                                  by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

pancreatic_mut_final <- pancreatic_mut_final[which(pancreatic_mut_final$Variant_Classification != "Silent"), ]
pancreatic_mut_final <- unique(pancreatic_mut_final[, 1:6])
pancreatic_mut_final <- droplevels(pancreatic_mut_final)

#calculate the second most frequent value
xx <- sort(table(pancreatic_mut_final$Hugo_Symbol),decreasing=TRUE)[1:300]
xx
#               

###########################################################export###########################################################
##### eliminate mutations that occur less than 4
pancreatic_mut_final <- pancreatic_mut_final[pancreatic_mut_final$Hugo_Symbol %in% names(which(table(pancreatic_mut_final$Hugo_Symbol) > 2)), ]
length(unique(pancreatic_mut_final$ccl)) # 18 cell lines 
length(unique(pancreatic_mut_final$Hugo_Symbol)) # 18 cell lines 
library(xlsx)
write.xlsx(pancreatic_mut_final, "pancreatic_AUC_mut.xlsx")
###########################################################export###########################################################


pancreatic_CASK <- pancreatic_mut_final[which(pancreatic_mut_final$Hugo_Symbol == "CASK"), ]
pancreatic_no_CASK <- pancreatic_mut_final[which(pancreatic_mut_final$Hugo_Symbol != "CASK"), ]

#find the most frequent variant type
calculate_mode(pancreatic_CASK$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with CASK_Missense
pancreatic_CASK_Mis <- pancreatic_CASK[which(pancreatic_CASK$Variant_Classification == "Missense_Mutation"), ]
pancreatic_CASK_No_Mis <- pancreatic_CASK[which(pancreatic_CASK$Variant_Classification != "Missense_Mutation"), ]
pancreatic_CASK <- droplevels(pancreatic_CASK)
#7 cell lines with CASK Missense mutation

#subset of all ccl
pancreatic_ccl <- distinct(pancreatic_mutation, ccl, .keep_all= TRUE )

y <- length(pancreatic_CASK_Mis$ID)
x <- length(pancreatic$ccl)
x <- x - y


final_table <- data.frame(mutation = c("CASK", "other"), counts = c(y, x))

bar_title <- paste("Mut: CASK Missense Mutation \n", y, "overlaps out of", x, "total")

bar <- ggplot(data = final_table, aes(x= mutation, y=counts, fill=mutation)) +
  ggtitle(bar_title) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("black", "darkgrey")) +
  theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
        panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                        size = 2, linetype = "solid"),
        legend.position = "none")


#######################################################################################################################
#(c) mut, lacks mut

#subset of all ccl
pancreatic_ccl <- distinct(pancreatic_mutation, ccl, .keep_all= TRUE )

pancreatic_ccl[match(pancreatic_CASK$ID, pancreatic_ccl$ID), ] <- pancreatic_CASK

pancreatic_ccl$Hugo_Symbol <- as.character(pancreatic_ccl$Hugo_Symbol)

pancreatic_ccl$Hugo_Symbol[pancreatic_ccl$Hugo_Symbol == "CASK"] <- "has Mut"
pancreatic_ccl$Hugo_Symbol[pancreatic_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =pancreatic_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(pancreatic_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
  ggtitle(box_title) +
  geom_boxplot(fill = "white", colour = "black", outlier.colour = "red")+
  theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
        axis.title.x=element_blank(),
        panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                        size = 2, linetype = "solid"))


library(gridExtra)

grid.arrange(grobs= list(enrich, bar, box), 
             layout_matrix = rbind(c(1, 1, 1, 1, 1, 1, 1, 1), 
                                   c(1, 1, 1, 1, 1, 1, 1, 1),
                                   c(1, 1, 1, 1, 1, 1, 1, 1),
                                   c(1, 1, 1, 1, 1, 1, 1, 1),
                                   c(1, 1, 1, 1, 1, 1, 1, 1),
                                   c(1, 1, 1, 1, 1, 1, 1, 1),
                                   c(1, 1, 1, 1, 1, 1, 1, 1),
                                   c(1, 1, 1, 1, 1, 1, 1, 1),
                                   c(1, 1, 1, 1, 1, 1, 1, 1),
                                   c(NA, NA, NA, NA, NA, NA, NA, NA),
                                   c(NA, 2, 2, NA, NA, 3, 3, NA),
                                   c(NA, 2, 2, NA, NA, 3, 3, NA),
                                   c(NA, 2, 2, NA, NA, 3, 3, NA),
                                   c(NA, 2, 2, NA, NA, 3, 3, NA),
                                   c(NA, 2, 2, NA, NA, 3, 3, NA),
                                   c(NA, 2, 2, NA, NA, 3, 3, NA),
                                   c(NA, 2, 2, NA, NA, 3, 3, NA),
                                   c(NA, 2, 2, NA, NA, 3, 3, NA),
                                   c(NA, 2, 2, NA, NA, 3, 3, NA),
                                   c(NA, 2, 2, NA, NA, 3, 3, NA),
                                   c(NA, 2, 2, NA, NA, 3, 3, NA)
             ))