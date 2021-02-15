#blood 

blood <- UBA2_merged[which(UBA2_merged$lineage == "blood"),]
blood <- blood[order(blood$UBA2_gene_effect), ]
blood_scores <-blood$UBA2_gene_effect
blood_average <- mean(blood_scores)
blood_sd <- sd(blood_scores)
min(blood_scores)
max(blood_scores)

length(blood$UBA2_gene_effect)
length(which(blood$UBA2_gene_effect < blood_average))

enrich_title <- paste("UBA2 gene effect in", blood$lineage ,"cancer cells \n (",
                      length(which(blood$UBA2_gene_effect < blood_average)),
                      "cell lines out of ",length(blood$UBA2_gene_effect)," with UBA2 score <" ,
                      round(blood_average, digits = 3), ")")

library(ggplot2)
library("scales")
enrich <- ggplot(data = blood , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","red", "white","white","white","white"),
                       guide = "colorbar", limits=c(min(blood_scores),0)) +
  scale_x_continuous(limits=c(min(blood_scores),0))+
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
########example with blood cells
blood_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "blood"), ]

#calculate the number of unique ccl
length(unique(blood_mutation$ccl))
#turns out to be 36 ccls

#subset of blood cell lines with below threshold gene effect score
blood_valid_mut <- blood_mutation[which(blood_mutation$UBA2_gene_effect <= blood_average), ]
blood_valid_mut <- droplevels(blood_valid_mut)
#subset of cell lines with above threshold
blood_invalid_mut <- blood_mutation[which(blood_mutation$UBA2_gene_effect > blood_average), ]
blood_invalid_mut <- droplevels(blood_invalid_mut)

#eliminate mutations that also occur in invalid
#blood_valid_mut1 <- blood_valid_mut[!duplicated(blood_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
blood_mut_final <- anti_join(blood_valid_mut, blood_invalid_mut, 
                                  by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

#most frequent mutations
xx <- sort(table(blood_mut_final$Hugo_Symbol),decreasing=TRUE)[1:10]
xx
#             

blood_ATHL1 <- blood_mut_final[which(blood_mut_final$Hugo_Symbol == "ATHL1"), ]
blood_no_ATHL1 <- blood_mut_final[which(blood_mut_final$Hugo_Symbol != "ATHL1"), ]

#find the most frequent variant type
calculate_mode(blood_ATHL1$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with ATHL1_Missense
blood_ATHL1_Mis <- blood_ATHL1[which(blood_ATHL1$Variant_Classification == "Missense_Mutation"), ]
blood_ATHL1_No_Mis <- blood_ATHL1[which(blood_ATHL1$Variant_Classification != "Missense_Mutation"), ]
blood_ATHL1 <- droplevels(blood_ATHL1)
#7 cell lines with ATHL1 Missense mutation

#subset of all ccl
blood_ccl <- distinct(blood_mutation, ccl, .keep_all= TRUE )

y <- length(blood_ATHL1_Mis$ID)
x <- length(blood$ccl)
x <- x - y


final_table <- data.frame(mutation = c("ATHL1", "other"), counts = c(y, x))

bar_title <- paste("Mut: ATHL1 Missense Mutation \n", y, "overlaps out of", x, "total")

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
blood_ccl <- distinct(blood_mutation, ccl, .keep_all= TRUE )

blood_ccl[match(blood_ATHL1$ID, blood_ccl$ID), ] <- blood_ATHL1

blood_ccl$Hugo_Symbol <- as.character(blood_ccl$Hugo_Symbol)

blood_ccl$Hugo_Symbol[blood_ccl$Hugo_Symbol == "ATHL1"] <- "has Mut"
blood_ccl$Hugo_Symbol[blood_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =blood_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(blood_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
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

