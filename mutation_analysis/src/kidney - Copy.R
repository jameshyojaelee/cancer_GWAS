#kidney 

kidney <- UBA2_merged[which(UBA2_merged$lineage == "kidney"),]
kidney <- kidney[order(kidney$UBA2_gene_effect), ]
kidney_scores <-kidney$UBA2_gene_effect
kidney_average <- mean(kidney_scores)
kidney_sd <- sd(kidney_scores)
min(kidney_scores)
max(kidney_scores)

length(kidney$UBA2_gene_effect)
length(which(kidney$UBA2_gene_effect < kidney_average))

enrich_title <- paste("UBA2 gene effect in", kidney$lineage ,"cancer cells \n (",
                      length(which(kidney$UBA2_gene_effect < kidney_average)),
                      "cell lines out of ",length(kidney$UBA2_gene_effect)," with UBA2 score <" ,
                      round(kidney_average, digits = 3), ")")

library(ggplot2)
library("scales")
enrich <- ggplot(data = kidney , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","white","white","white"),
                       guide = "colorbar", limits=c(min(kidney_scores),0)) +
  scale_x_continuous(limits=c(min(kidney_scores),0))+
  theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=11, angle = 90, vjust = 3, hjust= 1),
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())

#################################################################################################
########example with kidney cells
kidney_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "kidney"), ]

#calculate the number of unique ccl
length(unique(kidney_mutation$ccl))
#turns out to be 36 ccls

#subset of kidney cell lines with below threshold gene effect score
kidney_valid_mut <- kidney_mutation[which(kidney_mutation$UBA2_gene_effect <= kidney_average), ]
kidney_valid_mut <- droplevels(kidney_valid_mut)
#subset of cell lines with above threshold
kidney_invalid_mut <- kidney_mutation[which(kidney_mutation$UBA2_gene_effect > kidney_average), ]
kidney_invalid_mut <- droplevels(kidney_invalid_mut)

#eliminate mutations that also occur in invalid
#kidney_valid_mut1 <- kidney_valid_mut[!duplicated(kidney_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
kidney_mut_final <- anti_join(kidney_valid_mut, kidney_invalid_mut, 
                            by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

#calculate the second most frequent value
xx <- sort(table(kidney_mut_final$Hugo_Symbol),decreasing=TRUE)[1:10]
xx
#                                 

kidney_GPR174 <- kidney_mut_final[which(kidney_mut_final$Hugo_Symbol == "GPR174"), ]
kidney_no_GPR174 <- kidney_mut_final[which(kidney_mut_final$Hugo_Symbol != "GPR174"), ]

#find the most frequent variant type
calculate_mode(kidney_GPR174$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with GPR174_Missense
kidney_GPR174_Mis <- kidney_GPR174[which(kidney_GPR174$Variant_Classification == "Missense_Mutation"), ]
kidney_GPR174_No_Mis <- kidney_GPR174[which(kidney_GPR174$Variant_Classification != "Missense_Mutation"), ]
kidney_GPR174 <- droplevels(kidney_GPR174)
#7 cell lines with GPR174 Missense_Mutation

#subset of all ccl
kidney_ccl <- distinct(kidney_mutation, ccl, .keep_all= TRUE )

y <- length(kidney_GPR174_Mis$ID)
x <- length(kidney$ccl)
x <- x - y


final_table <- data.frame(mutation = c("GPR174", "other"), counts = c(y, x))

bar_title <- paste("Mut: GPR174 Missense_Mutation \n", y, "overlaps out of", x, "total")

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
kidney_ccl <- distinct(kidney_mutation, ccl, .keep_all= TRUE )

kidney_ccl[match(kidney_GPR174$ID, kidney_ccl$ID), ] <- kidney_GPR174

kidney_ccl$Hugo_Symbol <- as.character(kidney_ccl$Hugo_Symbol)

kidney_ccl$Hugo_Symbol[kidney_ccl$Hugo_Symbol == "GPR174"] <- "has Mut"
kidney_ccl$Hugo_Symbol[kidney_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =kidney_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(kidney_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
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