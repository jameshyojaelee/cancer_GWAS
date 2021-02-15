#skin 

skin <- UBA2_merged[which(UBA2_merged$lineage == "skin"),]
skin <- skin[order(skin$UBA2_gene_effect), ]
skin_scores <-skin$UBA2_gene_effect
skin_average <- mean(skin_scores)
skin_sd <- sd(skin_scores)
min(skin_scores)
max(skin_scores)

length(skin$UBA2_gene_effect)
length(which(skin$UBA2_gene_effect < skin_average))

enrich_title <- paste("UBA2 gene effect in", skin$lineage ,"cancer cells \n (",
                      length(which(skin$UBA2_gene_effect < skin_average)),
                      "cell lines out of ",length(skin$UBA2_gene_effect)," with UBA2 score <" ,
                      round(skin_average, digits = 3), ")")

library(ggplot2)
library("scales")
enrich <- ggplot(data = skin , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","white","white","white"),
                       guide = "colorbar", limits=c(min(skin_scores),0)) +
  scale_x_continuous(limits=c(min(skin_scores),0))+
  theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=11, angle = 90, vjust = 3, hjust= 1),
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())

#################################################################################################
########example with skin cells
skin_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "skin"), ]

#calculate the number of unique ccl
length(unique(skin_mutation$ccl))
#turns out to be 36 ccls

#subset of skin cell lines with below threshold gene effect score
skin_valid_mut <- skin_mutation[which(skin_mutation$UBA2_gene_effect <= skin_average), ]
skin_valid_mut <- droplevels(skin_valid_mut)
#subset of cell lines with above threshold
skin_invalid_mut <- skin_mutation[which(skin_mutation$UBA2_gene_effect > skin_average), ]
skin_invalid_mut <- droplevels(skin_invalid_mut)

#eliminate mutations that also occur in invalid
#skin_valid_mut1 <- skin_valid_mut[!duplicated(skin_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
skin_mut_final <- anti_join(skin_valid_mut, skin_invalid_mut, 
                             by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

#calculate the second most frequent value
xx <- sort(table(skin_mut_final$Hugo_Symbol),decreasing=TRUE)[1:10]
xx
#                    


skin_ACACA <- skin_mut_final[which(skin_mut_final$Hugo_Symbol == "ACACA"), ]
skin_no_ACACA <- skin_mut_final[which(skin_mut_final$Hugo_Symbol != "ACACA"), ]

#find the most frequent variant type
calculate_mode(skin_ACACA$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with ACACA_Missense
skin_ACACA_Mis <- skin_ACACA[which(skin_ACACA$Variant_Classification == "Missense_Mutation"), ]
skin_ACACA_No_Mis <- skin_ACACA[which(skin_ACACA$Variant_Classification != "Missense_Mutation"), ]
skin_ACACA <- droplevels(skin_ACACA)
#7 cell lines with ACACA Missense mutation

#subset of all ccl
skin_ccl <- distinct(skin_mutation, ccl, .keep_all= TRUE )

y <- length(skin_ACACA_Mis$ID)
x <- length(skin$ccl)
x <- x - y


final_table <- data.frame(mutation = c("ACACA", "other"), counts = c(y, x))

bar_title <- paste("Mut: ACACA Missense Mutation \n", y, "overlaps out of", x, "total")

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
skin_ccl <- distinct(skin_mutation, ccl, .keep_all= TRUE )

skin_ccl[match(skin_ACACA$ID, skin_ccl$ID), ] <- skin_ACACA

skin_ccl$Hugo_Symbol <- as.character(skin_ccl$Hugo_Symbol)

skin_ccl$Hugo_Symbol[skin_ccl$Hugo_Symbol == "ACACA"] <- "has Mut"
skin_ccl$Hugo_Symbol[skin_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =skin_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(skin_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
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
