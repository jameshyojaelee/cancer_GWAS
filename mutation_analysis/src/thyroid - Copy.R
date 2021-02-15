#thyroid 

thyroid <- UBA2_merged[which(UBA2_merged$lineage == "thyroid"),]
thyroid <- thyroid[order(thyroid$UBA2_gene_effect), ]
thyroid_scores <-thyroid$UBA2_gene_effect
thyroid_average <- mean(thyroid_scores)
thyroid_sd <- sd(thyroid_scores)
min(thyroid_scores)
max(thyroid_scores)

length(thyroid$UBA2_gene_effect)
length(which(thyroid$UBA2_gene_effect < thyroid_average))

enrich_title <- paste("UBA2 gene effect in", thyroid$lineage ,"cancer cells \n (",
                      length(which(thyroid$UBA2_gene_effect < thyroid_average)),
                      "cell lines out of ",length(thyroid$UBA2_gene_effect)," with UBA2 score <" ,
                      round(thyroid_average, digits = 3), ")")

library(ggplot2)
library("scales")
enrich <- ggplot(data = thyroid , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","white","white","white"),
                       guide = "colorbar", limits=c(min(thyroid_scores),0)) +
  scale_x_continuous(limits=c(min(thyroid_scores),0))+
  theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=11, angle = 90, vjust = 3, hjust= 1),
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())

#################################################################################################
########example with thyroid cells
thyroid_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "thyroid"), ]

#calculate the number of unique ccl
length(unique(thyroid_mutation$ccl))
#turns out to be 36 ccls

#subset of thyroid cell lines with below threshold gene effect score
thyroid_valid_mut <- thyroid_mutation[which(thyroid_mutation$UBA2_gene_effect <= thyroid_average), ]
thyroid_valid_mut <- droplevels(thyroid_valid_mut)
#subset of cell lines with above threshold
thyroid_invalid_mut <- thyroid_mutation[which(thyroid_mutation$UBA2_gene_effect > thyroid_average), ]
thyroid_invalid_mut <- droplevels(thyroid_invalid_mut)

#eliminate mutations that also occur in invalid
#thyroid_valid_mut1 <- thyroid_valid_mut[!duplicated(thyroid_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
thyroid_mut_final <- anti_join(thyroid_valid_mut, thyroid_invalid_mut, 
                             by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

#calculate the second most frequent value
xx <- sort(table(thyroid_mut_final$Hugo_Symbol),decreasing=TRUE)[1:10]
xx
#                     

thyroid_DYNC2H1 <- thyroid_mut_final[which(thyroid_mut_final$Hugo_Symbol == "DYNC2H1"), ]
thyroid_no_DYNC2H1 <- thyroid_mut_final[which(thyroid_mut_final$Hugo_Symbol != "DYNC2H1"), ]

#find the most frequent variant type
calculate_mode(thyroid_DYNC2H1$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with DYNC2H1_Missense
thyroid_DYNC2H1_Mis <- thyroid_DYNC2H1[which(thyroid_DYNC2H1$Variant_Classification == "Missense_Mutation"), ]
thyroid_DYNC2H1_No_Mis <- thyroid_DYNC2H1[which(thyroid_DYNC2H1$Variant_Classification != "Missense_Mutation"), ]
thyroid_DYNC2H1 <- droplevels(thyroid_DYNC2H1)
#7 cell lines with DYNC2H1 Missense mutation

#subset of all ccl
thyroid_ccl <- distinct(thyroid_mutation, ccl, .keep_all= TRUE )

y <- length(thyroid_DYNC2H1_Mis$ID)
x <- length(thyroid$ccl)
x <- x - y


final_table <- data.frame(mutation = c("DYNC2H1", "other"), counts = c(y, x))

bar_title <- paste("Mut: DYNC2H1 Missense Mutation \n", y, "overlaps out of", x, "total")

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
thyroid_ccl <- distinct(thyroid_mutation, ccl, .keep_all= TRUE )

thyroid_ccl[match(thyroid_DYNC2H1$ID, thyroid_ccl$ID), ] <- thyroid_DYNC2H1

thyroid_ccl$Hugo_Symbol <- as.character(thyroid_ccl$Hugo_Symbol)

thyroid_ccl$Hugo_Symbol[thyroid_ccl$Hugo_Symbol == "DYNC2H1"] <- "has Mut"
thyroid_ccl$Hugo_Symbol[thyroid_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =thyroid_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(thyroid_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
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
