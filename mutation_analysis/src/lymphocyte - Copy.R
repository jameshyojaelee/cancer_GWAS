#lymphocyte 

lymphocyte <- UBA2_merged[which(UBA2_merged$lineage == "lymphocyte"),]
lymphocyte <- lymphocyte[order(lymphocyte$UBA2_gene_effect), ]
lymphocyte_scores <-lymphocyte$UBA2_gene_effect
lymphocyte_average <- mean(lymphocyte_scores)
lymphocyte_sd <- sd(lymphocyte_scores)
min(lymphocyte_scores)
max(lymphocyte_scores)

length(lymphocyte$UBA2_gene_effect)
length(which(lymphocyte$UBA2_gene_effect < lymphocyte_average))

enrich_title <- paste("UBA2 gene effect in", lymphocyte$lineage ,"cancer cells \n (",
                      length(which(lymphocyte$UBA2_gene_effect < lymphocyte_average)),
                      "cell lines out of ",length(lymphocyte$UBA2_gene_effect)," with UBA2 score <" ,
                      round(lymphocyte_average, digits = 3), ")")

library(ggplot2)
library("scales")
enrich <- ggplot(data = lymphocyte , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","white","white","white"),
                       guide = "colorbar", limits=c(min(lymphocyte_scores),0)) +
  scale_x_continuous(limits=c(min(lymphocyte_scores),0))+
  theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=11, angle = 90, vjust = 3, hjust= 1),
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())

#################################################################################################
########example with lymphocyte cells
lymphocyte_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "lymphocyte"), ]

#calculate the number of unique ccl
length(unique(lymphocyte_mutation$ccl))
#turns out to be 36 ccls

#subset of lymphocyte cell lines with below threshold gene effect score
lymphocyte_valid_mut <- lymphocyte_mutation[which(lymphocyte_mutation$UBA2_gene_effect <= lymphocyte_average), ]
lymphocyte_valid_mut <- droplevels(lymphocyte_valid_mut)
#subset of cell lines with above threshold
lymphocyte_invalid_mut <- lymphocyte_mutation[which(lymphocyte_mutation$UBA2_gene_effect > lymphocyte_average), ]
lymphocyte_invalid_mut <- droplevels(lymphocyte_invalid_mut)

#eliminate mutations that also occur in invalid
#lymphocyte_valid_mut1 <- lymphocyte_valid_mut[!duplicated(lymphocyte_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
lymphocyte_mut_final <- anti_join(lymphocyte_valid_mut, lymphocyte_invalid_mut, 
                             by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

#calculate the second most frequent value
xx <- sort(table(lymphocyte_mut_final$Hugo_Symbol),decreasing=TRUE)[1:10]
xx
#                        

lymphocyte_PCDHA7 <- lymphocyte_mut_final[which(lymphocyte_mut_final$Hugo_Symbol == "PCDHA7"), ]
lymphocyte_no_PCDHA7 <- lymphocyte_mut_final[which(lymphocyte_mut_final$Hugo_Symbol != "PCDHA7"), ]

#find the most frequent variant type
calculate_mode(lymphocyte_PCDHA7$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with PCDHA7_Missense
lymphocyte_PCDHA7_Mis <- lymphocyte_PCDHA7[which(lymphocyte_PCDHA7$Variant_Classification == "Missense_Mutation"), ]
lymphocyte_PCDHA7_No_Mis <- lymphocyte_PCDHA7[which(lymphocyte_PCDHA7$Variant_Classification != "Missense_Mutation"), ]
lymphocyte_PCDHA7 <- droplevels(lymphocyte_PCDHA7)
#7 cell lines with PCDHA7 Missense mutation

#subset of all ccl
lymphocyte_ccl <- distinct(lymphocyte_mutation, ccl, .keep_all= TRUE )

y <- length(lymphocyte_PCDHA7_Mis$ID)
x <- length(lymphocyte$ccl)
x <- x - y


final_table <- data.frame(mutation = c("PCDHA7", "other"), counts = c(y, x))

bar_title <- paste("Mut: PCDHA7 Missense Mutation \n", y, "overlaps out of", x, "total")

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
lymphocyte_ccl <- distinct(lymphocyte_mutation, ccl, .keep_all= TRUE )

lymphocyte_ccl[match(lymphocyte_PCDHA7$ID, lymphocyte_ccl$ID), ] <- lymphocyte_PCDHA7

lymphocyte_ccl$Hugo_Symbol <- as.character(lymphocyte_ccl$Hugo_Symbol)

lymphocyte_ccl$Hugo_Symbol[lymphocyte_ccl$Hugo_Symbol == "PCDHA7"] <- "has Mut"
lymphocyte_ccl$Hugo_Symbol[lymphocyte_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =lymphocyte_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(lymphocyte_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
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
