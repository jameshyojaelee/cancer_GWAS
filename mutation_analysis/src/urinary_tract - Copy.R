#urinary_tract 

urinary_tract <- UBA2_merged[which(UBA2_merged$lineage == "urinary_tract"),]
urinary_tract <- urinary_tract[order(urinary_tract$UBA2_gene_effect), ]
urinary_tract_scores <-urinary_tract$UBA2_gene_effect
urinary_tract_average <- mean(urinary_tract_scores)
urinary_tract_sd <- sd(urinary_tract_scores)
min(urinary_tract_scores)
max(urinary_tract_scores)

length(urinary_tract$UBA2_gene_effect)
length(which(urinary_tract$UBA2_gene_effect < urinary_tract_average))

enrich_title <- paste("UBA2 gene effect in", urinary_tract$lineage ,"cancer cells \n (",
                      length(which(urinary_tract$UBA2_gene_effect < urinary_tract_average)),
                      "cell lines out of ",length(urinary_tract$UBA2_gene_effect)," with UBA2 score <" ,
                      round(urinary_tract_average, digits = 3), ")")

library(ggplot2)
library("scales")
enrich <- ggplot(data = urinary_tract , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","white","white","white"),
                       guide = "colorbar", limits=c(min(urinary_tract_scores),0)) +
  scale_x_continuous(limits=c(min(urinary_tract_scores),0))+
  theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=11, angle = 90, vjust = 3, hjust= 1),
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())

#################################################################################################
########example with urinary_tract cells
urinary_tract_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "urinary_tract"), ]

#calculate the number of unique ccl
length(unique(urinary_tract_mutation$ccl))
#turns out to be 36 ccls

#subset of urinary_tract cell lines with below threshold gene effect score
urinary_tract_valid_mut <- urinary_tract_mutation[which(urinary_tract_mutation$UBA2_gene_effect <= urinary_tract_average), ]
urinary_tract_valid_mut <- droplevels(urinary_tract_valid_mut)
#subset of cell lines with above threshold
urinary_tract_invalid_mut <- urinary_tract_mutation[which(urinary_tract_mutation$UBA2_gene_effect > urinary_tract_average), ]
urinary_tract_invalid_mut <- droplevels(urinary_tract_invalid_mut)

#eliminate mutations that also occur in invalid
#urinary_tract_valid_mut1 <- urinary_tract_valid_mut[!duplicated(urinary_tract_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
urinary_tract_mut_final <- anti_join(urinary_tract_valid_mut, urinary_tract_invalid_mut, 
                             by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

#calculate the second most frequent value
xx <- sort(table(urinary_tract_mut_final$Hugo_Symbol),decreasing=TRUE)[1:10]
xx
#                         

urinary_tract_MYCBP2 <- urinary_tract_mut_final[which(urinary_tract_mut_final$Hugo_Symbol == "MYCBP2"), ]
urinary_tract_no_MYCBP2 <- urinary_tract_mut_final[which(urinary_tract_mut_final$Hugo_Symbol != "MYCBP2"), ]

#find the most frequent variant type
calculate_mode(urinary_tract_MYCBP2$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with MYCBP2_Missense
urinary_tract_MYCBP2_Mis <- urinary_tract_MYCBP2[which(urinary_tract_MYCBP2$Variant_Classification == "Missense_Mutation"), ]
urinary_tract_MYCBP2_No_Mis <- urinary_tract_MYCBP2[which(urinary_tract_MYCBP2$Variant_Classification != "Missense_Mutation"), ]
urinary_tract_MYCBP2 <- droplevels(urinary_tract_MYCBP2)
#7 cell lines with MYCBP2 Missense mutation

#subset of all ccl
urinary_tract_ccl <- distinct(urinary_tract_mutation, ccl, .keep_all= TRUE )

y <- length(urinary_tract_MYCBP2_Mis$ID)
x <- length(urinary_tract$ccl)
x <- x - y


final_table <- data.frame(mutation = c("MYCBP2", "other"), counts = c(y, x))

bar_title <- paste("Mut: MYCBP2 Missense Mutation \n", y, "overlaps out of", x, "total")

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
urinary_tract_ccl <- distinct(urinary_tract_mutation, ccl, .keep_all= TRUE )

urinary_tract_ccl[match(urinary_tract_MYCBP2$ID, urinary_tract_ccl$ID), ] <- urinary_tract_MYCBP2

urinary_tract_ccl$Hugo_Symbol <- as.character(urinary_tract_ccl$Hugo_Symbol)

urinary_tract_ccl$Hugo_Symbol[urinary_tract_ccl$Hugo_Symbol == "MYCBP2"] <- "has Mut"
urinary_tract_ccl$Hugo_Symbol[urinary_tract_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =urinary_tract_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(urinary_tract_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
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

