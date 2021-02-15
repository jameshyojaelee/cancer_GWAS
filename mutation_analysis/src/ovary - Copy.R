#ovary 

ovary <- UBA2_merged[which(UBA2_merged$lineage == "ovary"),]
ovary <- ovary[order(ovary$UBA2_gene_effect), ]
ovary_scores <-ovary$UBA2_gene_effect
ovary_average <- mean(ovary_scores)
ovary_sd <- sd(ovary_scores)
min(ovary_scores)
max(ovary_scores)

length(ovary$UBA2_gene_effect)
length(which(ovary$UBA2_gene_effect < ovary_average))

enrich_title <- paste("UBA2 gene effect in", ovary$lineage ,"cancer cells \n (",
                      length(which(ovary$UBA2_gene_effect < ovary_average)),
                      "cell liREV3L out of ",length(ovary$UBA2_gene_effect)," with UBA2 score <" ,
                      round(ovary_average, digits = 3), ")")

library(ggplot2)
library("scales")
enrich <- ggplot(data = ovary , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","white","white","white"),
                       guide = "colorbar", limits=c(min(ovary_scores),0)) +
  scale_x_continuous(limits=c(min(ovary_scores),0))+
  theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=11, angle = 90, vjust = 3, hjust= 1),
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())

#################################################################################################
########example with ovary cells
ovary_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "ovary"), ]

#calculate the number of unique ccl
length(unique(ovary_mutation$ccl))
#turns out to be 36 ccls

#subset of ovary cell liREV3L with below threshold gene effect score
ovary_valid_mut <- ovary_mutation[which(ovary_mutation$UBA2_gene_effect <= ovary_average), ]
ovary_valid_mut <- droplevels(ovary_valid_mut)
#subset of cell liREV3L with above threshold
ovary_invalid_mut <- ovary_mutation[which(ovary_mutation$UBA2_gene_effect > ovary_average), ]
ovary_invalid_mut <- droplevels(ovary_invalid_mut)

#eliminate mutations that also occur in invalid
#ovary_valid_mut1 <- ovary_valid_mut[!duplicated(ovary_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
ovary_mut_final <- anti_join(ovary_valid_mut, ovary_invalid_mut, 
                             by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

#calculate the second most frequent value
xx <- sort(table(ovary_mut_final$Hugo_Symbol),decreasing=TRUE)[1:10]
xx
#                                   

ovary_REV3L <- ovary_mut_final[which(ovary_mut_final$Hugo_Symbol == "REV3L"), ]
ovary_no_REV3L <- ovary_mut_final[which(ovary_mut_final$Hugo_Symbol != "REV3L"), ]

#find the most frequent variant type
calculate_mode(ovary_REV3L$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with REV3L_Missense
ovary_REV3L_Mis <- ovary_REV3L[which(ovary_REV3L$Variant_Classification == "Missense_Mutation"), ]
ovary_REV3L_No_Mis <- ovary_REV3L[which(ovary_REV3L$Variant_Classification != "Missense_Mutation"), ]
ovary_REV3L <- droplevels(ovary_REV3L)
#7 cell liREV3L with REV3L Missense_Mutation

#subset of all ccl
ovary_ccl <- distinct(ovary_mutation, ccl, .keep_all= TRUE )

y <- length(ovary_REV3L_Mis$ID)
x <- length(ovary$ccl)
x <- x - y


final_table <- data.frame(mutation = c("REV3L", "other"), counts = c(y, x))

bar_title <- paste("Mut: REV3L Missense_Mutation \n", y, "overlaps out of", x, "total")

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
ovary_ccl <- distinct(ovary_mutation, ccl, .keep_all= TRUE )

ovary_ccl[match(ovary_REV3L$ID, ovary_ccl$ID), ] <- ovary_REV3L

ovary_ccl$Hugo_Symbol <- as.character(ovary_ccl$Hugo_Symbol)

ovary_ccl$Hugo_Symbol[ovary_ccl$Hugo_Symbol == "REV3L"] <- "has Mut"
ovary_ccl$Hugo_Symbol[ovary_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =ovary_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(ovary_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
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

