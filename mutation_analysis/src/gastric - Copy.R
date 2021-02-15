#gastric 

gastric <- UBA2_merged[which(UBA2_merged$lineage == "gastric"),]
gastric <- gastric[order(gastric$UBA2_gene_effect), ]
gastric_scores <-gastric$UBA2_gene_effect
gastric_average <- mean(gastric_scores)
gastric_sd <- sd(gastric_scores)
min(gastric_scores)
max(gastric_scores)

length(gastric$UBA2_gene_effect)
length(which(gastric$UBA2_gene_effect < gastric_average))

enrich_title <- paste("UBA2 gene effect in", gastric$lineage ,"cancer cells \n (",
                      length(which(gastric$UBA2_gene_effect < gastric_average)),
                      "cell lines out of ",length(gastric$UBA2_gene_effect)," with UBA2 score <" ,
                      round(gastric_average, digits = 3), ")")

library(ggplot2)
library("scales")
enrich <- ggplot(data = gastric , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","white","white","white"),
                       guide = "colorbar", limits=c(min(gastric_scores),0)) +
  scale_x_continuous(limits=c(min(gastric_scores),0))+
  theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=11, angle = 90, vjust = 3, hjust= 1),
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())

#################################################################################################
########example with gastric cells
gastric_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "gastric"), ]

#calculate the number of unique ccl
length(unique(gastric_mutation$ccl))
#turns out to be 36 ccls

#subset of gastric cell lines with below threshold gene effect score
gastric_valid_mut <- gastric_mutation[which(gastric_mutation$UBA2_gene_effect <= gastric_average), ]
gastric_valid_mut <- droplevels(gastric_valid_mut)
#subset of cell lines with above threshold
gastric_invalid_mut <- gastric_mutation[which(gastric_mutation$UBA2_gene_effect > gastric_average), ]
gastric_invalid_mut <- droplevels(gastric_invalid_mut)

#eliminate mutations that also occur in invalid
#gastric_valid_mut1 <- gastric_valid_mut[!duplicated(gastric_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
gastric_mut_final <- anti_join(gastric_valid_mut, gastric_invalid_mut, 
                            by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

#calculate the second most frequent value
xx <- sort(table(gastric_mut_final$Hugo_Symbol),decreasing=TRUE)[1:10]
xx
#                                  

gastric_MAGI3 <- gastric_mut_final[which(gastric_mut_final$Hugo_Symbol == "MAGI3"), ]
gastric_no_MAGI3 <- gastric_mut_final[which(gastric_mut_final$Hugo_Symbol != "MAGI3"), ]

#find the most frequent variant type
calculate_mode(gastric_MAGI3$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with MAGI3_Missense
gastric_MAGI3_Mis <- gastric_MAGI3[which(gastric_MAGI3$Variant_Classification == "Missense_Mutation"), ]
gastric_MAGI3_No_Mis <- gastric_MAGI3[which(gastric_MAGI3$Variant_Classification != "Missense_Mutation"), ]
gastric_MAGI3 <- droplevels(gastric_MAGI3)
#7 cell lines with MAGI3 Missense mutation

#subset of all ccl
gastric_ccl <- distinct(gastric_mutation, ccl, .keep_all= TRUE )

y <- length(gastric_MAGI3_Mis$ID)
x <- length(gastric$ccl)
x <- x - y


final_table <- data.frame(mutation = c("MAGI3", "other"), counts = c(y, x))

bar_title <- paste("Mut: MAGI3 Missense Mutation \n", y, "overlaps out of", x, "total")

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
gastric_ccl <- distinct(gastric_mutation, ccl, .keep_all= TRUE )

gastric_ccl[match(gastric_MAGI3$ID, gastric_ccl$ID), ] <- gastric_MAGI3

gastric_ccl$Hugo_Symbol <- as.character(gastric_ccl$Hugo_Symbol)

gastric_ccl$Hugo_Symbol[gastric_ccl$Hugo_Symbol == "MAGI3"] <- "has Mut"
gastric_ccl$Hugo_Symbol[gastric_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =gastric_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(gastric_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
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
