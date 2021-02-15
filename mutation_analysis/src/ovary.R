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
                      "cell lines out of ",length(ovary$UBA2_gene_effect)," with UBA2 score <" ,
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

#subset of ovary cell lines with below threshold gene effect score
ovary_valid_mut <- ovary_mutation[which(ovary_mutation$UBA2_gene_effect <= ovary_average), ]
ovary_valid_mut <- droplevels(ovary_valid_mut)
#subset of cell lines with above threshold
ovary_invalid_mut <- ovary_mutation[which(ovary_mutation$UBA2_gene_effect > ovary_average), ]
ovary_invalid_mut <- droplevels(ovary_invalid_mut)

#eliminate mutations that also occur in invalid
#ovary_valid_mut1 <- ovary_valid_mut[!duplicated(ovary_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
ovary_mut_final <- anti_join(ovary_valid_mut, ovary_invalid_mut, 
                             by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

#calculate the mode (most frequent mutation)
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}
calculate_mode(ovary_mut_final$Hugo_Symbol)
#turns out to be COL6A6!


ovary_COL6A6 <- ovary_mut_final[which(ovary_mut_final$Hugo_Symbol == "COL6A6"), ]
ovary_no_COL6A6 <- ovary_mut_final[which(ovary_mut_final$Hugo_Symbol != "COL6A6"), ]

#find the most frequent variant type
calculate_mode(ovary_COL6A6$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with COL6A6_Missense
ovary_COL6A6_Mis <- ovary_COL6A6[which(ovary_COL6A6$Variant_Classification == "Missense_Mutation"), ]
ovary_COL6A6_No_Mis <- ovary_COL6A6[which(ovary_COL6A6$Variant_Classification != "Missense_Mutation"), ]
ovary_COL6A6 <- droplevels(ovary_COL6A6)
#7 cell lines with COL6A6 Missense mutation

#subset of all ccl
ovary_ccl <- distinct(ovary_mutation, ccl, .keep_all= TRUE )

y <- length(ovary_COL6A6_Mis$ID)
x <- length(ovary$ccl)
x <- x - y


final_table <- data.frame(mutation = c("COL6A6", "other"), counts = c(y, x))

bar_title <- paste("Mut: COL6A6 Missense Mutation \n", y, "overlaps out of", x, "total")

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

ovary_ccl[match(ovary_COL6A6$ID, ovary_ccl$ID), ] <- ovary_COL6A6

ovary_ccl$Hugo_Symbol <- as.character(ovary_ccl$Hugo_Symbol)

ovary_ccl$Hugo_Symbol[ovary_ccl$Hugo_Symbol == "COL6A6"] <- "has Mut"
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

