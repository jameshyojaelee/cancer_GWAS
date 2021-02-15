#breast 

breast <- UBA2_merged[which(UBA2_merged$lineage == "breast"),]
breast <- breast[order(breast$UBA2_gene_effect), ]
breast_scores <-breast$UBA2_gene_effect
breast_average <- mean(breast_scores)
breast_sd <- sd(breast_scores)
min(breast_scores)
max(breast_scores)

length(breast$UBA2_gene_effect)
length(which(breast$UBA2_gene_effect < breast_average))

enrich_title <- paste("UBA2 gene effect in", breast$lineage ,"cancer cells \n (",
                      length(which(breast$UBA2_gene_effect < breast_average)),
                      "cell lines out of ",length(breast$UBA2_gene_effect)," with UBA2 score <" ,
                      round(breast_average, digits = 3), ")")

library(ggplot2)
library("scales")
enrich <- ggplot(data = breast , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","white","white","white"),
                       guide = "colorbar", limits=c(min(breast_scores),0)) +
  scale_x_continuous(limits=c(min(breast_scores),0))+
  theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=11, angle = 90, vjust = 3, hjust= 1),
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())

#################################################################################################
########example with breast cells
breast_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "breast"), ]

#calculate the number of unique ccl
length(unique(breast_mutation$ccl))
#turns out to be 36 ccls

#subset of breast cell lines with below threshold gene effect score
breast_valid_mut <- breast_mutation[which(breast_mutation$UBA2_gene_effect <= breast_average), ]
breast_valid_mut <- droplevels(breast_valid_mut)
#subset of cell lines with above threshold
breast_invalid_mut <- breast_mutation[which(breast_mutation$UBA2_gene_effect > breast_average), ]
breast_invalid_mut <- droplevels(breast_invalid_mut)

#eliminate mutations that also occur in invalid
#breast_valid_mut1 <- breast_valid_mut[!duplicated(breast_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
breast_mut_final <- anti_join(breast_valid_mut, breast_invalid_mut, 
                            by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

#calculate the mode (most frequent mutation)
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}
calculate_mode(breast_mut_final$Hugo_Symbol)
#turns out to be BAZ2A!


breast_BAZ2A <- breast_mut_final[which(breast_mut_final$Hugo_Symbol == "BAZ2A"), ]
breast_no_BAZ2A <- breast_mut_final[which(breast_mut_final$Hugo_Symbol != "BAZ2A"), ]

#find the most frequent variant type
calculate_mode(breast_BAZ2A$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with BAZ2A_Missense
breast_BAZ2A_Mis <- breast_BAZ2A[which(breast_BAZ2A$Variant_Classification == "Missense_Mutation"), ]
breast_BAZ2A_No_Mis <- breast_BAZ2A[which(breast_BAZ2A$Variant_Classification != "Missense_Mutation"), ]
breast_BAZ2A <- droplevels(breast_BAZ2A)
#7 cell lines with BAZ2A Missense mutation

#subset of all ccl
breast_ccl <- distinct(breast_mutation, ccl, .keep_all= TRUE )

y <- length(breast_BAZ2A_Mis$ID)
x <- length(breast$ccl)
x <- x - y


final_table <- data.frame(mutation = c("BAZ2A", "other"), counts = c(y, x))

bar_title <- paste("Mut: BAZ2A Missense Mutation \n", y, "overlaps out of", x, "total")

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
breast_ccl <- distinct(breast_mutation, ccl, .keep_all= TRUE )

breast_ccl[match(breast_BAZ2A$ID, breast_ccl$ID), ] <- breast_BAZ2A

breast_ccl$Hugo_Symbol <- as.character(breast_ccl$Hugo_Symbol)

breast_ccl$Hugo_Symbol[breast_ccl$Hugo_Symbol == "BAZ2A"] <- "has Mut"
breast_ccl$Hugo_Symbol[breast_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =breast_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(breast_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
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