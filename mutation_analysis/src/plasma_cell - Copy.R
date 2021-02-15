#plasma_cell 

plasma_cell <- UBA2_merged[which(UBA2_merged$lineage == "plasma_cell"),]
plasma_cell <- plasma_cell[order(plasma_cell$UBA2_gene_effect), ]
plasma_cell_scores <-plasma_cell$UBA2_gene_effect
plasma_cell_average <- mean(plasma_cell_scores)
plasma_cell_sd <- sd(plasma_cell_scores)
min(plasma_cell_scores)
max(plasma_cell_scores)

length(plasma_cell$UBA2_gene_effect)
length(which(plasma_cell$UBA2_gene_effect < plasma_cell_average))

enrich_title <- paste("UBA2 gene effect in", plasma_cell$lineage ,"cancer cells \n (",
                      length(which(plasma_cell$UBA2_gene_effect < plasma_cell_average)),
                      "cell lines out of ",length(plasma_cell$UBA2_gene_effect)," with UBA2 score <" ,
                      round(plasma_cell_average, digits = 3), ")")

library(ggplot2)
library("scales")
enrich <- ggplot(data = plasma_cell , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","white","white","white"),
                       guide = "colorbar", limits=c(min(plasma_cell_scores),0)) +
  scale_x_continuous(limits=c(min(plasma_cell_scores),0))+
  theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=11, angle = 90, vjust = 3, hjust= 1),
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())

#################################################################################################
########example with plasma_cell cells
plasma_cell_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "plasma_cell"), ]

#calculate the number of unique ccl
length(unique(plasma_cell_mutation$ccl))
#turns out to be 36 ccls

#subset of plasma_cell cell lines with below threshold gene effect score
plasma_cell_valid_mut <- plasma_cell_mutation[which(plasma_cell_mutation$UBA2_gene_effect <= plasma_cell_average), ]
plasma_cell_valid_mut <- droplevels(plasma_cell_valid_mut)
#subset of cell lines with above threshold
plasma_cell_invalid_mut <- plasma_cell_mutation[which(plasma_cell_mutation$UBA2_gene_effect > plasma_cell_average), ]
plasma_cell_invalid_mut <- droplevels(plasma_cell_invalid_mut)

#eliminate mutations that also occur in invalid
#plasma_cell_valid_mut1 <- plasma_cell_valid_mut[!duplicated(plasma_cell_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
plasma_cell_mut_final <- anti_join(plasma_cell_valid_mut, plasma_cell_invalid_mut, 
                             by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

xx <- sort(table(plasma_cell_mut_final$Hugo_Symbol),decreasing=TRUE)[1:10]
xx
#                    

plasma_cell_FBXL20 <- plasma_cell_mut_final[which(plasma_cell_mut_final$Hugo_Symbol == "FBXL20"), ]
plasma_cell_no_FBXL20 <- plasma_cell_mut_final[which(plasma_cell_mut_final$Hugo_Symbol != "FBXL20"), ]

#find the most frequent variant type
calculate_mode(plasma_cell_FBXL20$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with FBXL20_Missense
plasma_cell_FBXL20_Mis <- plasma_cell_FBXL20[which(plasma_cell_FBXL20$Variant_Classification == "Missense_Mutation"), ]
plasma_cell_FBXL20_No_Mis <- plasma_cell_FBXL20[which(plasma_cell_FBXL20$Variant_Classification != "Missense_Mutation"), ]
plasma_cell_FBXL20 <- droplevels(plasma_cell_FBXL20)
#7 cell lines with FBXL20 Missense mutation

#subset of all ccl
plasma_cell_ccl <- distinct(plasma_cell_mutation, ccl, .keep_all= TRUE )

y <- length(plasma_cell_FBXL20_Mis$ID)
x <- length(plasma_cell$ccl)
x <- x - y


final_table <- data.frame(mutation = c("FBXL20", "other"), counts = c(y, x))

bar_title <- paste("Mut: FBXL20 Missense Mutation \n", y, "overlaps out of", x, "total")

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
plasma_cell_ccl <- distinct(plasma_cell_mutation, ccl, .keep_all= TRUE )

plasma_cell_ccl[match(plasma_cell_FBXL20$ID, plasma_cell_ccl$ID), ] <- plasma_cell_FBXL20

plasma_cell_ccl$Hugo_Symbol <- as.character(plasma_cell_ccl$Hugo_Symbol)

plasma_cell_ccl$Hugo_Symbol[plasma_cell_ccl$Hugo_Symbol == "FBXL20"] <- "has Mut"
plasma_cell_ccl$Hugo_Symbol[plasma_cell_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =plasma_cell_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(plasma_cell_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
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
