#lung 

lung <- UBA2_merged[which(UBA2_merged$lineage == "lung"),]
lung <- lung[order(lung$UBA2_gene_effect), ]
lung_scores <-lung$UBA2_gene_effect
lung_average <- mean(lung_scores)
lung_sd <- sd(lung_scores)
min(lung_scores)
max(lung_scores)

length(lung$UBA2_gene_effect)
length(which(lung$UBA2_gene_effect < lung_average))

enrich_title <- paste("UBA2 gene effect in", lung$lineage ,"cancer cells \n (",
                      length(which(lung$UBA2_gene_effect < lung_average)),
                      "cell lines out of ",length(lung$UBA2_gene_effect)," with UBA2 score <" ,
                      round(lung_average, digits = 3), ")")

library(ggplot2)
library("scales")
enrich <- ggplot(data = lung , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","white","white","white"),
                       guide = "colorbar", limits=c(min(lung_scores),0)) +
  scale_x_continuous(limits=c(min(lung_scores),0))+
  theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=11, angle = 90, vjust = 3, hjust= 1),
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())

#################################################################################################
########example with lung cells
lung_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "lung"), ]

#calculate the number of unique ccl
length(unique(lung_mutation$ccl))
#turns out to be 36 ccls

#subset of lung cell lines with below threshold gene effect score
lung_valid_mut <- lung_mutation[which(lung_mutation$UBA2_gene_effect <= lung_average), ]
lung_valid_mut <- droplevels(lung_valid_mut)
#subset of cell lines with above threshold
lung_invalid_mut <- lung_mutation[which(lung_mutation$UBA2_gene_effect > lung_average), ]
lung_invalid_mut <- droplevels(lung_invalid_mut)

#eliminate mutations that also occur in invalid
#lung_valid_mut1 <- lung_valid_mut[!duplicated(lung_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
lung_mut_final <- anti_join(lung_valid_mut, lung_invalid_mut, 
                            by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

#calculate the second most frequent value
xx <- sort(table(lung_mut_final$Hugo_Symbol),decreasing=TRUE)[1:10]
xx
#                                                          

lung_PHKB <- lung_mut_final[which(lung_mut_final$Hugo_Symbol == "PHKB"), ]
lung_no_PHKB <- lung_mut_final[which(lung_mut_final$Hugo_Symbol != "PHKB"), ]

#find the most frequent variant type
calculate_mode(lung_PHKB$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with PHKB_Missense
lung_PHKB_Mis <- lung_PHKB[which(lung_PHKB$Variant_Classification == "Missense_Mutation"), ]
lung_PHKB_No_Mis <- lung_PHKB[which(lung_PHKB$Variant_Classification != "Missense_Mutation"), ]
lung_PHKB <- droplevels(lung_PHKB)
#7 cell lines with PHKB Missense_Mutation

#subset of all ccl
lung_ccl <- distinct(lung_mutation, ccl, .keep_all= TRUE )

y <- length(lung_PHKB_Mis$ID)
x <- length(lung$ccl)
x <- x - y


final_table <- data.frame(mutation = c("PHKB", "other"), counts = c(y, x))

bar_title <- paste("Mut: PHKB Missense_Mutation \n", y, "overlaps out of", x, "total")

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
lung_ccl <- distinct(lung_mutation, ccl, .keep_all= TRUE )

lung_ccl[match(lung_PHKB$ID, lung_ccl$ID), ] <- lung_PHKB

lung_ccl$Hugo_Symbol <- as.character(lung_ccl$Hugo_Symbol)

lung_ccl$Hugo_Symbol[lung_ccl$Hugo_Symbol == "PHKB"] <- "has Mut"
lung_ccl$Hugo_Symbol[lung_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =lung_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(lung_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
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