#colorectal 

colorectal <- UBA2_merged[which(UBA2_merged$lineage == "colorectal"),]
colorectal <- colorectal[order(colorectal$UBA2_gene_effect), ]
colorectal_scores <-colorectal$UBA2_gene_effect
colorectal_average <- mean(colorectal_scores)
colorectal_sd <- sd(colorectal_scores)
min(colorectal_scores)
max(colorectal_scores)
length(colorectal$UBA2_gene_effect)
length(which(colorectal$UBA2_gene_effect < colorectal_average))

enrich_title <- paste("UBA2 gene effect in", colorectal$lineage ,"cancer cells \n (",
                      length(which(colorectal$UBA2_gene_effect < colorectal_average)),
                      "cell lines out of ",length(colorectal$UBA2_gene_effect)," with UBA2 score <" ,
                      round(colorectal_average, digits = 3), ")")

library(ggplot2)
library("scales")
enrich <- ggplot(data = colorectal , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","white","white","white"),
                       guide = "colorbar", limits=c(min(colorectal_scores),0)) +
  scale_x_continuous(limits=c(min(colorectal_scores),0))+
  theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=11, angle = 90, vjust = 3, hjust= 1),
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())

#################################################################################################
########example with colorectal cells
colorectal_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "colorectal"), ]

#calculate the number of unique ccl
length(unique(colorectal_mutation$ccl))
#turns out to be 36 ccls

#subset of colorectal cell lines with below threshold gene effect score
colorectal_valid_mut <- colorectal_mutation[which(colorectal_mutation$UBA2_gene_effect <= colorectal_average), ]
colorectal_valid_mut <- droplevels(colorectal_valid_mut)
#subset of cell lines with above threshold
colorectal_invalid_mut <- colorectal_mutation[which(colorectal_mutation$UBA2_gene_effect > colorectal_average), ]
colorectal_invalid_mut <- droplevels(colorectal_invalid_mut)

#eliminate mutations that also occur in invalid
#colorectal_valid_mut1 <- colorectal_valid_mut[!duplicated(colorectal_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
colorectal_mut_final <- anti_join(colorectal_valid_mut, colorectal_invalid_mut, 
                              by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)
pan_unique <- unique(colorectal_mut_final$ccl)

#calculate the second most frequent value
xx <- sort(table(colorectal_mut_final$Hugo_Symbol),decreasing=TRUE)[1:300]
xx
#               

###########################################################export###########################################################
##### eliminate mutations that occur less than 6
colorectal_mut_final <- colorectal_mut_final[colorectal_mut_final$Hugo_Symbol %in% names(which(table(colorectal_mut_final$Hugo_Symbol) > 5)), ]
length(unique(colorectal_mut_final$ccl)) # 18 cell lines 
length(unique(colorectal_mut_final$Hugo_Symbol)) # 18 cell lines 
library(xlsx)
write.xlsx(colorectal_mut_final, "colorectal_mut.xlsx")
###########################################################export###########################################################


colorectal_CASK <- colorectal_mut_final[which(colorectal_mut_final$Hugo_Symbol == "CASK"), ]
colorectal_no_CASK <- colorectal_mut_final[which(colorectal_mut_final$Hugo_Symbol != "CASK"), ]

#find the most frequent variant type
calculate_mode(colorectal_CASK$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with CASK_Missense
colorectal_CASK_Mis <- colorectal_CASK[which(colorectal_CASK$Variant_Classification == "Missense_Mutation"), ]
colorectal_CASK_No_Mis <- colorectal_CASK[which(colorectal_CASK$Variant_Classification != "Missense_Mutation"), ]
colorectal_CASK <- droplevels(colorectal_CASK)
#7 cell lines with CASK Missense mutation

#subset of all ccl
colorectal_ccl <- distinct(colorectal_mutation, ccl, .keep_all= TRUE )

y <- length(colorectal_CASK_Mis$ID)
x <- length(colorectal$ccl)
x <- x - y


final_table <- data.frame(mutation = c("CASK", "other"), counts = c(y, x))

bar_title <- paste("Mut: CASK Missense Mutation \n", y, "overlaps out of", x, "total")

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
colorectal_ccl <- distinct(colorectal_mutation, ccl, .keep_all= TRUE )

colorectal_ccl[match(colorectal_CASK$ID, colorectal_ccl$ID), ] <- colorectal_CASK

colorectal_ccl$Hugo_Symbol <- as.character(colorectal_ccl$Hugo_Symbol)

colorectal_ccl$Hugo_Symbol[colorectal_ccl$Hugo_Symbol == "CASK"] <- "has Mut"
colorectal_ccl$Hugo_Symbol[colorectal_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =colorectal_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(colorectal_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
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
