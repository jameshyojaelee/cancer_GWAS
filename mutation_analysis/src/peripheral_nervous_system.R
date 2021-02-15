#peripheral_nervous_system 

peripheral_nervous_system <- UBA2_merged[which(UBA2_merged$lineage == "peripheral_nervous_system"),]
peripheral_nervous_system <- peripheral_nervous_system[order(peripheral_nervous_system$UBA2_gene_effect), ]
peripheral_nervous_system_scores <-peripheral_nervous_system$UBA2_gene_effect
peripheral_nervous_system_average <- mean(peripheral_nervous_system_scores)
peripheral_nervous_system_sd <- sd(peripheral_nervous_system_scores)
min(peripheral_nervous_system_scores)
max(peripheral_nervous_system_scores)

length(peripheral_nervous_system$UBA2_gene_effect)
length(which(peripheral_nervous_system$UBA2_gene_effect < peripheral_nervous_system_average))

enrich_title <- paste("UBA2 gene effect in", peripheral_nervous_system$lineage ,"cancer cells \n (",
                      length(which(peripheral_nervous_system$UBA2_gene_effect < peripheral_nervous_system_average)),
                      "cell lines out of ",length(peripheral_nervous_system$UBA2_gene_effect)," with UBA2 score <" ,
                      round(peripheral_nervous_system_average, digits = 3), ")")

library(ggplot2)
library("scales")
enrich <- ggplot(data = peripheral_nervous_system , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","white","white","white"),
                       guide = "colorbar", limits=c(min(peripheral_nervous_system_scores),0)) +
  scale_x_continuous(limits=c(min(peripheral_nervous_system_scores),0))+
  theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=11, angle = 90, vjust = 3, hjust= 1),
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())

#################################################################################################
########example with peripheral_nervous_system cells
peripheral_nervous_system_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "peripheral_nervous_system"), ]

#calculate the number of unique ccl
length(unique(peripheral_nervous_system_mutation$ccl))
#turns out to be 36 ccls

#subset of peripheral_nervous_system cell lines with below threshold gene effect score
peripheral_nervous_system_valid_mut <- peripheral_nervous_system_mutation[which(peripheral_nervous_system_mutation$UBA2_gene_effect <= peripheral_nervous_system_average), ]
peripheral_nervous_system_valid_mut <- droplevels(peripheral_nervous_system_valid_mut)
#subset of cell lines with above threshold
peripheral_nervous_system_invalid_mut <- peripheral_nervous_system_mutation[which(peripheral_nervous_system_mutation$UBA2_gene_effect > peripheral_nervous_system_average), ]
peripheral_nervous_system_invalid_mut <- droplevels(peripheral_nervous_system_invalid_mut)

#eliminate mutations that also occur in invalid
#peripheral_nervous_system_valid_mut1 <- peripheral_nervous_system_valid_mut[!duplicated(peripheral_nervous_system_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
peripheral_nervous_system_mut_final <- anti_join(peripheral_nervous_system_valid_mut, peripheral_nervous_system_invalid_mut, 
                             by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

#calculate the mode (most frequent mutation)
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}
calculate_mode(peripheral_nervous_system_mut_final$Hugo_Symbol)
#turns out to be MYH13!


peripheral_nervous_system_MYH13 <- peripheral_nervous_system_mut_final[which(peripheral_nervous_system_mut_final$Hugo_Symbol == "MYH13"), ]
peripheral_nervous_system_no_MYH13 <- peripheral_nervous_system_mut_final[which(peripheral_nervous_system_mut_final$Hugo_Symbol != "MYH13"), ]

#find the most frequent variant type
calculate_mode(peripheral_nervous_system_MYH13$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with MYH13_Missense
peripheral_nervous_system_MYH13_Mis <- peripheral_nervous_system_MYH13[which(peripheral_nervous_system_MYH13$Variant_Classification == "Missense_Mutation"), ]
peripheral_nervous_system_MYH13_No_Mis <- peripheral_nervous_system_MYH13[which(peripheral_nervous_system_MYH13$Variant_Classification != "Missense_Mutation"), ]
peripheral_nervous_system_MYH13 <- droplevels(peripheral_nervous_system_MYH13)
#7 cell lines with MYH13 Missense mutation

#subset of all ccl
peripheral_nervous_system_ccl <- distinct(peripheral_nervous_system_mutation, ccl, .keep_all= TRUE )

y <- length(peripheral_nervous_system_MYH13_Mis$ID)
x <- length(peripheral_nervous_system$ccl)
x <- x - y


final_table <- data.frame(mutation = c("MYH13", "other"), counts = c(y, x))

bar_title <- paste("Mut: MYH13 Missense Mutation \n", y, "overlaps out of", x, "total")

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
peripheral_nervous_system_ccl <- distinct(peripheral_nervous_system_mutation, ccl, .keep_all= TRUE )

peripheral_nervous_system_ccl[match(peripheral_nervous_system_MYH13$ID, peripheral_nervous_system_ccl$ID), ] <- peripheral_nervous_system_MYH13

peripheral_nervous_system_ccl$Hugo_Symbol <- as.character(peripheral_nervous_system_ccl$Hugo_Symbol)

peripheral_nervous_system_ccl$Hugo_Symbol[peripheral_nervous_system_ccl$Hugo_Symbol == "MYH13"] <- "has Mut"
peripheral_nervous_system_ccl$Hugo_Symbol[peripheral_nervous_system_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =peripheral_nervous_system_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(peripheral_nervous_system_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
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
