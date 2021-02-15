#upper_aerodigestive 

upper_aerodigestive <- UBA2_merged[which(UBA2_merged$lineage == "upper_aerodigestive"),]
upper_aerodigestive <- upper_aerodigestive[order(upper_aerodigestive$UBA2_gene_effect), ]
upper_aerodigestive_scores <-upper_aerodigestive$UBA2_gene_effect
upper_aerodigestive_average <- mean(upper_aerodigestive_scores)
upper_aerodigestive_sd <- sd(upper_aerodigestive_scores)
min(upper_aerodigestive_scores)
max(upper_aerodigestive_scores)

length(upper_aerodigestive$UBA2_gene_effect)
length(which(upper_aerodigestive$UBA2_gene_effect < upper_aerodigestive_average))

enrich_title <- paste("UBA2 gene effect in", upper_aerodigestive$lineage ,"cancer cells \n (",
                      length(which(upper_aerodigestive$UBA2_gene_effect < upper_aerodigestive_average)),
                      "cell lines out of ",length(upper_aerodigestive$UBA2_gene_effect)," with UBA2 score <" ,
                      round(upper_aerodigestive_average, digits = 3), ")")

library(ggplot2)
library("scales")
enrich <- ggplot(data = upper_aerodigestive , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","white","white","white"),
                       guide = "colorbar", limits=c(min(upper_aerodigestive_scores),0)) +
  scale_x_continuous(limits=c(min(upper_aerodigestive_scores),0))+
  theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=11, angle = 90, vjust = 3, hjust= 1),
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())

#################################################################################################
########example with upper_aerodigestive cells
upper_aerodigestive_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "upper_aerodigestive"), ]

#calculate the number of unique ccl
length(unique(upper_aerodigestive_mutation$ccl))
#turns out to be 36 ccls

#subset of upper_aerodigestive cell lines with below threshold gene effect score
upper_aerodigestive_valid_mut <- upper_aerodigestive_mutation[which(upper_aerodigestive_mutation$UBA2_gene_effect <= upper_aerodigestive_average), ]
upper_aerodigestive_valid_mut <- droplevels(upper_aerodigestive_valid_mut)
#subset of cell lines with above threshold
upper_aerodigestive_invalid_mut <- upper_aerodigestive_mutation[which(upper_aerodigestive_mutation$UBA2_gene_effect > upper_aerodigestive_average), ]
upper_aerodigestive_invalid_mut <- droplevels(upper_aerodigestive_invalid_mut)

#eliminate mutations that also occur in invalid
#upper_aerodigestive_valid_mut1 <- upper_aerodigestive_valid_mut[!duplicated(upper_aerodigestive_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
upper_aerodigestive_mut_final <- anti_join(upper_aerodigestive_valid_mut, upper_aerodigestive_invalid_mut, 
                             by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

#calculate the mode (most frequent mutation)
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}
calculate_mode(upper_aerodigestive_mut_final$Hugo_Symbol)
#turns out to be TP53!


upper_aerodigestive_TP53 <- upper_aerodigestive_mut_final[which(upper_aerodigestive_mut_final$Hugo_Symbol == "TP53"), ]
upper_aerodigestive_no_TP53 <- upper_aerodigestive_mut_final[which(upper_aerodigestive_mut_final$Hugo_Symbol != "TP53"), ]

#find the most frequent variant type
calculate_mode(upper_aerodigestive_TP53$Variant_Classification)
#turns out to be Nonsense_Mutation

#create subset with TP53_Missense
upper_aerodigestive_TP53_Mis <- upper_aerodigestive_TP53[which(upper_aerodigestive_TP53$Variant_Classification == "Nonsense_Mutation"), ]
upper_aerodigestive_TP53_No_Mis <- upper_aerodigestive_TP53[which(upper_aerodigestive_TP53$Variant_Classification != "Nonsense_Mutation"), ]
upper_aerodigestive_TP53 <- droplevels(upper_aerodigestive_TP53)
#7 cell lines with TP53 Nonsense_Mutation

#subset of all ccl
upper_aerodigestive_ccl <- distinct(upper_aerodigestive_mutation, ccl, .keep_all= TRUE )

y <- length(upper_aerodigestive_TP53_Mis$ID)
x <- length(upper_aerodigestive$ccl)
x <- x - y


final_table <- data.frame(mutation = c("TP53", "other"), counts = c(y, x))

bar_title <- paste("Mut: TP53 Nonsense_Mutation \n", y, "overlaps out of", x, "total")

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
upper_aerodigestive_ccl <- distinct(upper_aerodigestive_mutation, ccl, .keep_all= TRUE )

upper_aerodigestive_ccl[match(upper_aerodigestive_TP53$ID, upper_aerodigestive_ccl$ID), ] <- upper_aerodigestive_TP53

upper_aerodigestive_ccl$Hugo_Symbol <- as.character(upper_aerodigestive_ccl$Hugo_Symbol)

upper_aerodigestive_ccl$Hugo_Symbol[upper_aerodigestive_ccl$Hugo_Symbol == "TP53"] <- "has Mut"
upper_aerodigestive_ccl$Hugo_Symbol[upper_aerodigestive_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =upper_aerodigestive_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(upper_aerodigestive_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
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
