#cervix 

cervix <- UBA2_merged[which(UBA2_merged$lineage == "cervix"),]
cervix <- cervix[order(cervix$UBA2_gene_effect), ]
cervix_scores <-cervix$UBA2_gene_effect
cervix_average <- mean(cervix_scores)
cervix_sd <- sd(cervix_scores)
min(cervix_scores)
max(cervix_scores)

length(cervix$UBA2_gene_effect)
length(which(cervix$UBA2_gene_effect < cervix_average))

enrich_title <- paste("UBA2 gene effect in", cervix$lineage ,"cancer cells \n (",
                      length(which(cervix$UBA2_gene_effect < cervix_average)),
                      "cell lines out of ",length(cervix$UBA2_gene_effect)," with UBA2 score <" ,
                      round(cervix_average, digits = 3), ")")

library(ggplot2)
library("scales")
enrich <- ggplot(data = cervix , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","white","white","white"),
                       guide = "colorbar", limits=c(min(cervix_scores),0)) +
  scale_x_continuous(limits=c(min(cervix_scores),0))+
  theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=11, angle = 90, vjust = 3, hjust= 1),
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())

#################################################################################################
########example with cervix cells
cervix_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "cervix"), ]

#calculate the number of unique ccl
length(unique(cervix_mutation$ccl))
#turns out to be 36 ccls

#subset of cervix cell lines with below threshold gene effect score
cervix_valid_mut <- cervix_mutation[which(cervix_mutation$UBA2_gene_effect <= cervix_average), ]
cervix_valid_mut <- droplevels(cervix_valid_mut)
#subset of cell lines with above threshold
cervix_invalid_mut <- cervix_mutation[which(cervix_mutation$UBA2_gene_effect > cervix_average), ]
cervix_invalid_mut <- droplevels(cervix_invalid_mut)

#eliminate mutations that also occur in invalid
#cervix_valid_mut1 <- cervix_valid_mut[!duplicated(cervix_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
cervix_mut_final <- anti_join(cervix_valid_mut, cervix_invalid_mut, 
                             by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

#calculate the mode (most frequent mutation)
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}
calculate_mode(cervix_mut_final$Hugo_Symbol)
#turns out to be TENM4!


cervix_TENM4 <- cervix_mut_final[which(cervix_mut_final$Hugo_Symbol == "TENM4"), ]
cervix_no_TENM4 <- cervix_mut_final[which(cervix_mut_final$Hugo_Symbol != "TENM4"), ]

#find the most frequent variant type
calculate_mode(cervix_TENM4$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with TENM4_Missense
cervix_TENM4_Mis <- cervix_TENM4[which(cervix_TENM4$Variant_Classification == "Missense_Mutation"), ]
cervix_TENM4_No_Mis <- cervix_TENM4[which(cervix_TENM4$Variant_Classification != "Missense_Mutation"), ]
cervix_TENM4 <- droplevels(cervix_TENM4)
#7 cell lines with TENM4 Missense mutation

#subset of all ccl
cervix_ccl <- distinct(cervix_mutation, ccl, .keep_all= TRUE )

y <- length(cervix_TENM4_Mis$ID)
x <- length(cervix$ccl)
x <- x - y


final_table <- data.frame(mutation = c("TENM4", "other"), counts = c(y, x))

bar_title <- paste("Mut: TENM4 Missense Mutation \n", y, "overlaps out of", x, "total")

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
cervix_ccl <- distinct(cervix_mutation, ccl, .keep_all= TRUE )

cervix_ccl[match(cervix_TENM4$ID, cervix_ccl$ID), ] <- cervix_TENM4

cervix_ccl$Hugo_Symbol <- as.character(cervix_ccl$Hugo_Symbol)

cervix_ccl$Hugo_Symbol[cervix_ccl$Hugo_Symbol == "TENM4"] <- "has Mut"
cervix_ccl$Hugo_Symbol[cervix_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =cervix_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(cervix_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
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
