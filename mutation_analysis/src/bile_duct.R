#bile_duct 

bile_duct <- UBA2_merged[which(UBA2_merged$lineage == "bile_duct"),]
bile_duct <- bile_duct[order(bile_duct$UBA2_gene_effect), ]
bile_duct_scores <-bile_duct$UBA2_gene_effect
bile_duct_average <- mean(bile_duct_scores)
bile_duct_sd <- sd(bile_duct_scores)
min(bile_duct_scores)
max(bile_duct_scores)

length(bile_duct$UBA2_gene_effect)
length(which(bile_duct$UBA2_gene_effect < bile_duct_average))

enrich_title <- paste("UBA2 gene effect in", bile_duct$lineage ,"cancer cells \n (",
                      length(which(bile_duct$UBA2_gene_effect < bile_duct_average)),
                      "cell lines out of ",length(bile_duct$UBA2_gene_effect)," with UBA2 score <" ,
                      round(bile_duct_average, digits = 3), ")")

library(ggplot2)
library("scales")
enrich <- ggplot(data = bile_duct , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","white","white","white"),
                       guide = "colorbar", limits=c(min(bile_duct_scores),0)) +
  scale_x_continuous(limits=c(min(bile_duct_scores),0))+
  theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=11, angle = 90, vjust = 3, hjust= 1),
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())

#################################################################################################
########example with bile_duct cells
bile_duct_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "bile_duct"), ]

#calculate the number of unique ccl
length(unique(bile_duct_mutation$ccl))
#turns out to be 36 ccls

#subset of bile_duct cell lines with below threshold gene effect score
bile_duct_valid_mut <- bile_duct_mutation[which(bile_duct_mutation$UBA2_gene_effect <= bile_duct_average), ]
bile_duct_valid_mut <- droplevels(bile_duct_valid_mut)
#subset of cell lines with above threshold
bile_duct_invalid_mut <- bile_duct_mutation[which(bile_duct_mutation$UBA2_gene_effect > bile_duct_average), ]
bile_duct_invalid_mut <- droplevels(bile_duct_invalid_mut)

#eliminate mutations that also occur in invalid
#bile_duct_valid_mut1 <- bile_duct_valid_mut[!duplicated(bile_duct_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
bile_duct_mut_final <- anti_join(bile_duct_valid_mut, bile_duct_invalid_mut, 
                             by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

#calculate the mode (most frequent mutation)
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}
calculate_mode(bile_duct_mut_final$Hugo_Symbol)
#turns out to be DACH1!


bile_duct_DACH1 <- bile_duct_mut_final[which(bile_duct_mut_final$Hugo_Symbol == "DACH1"), ]
bile_duct_no_DACH1 <- bile_duct_mut_final[which(bile_duct_mut_final$Hugo_Symbol != "DACH1"), ]

#find the most frequent variant type
calculate_mode(bile_duct_DACH1$Variant_Classification)
#turns out to be In_Frame_Ins

#create subset with DACH1_Missense
bile_duct_DACH1_Mis <- bile_duct_DACH1[which(bile_duct_DACH1$Variant_Classification == "In_Frame_Ins"), ]
bile_duct_DACH1_No_Mis <- bile_duct_DACH1[which(bile_duct_DACH1$Variant_Classification != "In_Frame_Ins"), ]
bile_duct_DACH1 <- droplevels(bile_duct_DACH1)
#7 cell lines with DACH1 In_Frame_Ins

#subset of all ccl
bile_duct_ccl <- distinct(bile_duct_mutation, ccl, .keep_all= TRUE )

y <- length(bile_duct_DACH1_Mis$ID)
x <- length(bile_duct$ccl)
x <- x - y


final_table <- data.frame(mutation = c("DACH1", "other"), counts = c(y, x))

bar_title <- paste("Mut: DACH1 In_Frame_Ins \n", y, "overlaps out of", x, "total")

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
bile_duct_ccl <- distinct(bile_duct_mutation, ccl, .keep_all= TRUE )

bile_duct_ccl[match(bile_duct_DACH1$ID, bile_duct_ccl$ID), ] <- bile_duct_DACH1

bile_duct_ccl$Hugo_Symbol <- as.character(bile_duct_ccl$Hugo_Symbol)

bile_duct_ccl$Hugo_Symbol[bile_duct_ccl$Hugo_Symbol == "DACH1"] <- "has Mut"
bile_duct_ccl$Hugo_Symbol[bile_duct_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =bile_duct_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(bile_duct_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
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
