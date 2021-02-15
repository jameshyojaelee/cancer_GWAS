#blood 

blood <- UBA2_merged[which(UBA2_merged$lineage == "blood"),]
blood <- blood[order(blood$UBA2_gene_effect), ]
blood_scores <-blood$UBA2_gene_effect
blood_average <- mean(blood_scores)
blood_sd <- sd(blood_scores)
min(blood_scores)
max(blood_scores)

length(blood$UBA2_gene_effect)
length(which(blood$UBA2_gene_effect < blood_average))

enrich_title <- paste("UBA2 gene effect in", blood$lineage ,"cancer cells \n (",
                      length(which(blood$UBA2_gene_effect < blood_average)),
                      "cell lines out of ",length(blood$UBA2_gene_effect)," with UBA2 score <" ,
                      round(blood_average, digits = 3), ")")

library(ggplot2)
library("scales")
enrich <- ggplot(data = blood , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","red", "white","white","white","white"),
                       guide = "colorbar", limits=c(min(blood_scores),0)) +
  scale_x_continuous(limits=c(min(blood_scores),0))+
  theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=11, angle = 90, vjust = 3, hjust= 1),
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())
enrich
#################################################################################################
########example with blood cells
blood_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "blood"), ]

#calculate the number of unique ccl
length(unique(blood_mutation$ccl))
#turns out to be 36 ccls

#subset of blood cell lines with below threshold gene effect score
blood_valid_mut <- blood_mutation[which(blood_mutation$UBA2_gene_effect <= blood_average), ]
blood_valid_mut <- droplevels(blood_valid_mut)
#subset of cell lines with above threshold
blood_invalid_mut <- blood_mutation[which(blood_mutation$UBA2_gene_effect > blood_average), ]
blood_invalid_mut <- droplevels(blood_invalid_mut)

#eliminate mutations that also occur in invalid
#blood_valid_mut1 <- blood_valid_mut[!duplicated(blood_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
blood_mut_final <- anti_join(blood_valid_mut, blood_invalid_mut, 
                                  by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

#calculate the mode (most frequent mutation)
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}
calculate_mode(blood_mut_final$Hugo_Symbol)
#turns out to be WDR87!


blood_WDR87 <- blood_mut_final[which(blood_mut_final$Hugo_Symbol == "WDR87"), ]
blood_no_WDR87 <- blood_mut_final[which(blood_mut_final$Hugo_Symbol != "WDR87"), ]

#find the most frequent variant type
calculate_mode(blood_WDR87$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with WDR87_Missense
blood_WDR87_Mis <- blood_WDR87[which(blood_WDR87$Variant_Classification == "Missense_Mutation"), ]
blood_WDR87_No_Mis <- blood_WDR87[which(blood_WDR87$Variant_Classification != "Missense_Mutation"), ]
blood_WDR87 <- droplevels(blood_WDR87)
#7 cell lines with WDR87 Missense mutation

#subset of all ccl
blood_ccl <- distinct(blood_mutation, ccl, .keep_all= TRUE )

y <- length(blood_WDR87_Mis$ID)
x <- length(blood$ccl)
x <- x - y


final_table <- data.frame(mutation = c("WDR87", "other"), counts = c(y, x))

bar_title <- paste("Mut: WDR87 Missense Mutation \n", y, "overlaps out of", x, "total")

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
blood_ccl <- distinct(blood_mutation, ccl, .keep_all= TRUE )

blood_ccl[match(blood_WDR87$ID, blood_ccl$ID), ] <- blood_WDR87

blood_ccl$Hugo_Symbol <- as.character(blood_ccl$Hugo_Symbol)

blood_ccl$Hugo_Symbol[blood_ccl$Hugo_Symbol == "WDR87"] <- "has Mut"
blood_ccl$Hugo_Symbol[blood_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =blood_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(blood_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
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

