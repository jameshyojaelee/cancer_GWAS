#bone 

bone <- UBA2_merged[which(UBA2_merged$lineage == "bone"),]
bone <- bone[order(bone$UBA2_gene_effect), ]
bone_scores <-bone$UBA2_gene_effect
bone_average <- mean(bone_scores)
bone_sd <- sd(bone_scores)
min(bone_scores)
max(bone_scores)

length(bone$UBA2_gene_effect)
length(which(bone$UBA2_gene_effect < bone_average))

enrich_title <- paste("UBA2 gene effect in", bone$lineage ,"cancer cells \n (",
                      length(which(bone$UBA2_gene_effect < bone_average)),
                      "cell lines out of ",length(bone$UBA2_gene_effect)," with UBA2 score <" ,
                      round(bone_average, digits = 3), ")")

library(ggplot2)
library("scales")
enrich <- ggplot(data = bone , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","white","white","white"),
                       guide = "colorbar", limits=c(min(bone_scores),0)) +
  scale_x_continuous(limits=c(min(bone_scores),0))+
  theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=11, angle = 90, vjust = 3, hjust= 1),
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())

#################################################################################################
########example with bone cells
bone_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "bone"), ]

#calculate the number of unique ccl
length(unique(bone_mutation$ccl))
#turns out to be 36 ccls

#subset of bone cell lines with below threshold gene effect score
bone_valid_mut <- bone_mutation[which(bone_mutation$UBA2_gene_effect <= bone_average), ]
bone_valid_mut <- droplevels(bone_valid_mut)
#subset of cell lines with above threshold
bone_invalid_mut <- bone_mutation[which(bone_mutation$UBA2_gene_effect > bone_average), ]
bone_invalid_mut <- droplevels(bone_invalid_mut)

#eliminate mutations that also occur in invalid
#bone_valid_mut1 <- bone_valid_mut[!duplicated(bone_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
bone_mut_final <- anti_join(bone_valid_mut, bone_invalid_mut, 
                                  by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

#calculate the mode (most frequent mutation)
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}
calculate_mode(bone_mut_final$Hugo_Symbol)
#turns out to be KMT2A!


bone_KMT2A <- bone_mut_final[which(bone_mut_final$Hugo_Symbol == "KMT2A"), ]
bone_no_KMT2A <- bone_mut_final[which(bone_mut_final$Hugo_Symbol != "KMT2A"), ]

#find the most frequent variant type
calculate_mode(bone_KMT2A$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with KMT2A_Missense
bone_KMT2A_Mis <- bone_KMT2A[which(bone_KMT2A$Variant_Classification == "Missense_Mutation"), ]
bone_KMT2A_No_Mis <- bone_KMT2A[which(bone_KMT2A$Variant_Classification != "Missense_Mutation"), ]
bone_KMT2A <- droplevels(bone_KMT2A)
#7 cell lines with KMT2A Missense mutation

#subset of all ccl
bone_ccl <- distinct(bone_mutation, ccl, .keep_all= TRUE )

y <- length(bone_KMT2A_Mis$ID)
x <- length(bone$ccl)
x <- x - y


final_table <- data.frame(mutation = c("KMT2A", "other"), counts = c(y, x))

bar_title <- paste("Mut: KMT2A Missense Mutation \n", y, "overlaps out of", x, "total")

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
bone_ccl <- distinct(bone_mutation, ccl, .keep_all= TRUE )

bone_ccl[match(bone_KMT2A$ID, bone_ccl$ID), ] <- bone_KMT2A

bone_ccl$Hugo_Symbol <- as.character(bone_ccl$Hugo_Symbol)

bone_ccl$Hugo_Symbol[bone_ccl$Hugo_Symbol == "KMT2A"] <- "has Mut"
bone_ccl$Hugo_Symbol[bone_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =bone_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(bone_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
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
