#soft_tissue 

soft_tissue <- UBA2_merged[which(UBA2_merged$lineage == "soft_tissue"),]
soft_tissue <- soft_tissue[order(soft_tissue$UBA2_gene_effect), ]
soft_tissue_scores <-soft_tissue$UBA2_gene_effect
soft_tissue_average <- mean(soft_tissue_scores)
soft_tissue_sd <- sd(soft_tissue_scores)
min(soft_tissue_scores)
max(soft_tissue_scores)

length(soft_tissue$UBA2_gene_effect)
length(which(soft_tissue$UBA2_gene_effect < soft_tissue_average))

enrich_title <- paste("UBA2 gene effect in", soft_tissue$lineage ,"cancer cells \n (",
                      length(which(soft_tissue$UBA2_gene_effect < soft_tissue_average)),
                      "cell lines out of ",length(soft_tissue$UBA2_gene_effect)," with UBA2 score <" ,
                      round(soft_tissue_average, digits = 3), ")")

library(ggplot2)
library("scales")
enrich <- ggplot(data = soft_tissue , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","white","white","white"),
                       guide = "colorbar", limits=c(min(soft_tissue_scores),0)) +
  scale_x_continuous(limits=c(min(soft_tissue_scores),0))+
  theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=11, angle = 90, vjust = 3, hjust= 1),
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())

#################################################################################################
########example with soft_tissue cells
soft_tissue_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "soft_tissue"), ]

#calculate the number of unique ccl
length(unique(soft_tissue_mutation$ccl))
#turns out to be 36 ccls

#subset of soft_tissue cell lines with below threshold gene effect score
soft_tissue_valid_mut <- soft_tissue_mutation[which(soft_tissue_mutation$UBA2_gene_effect <= soft_tissue_average), ]
soft_tissue_valid_mut <- droplevels(soft_tissue_valid_mut)
#subset of cell lines with above threshold
soft_tissue_invalid_mut <- soft_tissue_mutation[which(soft_tissue_mutation$UBA2_gene_effect > soft_tissue_average), ]
soft_tissue_invalid_mut <- droplevels(soft_tissue_invalid_mut)

#eliminate mutations that also occur in invalid
#soft_tissue_valid_mut1 <- soft_tissue_valid_mut[!duplicated(soft_tissue_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
soft_tissue_mut_final <- anti_join(soft_tissue_valid_mut, soft_tissue_invalid_mut, 
                             by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)


#calculate the second most frequent value
xx <- sort(table(soft_tissue_mut_final$Hugo_Symbol),decreasing=TRUE)[1:10]
xx
#                                  

soft_tissue_DSCAM <- soft_tissue_mut_final[which(soft_tissue_mut_final$Hugo_Symbol == "DSCAM"), ]
soft_tissue_no_DSCAM <- soft_tissue_mut_final[which(soft_tissue_mut_final$Hugo_Symbol != "DSCAM"), ]

#find the most frequent variant type
calculate_mode(soft_tissue_DSCAM$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with DSCAM_Missense
soft_tissue_DSCAM_Mis <- soft_tissue_DSCAM[which(soft_tissue_DSCAM$Variant_Classification == "Missense_Mutation"), ]
soft_tissue_DSCAM_No_Mis <- soft_tissue_DSCAM[which(soft_tissue_DSCAM$Variant_Classification != "Missense_Mutation"), ]
soft_tissue_DSCAM <- droplevels(soft_tissue_DSCAM)
#7 cell lines with DSCAM Missense mutation

#subset of all ccl
soft_tissue_ccl <- distinct(soft_tissue_mutation, ccl, .keep_all= TRUE )

y <- length(soft_tissue_DSCAM_Mis$ID)
x <- length(soft_tissue$ccl)
x <- x - y


final_table <- data.frame(mutation = c("DSCAM", "other"), counts = c(y, x))

bar_title <- paste("Mut: DSCAM Missense Mutation \n", y, "overlaps out of", x, "total")

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
soft_tissue_ccl <- distinct(soft_tissue_mutation, ccl, .keep_all= TRUE )

soft_tissue_ccl[match(soft_tissue_DSCAM$ID, soft_tissue_ccl$ID), ] <- soft_tissue_DSCAM

soft_tissue_ccl$Hugo_Symbol <- as.character(soft_tissue_ccl$Hugo_Symbol)

soft_tissue_ccl$Hugo_Symbol[soft_tissue_ccl$Hugo_Symbol == "DSCAM"] <- "has Mut"
soft_tissue_ccl$Hugo_Symbol[soft_tissue_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =soft_tissue_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(soft_tissue_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
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

