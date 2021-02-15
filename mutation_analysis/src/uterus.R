#uterus 

uterus <- UBA2_merged[which(UBA2_merged$lineage == "uterus"),]
uterus <- uterus[order(uterus$UBA2_gene_effect), ]
uterus_scores <-uterus$UBA2_gene_effect
uterus_average <- mean(uterus_scores)
uterus_sd <- sd(uterus_scores)
min(uterus_scores)
max(uterus_scores)

length(uterus$UBA2_gene_effect)
length(which(uterus$UBA2_gene_effect < uterus_average))

enrich_title <- paste("UBA2 gene effect in", uterus$lineage ,"cancer cells \n (",
                      length(which(uterus$UBA2_gene_effect < uterus_average)),
                      "cell lines out of ",length(uterus$UBA2_gene_effect)," with UBA2 score <" ,
                      round(uterus_average, digits = 3), ")")

library(ggplot2)
library("scales")
enrich <- ggplot(data = uterus , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","white","white"),
                       guide = "colorbar", limits=c(min(uterus_scores),0)) +
  scale_x_continuous(limits=c(min(uterus_scores),0))+
  theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=11, angle = 90, vjust = 3, hjust= 1),
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())

#################################################################################################
########example with uterus cells
uterus_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "uterus"), ]

#calculate the number of unique ccl
length(unique(uterus_mutation$ccl))
#turns out to be 36 ccls

#subset of uterus cell lines with below threshold gene effect score
uterus_valid_mut <- uterus_mutation[which(uterus_mutation$UBA2_gene_effect <= uterus_average), ]
uterus_valid_mut <- droplevels(uterus_valid_mut)
#subset of cell lines with above threshold
uterus_invalid_mut <- uterus_mutation[which(uterus_mutation$UBA2_gene_effect > uterus_average), ]
uterus_invalid_mut <- droplevels(uterus_invalid_mut)

#eliminate mutations that also occur in invalid
#uterus_valid_mut1 <- uterus_valid_mut[!duplicated(uterus_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
uterus_mut_final <- anti_join(uterus_valid_mut, uterus_invalid_mut, 
                             by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

#calculate the mode (most frequent mutation)
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}
calculate_mode(uterus_mut_final$Hugo_Symbol)
#turns out to be ZNF318!


uterus_ZNF318 <- uterus_mut_final[which(uterus_mut_final$Hugo_Symbol == "ZNF318"), ]
uterus_no_ZNF318 <- uterus_mut_final[which(uterus_mut_final$Hugo_Symbol != "ZNF318"), ]

#find the most frequent variant type
calculate_mode(uterus_ZNF318$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with ZNF318_Missense
uterus_ZNF318_Mis <- uterus_ZNF318[which(uterus_ZNF318$Variant_Classification == "Missense_Mutation"), ]
uterus_ZNF318_No_Mis <- uterus_ZNF318[which(uterus_ZNF318$Variant_Classification != "Missense_Mutation"), ]
uterus_ZNF318 <- droplevels(uterus_ZNF318)
#7 cell lines with ZNF318 Missense mutation

#subset of all ccl
uterus_ccl <- distinct(uterus_mutation, ccl, .keep_all= TRUE )

y <- length(uterus_ZNF318_Mis$ID)
x <- length(uterus$ccl)
x <- x - y


final_table <- data.frame(mutation = c("ZNF318", "other"), counts = c(y, x))

bar_title <- paste("Mut: ZNF318 Missense Mutation \n", y, "overlaps out of", x, "total")

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
uterus_ccl <- distinct(uterus_mutation, ccl, .keep_all= TRUE )

uterus_ccl[match(uterus_ZNF318$ID, uterus_ccl$ID), ] <- uterus_ZNF318

uterus_ccl$Hugo_Symbol <- as.character(uterus_ccl$Hugo_Symbol)

uterus_ccl$Hugo_Symbol[uterus_ccl$Hugo_Symbol == "ZNF318"] <- "has Mut"
uterus_ccl$Hugo_Symbol[uterus_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =uterus_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(uterus_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
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
