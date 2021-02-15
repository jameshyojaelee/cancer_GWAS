#central_nervous_system 

central_nervous_system <- UBA2_merged[which(UBA2_merged$lineage == "central_nervous_system"),]
central_nervous_system <- central_nervous_system[order(central_nervous_system$UBA2_gene_effect), ]
central_nervous_system_scores <-central_nervous_system$UBA2_gene_effect
central_nervous_system_average <- mean(central_nervous_system_scores)
central_nervous_system_sd <- sd(central_nervous_system_scores)
min(central_nervous_system_scores)
max(central_nervous_system_scores)

length(central_nervous_system$UBA2_gene_effect)
length(which(central_nervous_system$UBA2_gene_effect < central_nervous_system_average))

enrich_title <- paste("UBA2 gene effect in", central_nervous_system$lineage ,"cancer cells \n (",
                      length(which(central_nervous_system$UBA2_gene_effect < central_nervous_system_average)),
                      "cell lines out of ",length(central_nervous_system$UBA2_gene_effect)," with UBA2 score <" ,
                      round(central_nervous_system_average, digits = 3), ")")

library(ggplot2)
library("scales")
enrich <- ggplot(data = central_nervous_system , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","white","white","white"),
                       guide = "colorbar", limits=c(min(central_nervous_system_scores),0)) +
  scale_x_continuous(limits=c(min(central_nervous_system_scores),0))+
  theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=11, angle = 90, vjust = 3, hjust= 1),
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())

#################################################################################################
########example with central_nervous_system cells
central_nervous_system_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "central_nervous_system"), ]

#calculate the number of unique ccl
length(unique(central_nervous_system_mutation$ccl))
#turns out to be 36 ccls

#subset of central_nervous_system cell lines with below threshold gene effect score
central_nervous_system_valid_mut <- central_nervous_system_mutation[which(central_nervous_system_mutation$UBA2_gene_effect <= central_nervous_system_average), ]
central_nervous_system_valid_mut <- droplevels(central_nervous_system_valid_mut)
#subset of cell lines with above threshold
central_nervous_system_invalid_mut <- central_nervous_system_mutation[which(central_nervous_system_mutation$UBA2_gene_effect > central_nervous_system_average), ]
central_nervous_system_invalid_mut <- droplevels(central_nervous_system_invalid_mut)

#eliminate mutations that also occur in invalid
#central_nervous_system_valid_mut1 <- central_nervous_system_valid_mut[!duplicated(central_nervous_system_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
central_nervous_system_mut_final <- anti_join(central_nervous_system_valid_mut, central_nervous_system_invalid_mut, 
                                  by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

#calculate the mode (most frequent mutation)
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}
calculate_mode(central_nervous_system_mut_final$Hugo_Symbol)
#turns out to be SETX!


central_nervous_system_SETX <- central_nervous_system_mut_final[which(central_nervous_system_mut_final$Hugo_Symbol == "SETX"), ]
central_nervous_system_no_SETX <- central_nervous_system_mut_final[which(central_nervous_system_mut_final$Hugo_Symbol != "SETX"), ]

#find the most frequent variant type
calculate_mode(central_nervous_system_SETX$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with SETX_Missense
central_nervous_system_SETX_Mis <- central_nervous_system_SETX[which(central_nervous_system_SETX$Variant_Classification == "Missense_Mutation"), ]
central_nervous_system_SETX_No_Mis <- central_nervous_system_SETX[which(central_nervous_system_SETX$Variant_Classification != "Missense_Mutation"), ]
central_nervous_system_SETX <- droplevels(central_nervous_system_SETX)
#7 cell lines with SETX Missense mutation

#subset of all ccl
central_nervous_system_ccl <- distinct(central_nervous_system_mutation, ccl, .keep_all= TRUE )

y <- length(central_nervous_system_SETX_Mis$ID)
x <- length(central_nervous_system$ccl)
x <- x - y


final_table <- data.frame(mutation = c("SETX", "other"), counts = c(y, x))

bar_title <- paste("Mut: SETX Missense Mutation \n", y, "overlaps out of", x, "total")

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
central_nervous_system_ccl <- distinct(central_nervous_system_mutation, ccl, .keep_all= TRUE )

central_nervous_system_ccl[match(central_nervous_system_SETX$ID, central_nervous_system_ccl$ID), ] <- central_nervous_system_SETX

central_nervous_system_ccl$Hugo_Symbol <- as.character(central_nervous_system_ccl$Hugo_Symbol)

central_nervous_system_ccl$Hugo_Symbol[central_nervous_system_ccl$Hugo_Symbol == "SETX"] <- "has Mut"
central_nervous_system_ccl$Hugo_Symbol[central_nervous_system_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =central_nervous_system_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(central_nervous_system_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
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
