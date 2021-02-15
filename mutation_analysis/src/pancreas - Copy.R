#pancreas 

pancreas <- UBA2_merged[which(UBA2_merged$lineage == "pancreas"),]
pancreas <- pancreas[order(pancreas$UBA2_gene_effect), ]
pancreas_scores <-pancreas$UBA2_gene_effect
pancreas_average <- mean(pancreas_scores)
pancreas_sd <- sd(pancreas_scores)
min(pancreas_scores)
max(pancreas_scores)

length(pancreas$UBA2_gene_effect)
length(which(pancreas$UBA2_gene_effect < pancreas_average))

enrich_title <- paste("UBA2 gene effect in", pancreas$lineage ,"cancer cells \n (",
                      length(which(pancreas$UBA2_gene_effect < pancreas_average)),
                      "cell lines out of ",length(pancreas$UBA2_gene_effect)," with UBA2 score <" ,
                      round(pancreas_average, digits = 3), ")")

enrich <- ggplot(data = pancreas , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","white","white","white"),
                       guide = "colorbar", limits=c(min(pancreas_scores),0)) +
  scale_x_continuous(limits=c(min(pancreas_scores),0))+
  theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=11, angle = 90, vjust = 3, hjust= 1),
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())

#################################################################################################
########example with pancreas cells
pancreas_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "pancreas"), ]

#calculate the number of unique ccl
length(unique(pancreas_mutation$ccl))
#turns out to be 36 ccls

#subset of pancreas cell lines with below threshold gene effect score
pancreas_valid_mut <- pancreas_mutation[which(pancreas_mutation$UBA2_gene_effect <= pancreas_average), ]
pancreas_valid_mut <- droplevels(pancreas_valid_mut)
#subset of cell lines with above threshold
pancreas_invalid_mut <- pancreas_mutation[which(pancreas_mutation$UBA2_gene_effect > pancreas_average), ]
pancreas_invalid_mut <- droplevels(pancreas_invalid_mut)

#eliminate mutations that also occur in invalid
#pancreas_valid_mut1 <- pancreas_valid_mut[!duplicated(pancreas_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
pancreas_mut_final <- anti_join(pancreas_valid_mut, pancreas_invalid_mut, 
                             by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

#calculate the second most frequent value
xx <- sort(table(pancreas_mut_final$Hugo_Symbol),decreasing=TRUE)[1:300]
xx
# GPATCH8 HLA-DQB2 KIAA1109    KIF4B    MAPK6  METTL22    MROH1     MYH1    MYH15    NACAD  NCKAP1L    NFXL1   PABPC1     PCNX   PLXNA4               



###########################################################export###########################################################
##### eliminate mutations that occur less than 3
pancreas_mut_final <- pancreas_mut_final[pancreas_mut_final$Hugo_Symbol %in% names(which(table(pancreas_mut_final$Hugo_Symbol) > 2)), ]
length(unique(pancreas_mut_final$ccl)) # 17 cell lines 
length(unique(pancreas_mut_final$Hugo_Symbol)) # 39 muts 
library(xlsx)
write.xlsx(pancreas_mut_final, "pancreatic_mut.xlsx")
###########################################################export###########################################################




pancreas_GLTSCR1L <- pancreas_mut_final[which(pancreas_mut_final$Hugo_Symbol == "GLTSCR1L"), ]
pancreas_no_GLTSCR1L <- pancreas_mut_final[which(pancreas_mut_final$Hugo_Symbol != "GLTSCR1L"), ]

#find the most frequent variant type
calculate_mode(pancreas_GLTSCR1L$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with GLTSCR1L_Missense
pancreas_GLTSCR1L_Mis <- pancreas_GLTSCR1L[which(pancreas_GLTSCR1L$Variant_Classification == "Missense_Mutation"), ]
pancreas_GLTSCR1L_No_Mis <- pancreas_GLTSCR1L[which(pancreas_GLTSCR1L$Variant_Classification != "Missense_Mutation"), ]
pancreas_GLTSCR1L <- droplevels(pancreas_GLTSCR1L)
#7 cell lines with GLTSCR1L Missense_Mutation

#subset of all ccl
pancreas_ccl <- distinct(pancreas_mutation, ccl, .keep_all= TRUE )

y <- length(pancreas_GLTSCR1L_Mis$ID)
x <- length(pancreas$ccl)
x <- x - y


final_table <- data.frame(mutation = c("GLTSCR1L", "other"), counts = c(y, x))

bar_title <- paste("Mut: GLTSCR1L Missense_Mutation \n", y, "overlaps out of", x, "total")

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
pancreas_ccl <- distinct(pancreas_mutation, ccl, .keep_all= TRUE )

pancreas_ccl[match(pancreas_GLTSCR1L$ID, pancreas_ccl$ID), ] <- pancreas_GLTSCR1L

pancreas_ccl$Hugo_Symbol <- as.character(pancreas_ccl$Hugo_Symbol)

pancreas_ccl$Hugo_Symbol[pancreas_ccl$Hugo_Symbol == "GLTSCR1L"] <- "has Mut"
pancreas_ccl$Hugo_Symbol[pancreas_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =pancreas_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(pancreas_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
  ggtitle(box_title) +
  geom_boxplot(fill = "white", colour = "black", outlier.colour = "red")+
  theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
        axis.title.x=element_blank(),
        panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                        size = 2, linetype = "solid"))


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

