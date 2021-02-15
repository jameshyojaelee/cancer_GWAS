library(devtools)
library(ggplot2)
library("scales")
library(dplyr)
library(grid)
library(gridExtra)


#read necessary csv files
setwd("/Users/james/Desktop/Sumoylation_Analysis/project")
gene_effect <- read.csv("C:/Users/james/Desktop/Sumoylation_Analysis/data/Achilles_gene_effect.csv") #CRISPR(Avana) data
cclID <- read.csv("C:/Users/james/Desktop/Sumoylation_Analysis/data/sample_info.csv") #list of all cell lines, ID and lineage
ccl_mut <- read.csv("C:/Users/james/Desktop/Sumoylation_Analysis/data/CCLE_mutations.csv") #import mutation info
exp <- read.csv("C:/Users/james/Desktop/Sumoylation_Analysis/data/CCLE_expression.csv", fileEncoding="UTF-8-BOM") #gene expression data for just protein coding genes

#dataframe with the cancer cell line names
ccl_df <- data.frame(ID = cclID$DepMap_ID, ccl = cclID$stripped_cell_line_name, lineage = cclID$lineage)

#dataframe with the gene effect score (example with UBA2 Gene ONLY)
UBA2_df <- data.frame(ID = gene_effect$DepMap_ID, UBA2_gene_effect = gene_effect$UBA2)

#Merge cell line info and gene effect score
UBA2_merged <- merge(x = ccl_df, 
                  y = UBA2_df,
                  by = "ID")
UBA2_merged <- droplevels(UBA2_merged)

################################################################################################################
#part 1 - UBA2 gene effect per lineage. below -1 is a viable threshold

#boxplot with different cell types and their gene effects
boxplot(UBA2_gene_effect ~ lineage,
        data=UBA2_merged,
        main="UBA2 gene KO effect per lineage",
        xlab="Lineage",
        ylab="Gene effect",
        col="red",
        border="black"
)

#(optional) individual dots added
stripchart(UBA2_gene_effect ~ lineage, vertical = TRUE, data = UBA2_merged, 
           method = "jitter", add = TRUE, pch = 21, col = 'black')


################################################################################################################
#part2 - Enrichment Analysis

#(a) Explore the UBA2 effect scores on each lineage

#example with colorectal 

colorectal <- UBA2_merged[which(UBA2_merged$lineage == "colorectal"),]
colorectal <- colorectal[order(colorectal$UBA2_gene_effect), ]
colorectal_scores <-colorectal$UBA2_gene_effect
colorectal_average <- mean(colorectal_scores)
colorectal_sd <- sd(colorectal_scores)
min(colorectal_scores)
max(colorectal_scores)

length(colorectal$UBA2_gene_effect)
length(which(colorectal$UBA2_gene_effect < colorectal_average))

enrich_title <- paste("UBA2 gene effect in", colorectal$lineage ,"cancer cells \n (",
                      length(which(colorectal$UBA2_gene_effect < colorectal_average)),
                      "cell lines out of ",length(colorectal$UBA2_gene_effect)," with UBA2 score <" ,
                      round(colorectal_average, digits = 3), ")")


enrich <- ggplot(data = colorectal , aes(x = UBA2_gene_effect, y = lineage, width= 0.05)) +
  ggtitle(enrich_title) +
  geom_tile(aes(fill = UBA2_gene_effect)) +
  scale_fill_gradientn(colours = c("red","white","white","white"),
                       guide = "colorbar", limits=c(min(colorectal_scores),0)) +
  scale_x_continuous(limits=c(min(colorectal_scores),0))+
  theme_bw()  + theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"), 
                      panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                      size = 2, linetype = "solid"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())
enrich




#Build a function to repeat the same enrichment analysis with other cell lines

enrichment <- function(cell_line){
  x <- UBA2_merged[which(UBA2_merged$lineage == cell_line),]
  x <- x[order(x$UBA2_gene_effect), ]
  scores <-x$UBA2_gene_effect
  average <- mean(scores)
  sdev <- sd(scores)

  enrich_title <- paste("UBA2 gene knockout effect in", x$lineage ,"cancer cells \n (",
                        length(which(x$UBA2_gene_effect < average)),
                        "cell lines out of ",length(x$UBA2_gene_effect)," with UBA2 score <" ,
                        round(average, digits = 3), ")")
  
  enrich <- ggplot(data = x , aes(x = UBA2_gene_effect, y = lineage, width= 0.1)) +
    ggtitle(enrich_title) +
    geom_tile(aes(fill = UBA2_gene_effect)) +
    scale_fill_gradientn(colours = c("red","white","white"),
                         guide = "colorbar", limits=c(min(scores, na.rm=TRUE),0)) +
    scale_x_continuous(limits=c(min(scores),0))+
    theme_bw()  + theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
                        panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                                        size = 2, linetype = "solid"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank())
  
  return (enrich)
}
enrichment("bile_duct")
#lineage <- unique(UBA2_merged$lineage)
#lineage <- as.list(lineage)
#all_enrich <- list()

#for (i in length(lineage)){
#  assign(all_enrich[i], enrichment(lineage[i]))
#  i <- i + 1
#}


################################################################################################################
#(b) Common mutation

colnames(ccl_mut)[which(names(ccl_mut) == "DepMap_ID")] <- "ID"
ccl_mut <- subset(ccl_mut, select = c("ID", "Hugo_Symbol", "Variant_Classification"))

#import the auc data to manipulate
gene_effect_mutation <- merge(x = UBA2_merged, 
                            y = ccl_mut,
                            by = "ID")

#eliminate silent mutations
gene_effect_mutation <- gene_effect_mutation[which(gene_effect_mutation$Variant_Classification != "Silent"),]

#eliminate duplicated rows
gene_effect_mutation <- unique(gene_effect_mutation)



#############################example with colorectal cells############################
colorectal_mutation <- gene_effect_mutation[which(gene_effect_mutation$lineage == "colorectal"), ]

#calculate the number of unique ccl
length(unique(colorectal_mutation$ccl))
#turns out to be 36 ccls

#subset of colorectal cell lines with below threshold gene effect score
colorectal_valid_mut <- colorectal_mutation[which(colorectal_mutation$UBA2_gene_effect <= colorectal_average), ]
colorectal_valid_mut <- droplevels(colorectal_valid_mut)
#subset of cell lines with above threshold
colorectal_invalid_mut <- colorectal_mutation[which(colorectal_mutation$UBA2_gene_effect > colorectal_average), ]
colorectal_invalid_mut <- droplevels(colorectal_invalid_mut)

#eliminate mutations that also occur in invalid
#colorectal_valid_mut1 <- colorectal_valid_mut[!duplicated(colorectal_invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
colorectal_mut_final <- anti_join(colorectal_valid_mut, colorectal_invalid_mut, 
                       by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)

#calculate the mode (most frequent mutation)
calculate_mode <- function(x) {
        uniqx <- unique(na.omit(x))
        uniqx[which.max(tabulate(match(x, uniqx)))]
}
calculate_mode(colorectal_mut_final$Hugo_Symbol)
#turns out to be HUWE1

colorectal_HUWE1 <- colorectal_mut_final[which(colorectal_mut_final$Hugo_Symbol == "HUWE1"), ]

#find the most frequent variant type
calculate_mode(colorectal_HUWE1$Variant_Classification)
#turns out to be Missense_Mutation

#create subset with HUWE1_Missense
colorectal_HUWE1_Mis <- colorectal_HUWE1[which(colorectal_HUWE1$Variant_Classification == "Missense_Mutation"), ]
colorectal_HUWE1 <- droplevels(colorectal_HUWE1)
colorectal_HUWE1_Mis <-droplevels(colorectal_HUWE1_Mis)
#7 cell lines with HUWE1 Missense mutation

y <- length(colorectal_HUWE1_Mis$ccl)
x <- length(colorectal$ccl)
x <- x - y

final_table <- data.frame(mutation = c("HUWE1", "other"), counts = c(y, x))

bar_title <- paste("Mut: HUWE1 Missense Mutation \n", y, "overlaps out of", x, "total")

bar <- ggplot(data = final_table, aes(x= mutation, y=counts)) +
        ggtitle(bar_title) +
        geom_bar(stat="identity", fill = "white", colour = "black")+
        theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"), 
              panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                              size = 2, linetype = "solid"))

#(c) mut, lacks mut

#example with colorectal
colorectal_ccl <- distinct(colorectal_mutation, ccl, .keep_all= TRUE )

colorectal_ccl[match(colorectal_HUWE1$ID, colorectal_ccl$ID), ] <- colorectal_HUWE1

colorectal_ccl$Hugo_Symbol <- as.character(colorectal_ccl$Hugo_Symbol)

colorectal_ccl$Hugo_Symbol[colorectal_ccl$Hugo_Symbol == "HUWE1"] <- "has Mut"
colorectal_ccl$Hugo_Symbol[colorectal_ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"

t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data =colorectal_ccl) 

box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))

box <- ggplot(colorectal_ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
  ggtitle(box_title) +
  geom_boxplot(fill = "white", colour = "black", outlier.colour = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"), 
        panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                        size = 2, linetype = "solid"))



##create function to repeat the process (b) and (c)  

mutation_analysis <- function(cell_line){
  enrich <- enrichment(cell_line)
  mut <- gene_effect_mutation[which(gene_effect_mutation$lineage == cell_line), ]

  
  #subset of colorectal cell lines with below threshold gene effect score
  valid_mut <- mut[which(mut$UBA2_gene_effect <= colorectal_average), ]
  valid_mut <- droplevels(valid_mut)
  #subset of cell lines with above threshold
  invalid_mut <- mut[which(mut$UBA2_gene_effect > colorectal_average), ]
  invalid_mut <- droplevels(invalid_mut)
  
  #eliminate mutations that also occur in invalid
  #colorectal_valid_mut1 <- valid_mut[!duplicated(invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
  mut_final <- anti_join(valid_mut, invalid_mut, 
                                    by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)
  
  #calculate the mode (most frequent mutation)
  mutated_gene <- calculate_mode(mut_final$Hugo_Symbol)

  #info with the mutated_gene above
  mut_mutated_gene <- mut_final[which(mut_final$Hugo_Symbol == mutated_gene), ]
 
  #find the most frequent variant type
  variant <- calculate_mode(mut_mutated_gene$Variant_Classification)
  
  #create subset with the specific variant from above
  mut_mutated_gene_var <- mut_mutated_gene[which(mut_mutated_gene$Variant_Classification == variant), ]
  
  mut_mutated_gene <- droplevels(mut_mutated_gene)
  #7 cell lines with HUWE1 Missense mutation
  
  y <- length(mut_mutated_gene_var$ccl)
  x <- length(unique(mut$ccl))
  x <- x - y
  
  final_table <- data.frame(mutation = c(as.character(variant), "other"), counts = c(y, x))
  
  bar_title <- paste("Mut:", as.character(mutated_gene), as.character(variant),
                     "Mutation \n", y, "overlaps out of", x, "total")
  
  bar <- ggplot(data = final_table, aes(x= mutation, y=counts)) +
    ggtitle(bar_title) +
    geom_bar(stat="identity", fill = "white", colour = "black")+
    theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
          panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                          size = 2, linetype = "solid"))
  
  
  ##### Statistical analysis part 
  ccl <- distinct(mut, mut$ccl, .keep_all= TRUE )
  
  ccl[match(mut_mutated_gene$ID, ccl$ID), ] <- mut_mutated_gene
  
  ccl$Hugo_Symbol <- as.character(ccl$Hugo_Symbol)
  
  ccl$Hugo_Symbol[ccl$Hugo_Symbol == mutated_gene] <- "has Mut"
  ccl$Hugo_Symbol[ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"
  
  t <- t.test(UBA2_gene_effect ~ Hugo_Symbol, data = ccl) 
  
  box_title <- paste("UBA2 gene effect \n", "T-test p-value =", round(t$p.value, digits = 4))
  
  box <- ggplot(ccl, aes(x= Hugo_Symbol, y=UBA2_gene_effect)) +
    ggtitle(box_title) +
    geom_boxplot(fill = "white", colour = "black", outlier.colour = "red")+
    theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
          panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                          size = 2, linetype = "solid"))
  
  #arrange all 3 plots from (a), (b), and (c)
  # viewport() function controls the margin of the graph. (prevent graphs from getting cut off)
  final_analysis <- grid.arrange(grobs= list(enrich, bar, box), vp=viewport(width=0.8, height=1),
                                 layout_matrix = rbind(c(1, 1), c(1, 1), c(2, 3), c(2, 3), c(2, 3)))
  
  return (final_analysis)
}


#must type in argument as a string (wrapped with "")
mutation_analysis("uterus")



