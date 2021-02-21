library(devtools)
library(ggplot2)
library("scales")
library(dplyr)
library(grid)
library(gridExtra)
library(forcats) # for ordering x variables in ggplot


#read necessary csv files
setwd("/Users/james/Desktop/Sumoylation_Analysis/project")
gene_effect <- read.csv("C:/Users/james/Desktop/Sumoylation_Analysis/data/Achilles_gene_effect.csv") #CRISPR(Avana) data
cclID <- read.csv("C:/Users/james/Desktop/Sumoylation_Analysis/data/sample_info.csv") #list of all cell lines, ID and lineage
ccl_mut <- read.csv("C:/Users/james/Desktop/Sumoylation_Analysis/data/CCLE_mutations.csv") #import mutation info
exp <- read.csv("C:/Users/james/Desktop/Sumoylation_Analysis/data/CCLE_expression.csv", fileEncoding="UTF-8-BOM") #gene expression data for just protein coding genes

#Enter the gene of interest
Gene <- "UBR1"

#dataframe with the cancer cell line names
ccl_df <- data.frame(ID = cclID$DepMap_ID, ccl = cclID$stripped_cell_line_name, lineage = cclID$lineage)

# delete all the suffixes in col names
names(gene_effect) <- sub("\\..*", "", names(gene_effect))

#dataframe with the gene effect score (example with MSH2 Gene ONLY)
Gene_df <- data.frame(ID = gene_effect$DepMap_ID, KO_effect = gene_effect[, c(Gene)])


#Merge cell line info and gene effect score
Gene_merged <- merge(x = ccl_df, 
                  y = Gene_df,
                  by = "ID")
Gene_merged <- droplevels(Gene_merged)

#Merge with mutation dataset
colnames(ccl_mut)[which(names(ccl_mut) == "DepMap_ID")] <- "ID"
ccl_mut <- subset(ccl_mut, select = c("ID", "Hugo_Symbol", "Variant_Classification"))

#import the auc data to manipulate
gene_effect_mutation <- merge(x = Gene_merged, 
                              y = ccl_mut,
                              by = "ID")

#eliminate silent mutations
gene_effect_mutation <- gene_effect_mutation[which(gene_effect_mutation$Variant_Classification != "Silent"),]

#eliminate duplicated rows
gene_effect_mutation <- unique(gene_effect_mutation)


lineage <- unique(gene_effect_mutation$lineage)
lineage <- sort(lineage)



################################################################################################################
#part 1 - Gene gene effect per lineage. below -1 is a viable threshold

#boxplot with different cell types and their gene effects
png(filename="gene_knockout_analysis/UBR1_KO_effect.png", width=850, bg="white")
par(mar=c(12,5,3,1))
boxplot(KO_effect ~ lineage,
        data=Gene_merged,
        main=paste(Gene, "KO effect by lineage"),
        xlab = "", 
        xaxt = "n", ## delete all x axis labels
        ylab = "Knockout Effect",
        col="red",
        border="black"
)
#axis(1, labels = FALSE)
text(x = 1:length(lineage),
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.05,
     ## Use names from the data list.
     labels = lineage,
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 60,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965,
     cex = 1.2)
dev.off()

#(optional) individual dots added
stripchart(KO_effect ~ lineage, vertical = TRUE, data = Gene_merged, 
           method = "jitter", add = TRUE, pch = 21, col = 'black')


################################################################################################################
#part2 - Enrichment Analysis

#Build a function to repeat the same enrichment analysis with other cell lines
enrichment <- function(cell_line){
  x <- gene_effect_mutation[which(gene_effect_mutation$lineage == cell_line),]
  #separate df with unique KO effect score for each ccl
  x_average_KO <- x %>%
                    group_by(ccl) %>%
                    summarize(KO_effect = mean(KO_effect))
  
  average <- mean(x_average_KO$KO_effect)
  #extract the ccl with KO_effect below average and the max value among them (for drawing the threshold line on the graph)
  significant_ccl <- x_average_KO[x_average_KO$KO_effect <= average, ]
  max_ccl <- significant_ccl$ccl[which.max(significant_ccl$KO_effect)]
  
  enrich_title <- paste(Gene, "knockout effect in", x$lineage ,"cancer cells \n (",
                        length(unique(x$ccl[which(x$KO_effect < average)])),
                        "cell lines out of ",length(unique(x$ccl)),"with KO score <" ,
                        round(average, digits = 3), ")")
  
  enrich <- ggplot(data = x , aes(x = fct_reorder(ccl, KO_effect), y = lineage, width= 1)) +
    ggtitle(enrich_title) +
    geom_tile(aes(fill = KO_effect)) +
    geom_vline(xintercept = max_ccl) + 
    labs(caption = "*Black line indicates the inclusive threshold 
         (all cell lines on the left including the one on the line have KO score below average)") + 
    # limit of the fill colors (color limits at minimum and maximum KO effect values)
    scale_fill_gradient(low = "red", high = "white", 
                        guide = "colorbar", limits=c(min(gene_effect_mutation$KO_effect),max(gene_effect_mutation$KO_effect))) +
    theme(plot.title = element_text(size=12, hjust = 0.5, face = "bold"), 
          panel.background = element_rect(fill = "white", colour = "lightgrey", 
                                          size = 1, linetype = "solid"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), 
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          #margin(top, right, bottom, left)
          plot.margin = unit(c(1,1,1,1), "cm")) 
  
  return (enrich)
}

enrichment("colorectal")


################################################################################################################
#(b) Common mutation


calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}


### IMPORTANT 
### the variable "mutated_gene" can be set to any gene. 
### Refer to the list printed by the function and set the gene accordingly

##create function to repeat the process (b) and (c)  
mutation_analysis <- function(cell_line, index){
  enrich <- enrichment(cell_line)
  mut <- gene_effect_mutation[which(gene_effect_mutation$lineage == cell_line), ]
  average <- mean(mut$KO_effect)
  
  #subset of target cell lines with below threshold gene effect score (average)
  valid_mut <- mut[which(mut$KO_effect <= average), ]
  valid_mut <- droplevels(valid_mut)
  #subset of cell lines with above threshold
  invalid_mut <- mut[which(mut$KO_effect > average), ]
  invalid_mut <- droplevels(invalid_mut)
  
  #eliminate mutations that also occur in invalid
  #colorectal_valid_mut1 <- valid_mut[!duplicated(invalid_mut[c('Hugo_Symbol', 'Variant_Classification')]),] 
  mut_final <- anti_join(valid_mut, invalid_mut, 
                                    by=c("Hugo_Symbol", "Variant_Classification"), copy = TRUE)
  
  frequent_mut <- sort(table(mut_final$Hugo_Symbol),decreasing=TRUE)[1:10]
  print(frequent_mut)
  #calculate the mode (most frequent mutation)
  #mutated_gene <- calculate_mode(mut_final$Hugo_Symbol)
  mutated_gene <- names(frequent_mut)[index]
  #info with the mutated_gene above
  mut_mutated_gene <- mut_final[which(mut_final$Hugo_Symbol == mutated_gene), ]
 
  #find the most frequent variant type
  variant <- calculate_mode(mut_mutated_gene$Variant_Classification)
  
  #create subset with the specific variant from above
  mut_mutated_gene_var <- mut_mutated_gene[which(mut_mutated_gene$Variant_Classification == variant), ]
  
  mut_mutated_gene <- droplevels(mut_mutated_gene)

  y <- length(mut_mutated_gene_var$ccl)
  x <- length(unique(mut$ccl))
  x <- x - y
  
  final_table <- data.frame(mutation = c(paste(as.character(mutated_gene),"mutation"), "other"), counts = c(y, x))
  
  bar_title <- paste("Mut:", as.character(mutated_gene), as.character(variant),
                     "\n", y, "overlaps out of", x, "total")
  
  bar <- ggplot(data = final_table, aes(x= mutation, y=counts)) +
    ggtitle(bar_title) +
    geom_bar(stat="identity", fill = "white", colour = "black")+
    theme(plot.title = element_text(size=10, hjust = 0.5, face = "bold"), 
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, hjust=0.8),
          panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                          size = 2, linetype = "solid"),
          #margin(top, right, bottom, left)
          plot.margin=unit(c(1,3,1,4),"cm"))

  
  ##### Statistical analysis part 
  ccl <- distinct(mut, mut$ccl, .keep_all= TRUE )
  
  ccl[match(mut_mutated_gene$ID, ccl$ID), ] <- mut_mutated_gene
  
  ccl$Hugo_Symbol <- as.character(ccl$Hugo_Symbol)
  
  ccl$Hugo_Symbol[ccl$Hugo_Symbol == mutated_gene] <- "has Mut"
  ccl$Hugo_Symbol[ccl$Hugo_Symbol != "has Mut"] <- "lacks Mut"
  
  t <- t.test(KO_effect ~ Hugo_Symbol, data = ccl) 
  
  box_title <- paste("Gene gene effect \n", "T-test p-value =", round(t$p.value, digits = 5))
  
  box <- ggplot(ccl, aes(x= Hugo_Symbol, y=KO_effect)) +
    ggtitle(box_title) +
    geom_boxplot(fill = "white", colour = "black", outlier.colour = "red")+
    theme(plot.title = element_text(size=10, hjust = 0.5, face = "bold"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, hjust=0.8),
          panel.background = element_rect(fill = "lightgrey", colour = "lightgrey",
                                          size = 2, linetype = "solid"), 
          #margin(top, right, bottom, left)
          plot.margin=unit(c(1,4,1,3),"cm"))
  
  #arrange all 3 plots from (a), (b), and (c)
  # viewport() function controls the margin of the graph. (prevent graphs from getting cut off)
  final_analysis <- grid.arrange(grobs= list(enrich, bar, box), 
                                 layout_matrix = rbind(c(1, 1), 
                                                       c(1, 1),
                                                       c(1, 1),
                                                       c(2, 3), 
                                                       c(2, 3), 
                                                       c(2, 3)))
  
  return (final_analysis)
}


#must type in argument as a string (wrapped with "")
# second parameter indicates which common mutation to analyze. 1 being the most common mutation, 2 being the second most common, and so on. 
mutation_analysis("pancreas", 3)



