
########################################################################################################################
#Part 3 - correlation analysis

#dataframe with gene expression
exp_df <- data.frame(ID = exp$X, UBA2_expression_level = exp$UBA2)


#Merge cell line info and gene expression data
final_df <- merge(x = UBA2_merged,
                  y = exp_df,
                  by = "ID")

#eliminate empty data columns and rows
final_df <- droplevels(final_df)


library(ggpmisc)

#example with colorectal cell lines
colorectal_final <- final_df[which(final_df$lineage == "colorectal"), ]

colorectal_plot <- ggplot(colorectal_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("colorectal") + 
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()


blood_final <- final_df[which(final_df$lineage == "blood"), ]

blood_plot <- ggplot(blood_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("blood") +
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

urinary_tract_final <- final_df[which(final_df$lineage == "urinary_tract"), ]

urinary_tract_plot <- ggplot(urinary_tract_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("urinary_tract") +
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

lung_final <- final_df[which(final_df$lineage == "lung"), ]

lung_plot <- ggplot(lung_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("lung") +
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

ovary_final <- final_df[which(final_df$lineage == "ovary"), ]

ovary_plot <- ggplot(ovary_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("ovary") +
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

skin_final <- final_df[which(final_df$lineage == "skin"), ]

skin_plot <- ggplot(skin_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("skin") +
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

breast_final <- final_df[which(final_df$lineage == "breast"), ]

breast_plot <- ggplot(breast_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("breast") +
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

pancreas_final <- final_df[which(final_df$lineage == "pancreas"), ]

pancreas_plot <- ggplot(pancreas_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("pancreas") +
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

plasma_cell_final <- final_df[which(final_df$lineage == "plasma_cell"), ]

plasma_cell_plot <- ggplot(plasma_cell_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("plasma_cell") +
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

central_nervous_system_final <- final_df[which(final_df$lineage == "central_nervous_system"), ]

central_nervous_system_plot <- ggplot(central_nervous_system_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("central_nervous_system") +
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

soft_tissue_final <- final_df[which(final_df$lineage == "soft_tissue"), ]

soft_tissue_plot <- ggplot(soft_tissue_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("soft_tissue") +
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

bone_final <- final_df[which(final_df$lineage == "bone"), ]

bone_plot <- ggplot(bone_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("bone") +
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

gastric_final <- final_df[which(final_df$lineage == "gastric"), ]

formula <- y ~ x

library(ggpmisc)
gastric_plot <- ggplot(gastric_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("gastric") + 
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.15,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+ 
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

lymphocyte_final <- final_df[which(final_df$lineage == "lymphocyte"), ]

lymphocyte_plot <- ggplot(lymphocyte_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("lymphocyte") +
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

peripheral_nervous_system_final <- final_df[which(final_df$lineage == "peripheral_nervous_system"), ]

peripheral_nervous_system_plot <- ggplot(peripheral_nervous_system_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("peripheral_nervous_system") +
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

kidney_final <- final_df[which(final_df$lineage == "kidney"), ]

kidney_plot <- ggplot(kidney_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("kidney") +
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

prostate_final <- final_df[which(final_df$lineage == "prostate"), ]

prostate_plot <- ggplot(prostate_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("prostate") +
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

fibroblast_final <- final_df[which(final_df$lineage == "fibroblast"), ]

fibroblast_plot <- ggplot(fibroblast_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("fibroblast") +
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

bile_duct_final <- final_df[which(final_df$lineage == "bile_duct"), ]

bile_duct_plot <- ggplot(bile_duct_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("bile_duct") +
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

thyroid_final <- final_df[which(final_df$lineage == "thyroid"), ]

thyroid_plot <- ggplot(thyroid_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("thyroid") +
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

upper_aerodigestive_final <- final_df[which(final_df$lineage == "upper_aerodigestive"), ]

upper_aerodigestive_plot <- ggplot(upper_aerodigestive_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("upper_aerodigestive") +
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

liver_final <- final_df[which(final_df$lineage == "liver"), ]

liver_plot <- ggplot(liver_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("liver") +
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

uterus_final <- final_df[which(final_df$lineage == "uterus"), ]

uterus_plot <- ggplot(uterus_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("uterus") +
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

esophagus_final <- final_df[which(final_df$lineage == "esophagus"), ]


esophagus_plot <- ggplot(esophagus_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("esophagus") +
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

cervix_final <- final_df[which(final_df$lineage == "cervix"), ]

cervix_plot <- ggplot(cervix_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("cervix") +
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()

eye_final <- final_df[which(final_df$lineage == "eye"), ]

eye_plot <- ggplot(eye_final, aes(x= UBA2_gene_effect, y =UBA2_expression_level)) + 
  ggtitle("eye") +
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.7,
               formula = formula, parse = TRUE, size = 5)+
  geom_point(color = "red")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold")) +
  theme_bw()



grid.arrange(eye_plot, esophagus_plot, gastric_plot, bone_plot, pancreas_plot, kidney_plot, peripheral_nervous_system_plot,
             skin_plot, plasma_cell_plot, lymphocyte_plot, uterus_plot, thyroid_plot, central_nervous_system_plot,  bile_duct_plot,
             breast_plot, upper_aerodigestive_plot, blood_plot, lung_plot, liver_plot, urinary_tract_plot, colorectal_plot, 
             ovary_plot,  cervix_plot, soft_tissue_plot, prostate_plot, fibroblast_plot
             ,ncol =5, nrow =6)

