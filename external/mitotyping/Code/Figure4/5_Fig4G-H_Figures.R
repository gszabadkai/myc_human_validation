rm(list = ls())

library(tidyverse)
library(ComplexHeatmap)
library(ggfortify)
library(fmsb)
library(broom)
# Setup -------------------------------------------------------------------

color_conditions = c("Control" ="#929292", 
                     "DEX"= "#E7298A", 
                     "Oligomycin"= "#1E78B4", 
                     "Oligomycin_DEX"= "#A5CEE3", 
                     "SURF1"= "#E3211D", 
                     "SURF1_DEX"= "#FC9A99", 
                     "MitoNUITs"="#33A12C", 
                     "MitoNUITs_DEX"= "#B2DF8A", 
                     "Galactose" = "#D96004", 
                     "2DG" = "#FFDA2E", 
                     "BHB"= "#B05928",
                     "Hypoxia"= "#F681C0",
                     "Contact_Inhibition" = "#6A3E9A",  
                     "Contact_Inhibition_Hypoxia"= "#CBB2D5")

# Figure 6H Heatmap -------------------------------------------------------

mitoPPS <- read_csv(here::here("Data", "Fibroblasts", 
                                        "ProcessedData", "mitoPPS_All.csv"))
data_heatmap <- mitoPPS %>%
  column_to_rownames("RNAseq_sampleID") %>%
  arrange(fct_relevel(Condition, c("Control",
                                   
                                   "Hypoxia","Contact_Inhibition", "Contact_Inhibition_Hypoxia",
                                   
                                   "2DG", "Galactose", "BHB",
                                   "SURF1","Oligomycin","MitoNUITs", 
                                   "DEX", "MitoNUITs_DEX", "Oligomycin_DEX", "SURF1_DEX")), Passage)

exprs <- log10(data_heatmap[,8:ncol(data_heatmap)])
anno_df <- data_heatmap[,1:7]

#Elbow Method for finding the optimal number of clusters
set.seed(123)
# Compute and plot wss for k = 2 to k = 15.
k.max <- 15
wss <- sapply(1:k.max, 
              function(k){kmeans(exprs, k, nstart=50,iter.max = 15 )$tot.withinss})
wss
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")


col_anno <- columnAnnotation(Passage = anno_barplot(anno_df$Passage),
                             Condition = anno_df$Condition,
                             col=list(#Treatment  = color_conditions,
                               Condition = color_conditions
                             ))

set.seed(123)
hm <- ComplexHeatmap::Heatmap(t(exprs),
                              column_title = "",
                              column_title_side = c("top"),
                              cluster_columns = F,
                              cluster_rows = T,
                              top_annotation = col_anno,
                              show_column_names = F,
                              row_names_gp = gpar(fontsize = 6), row_km = 4)

hm <- draw(hm)


row_order <- reshape2::melt(row_order(hm)) %>%
  dplyr::rename(Position = value, Cluster = L1) %>%
  dplyr::mutate(row_order_in_heatmap = seq_len(149))
clusters <- data.frame(Pathway = colnames(exprs), Position = seq(1:149)) %>%
  inner_join(row_order) %>%
  arrange()

hm

pdf(here::here("Figures","Figure4","Fig4G_Heatmap_mitoPPS_noLabels_small.pdf"), height = 6, width = 10, "in")
set.seed(123)
hm <- ComplexHeatmap::Heatmap(t(exprs),
                              column_title = "",
                              column_title_side = c("top"),
                              cluster_columns = F,
                              cluster_rows = T,
                              top_annotation = col_anno,
                              show_column_names = F,
                              show_row_names = F, row_km = 4)
draw(hm)
dev.off()



# Suppl FigS9 ------------------------------------------------------------


## B Heatmap ---------------------------------------------------------------

pdf(here::here("Figures","Figure4","Suppl_FigS9B_Heatmap_mitoPPS_withLabels_large.pdf"), height = 12, width = 10, "in")
set.seed(123)
hm <- ComplexHeatmap::Heatmap(t(exprs),
                              column_title = "",
                              column_title_side = c("top"),
                              cluster_columns = F,
                              cluster_rows = T,
                              top_annotation = col_anno,
                              show_column_names = F,
                              show_row_names = T, row_km = 4,
                              row_names_gp = gpar(fontsize = 6))

draw(hm)
dev.off()

## C PCA - combine datasets ---------------------------------------------------

mitocarta_sheet4 <- readxl::read_xls(here::here("Data","MitoCarta", "OriginalData","HumanMitoCarta3_0.xls"), sheet = 4) %>%
  dplyr::select(MitoPathway, Genes) %>%
  na.omit()
gene_to_pathway <- splitstackshape::cSplit(mitocarta_sheet4, 'Genes', ',') %>%
  column_to_rownames("MitoPathway") %>%
  t() %>%
  as.data.frame()
gene_to_pathway <- gene_to_pathway %>%
  pivot_longer(cols = colnames(gene_to_pathway), names_to = "Pathway", values_to = "Gene") %>%
  na.omit()

all_pathways <- unique(gene_to_pathway$Pathway)
mitoPPS <- read_csv(here::here("Data", "Fibroblasts", 
                                        "ProcessedData", "mitoPPS_All.csv"))
tissue_mitoPPS <- read_csv(here::here("Data", "HumanProteinAtlas", "ProcessedData", "mitoPPS_all_tissues.csv"))

combined <- mitoPPS %>%
  mutate(Group = "Fibroblasts") %>%
  mutate(Tissue = "Control") %>%
  #  inner_join(meta) %>%
  select(RNAseq_sampleID, Line, Study_part, Passage, Days_grown_Udays, Tissue, Group, all_of(all_pathways)) %>%
  pivot_longer(cols = -c(RNAseq_sampleID, Line, Study_part, Passage, Days_grown_Udays, Tissue, Group), names_to = "Pathway") %>%
  full_join(tissue_mitoPPS %>%pivot_longer(cols = -c(Tissue, Group) , names_to = "Pathway")%>%
              mutate(Passage =NA) ) %>%
  pivot_wider(names_from = "Pathway", values_from = "value") %>%
  arrange(desc(Group))

theme_pca <- theme(axis.text = element_text(size = 6),
                   axis.title = element_text(size = 7),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.position = "none"
)


pca <- prcomp(combined[, -c(1:7)], scale = T)
loadings <- pca$rotation
write_csv(as.data.frame(loadings) %>% rownames_to_column("Pathway"), 
          here::here("Supplemental/Suppl_Fig8C_loadings.csv"))

summary(pca)
autoplot(pca, data = combined , colour = 'Passage', alpha = 0.6, scale = 1,
         fill = 'Passage', size = 1, shape = 21, lwd = 0.5, )+
  viridis::scale_color_viridis(direction = -1, na.value = "grey") +
  viridis::scale_fill_viridis(direction = -1, na.value = "grey") +
  theme_minimal() +
  scale_x_reverse() +
  scale_y_reverse() +
  theme_pca## Note the switch in sign in PCA, axis reversed for visualization purposes. 

ggsave(here::here("Figures", "Figure4", "Suppl_FigS9C_PCA_Tissue_Fibroblasts_AllTreatments_AllPassages.png"), height = 1.35, width = 1.5,unit = "in", dpi = 1200)


pca$x <- pca$x[which(combined$Group == "Fibroblasts"), ]
autoplot(pca, data = combined %>%
           filter(Group %in% "Fibroblasts"), colour = 'Passage', alpha = 0.6, scale = 1,
         fill = 'Passage', size = 2, shape = 21, lwd = 0.5)+
  viridis::scale_color_viridis(direction = -1) +
  viridis::scale_fill_viridis(direction = -1) +
  theme_minimal() +
  theme_pca + 
  scale_x_reverse()+
  scale_y_reverse() 
ggsave(here::here("Figures", "Figure4", "Suppl_FigS9C_PCA_Tissue_Fibroblasts_AllTreatments_AllPassages_subset.png"), 
       height = 1.35, width = 1.5,unit = "in", dpi = 1200)

color_conditions = c("Control" ="#929292", 
                     "Tissue" = "gray",
                     "DEX"= "#E7298A", 
                     "Oligomycin"= "#1E78B4", 
                     "Oligomycin_DEX"= "#A5CEE3", 
                     "SURF1"= "#E3211D", 
                     "SURF1_DEX"= "#FC9A99", 
                     "MitoNUITs"="#33A12C", 
                     "MitoNUITs_DEX"= "#B2DF8A", 
                     "Galactose" = "#D96004", 
                     "2DG" = "#FFDA2E", 
                     "BHB"= "#B05928",
                     "Hypoxia"= "#F681C0",
                     "Contact_Inhibition" = "#6A3E9A",  
                     "Contact_Inhibition_Hypoxia"= "#CBB2D5")



## D PCA All passages averaged ------------------------------------------------

combined <- mitoPPS %>%
  mutate(Group = "Fibroblasts") %>%
  mutate(Tissue = "Control") %>%
  select(RNAseq_sampleID, Line, Study_part, Passage, Days_grown_Udays, 
         Tissue, Group,Condition, Treatments, Days_grown_Udays, all_of(all_pathways)) %>%
  pivot_longer(cols = -c(RNAseq_sampleID, Line, Study_part, 
                         Passage, Days_grown_Udays, Tissue, 
                         Group, Condition, Treatments, Days_grown_Udays), names_to = "Pathway") %>%
  unique() %>%
  group_by(Line, Condition, Treatments, Pathway ) %>%
  mutate(mean = mean(value)) %>% 
  select(-c(value, Passage, Days_grown_Udays, RNAseq_sampleID, Study_part)) %>%
  unique() %>%
  full_join(tissue_mitoPPS %>%pivot_longer(cols = -c(Tissue, Group) , names_to = "Pathway", values_to = "mean")%>%
              mutate(Passage =NA) %>%
              mutate(Condition = "Tissue")) %>%
  pivot_wider(names_from = "Pathway", values_from = "mean") %>%
  arrange(desc(Group))

pca <- prcomp(combined[, -c(1:7)], scale = T)
summary(pca)
loadings <- pca$rotation
write_csv(as.data.frame(loadings) %>% rownames_to_column("Pathway"), 
          here::here("Supplemental/Suppl_Fig8D_loadings.csv"))

autoplot(pca, data = combined , colour = 'Condition', alpha = 0.6, scale = 1,
         fill = 'Condition', size = 1, shape = 21, lwd = 0.5, )+
  scale_color_manual(values = color_conditions  ) +
  scale_fill_manual(values = color_conditions  ) +
  theme_minimal() +
  theme_pca +
  scale_y_reverse()
ggsave(here::here("Figures", "Figure4", "Suppl_FigS9D_PCA_Tissue_Fibroblasts_AllTreatments_AllPassagesAverage.png"),
       height = 1.35, width = 1.5,unit = "in", dpi = 1200)


pca$x <- pca$x[which(combined$Group == "Fibroblasts"), ]
autoplot(pca, data = combined %>%
           filter(Group %in% "Fibroblasts"), colour = 'Condition', alpha = 0.6, scale = 1,
         fill = 'Condition', size = 2, shape = 21, lwd = 0.5, )+
  scale_color_manual(values = color_conditions  ) + 
  scale_fill_manual(values = color_conditions  ) + 
  theme_minimal() +
  theme_pca + 
  scale_y_reverse() 
ggsave(here::here("Figures", "Figure4", "Suppl_FigS9D_PCA_Tissue_Fibroblasts_AllTreatments_AllPassagesAveraged_subset.png"), 
       height = 1.35, width = 1.5,unit = "in", dpi = 1200)


## PC3 vs 4
autoplot(pca, x = 3, y = 4, data = combined %>%
           filter(Group %in% "Fibroblasts"), colour = 'Condition', alpha = 0.6, scale = 1,
         fill = 'Condition', size = 2, shape = 21, lwd = 0.5, )+
  scale_color_manual(values = color_conditions  ) + 
  scale_fill_manual(values = color_conditions  ) + 
  theme_minimal() +
  theme_pca + 
  scale_y_reverse() 

## PC3 vs 4
autoplot(pca, x = 5, y = 6, data = combined %>%
           filter(Group %in% "Fibroblasts"), colour = 'Condition', alpha = 0.6, scale = 1,
         fill = 'Condition', size = 2, shape = 21, lwd = 0.5, )+
  scale_color_manual(values = color_conditions  ) + 
  scale_fill_manual(values = color_conditions  ) + 
  theme_minimal() +
  theme_pca + 
  scale_y_reverse() 


## E 3D plot -----------------------------------------------------------------
#pca <- prcomp(combined[, -c(1:7)], scale = T)

summary_pca <- summary(pca)
color_conditions = c("Control" ="#929292", 
                     #"Tissue" = "gray",
                     "DEX"= "#E7298A", 
                     "Oligomycin"= "#1E78B4", 
                     "Oligomycin_DEX"= "#A5CEE3", 
                     "SURF1"= "#E3211D", 
                     "SURF1_DEX"= "#FC9A99", 
                     "MitoNUITs"="#33A12C", 
                     "MitoNUITs_DEX"= "#B2DF8A", 
                     "Galactose" = "#D96004", 
                     "2DG" = "#FFDA2E", 
                     "BHB"= "#B05928",
                     "Hypoxia"= "#F681C0",
                     "Contact_Inhibition" = "#6A3E9A",  
                     "Contact_Inhibition_Hypoxia"= "#CBB2D5")

## Commented out because it requires developer tools to be installed
# var_1 <- round(summary_pca$importance[2,1]*100,2)
# var_2 <- round(summary_pca$importance[2,2]*100,2)
# var_3 <- round(summary_pca$importance[2,3]*100,2)
# group_color_df <- data.frame(Color = color_conditions) %>%
#   rownames_to_column("Condition")
# df <- pca$x
# conditions <- combined %>% filter(!Condition %in% "Tissue") %>%
#   pull(Condition)
# df <- data.frame(PC1=df[,1], PC2=df[,2], PC3=df[,3], Condition = as.factor(conditions)) %>%
#   full_join(group_color_df, by = "Condition") %>%
#   na.omit() 
# 
# rgl::par3d(cex=0.001)
# rgl::bg3d(color = "white")
# with(df, rgl::plot3d(PC1,PC2,PC3, col= Color, alpha = 0.6, size = 15, type = "p", 
#                      xlab = "", ylab = "", zlab = "",
#                      ann = FALSE, axes = FALSE
# ))
# rgl::box3d()


# Suppl FigS10 Radar charts ---------------------------------------------------

## Fibroblast data
fibs <- mitoPPS %>%
  filter(!Condition %in% "Control") %>%
  select(Condition, all_of(all_pathways)) %>%
  pivot_longer(cols = -c(Condition), names_to = "Pathway") %>%
  unique() %>%
  group_by(Condition, Pathway ) %>%
  mutate(mean = mean(value)) %>% 
  select(-value) %>%
  unique() 
## Axis label with percent changes
label_fibs <- fibs %>%
  group_by(Pathway) %>%
  mutate(mean_all = mean(mean)) %>%
  mutate(pct_change = round((mean - mean_all) / mean_all * 100)) %>% 
  mutate(pct = "%") %>%
  mutate(sign = case_when(
    pct_change <0 ~ "",
    TRUE ~ "+"
  )) %>%
  unite(label, sign, pct_change, pct, sep = "") %>%
  unite(label, Condition, label, sep = " ", remove= FALSE) %>%
  select(Condition, label, Pathway)

## Tissue data
tissues <- tissue_mitoPPS %>%
  pivot_longer(cols = -c(Tissue, Group) , names_to = "Pathway", values_to = "mitoPPS")%>%
  filter(Tissue %in% c("Liver", "Cerebral cortex", "Skeletal muscle", "Bone marrow", "Small intestine")) 
## Axis label with percent changes
label_tissues <- tissue_mitoPPS %>%
  pivot_longer(cols = -c(Tissue, Group) , names_to = "Pathway", values_to = "mitoPPS") %>%
  group_by(Pathway) %>%
  mutate(mean_all = mean(mitoPPS)) %>%
  mutate(pct_change = round((mitoPPS - mean_all) / mean_all * 100)) %>% 
  mutate(pct = "%") %>%
  mutate(sign = case_when(
    pct_change <0 ~ "",
    TRUE ~ "+"
  )) %>%
  unite(label, sign, pct_change, pct, sep = "") %>%
  unite(label, Tissue, label, sep = " ", remove= FALSE) %>%
  select(Tissue, label, Pathway) %>%
  filter(Tissue %in% c("Liver", "Cerebral cortex", "Skeletal muscle", "Bone marrow", "Small intestine")) 


## Complex I ---------------------------------------------------------------
tissues_sub <- tissues %>%
  full_join(label_tissues) %>%
  ungroup() %>%
  select(-c(Group, Tissue)) %>%
  dplyr::rename(Tissue = label) %>%
  filter(Pathway %in% "Complex I") 
fibs_sub <- fibs %>%
  filter(Pathway %in% "Complex I") #%>%
fibs_average <- mean(fibs_sub$mean)
tissues_average <- tissue_mitoPPS %>%
  pivot_longer(cols = -c(Tissue, Group) , names_to = "Pathway", values_to = "mitoPPS") %>%
  group_by(Pathway) %>%
  mutate(Grand_mean = mean(mitoPPS)) %>%
  filter(Pathway %in% "Complex I") %>%
  pull(Grand_mean) %>%
  unique()
min <- min(tissues_sub$mitoPPS, fibs_sub$mean)
min <- min + min/100
max <- max(tissues_sub$mitoPPS, fibs_sub$mean)
max <- max + max/100
data <- fibs_sub %>%
  pivot_wider(names_from = "Condition", values_from = "mean") %>%
  ungroup() %>%
  select(-Pathway) %>%
  select(Hypoxia, SURF1_DEX, Oligomycin_DEX,
         MitoNUITs_DEX, DEX, Oligomycin, MitoNUITs, 
         SURF1, BHB, Galactose, `2DG`, Contact_Inhibition_Hypoxia,
         Contact_Inhibition, Hypoxia) %>%
  as.data.frame()
new_colnames <- as.data.frame(colnames(data)) %>%
  dplyr::rename(Condition = `colnames(data)`) %>%
  full_join(label_fibs %>%
              filter(Pathway %in% "Complex I"))
colnames(data) <- new_colnames$label

data <- rbind(rep(max,13) , rep(min,13), rep(fibs_average, 13),data) 
rownames(data) <- c("Max", "Min", "Average", "Fibroblasts")
colors_border=c( 
  rgb(210/255, 210/255, 210/255,0.9), 
  rgb(73/255, 116/255, 159/255,0.9) 
  )
colors_in=c( 
  rgb(210/255, 210/255, 210/255,0.4), 
  rgb(73/255, 116/255, 159/255,0.4) 
) 
png("Figures/Figure4/Suppl_FigS10A_radar_Fibs_CI.png", width = 6.4, height = 6, units = "in", res = 1200)
radarchart(data  , axistype=0, pcol = colors_border, pfcol = colors_in,
           plwd=c(0.5,1.5),  pty=16,plty=c(1,3),
           cglcol="grey",cglwd = 1, cglty=3, axislabcol="grey",
           vlcex=0.8, title="Complex I in Fibroblasts", calcex=0.8, seg = 3
) 
dev.off()

# Tissues
data <- tissues_sub %>%
  pivot_wider(names_from = "Tissue", values_from = "mitoPPS") %>%
  ungroup() %>%
  select(-Pathway) %>%
  as.data.frame() 
data <- rbind(rep(max,5) , rep(min,5), rep(tissues_average, 5), data ) 
rownames(data) <- c("Max", "Min", "Average", "Tissues")
colors_border=c(
  rgb(210/255, 210/255, 210/255,0.9), 
  rgb(0,0,0,0.6) 
)
colors_in=c(
  rgb(210/255, 210/255, 210/255,0.4), 
  rgb(0,0,0,0.2) 
) 
png("Figures/Figure4/Suppl_FigS10A_radar_Tissues_CI.png", width = 6.4, height = 6, units = "in", res = 1200)

radarchart(data  , axistype=0, pcol = colors_border, pfcol = colors_in,
           plwd=c(0.5,1.5),  pty=16,plty=c(1,3),
           cglcol="grey",cglwd = 1, cglty=3, axislabcol="grey",
           vlcex=0.8, title="Complex I in Tissues", calcex=0.8, seg = 3
) 
dev.off()

## FAO ---------------------------------------------------------------
tissues_sub <- tissues %>%
  full_join(label_tissues) %>%
  ungroup() %>%
  select(-c(Group, Tissue)) %>%
  dplyr::rename(Tissue = label) %>%
  filter(Pathway %in% "Fatty acid oxidation") 
fibs_sub <- fibs %>%
  filter(Pathway %in% "Fatty acid oxidation") #%>%
fibs_average <- mean(fibs_sub$mean)
tissues_average <- tissue_mitoPPS %>%
  pivot_longer(cols = -c(Tissue, Group) , names_to = "Pathway", values_to = "mitoPPS") %>%
  group_by(Pathway) %>%
  mutate(Grand_mean = mean(mitoPPS)) %>%
  filter(Pathway %in% "Fatty acid oxidation") %>%
  pull(Grand_mean) %>%
  unique()
min <- min(tissues_sub$mitoPPS, fibs_sub$mean)
min <- min + min/100
max <- max(tissues_sub$mitoPPS, fibs_sub$mean)
max <- max + max/100
data <- fibs_sub %>%
  pivot_wider(names_from = "Condition", values_from = "mean") %>%
  ungroup() %>%
  select(-Pathway) %>%
  select(Hypoxia, SURF1_DEX, Oligomycin_DEX,
         MitoNUITs_DEX, DEX, Oligomycin, MitoNUITs, 
         SURF1, BHB, Galactose, `2DG`, Contact_Inhibition_Hypoxia,
         Contact_Inhibition, Hypoxia) %>%
  as.data.frame()
new_colnames <- as.data.frame(colnames(data)) %>%
  dplyr::rename(Condition = `colnames(data)`) %>%
  full_join(label_fibs %>%
              filter(Pathway %in% "Fatty acid oxidation"))
colnames(data) <- new_colnames$label

data <- rbind(rep(max,13) , rep(min,13), rep(fibs_average, 13),data) 
rownames(data) <- c("Max", "Min", "Average", "Fibroblasts")
colors_border=c( 
  rgb(210/255, 210/255, 210/255,0.9), 
  rgb(73/255, 116/255, 159/255,0.9) 
)
colors_in=c( 
  rgb(210/255, 210/255, 210/255,0.4), 
  rgb(73/255, 116/255, 159/255,0.4) 
) 
png("Figures/Figure4/Suppl_FigS10B_radar_Fibs_FAO.png", width = 6.4, height = 6, units = "in", res = 1200)
radarchart(data  , axistype=0, pcol = colors_border, pfcol = colors_in,
           plwd=c(0.5,1.5),  pty=16,plty=c(1,3),
           cglcol="grey",cglwd = 1, cglty=3, axislabcol="grey",
           vlcex=0.8, title="Fatty acid oxidation in Fibroblasts", calcex=0.8, seg = 3
) 
dev.off()

# Tissues
data <- tissues_sub %>%
  pivot_wider(names_from = "Tissue", values_from = "mitoPPS") %>%
  ungroup() %>%
  select(-Pathway) %>%
  as.data.frame() 
data <- rbind(rep(max,5) , rep(min,5), rep(tissues_average, 5), data ) 
rownames(data) <- c("Max", "Min", "Average", "Tissues")
colors_border=c(
  rgb(210/255, 210/255, 210/255,0.9), 
  rgb(0,0,0,0.6) 
)
colors_in=c(
  rgb(210/255, 210/255, 210/255,0.4), 
  rgb(0,0,0,0.2) 
) 
png("Figures/Figure4/Suppl_FigS10B_radar_Tissues_FAO.png", width = 6.4, height = 6, units = "in", res = 1200)

radarchart(data  , axistype=0, pcol = colors_border, pfcol = colors_in,
           plwd=c(0.5,1.5),  pty=16,plty=c(1,3),
           cglcol="grey",cglwd = 1, cglty=3, axislabcol="grey",
           vlcex=0.8, title="Fatty acid oxidation in Tissues", calcex=0.8, seg = 3
) 
dev.off()


## mtDNA repair ---------------------------------------------------------------
## Fibroblast data only controls!!

fibs <- mitoPPS %>%
  filter(Condition %in% "Control") %>%
  filter(Study_part %in% 2) %>%
  select(Line, Passage, all_of(all_pathways)) %>%
  pivot_longer(cols = -c(Line, Passage), names_to = "Pathway") %>%
  unique() %>%
  filter(Pathway %in% "mtDNA repair") 

## scatterplot
fibs %>% 
  ggplot(aes(x = Passage, y = value)) +
  geom_point(color ="darkmagenta") +
  facet_wrap(~Line) +
  theme_bw() +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        strip.text = element_text(size = 8))
ggsave("Figures/Figure4/Suppl_FigS10C_mtDNA_repair_scatter.png", width = 3.7, height = 2, units = "in", dpi = 1200)


label_tissues <- tissue_mitoPPS %>%
  pivot_longer(cols = -c(Tissue, Group) , names_to = "Pathway", values_to = "mitoPPS") %>%
  group_by(Pathway) %>%
  mutate(mean_all = mean(mitoPPS)) %>%
  mutate(pct_change = round((mitoPPS - mean_all) / mean_all * 100)) %>% 
  mutate(pct = "%") %>%
  mutate(sign = case_when(
    pct_change <0 ~ "",
    TRUE ~ "+"
  )) %>%
  unite(label, sign, pct_change, pct, sep = "") %>%
  unite(label, Tissue, label, sep = " ", remove= FALSE) %>%
  select(Tissue, label, Pathway) %>%
  filter(Tissue %in% c("Liver", "Cerebral cortex", "Skeletal muscle", "Bone marrow", "Small intestine")) 


### Tissue radar
tissues_sub <- tissues %>%
  full_join(label_tissues) %>%
  ungroup() %>%
  select(-c(Group, Tissue)) %>%
  dplyr::rename(Tissue = label) %>%
  filter(Pathway %in% "mtDNA repair") 

# Tissues
data <- tissues_sub %>%
  pivot_wider(names_from = "Tissue", values_from = "mitoPPS") %>%
  ungroup() %>%
  select(-Pathway) %>%
  as.data.frame() 
data <- rbind(rep(max,5) , rep(min,5), rep(tissues_average, 5), data ) 
rownames(data) <- c("Max", "Min", "Average", "Tissues")
colors_border=c(
  rgb(210/255, 210/255, 210/255,0.9), 
  rgb(0,0,0,0.6) 
)
colors_in=c(
  rgb(210/255, 210/255, 210/255,0.4), 
  rgb(0,0,0,0.2) 
) 
png("Figures/Figure4/Suppl_FigS10C_radar_Tissues_mtDNA.png", width = 6.4, height = 6, units = "in", res = 1200)

radarchart(data  , axistype=0, pcol = colors_border, pfcol = colors_in,
           plwd=c(0.5,1.5),  pty=16,plty=c(1,3),
           cglcol="grey",cglwd = 1, cglty=3, axislabcol="grey",
           vlcex=0.8, title="mtDNA repair in Tissues", calcex=0.8, seg = 3
) 
dev.off()





# Suppl FigS11 Dynamic range -----------------------------------------------------------
mitoPPS <- read_csv(here::here("Data", "Fibroblasts", 
                                        "ProcessedData", "mitoPPS_All.csv")) %>%
  pivot_longer(cols = -c(RNAseq_sampleID, Line,  Study_part, Passage,
                         Days_grown_Udays, Treatments, Clinical_condition, Condition), 
               names_to = "Pathway", values_to = "mitoPPS") 

mitoPPS_Tissues <- read_csv(here::here("Data", "HumanProteinAtlas","ProcessedData", "mitoPPS_all_tissues.csv" )) %>%
  pivot_longer(cols = -c(Tissue, Group), names_to = "Pathway", values_to = "mitoPPS")

delta_controls <- mitoPPS %>%
  filter(Condition %in% "Control") %>%
  filter(Study_part %in% 2) %>%
  group_by(Pathway) %>%
  mutate(delta_controls = max(mitoPPS) - min(mitoPPS)) %>%
  select(Pathway, delta_controls) %>%
  unique()

delta_treatments <- mitoPPS %>%
  filter(!Condition %in% "Control") %>%
  group_by(Condition, Pathway, Line) %>%
  mutate(mitoPPS = mean(mitoPPS)) %>%
  select(Condition, Pathway, Line, mitoPPS) %>%
  unique() %>%
  group_by(Pathway) %>%
  mutate(delta_treatments = max(mitoPPS) - min(mitoPPS)) %>%
  select(Pathway, delta_treatments) %>%
  unique()

delta_tissues <- mitoPPS_Tissues %>%
  group_by(Pathway) %>%
  mutate(delta_tissues = max(mitoPPS) - min(mitoPPS)) %>%
  select(Pathway, delta_tissues) %>%
  unique()

combined <- delta_controls %>%
  full_join(delta_treatments) %>%
  full_join(delta_tissues) %>%
  dplyr::rename(`Treated Fibroblasts` = delta_treatments, Tissues = delta_tissues, `Control Fibroblasts` =delta_controls ) %>%
  pivot_longer(cols = -Pathway, names_to = "dataset", values_to = "delta") %>%
  unite(Group, Pathway, dataset, remove = F) %>%
  mutate(dataset = as.factor(dataset))
combined$dataset <- factor(combined$dataset, levels = c("Tissues", "Control Fibroblasts", "Treated Fibroblasts"))
order <- combined %>%
  ungroup() %>%
  filter(dataset %in% "Tissues") %>%
  arrange(desc(delta)) %>%
  mutate(order = seq(1:149)) %>%
  select(Pathway, order)

data1 <- combined %>%
  inner_join(order) %>%
  arrange(desc(order)) %>%
  mutate(Pathway = as.factor(Pathway)) %>%
  arrange(desc(order)) %>% 
  filter(dataset %in% "Tissues")
data2 <- combined %>%
  inner_join(order) %>%
  arrange(desc(order)) %>%
  mutate(Pathway = as.factor(Pathway)) %>%
  arrange(desc(order)) %>% 
  filter(dataset %in% "Treated Fibroblasts")
data3 <- combined %>%
  inner_join(order) %>%
  arrange(desc(order)) %>%
  mutate(Pathway = as.factor(Pathway)) %>%
  arrange(desc(order)) %>% 
  filter(dataset %in% "Control Fibroblasts")

ggplot() +
  geom_bar(data = data1, aes(x = Pathway, y = delta), stat = "identity", fill = "grey", alpha =0.7) + 
  geom_bar(data = data2, aes(x = Pathway, y = delta), stat = "identity", fill = "dodgerblue4", alpha =0.7) + # must include argument label "data"
  geom_bar(data = data3, aes(x = Pathway, y = delta), stat = "identity", fill = "darkmagenta", alpha =0.7) + 
  theme_classic() +
  scale_x_discrete(limits=data1$Pathway) +
  coord_flip() +
  labs(x = "MitoPathways",y = "Dynamic range of mitoPPS") + 
  theme(axis.text.y= element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 7))
ggsave(here::here("Figures", "Figure4", "Fig4H_dynamic_range.png"), 
       height = 5, width = 1.81,unit = "in", dpi = 1200)

ggplot() +
  geom_bar(data = data1, aes(x = Pathway, y = delta), stat = "identity", fill = "grey", alpha =0.7) + 
  geom_bar(data = data2, aes(x = Pathway, y = delta), stat = "identity", fill = "dodgerblue4", alpha =0.7) + # must include argument label "data"
  geom_bar(data = data3, aes(x = Pathway, y = delta), stat = "identity", fill = "darkmagenta", alpha =0.7) + 
  theme_bw() +
  scale_x_discrete(limits=data1$Pathway) +
  coord_flip() +
  labs(x = "MitoPathways",y = "Dynamic range of mitoPPS") + 
  theme(axis.text= element_text(size =5),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 5))
ggsave(here::here("Figures", "Figure4", "Suppl_FigS11A_dynamic_range_large.png"), 
       height = 10, width = 3.5,unit = "in", dpi = 1200)


## Bivariate mitoPPS ---------------------------------------------------------
color_conditions = c("Control" ="#929292", 
                     "DEX"= "#E7298A", 
                     "Oligomycin"= "#1E78B4", 
                     "Oligomycin_DEX"= "#A5CEE3", 
                     "SURF1"= "#E3211D", 
                     "SURF1_DEX"= "#FC9A99", 
                     "MitoNUITs"="#33A12C", 
                     "MitoNUITs_DEX"= "#B2DF8A", 
                     "Galactose" = "#D96004", 
                     "2DG" = "#FFDA2E", 
                     "BHB"= "#B05928",
                     "Hypoxia"= "#F681C0",
                     "Contact_Inhibition" = "#6A3E9A",  
                     "Contact_Inhibition_Hypoxia"= "#CBB2D5",
                     "Tissue" = "lightgrey")

mitoPPS <- read_csv(here::here("Data", "Fibroblasts", 
                                        "ProcessedData", "mitoPPS_All.csv")) %>%
  pivot_longer(cols = -c(RNAseq_sampleID, Line,  Study_part, Passage,
                         Days_grown_Udays, Treatments, Clinical_condition, Condition), 
               names_to = "Pathway", values_to = "mitoPPS") %>%
  select(RNAseq_sampleID, Condition, Pathway, Line, mitoPPS) %>%
  unique() 

mitoPPS_Tissues <- read_csv(here::here("Data", "HumanProteinAtlas","ProcessedData", "mitoPPS_all_tissues.csv" )) %>%
  pivot_longer(cols = -c(Tissue, Group), names_to = "Pathway", values_to = "mitoPPS") %>%
  mutate(Condition = "Tissue")

combined <- mitoPPS %>%
  full_join(mitoPPS_Tissues)

test <- combined %>%
  filter(Pathway %in% c("Phospholipid metabolism", "Catechol metabolism")) %>%
  pivot_wider(names_from = "Pathway", values_from = "mitoPPS") 

test %>%
  ggplot(aes(x = `Phospholipid metabolism`, y= `Catechol metabolism`, color = Condition)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "mitoPPS") + 
  scale_color_manual(values = color_conditions)  +
  theme_bw() +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        title = element_text(size = 8),
        legend.position = "none")
ggsave(here::here("Figures", "Figure4", "Suppl_FigS11B_Tissues_vs_fibs_mitoPPS_example1.png"), 
       height = 2.9, width = 2.61,unit = "in", dpi = 1200)

test <- combined %>%
  filter(Pathway %in% c("Cholesterol-associated", "Iron homeostasis")) %>%
  pivot_wider(names_from = "Pathway", values_from = "mitoPPS") 
test %>%
  ggplot(aes(x = `Cholesterol-associated`, y= `Iron homeostasis`, color = Condition)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "mitoPPS") + 
  scale_color_manual(values = color_conditions)  +
  theme_bw() +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        title = element_text(size = 8),
        legend.position = "none")
ggsave(here::here("Figures", "Figure4", "Suppl_FigS11C_Tissues_vs_fibs_mitoPPS_example2.png"), 
       height = 2.9, width = 2.61,unit = "in", dpi = 1200)




# Seahorse correlation ----------------------------------------------------

sh <- read_csv(here::here("Data", "Fibroblasts", "ProcessedData", "Lifespan_Study_data_seahorse.csv")) %>%
  dplyr::select(c(RNAseq_sampleID,
           ATPtotal_UpmolATP_per_min_per_20kcells,	
           ATPglyc_UpmolATP_per_min_per_20kcells,	
           ATPox_UpmolATP_per_min_per_20kcells,	
           ATPox_max_UpmolATP_per_min_per_20kcells,	
           ATPox_spare_UpmolATP_per_min_per_20kcells,	
           Proton_leak_respiration_UpmolO2per_min_per_20kcells,	
           Coupling_efficiency_Upercentage,	
           Non_mitochondrial_respiration_UpmolO2_per_min_per_20kcells)) %>%
  mutate(tmp = "Sample") %>%
  filter(!is.na(RNAseq_sampleID)) %>%
  unite(RNAseq_sampleID, tmp, RNAseq_sampleID) %>%
  dplyr::rename(ATPtotal = ATPtotal_UpmolATP_per_min_per_20kcells,	
         ATPglyc = ATPglyc_UpmolATP_per_min_per_20kcells,	
         ATPox = ATPox_UpmolATP_per_min_per_20kcells,	
         ATPox_max = ATPox_max_UpmolATP_per_min_per_20kcells,	
         ATPox_spare= ATPox_spare_UpmolATP_per_min_per_20kcells,	
         Proton_leak = Proton_leak_respiration_UpmolO2per_min_per_20kcells,	
         Coupling_efficiency = Coupling_efficiency_Upercentage	,
         Non_mito_resp = Non_mitochondrial_respiration_UpmolO2_per_min_per_20kcells) %>%
  pivot_longer(cols = -c(RNAseq_sampleID)) %>%
  group_by(name) %>%
  # remove extreme outliers based on IQR
  mutate(Q1 = quantile(value, 0.25, na.rm = T),
         Q3 = quantile(value, 0.75, na.rm = T),
         IQR = Q3 - Q1,
         lower_bound = Q1 - 3 * IQR,
         upper_bound = Q3 + 3 * IQR,
         is_outlier = ifelse(value < lower_bound | value > upper_bound, TRUE, FALSE)) %>%
  filter(!is_outlier) %>%
  dplyr::select(-c(Q1, Q3, IQR, lower_bound, upper_bound, is_outlier)) %>%
  pivot_wider(names_from = "name", values_from = "value") 

combined <- mitoPPS %>%
  pivot_wider(names_from = "Pathway", values_from = "mitoPPS") %>%
  inner_join(sh) %>%
  pivot_longer(cols = -c(RNAseq_sampleID, Condition, Line), names_to = "var", values_to = "val") %>%
  group_by(Condition, Line, var) %>%
  mutate(mean_val = mean(val)) %>%
  select(-c(RNAseq_sampleID, val)) %>%
  unique() %>%
  pivot_wider(names_from = "var", values_from = "mean_val")

all_pathways <- unique(mitoPPS$Pathway)

cor_res <- combined %>%
  pivot_longer(cols = c("ATPox","ATPglyc","ATPtotal","ATPox_max","ATPox_spare",
                        "Proton_leak","Coupling_efficiency","Non_mito_resp"),
               names_to = "SHvar", values_to= "value") %>%
  pivot_longer(cols = all_of(all_pathways), names_to = "Pathway", values_to = "mitoPPS") %>%
  group_by(Pathway, SHvar) %>%
  nest() %>%
  mutate(cor    = purrr::map(data, ~cor.test(.x$mitoPPS, .x$value, method = "spearman")),
         tidied = purrr::map(cor, broom::tidy)) %>%
  tidyr::unnest(tidied) %>%
  ungroup() %>%                             
  mutate(padj = p.adjust(p.value, method = "BH"))


## Heatmap -----------------------------------------------------------------


set.seed(123)
levels <- readxl::read_xls(here::here("Data/MitoCarta/OriginalData/HumanMitoCarta3_0.xls"), sheet = 4) %>%
  select(MitoPathway, `MitoPathways Hierarchy`) %>%
  separate(`MitoPathways Hierarchy`, into = c("Pathway_Level1", "Pathway_Level2", "Pathway_Level3"), sep = " > ") %>%
  mutate(Level = case_when(
    MitoPathway == Pathway_Level1 ~ "Pathway_Level1", 
    MitoPathway == Pathway_Level2 ~ "Pathway_Level2", 
    MitoPathway == Pathway_Level3 ~ "Pathway_Level3"
  )) %>%
  unique() %>%
  dplyr::rename(Pathway = MitoPathway)

levels_to_color <- levels %>%
  mutate(color = case_when(
    Pathway_Level1 == "Mitochondrial central dogma" ~ "#D95F02",
    Pathway_Level1 == "Mitochondrial dynamics and surveillance" ~ "#7570B3",
    Pathway_Level1 == "Metabolism" ~ "#1D9E77",
    Pathway_Level1 == "OXPHOS" ~ "#E7298A",
    Pathway_Level1 == "Protein import, sorting and homeostasis" ~ "#67A620",
    Pathway_Level1 == "Signaling" ~ "#E6AB00",
    Pathway_Level1 == "Small molecule transport" ~ "#C2A371",
    TRUE ~ "lightgrey"
  )) %>%
  filter(!is.na(Pathway ))
pathway_cols <- setNames(levels_to_color$color, levels_to_color$Pathway)
data_heatmap <- cor_res %>%
  select(SHvar, Pathway, estimate) %>%
  pivot_wider(names_from = "SHvar", values_from = "estimate") %>%
 # mutate
   #full_join(levels_to_color) %>%
  # filter(!is.na(Pathway)) %>%
  column_to_rownames("Pathway")

pathways <- rownames(data_heatmap)
ha_rows <- rowAnnotation(
  Pathway = pathways,                              # categorical vector (one per row)
  col = list(Pathway = pathway_cols),              # mapping: category -> color
  show_annotation_name = FALSE,                    # optional: hide label text
  width = unit(3, "mm"),
  show_legend = FALSE# optional: thickness of the strip
)

sig_mat <- cor_res %>%
  select(SHvar, Pathway, padj) %>%
  mutate(star = case_when(
    padj < 0.001 ~ "***",
    padj < 0.01  ~ "**",
    padj < 0.05  ~ "*",
    TRUE         ~ ""
  )) %>%
  select(Pathway, SHvar, star) %>%
  pivot_wider(names_from = "SHvar", values_from = "star") %>%
  tibble::column_to_rownames("Pathway") %>%
  .[rownames(data_heatmap), colnames(data_heatmap)] %>%
  as.matrix()
mat <- as.matrix(data_heatmap)
hm <- ComplexHeatmap::Heatmap(
  data_heatmap,
  name = "Spearman's rho",
  column_title = "",
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_column_names = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 7, angle = 45),
  column_dend_height = unit(4, "mm"),
  row_dend_width = unit(4, "mm"),
  right_annotation = ha_rows,
  cell_fun = function(j, i, x, y, w, h, fill) {
    lab <- sig_mat[i, j]
    if (!is.na(lab) && nzchar(lab)) {
      txt_col <- ifelse(abs(mat[i, j]) > 0.5, "white", "black")
      y_shifted = y - unit(2, "mm")
      
      grid.text(lab, x = x, y = y,
                gp = gpar(fontsize = 4, fontface = "bold", col = txt_col))
      
    }
  }
)
draw(hm)
pdf(here::here("Figures","Figure4","Suppl_FigS11F_Heatmap_Bioenergetics_cor_withLabels_large.pdf"), 
    height = 15, width = 5)
draw(hm)
dev.off()



# Suppl Fig S11 ATPox correlations ------------------------------------------------------


cor_res %>%
  filter(padj < 0.05) %>%
  filter(SHvar %in% "ATPox") %>%
  arrange(estimate) 

cor_res %>%
  filter(padj < 0.05) %>%
  filter(SHvar %in% "ATPox") %>%
  arrange(desc(estimate)) 

combined %>%
  ggplot(aes(x = `Mitochondrial central dogma`, y= ATPox, color = Condition)) +
  geom_point(size = 1.5, alpha =0.7) +
  scale_color_manual(values = color_conditions) +
  theme_bw() +
  theme_pca 

ggsave(here::here("Figures", "Figure4", "SupplFig_S11D_CentralDogma_vs_ATPox_spare.png"), 
       height = 2.5, width = 2.65,unit = "in", dpi = 1200)


combined %>%
  ggplot(aes(x = `Cholesterol-associated`, y= ATPox, color = Condition)) +
  geom_point(size = 1.5, alpha =0.7) +
  scale_color_manual(values = color_conditions) +
  theme_bw() +
  theme_pca 
ggsave(here::here("Figures", "Figure4", "SupplFig_S11E_chol_assoc_vs_ATPox_spare.png"), 
       height = 2.5, width = 2.65,unit = "in", dpi = 1200)

# Figure 4I Bivariate -----------------------------------------------------


cor_res %>%
  filter(padj < 0.05) %>%
  filter(Pathway %in% "OXPHOS subunits") %>%
  arrange(desc(estimate))  

combined %>%
  ggplot(aes(x = `OXPHOS subunits`, y= ATPox_spare, color = Condition)) +
  geom_point(size = 1.5, alpha =0.7) +
  scale_color_manual(values = color_conditions) +
  theme_bw() +
  geom_smooth(method = "lm", color = "black", se = T, linewidth = 0.3, linetype = "dotted") +
  theme_pca 
ggsave(here::here("Figures", "Figure4", "Fig4I_OXPHOS_vs_ATPox_spare.png"), 
       height = 1.79, width = 1.89,unit = "in", dpi = 1200)




