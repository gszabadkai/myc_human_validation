
rm(list = ls())

# Setup -------------------------------------------------------------------

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggfortify)

color_groups <- c( "CNS"= "#1DB100",
                   "Contractile"="#F8BA00",
                   "Reproductive"="#EE220C",
                   "Digestive"="#265A8C",
                   "Anabolic"="#EF5FA7",
                   "Secretory"="#A8A8A8",
                   "Other"="khaki",
                   "Immune"="#00A89D"
)

tissue_to_group_ms <- function(x) {
  mutate(x, Group = case_when(
    Tissue == "Cerebrum" ~ "CNS",
    Tissue == "Cerebellum" ~ "CNS",
    Tissue == "Brainstem" ~ "CNS",
    Tissue == "Spinal cord" ~ "CNS",
    Tissue == "Kidney" ~ "Anabolic",
    Tissue == "Liver" ~ "Anabolic",
    Tissue == "Heart" ~ "Contractile",
    Tissue == "Skeletal muscle" ~ "Contractile",
    Tissue == "Adipose" ~ "Secretory",
    Tissue == "Small intestine" ~ "Digestive",
    Tissue == "Large intestine" ~ "Digestive",
    Tissue == "Stomach" ~ "Digestive",
    Tissue == "Placenta" ~ "Reproductive",
    Tissue == "Testis" ~ "Reproductive"
  ), .after = "Tissue")
}

# Data_preparation ms ---------------------------------------------------------------

## Read data
original_data <- readxl::read_xls(here::here("Data", "MitoCarta","OriginalData","MouseMitoCarta3_0.xls"), sheet = 2) %>%
  dplyr::rename(Gene = Symbol)
original_data <- original_data[,c(3, 34:47)]

## Find proteins with no expression in any tissue
no_exprs <- original_data %>%
  mutate(count_na = rowSums(is.na(.))) %>% 
  filter(count_na %in% 14) 
no_exprs <- no_exprs$Gene 
excluded_data <- readxl::read_xls(here::here("Data","MitoCarta","OriginalData","MouseMitoCarta3_0.xls"), sheet = 2) %>%
  dplyr::rename(Gene = Symbol) %>%
  filter(Gene %in% no_exprs)
write.csv(excluded_data, here::here("Data","MitoCarta","ProcessedData", "Mouse_MitoCarta_163_Excluded_Genes.csv"))


## Revert log10 transform, remove genes without expression in any tissue
### 1140 genes reduced to 977, 63 genes removed
filtered_data <- original_data %>% 
  pivot_longer(cols = colnames(original_data[,2:ncol(original_data)])) %>%
  mutate(value = 10^value) %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  filter(!Gene %in% no_exprs)

## NA extrapolation
gene_min <- filtered_data %>%
  pivot_longer(cols = colnames(filtered_data[,2:ncol(filtered_data)])) %>%
  na.omit() %>%
  group_by(Gene) %>%
  mutate(min = min(value)) %>%
  select(Gene, min) %>%
  unique()

NAextrapolated <- filtered_data %>%
  full_join(gene_min, by = "Gene") %>%
  mutate(NArep = min/2)

processed_data <- NAextrapolated %>%
  pivot_longer(cols = colnames(NAextrapolated[,2:15])) %>%
  select(-min) %>%
  mutate(NAextrapolated = case_when(
    is.na(value) ~ NArep,
    !is.na(value) ~ value
  )) %>%
  select(Gene, name, NAextrapolated) %>%
  na.omit() %>%
  pivot_wider(names_from = "name", values_from = "NAextrapolated") %>%
  dplyr::rename(Cerebrum = cerebrum_total_peak_intensity_log10,
         Cerebellum = cerebellum_total_peak_intensity_log10,    
         Brainstem = brainstem_total_peak_intensity_log10,
         `Spinal cord` = spinalcord_total_peak_intensity_log10,
         Kidney = kidney_total_peak_intensity_log10,
         Liver = liver_total_peak_intensity_log10,
         Heart = heart_total_peak_intensity_log10,
         `Skeletal muscle` = skeletalmuscle_total_peak_intensity_log10,
         Adipose = adipose_total_peak_intensity_log10,
         `Small intestine`= smallintestine_total_peak_intensity_log10,
         `Large intestine` = largeintestine_total_peak_intensity_log10,
         Stomach = stomach_total_peak_intensity_log10,
         Placenta = placenta_total_peak_intensity_log10,
         Testis = testis_total_peak_intensity_log10) 

write.csv(processed_data, here::here("Data","MitoCarta","ProcessedData", "Mouse_ProcessedData.csv"), row.names=F)

rm(list = setdiff(ls(), c("processed_data", "color_groups", "color_tissues", "tissue_to_group_ms")))

# Figure 2A Cluster Heatmap ms --------------------------------------------

## Read the data
data <- processed_data #read_csv(here::here("Data","MitoCarta", "ProcessedData", "Mouse_ProcessedData.csv"))
data <- processed_data %>%
  pivot_longer(cols = colnames(data[2:ncol(data)])) %>%
  pivot_wider(names_from = "Gene", values_from = "value") %>%
  column_to_rownames("name")

## Heatmap data preparation
exprs <- log10(data)
col_fun = colorRamp2(c(4, 8,  12),c("#fff899", "#ffa76b", "#ff4242"))

### Clustering Columns
distance    = dist(t(as.matrix(exprs)), method="euclidean" ) 
coldend     = hclust(distance, method="ward.D2")  
### Clustering Rows
row_dist    = dist(as.matrix(exprs), method="euclidean" ) 
rowdend     = hclust(row_dist, method="ward.D2")  

## Column annotation
row_anno_df <- data %>%
  as.data.frame() %>%
  rownames_to_column("Tissue") %>%
  tissue_to_group_ms() %>%
  select(Tissue, Group)

row_anno = rowAnnotation(
  `Tissue group`=row_anno_df$Group,
  col=list(`Tissue group`  = color_groups),
  show_annotation_name = F,
  show_legend =  F,
  simple_anno_size = unit(0.1, "in"),
  annotation_legend_param = list(nrow=6))

## Heatmap
HM <- Heatmap(exprs, 
              name = "Log(protein exprs)", 
              col=col_fun,
              column_title="977 mitochondrial proteins",
              show_column_dend=F,
              show_column_names = F,
              cluster_columns =coldend,
              row_title="",
              row_dend_side = "right",
              row_names_side = "left",
              show_row_dend=TRUE,
              show_row_names=T,
              cluster_rows =rowdend,
              row_names_gp = grid::gpar(fontsize = 6),
              column_title_gp = grid::gpar(fontsize = 6),
              row_dend_width = unit(0.15, "in"),
              right_annotation=row_anno,
              width = unit(0.0013, "in")*977,
              heatmap_legend_param = list(
                at = c(4, 8, 12),
                title = "Log10 exprs",
                grid_height = unit(0.6, "in"),
                grid_width = unit(0.1, "in"),
                title_gp = gpar(fontsize = 6),
                labels_gp = gpar(fontsize = 6)
              )
)

png(here::here("Figures","Figure2", "Fig2A_ClusterHeatmap_ms.png"),height=1.97,width=3,units="in",
    res=1200, 
    pointsize = 6)
draw(HM, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off() 

rm(list = setdiff(ls(), c("processed_data", "color_groups", "color_tissues", "tissue_to_group_ms")))

# Figure 2B PCA Mouse -----------------------------------------------------

data <- processed_data %>%
  pivot_longer(cols = -Gene, names_to = "Tissue") %>%
  tissue_to_group_ms() %>%
  mutate(value = log10(value)) %>%
  pivot_wider(names_from = "Gene", values_from = "value")
pca <- prcomp(data[,-c(1,2)], scale. = T)
summary(pca)
p <- autoplot(pca, data = data, colour = 'Group', size = 2.8, alpha = 0.7)+#,loadings = T, loadings.label = TRUE, loadings.label.size  = 1) +
  theme_bw() +
  scale_color_manual(values = c(c("CNS"= "#1DB100",
                                  "Contractile"="#F8BA00",
                                  "Reproductive"="#EE220C",
                                  "Digestive"="#265A8C",
                                  "Anabolic"="#EF5FA7",
                                  "Secretory"="#A8A8A8",
                                  "Other"="khaki",
                                  "Immune"="#00A89D"))) +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())
p
png(here::here( "Figures", "Figure2","Fig2B_PCA_ms.png"),width=2.22,height=2.22,units="in",res=1200)
print(p)
dev.off()

p <- autoplot(pca,x = 3, y= 4, data = data, colour = 'Group', size = 2.8, alpha = 0.7)+#,loadings = T, loadings.label = TRUE, loadings.label.size  = 1) +
  theme_bw() +
  scale_color_manual(values = c(c("CNS"= "#1DB100",
                                  "Contractile"="#F8BA00",
                                  "Reproductive"="#EE220C",
                                  "Digestive"="#265A8C",
                                  "Anabolic"="#EF5FA7",
                                  "Secretory"="#A8A8A8",
                                  "Other"="khaki",
                                  "Immune"="#00A89D"))) +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())
p
png(here::here( "Figures", "Figure2","Suppl_Fig2A_PCA2_ms.png"),width=2.22,height=2.22,units="in",res=1200)
print(p)
dev.off()
plotly::ggplotly(autoplot(pca,x = 3, y= 4,data = data, colour = 'Tissue',
                          fill = 'Group', size = 2.8, alpha = 0.6, shape = 21)+
                   scale_fill_manual(values = color_groups) +
                   theme_bw())
# Supplemental figures


# FAO in mouse ------------------------------------------------------------


FAO_genes <- readxl::read_xls(here::here("Data", "MitoCarta", "OriginalData",
                                         "MouseMitoCarta3_0.xls"), sheet = 4) %>%
  select(MitoPathway, Genes) %>%
  na.omit() %>%
  splitstackshape::cSplit('Genes', ',') %>%
  column_to_rownames("MitoPathway") %>%
  t() %>%
  as.data.frame() %>%
  select(`Fatty acid oxidation`)  %>%
  na.omit() %>%
  pull(`Fatty acid oxidation`) 

data_FAO <- processed_data %>%
  filter(Gene %in% FAO_genes) %>%
  pivot_longer(cols = -Gene, names_to = "Tissue") %>%
  tissue_to_group_ms() %>%
  group_by(Tissue) %>%
  mutate(average_FAO_score = mean(value)) %>%
  select(Tissue, Group, average_FAO_score) %>%
  unique() %>%
  arrange(desc(average_FAO_score)) 
data_FAO


p <- data_FAO  %>%
  ggplot(aes(fill = Group)) +
  geom_segment(aes(y=reorder(Tissue, average_FAO_score, decreasing = F),yend=reorder(Tissue, average_FAO_score, decreasing = F),x=1, xend=average_FAO_score,
                   color = Group), size=2.5, alpha =0.6) +
  scale_fill_manual(values = color_groups) +
  scale_color_manual(values = color_groups) +
  labs(x= "FAO Pathway score") +
  
  theme_bw() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(),
        axis.text.y = element_text(size =6),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 6),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )
print(p)
ggsave(here::here("Figures", "Figure2", "Suppl_FigS2A_FAO_ms.png"), width = 2, height = 1.84, units = "in")


## Fold change
data_FAO %>%
  mutate(Group = case_when(Group == "CNS" ~ "CNS",
                           TRUE ~ Tissue)) %>%
  group_by(Group) %>%
  mutate(Group_average = mean(average_FAO_score)) %>%
  select(Group, Group_average) %>%
  unique() %>%
  pivot_wider(names_from = "Group", values_from = "Group_average") %>% 
  pivot_longer(cols = -CNS, names_to = "Tissue") %>%
  mutate(fc = value / CNS)

  



