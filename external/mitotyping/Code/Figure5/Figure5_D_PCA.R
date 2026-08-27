library(tidyverse)
library(broom)
library(ComplexHeatmap)

color_celltypes =  c(
  "unspecific" = "grey",
  "Inhibitory" = "#FF7F00",# 
  "Excitatory" = "#007200",
  "OPC" = "#A6CEE3", 
  "Oligodendrocytes" = "#1F78B4",
  "Astrocytes" = "#E15A57",
  "Microglia" ="#FFFF99",
  "VLMCs" =  "#6A3D9A",
  "Pericytes" = "#D0BEDA",
  "Endothelial" = "#B789DA",
  "Brain"= "#1DB100"
  
)

theme_pca <-   theme(
  axis.title = element_text(size = 7),
  axis.text = element_text(size = 6),
  legend.text = element_text(size = 6),
  legend.title = element_blank(),
  axis.line = element_line(linewidth  = 0.25),
  axis.ticks = element_line(linewidth  = 0.25),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank())


# Averages ----------------------------------------------------------------


mitoPPS_MBM <- read_csv(here::here("Data", "MBM_ROSMAP", "ProcessedData", "MBM_mitoPPS.csv")) %>%
  filter(!Voxel %in% "WM") %>%
  group_by(CellType, Pathway1) %>%
  mutate(mitoPPS = mean(mitoPPS)) %>%
  select(-Voxel) %>%
  unique() %>%
  mutate(Dataset = "MBM")

mitoPPS_RM <- read_csv(here::here("Data", "MBM_ROSMAP", "ProcessedData", "ROSMAP_mitoPPS.csv")) %>%

  separate(Group, into = c("CellType", "ID"), sep = "_", extra = "merge") %>%
  group_by(CellType, Pathway1) %>%
  mutate(mitoPPS = mean(mitoPPS)) %>%
  select(-ID) %>%
  unique() %>%
  mutate(CellType = case_when(
    CellType == "Ast" ~ "Astrocytes",
    CellType == "End" ~ "Endothelial",
    CellType == "Exc" ~ "Excitatory",
    CellType == "Inh" ~ "Inhibitory",
    CellType == "Mic" ~ "Microglia",
    CellType == "Oli" ~ "Oligodendrocytes",
    CellType == "OPC" ~ "OPC",
    CellType == "Peri" ~ "Pericytes",
    TRUE ~ CellType)) %>%
  mutate(Dataset = "ROSMAP") 

combined <- bind_rows(mitoPPS_MBM, mitoPPS_RM) %>%
  unite(Group, Dataset, CellType,  sep = "_", remove = FALSE) %>%
  pivot_wider(names_from = Pathway1, values_from = mitoPPS) %>%
  column_to_rownames(var = "Group") 


# PCA ---------------------------------------------------------------------

library(ggfortify)
pca <- prcomp(combined[,-c(1:3)], center = T, scale. = T)
loadings <- pca$rotation
p1 <-autoplot(pca, data = combined, colour = 'CellType', shape = 'Dataset', label = F,
              size =1.9) +
  scale_color_manual(values = unlist(color_celltypes)) +
  theme_classic() +
  theme_pca

p1
ggsave(here::here("Figures", "Figure5",  "Figure_5D_PCA_1.pdf"), width = 3.25, height =1.9, units = "in", dpi = 1200)

p1 <-autoplot(pca,x=2, y=3, data = combined, colour = 'CellType', shape = 'Dataset', label = F,
              size =1.9) +
  scale_color_manual(values = unlist(color_celltypes)) +
  theme_classic() +
  theme_pca
p1
ggsave(here::here("Figures", "Figure5",  "Figure_5D_PCA_2.pdf"), width = 3.25, height =1.9, units = "in", dpi = 1200)




mitoPPS_MBM <- read_csv(here::here("Data", "MBM_ROSMAP", "ProcessedData", "MBM_mitoPPS.csv")) %>%
  filter(!Voxel %in% "WM") %>%
  dplyr::rename(ID = Voxel) %>%
  mutate(Dataset = "MBM")

mitoPPS_RM <- read_csv(here::here("Data", "MBM_ROSMAP", "ProcessedData", "ROSMAP_mitoPPS.csv")) %>%
  separate(Group, into = c("CellType", "ID"), sep = "_", extra = "merge") %>%
  mutate(CellType = case_when(
    CellType == "Ast" ~ "Astrocytes",
    CellType == "End" ~ "Endothelial",
    CellType == "Exc" ~ "Excitatory",
    CellType == "Inh" ~ "Inhibitory",
    CellType == "Mic" ~ "Microglia",
    CellType == "Oli" ~ "Oligodendrocytes",
    CellType == "OPC" ~ "OPC",
    CellType == "Peri" ~ "Pericytes",
    TRUE ~ CellType)) %>%
  mutate(Dataset = "ROSMAP")

combined <- bind_rows(mitoPPS_MBM, mitoPPS_RM) %>%
  unite(Group, Dataset, CellType, ID, sep = "_", remove = FALSE) %>%
  pivot_wider(names_from = Pathway1, values_from = mitoPPS) %>%
  column_to_rownames(var = "Group") 

pca <- prcomp(combined[,-c(1:3)], center = T, scale. = T)
loadings <- pca$rotation

p1 <-autoplot(pca, data = combined, colour = 'CellType', shape = 'Dataset',alpha = 'Dataset', size = 'Dataset', label = F) +
  scale_color_manual(values = unlist(color_celltypes)) +
  scale_alpha_manual(values = c('MBM' = 1, 'ROSMAP' = 0.2)) +
  scale_size_manual(values = c('MBM' = 1, 'ROSMAP' = 0.5)) +
  theme_classic() +
  theme_pca

p1

ggsave(here::here("Figures", "Figure5",  "Suppl_Figure_12B_PCA_1.pdf"), width = 3.25, height =1.9, units = "in", dpi = 1200)
plotly::ggplotly(p1)
p1 <-autoplot(pca,x=2, y=3, data = combined, colour = 'CellType', shape = 'Dataset',alpha = 'Dataset', size = 'Dataset', label = F) +
  scale_color_manual(values = unlist(color_celltypes)) +
  scale_alpha_manual(values = c('MBM' = 1, 'ROSMAP' = 0.2)) +
  scale_size_manual(values = c('MBM' = 1, 'ROSMAP' = 0.5)) +
  scale_color_manual(values = unlist(color_celltypes)) +
  theme_classic() +
  theme_pca
p1
ggsave(here::here("Figures", "Figure5",  "Suppl_Figure_12C_PCA_2.pdf"), width = 3.25, height =1.9, units = "in", dpi = 1200)

