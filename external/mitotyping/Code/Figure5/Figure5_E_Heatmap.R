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



color_dataset <- c(
  "MBM" = "#696969",
  "ROSMAP" = "#AFC1E2"
)

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


combined <- bind_rows(mitoPPS_MBM, mitoPPS_RM) %>%
  unite(Group, Dataset, CellType, ID, sep = "_", remove = FALSE) %>%
  pivot_wider(names_from = Pathway1, values_from = mitoPPS) %>%
  column_to_rownames(var = "Group") 



exprs <- scale(t(combined[,-c(1:3)])) ## scaled because different magnitude scales
## Heatmap annotation
anno_df <- as.data.frame(combined) %>%
  pivot_longer(cols = -c(ID, CellType,Dataset), names_to= "Pathway") %>%
  na.omit() %>%
  pivot_wider(names_from = "Pathway", values_from = "value") %>%
  dplyr::select(ID, CellType)
anno_df <- combined  %>%
  select(ID, CellType,Dataset)

col_anno = columnAnnotation(
  Celltype=anno_df$CellType,
  Dataset = anno_df$Dataset,
  col=list(
    Celltype  =unlist(color_celltypes),
    Dataset  =unlist(color_dataset)
    
  ),
  show_annotation_name = F,
  show_legend =  T,
  annotation_legend_param = list(nrow=12)) 

col_dist    = dist(t(as.matrix(exprs)), method="euclidean" )
coldend     = hclust(col_dist, method="ward.D2")
row_dist    = dist(as.matrix(exprs), method="euclidean" )
rowdend     = hclust(row_dist, method="ward.D2")

hm <- Heatmap(exprs,  
              name = "Rel. expression",
              show_column_dend=T,
              show_column_names =T,
              show_row_names =T,
              show_row_dend=T,
              row_title="",
              cluster_columns =coldend,
              cluster_rows =rowdend,
              column_names_rot = 45,
              row_names_gp = grid::gpar(fontsize = 4),
              row_title_gp = grid::gpar(fontsize = 10),
              column_names_gp = grid::gpar(fontsize = 0.1),
              column_title_gp = grid::gpar(fontsize = 10),
              top_annotation=col_anno,
              width = unit(13, "cm"),
              height = unit(20, "cm"))
pdf(here::here("Figures", "Figure5",  "Figure_5E_Heatmap.pdf"), width = 11, height = 10)
draw(hm)
dev.off()

labels <- hm@column_names_param[["labels"]]
order <- hm@column_dend_param[["obj"]][["order"]]
# reorder test based on the vector order
labels_ordered <- as.data.frame(labels[order])
labels_ordered$`labels[order]`[which(grepl("MBM", labels_ordered$`labels[order]`))]
