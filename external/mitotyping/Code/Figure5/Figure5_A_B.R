library(Seurat)
library(ggfortify)

color_celltypes =  c(
  "Inhibitory" = "#FF7F00",# 
  "Excitatory" = "#007200",
  "OPC" = "#A6CEE3", 
  "Oligodendrocytes" = "#1F78B4",
  "Astrocytes" = "#E15A57",
  "Microglia" ="#FFFF99",
  "VLMCs" =  "#6A3D9A",
  "Pericytes" = "#D0BEDA",
  "Endothelial" = "#B789DA"
  
)


# MBM UMAP ----------------------------------------------------------------


mitobrain_data <- readRDS(here::here("Data", "MBM_ROSMAP", "OriginalData",  "MitoBrainMap_SingleCellData_Downsampled_clean.rds"))
mitobrain_data <- subset(mitobrain_data, subset = orig.ident != "voxel4")
mitobrain_data <- subset(mitobrain_data, subset = CellType != "VLMCs")
mitobrain_data <- RunUMAP(mitobrain_data, dims = 1:30, reduction = "harmony")
Idents(mitobrain_data) <- mitobrain_data$CellType
DimPlot(mitobrain_data, cols = color_celltypes, label = F, pt.size = 0.0001) +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank())
ggsave(here::here("Figures", "Figure5", "Figure5A_MitoBrainMap_UMAP.png"), width = 4, height = 3.8, units = "in", dpi = 1200)



# ROSMAP PCA --------------------------------------------------------------


data <- read.csv(here::here("Data", "MBM_ROSMAP", "OriginalData", 
                            "unified-cpm.log2cpm.exclude-Broad.rename-gene-and-donor.csv")) 
data_pca <- data %>% 
  select(-gene_id) %>%
  pivot_longer(cols = -gene_symbol, names_to = "SampleID", values_to = "value") %>%
  filter(grepl("Inh|Exc|OPC|Oli|Ast|Mic|Peri|End", SampleID)) %>%
  group_by(gene_symbol, SampleID) %>%
  summarize(value = mean(value, na.rm= T),.groups = "drop") %>%
  pivot_wider(names_from = "gene_symbol", values_from = "value") %>%
  separate(SampleID, into = c("CellType", "ID")) 
pca <- prcomp(data_pca[,-c(1:2)])
autoplot(pca, data = data_pca, colour = "CellType", size = 0.5) +
  scale_color_manual(values = c(
    "Inh" = "#FF7F00",# 
    "Exc" = "#007200",
    "OPC" = "#A6CEE3", 
    "Oli" = "#1F78B4",
    "Ast" = "#E15A57",
    "Mic" ="#FFFF99",
    "Peri" = "#D0BEDA",
    "End" = "#B789DA")
  ) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank())
ggsave(here::here("Figures", "Figure5", "Figure5A_RM_PCA.png"), width = 4, height = 3.8, units = "in", dpi = 1200)

