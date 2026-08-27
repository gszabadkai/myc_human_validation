
library(tidyverse)
library(ComplexHeatmap)
library(Seurat)

scGSVA_list <- readRDS(here::here("Data/Guo_et_al", "ProcessedData", "scGSEA_from_docker.rds"))


scGSVA_list <- scGSVA_list[sapply(scGSVA_list, length) != 1]
scGSVA <- list()
names <- names(scGSVA_list)
for (i in 1:length(scGSVA_list)) {
  data_sub <- as.data.frame(scGSVA_list[[names[i]]][["scgsea"]][["x"]]) %>%
    rownames_to_column("Cell") %>%
    mutate(Run = names[i]) %>%
    pivot_longer(cols = -c(Cell, Run), names_to = "Pathway", values_to = "score")
  scGSVA[[i]] <- data_sub
}
scGSVA_tmp <- bind_rows(scGSVA) 

pathways_to_keep <- scGSVA_tmp %>%
  select(Run, Pathway) %>%
  unique() %>%
  group_by(Pathway) %>%
  mutate(n = n()) %>%
  arrange(n, Pathway) %>%
  select(Pathway, n) %>%
  unique() %>%
  filter(n>15) %>%
  pull(Pathway)
scGSVA <- scGSVA_tmp %>%
  #filter(Pathway %in% c("cAMP-PKA signaling", "CI subunits")) %>%
  filter(Pathway %in% pathways_to_keep) %>%
  group_by(Cell, Pathway) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    #sd_score   = sd(score, na.rm = TRUE),
    #cv         = sd_score / mean_score,
    .groups = "drop"
  ) %>%
  group_by(Pathway) %>%
  mutate(mean_score = scale(mean_score)) %>%
  #select(-c(cv, sd_score)) %>%
  pivot_wider(names_from = Pathway, values_from = mean_score) %>%
  column_to_rownames("Cell") %>%
  t()


umap_coordinates <- read_csv(here::here("Data/Guo_et_al", "ProcessedData", "seurat_umap.csv")) %>%
  mutate(Cell = as.character(Cell)) %>%
 # filter(Cell %in% cells_keep) %>%
  column_to_rownames("Cell") %>%
  as.matrix()

meta <- as.data.frame(readRDS(here::here("Data/Guo_et_al", "ProcessedData","seurat_meta.rds"))) %>%
  rownames_to_column("Cell") #%>%
  #filter(Cell %in% cells_keep)
meta2 <- meta %>%
  column_to_rownames("Cell")
scGSVA_seurat <- CreateSeuratObject(counts = scGSVA, data = scGSVA,meta.data = meta2)
scGSVA_seurat[['umap']] <- CreateDimReducObject(embeddings = umap_coordinates, key = "umap_", 
                                                global = T, assay = "RNA")
Idents(scGSVA_seurat) <- "CellSubtype"
DimPlot(scGSVA_seurat, label = T, group.by = "CellSubtype", reduction = "umap") +
  ggtitle("scGSVA UMAP")
#saveRDS(scGSVA_seurat,here::here("HPA_CellType", "Data", "Processed", "scGSVA_testis_seurat.rds"))

pathway_selection <- rownames(scGSVA_seurat@assays[["RNA"]])
all_expr_values <- FetchData(scGSVA_seurat, vars = pathway_selection)
for (i in 1:length(pathway_selection)) {
  pathway_this <- pathway_selection[i]
  print(pathway_this)
  FeaturePlot(scGSVA_seurat,
              features = pathway_this)+
    scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 6, name = "RdBu")),
                           limits = c(-1, 1), oob = scales::squish)
  ggsave(here::here("Figures", "Figure6", "scGSEA_FeaturePlots", paste0(pathway_this, "_FeaturePlot.png")), width = 5, height = 5)
}



# Main figure: select pathways --------------------------------------------

# pathway_selection <-    c("Cytochromes", "Heme-containing proteins",  "PentosePhosphatePathway",
#                           "Propanoate metabolism", "Vitamin B12 metabolism",
#                           "Glycolysis", "MICOS complex", #"Selenoproteins",
#                           "Fission", "Fusion"
# )
# 
#   
# for (i in 1:length(pathway_selection)) {
#   pathway_this <- pathway_selection[i]
#   print(pathway_this)
#   FeaturePlot(scGSVA_seurat,
#               features = pathway_this)+
#     scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 6, name = "RdBu")),
#                            limits = c(-0.5, 0.5), oob = scales::squish)
#   ggsave(here::here("Figures", "Figure6", "scGSEA_FeaturePlots_MainFigure", paste0(pathway_this, "_FeaturePlot.png")), width = 5, height = 5)
# }



# range(scGSVA["CV subunits", ])
# FeaturePlot(scGSVA_seurat,
#             features = "CV subunits")+
#   scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 6, name = "RdBu")),
#                          limits = c(-1.5, 1.5), oob = scales::squish)
# ggsave(here::here("Figures", "Figure6_sperm", "FeaturePlots_Final", "CV_subunits_FeaturePlot_manual_scale.png"), width = 5, height = 5)
# 
# range(scGSVA["CIV subunits", ])
# FeaturePlot(scGSVA_seurat,
#             features = "CIV subunits")+
#   scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 6, name = "RdBu")),
#                          limits = c(-1.2, 1.2), oob = scales::squish)
# ggsave(here::here("Figures", "Figure6_sperm", "FeaturePlots_Final", "CIV_subunits_FeaturePlot_manual_scale.png"), width = 5, height = 5)


# 
# seu <- CreateSeuratObject(counts = scGSVA, data = scGSVA,meta.data = meta2)
# seu[['umap']] <- CreateDimReducObject(embeddings = umap_coordinates, key = "umap_", 
#                                                 global = T, assay = "RNA")
# 
# 
# #seu <- NormalizeData(seurat)
# seu <- FindVariableFeatures(seu, selection.method = "vst")
# seu <- ScaleData(seu)
# seu <- RunPCA(seu, features = VariableFeatures(object = seu))
# seu <- RunUMAP(seu, dims = 1:30)
# seu <- FindNeighbors(seu)
# seu <- FindClusters(seu, resolution = 0.5)
# DimPlot(seu, group.by = "seurat_clusters", label = T)
# for (i in 1:length(pathway_selection)) {
#   pathway_this <- pathway_selection[i]
#   print(pathway_this)
#   FeaturePlot(seu,
#               features = pathway_this)+
#     scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 6, name = "RdBu")),
#                            limits = c(-3, 3), oob = scales::squish)
#   ggsave(here::here("Figures", "Figure6_sperm", "FeaturePlots_recluster", paste0(pathway_this, "_FeaturePlot_unscaled.png")), width = 5, height = 5)
# }
# 
# Idents(seurat) <- "CellSubtype"
# DimPlot(seurat, label = F)+
#   theme(legend.position = "none")
# ggsave(here::here("Figures", "Figure6_sperm", "seurat_subset_final.png"), width =5, height = 5)
# 
# 
