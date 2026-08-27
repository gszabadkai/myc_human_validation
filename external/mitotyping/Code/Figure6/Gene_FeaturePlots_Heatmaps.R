library(future)
workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = workers) 

library(SingleCellExperiment)
library(slingshot)
library(tidyverse)
library(ComplexHeatmap)
library(Seurat)
# library(SeuratDisk)
# library(SeuratData)
library(broom)
library(reshape2)
library(mgcv)


# Load the data -----------------------------------------------------------

gene_sets <- readRDS(here::here("Data/Guo_et_al", "ProcessedData", "gene_sets.rds"))
gene_to_pathway <- stack(gene_sets)
colnames(gene_to_pathway) <- c("Gene", "Pathway")

seurat <- readRDS(here::here("Data/Guo_et_al", "ProcessedData", "seurat_object.rds"))
# 
# library(patchwork)
# all_pathways <- unique(gene_to_pathway$Pathway)
# for (i in 1:length(all_pathways)) {
#   pathway_this <- all_pathways[i]
#   genes_this <- gene_to_pathway %>%
#     filter(Pathway %in% pathway_this) %>%
#     pull(Gene)                  
#   
#   if ((length(genes_this) >20 & length(genes_this) <40 )) {
#     width <- 20
#     height <- 20
#   } else if (length(genes_this) ==1){
#     width = 5
#     height = 5
#   } else if ((length(genes_this) >1 & length(genes_this) <10 )){
#     width = 10
#     height = 10
#   }else if ((length(genes_this) >10 & length(genes_this) <20 )){
#     width = 15
#     height = 15
#   }  else if (length(genes_this) >40 ){
#     width = 30
#     height = 30
#   } else {
#     width = 20
#     height = 20
#   }
#   
#   # 1) get a LIST of single‐feature FeaturePlots
#   plot_list <- FeaturePlot(
#     object   = seurat,
#     features = genes_this,
#     reduction = "umap",
#     pt.size  = 1,
#     combine  = FALSE
#   )
#   
#   p <- wrap_plots(plot_list) +
#     plot_annotation(title = pathway_this) 
#   p <- p & 
#     theme_minimal() & 
#     theme(legend.position = "none")
#   ggsave(here::here("Figures/Figure6", "Gene_FeaturePlots", paste0(pathway_this, "_FeaturePlot.png")), 
#          plot = p, width = width, height = height,   bg = "white" )
# }
# 
# 



mat <- as.matrix(GetAssayData(seurat, slot = "data"))
#mat <- mat[unique(gene_to_pathway$Gene), ]
mat <- as.data.frame(mat) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Cell", values_to = "exprs")

meta <- readRDS(here::here("Data/Guo_et_al", "ProcessedData", "seurat_meta.rds")) %>%
  rownames_to_column("Cell") 
cell_to_bin <- read_csv(here::here("Data/Guo_et_al", "ProcessedData", "cell_to_bin.csv")) %>%
  mutate(Cell = as.character(Cell))
df <- full_join(mat, cell_to_bin) %>%
  filter(Gene %in% unique(gene_to_pathway$Gene)) %>%
  group_by(bin, Gene) %>%
  mutate(exprs = mean(exprs, na.rm = T)) %>%
  select(bin, mean_pt, Gene, exprs) %>%
  unique() %>%
  full_join(read_csv(here::here("Data/Guo_et_al", "ProcessedData", "meta_binned.csv"))) %>%
  full_join(gene_to_pathway) %>%
  na.omit() %>%
  group_by(Gene) %>%
  mutate(exprs = scale(exprs))
head(df)
colors <- c(
  "SSC" = "#A98FFF",
  "Diff_spermatogonia" = "#F87169",
  "Early_prim_spermatocytes" = "#5AB711",
  "Late_prim_spermatocytes" = "#FB62D8",
  "Round_spermatids" = "#01B6EB",
  "Elongated_spermatids" = "#1AC49B",
  "Late_spermatids" = "#C49B02"
  
)
all_pathways <- unique(gene_to_pathway$Pathway)
for (i in 1:length(all_pathways)) {
  pathway_this <- all_pathways[i]
  genes_this <- gene_to_pathway %>%
    filter(Pathway %in% pathway_this) %>%
    pull(Gene)      
  df_sub <- df %>%
    filter(Pathway %in% pathway_this) %>%
    pivot_wider(names_from = "Gene", values_from = "exprs") %>%
    arrange(bin)
  exprs <- t(df_sub[,-c(1:4)])
  col_anno <- columnAnnotation(
    binned_pseudotime = anno_barplot(df_sub$mean_pt),
    celltype = df_sub$dominant_cell_type,
    col = list(celltype = colors)
    
  )
  hm <- try(ComplexHeatmap::Heatmap(exprs,
                                    column_title = pathway_this,
                          row_names_gp = gpar(fontsize = 5),
                          cluster_columns = F,
                          top_annotation = col_anno
  ))
  pdf(here::here("Figures/Figure6", "Gene_Heatmaps", paste0(pathway_this, "_Gene_Heatmap.pdf")), height = 10, width =10)
  try(draw(hm))
  dev.off()

}

library(viridisLite)
library(circlize)
library(ComplexHeatmap)
library(grid)  # for unit()

pathway_this <- "Fusion"
genes_this <- gene_to_pathway %>%
  filter(Pathway %in% pathway_this) %>%
  pull(Gene)      
df_sub <- df %>%
  filter(Pathway %in% pathway_this) %>%
  pivot_wider(names_from = "Gene", values_from = "exprs") %>%
  arrange(bin)
exprs <- t(df_sub[,-c(1:4)])
rng <- range(df_sub$mean_pt, na.rm = TRUE)
pal    <- rev(viridis(256))  # "yellow -> green -> teal -> blue -> purple"
breaks <- seq(rng[1], rng[2], length.out = length(pal))
col_fun <- circlize::colorRamp2(breaks, pal)

col_anno <- columnAnnotation(
  Pseudotime = anno_simple(df_sub$mean_pt, col = col_fun, border = FALSE),
  `Cell state`   = df_sub$dominant_cell_type,
  col        = list(`Cell state` = colors),
  annotation_height = unit.c(unit(6, "mm"), unit(6, "mm")),
  annotation_name_gp = gpar(fontsize = 14)
)

hm1 <- ComplexHeatmap::Heatmap(exprs,
                                  row_names_gp = gpar(fontsize = 10),
                                  cluster_columns = F,
                                  top_annotation = col_anno,
                               col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                        row_dend_width = unit(0.5, "cm")
)


pathway_this <- "Fission"
genes_this <- gene_to_pathway %>%
  filter(Pathway %in% pathway_this) %>%
  pull(Gene)      
df_sub <- df %>%
  filter(Pathway %in% pathway_this) %>%
  pivot_wider(names_from = "Gene", values_from = "exprs") %>%
  arrange(bin)
exprs <- t(df_sub[,-c(1:4)])
rng <- range(df_sub$mean_pt, na.rm = TRUE)
pal    <- rev(viridis(256))  # "yellow -> green -> teal -> blue -> purple"
breaks <- seq(rng[1], rng[2], length.out = length(pal))
col_fun <- circlize::colorRamp2(breaks, pal)
col_anno <- columnAnnotation(
  pseudotime = anno_simple(df_sub$mean_pt, col = col_fun, border = FALSE),
  celltype   = df_sub$dominant_cell_type,
  col        = list(celltype = colors)
)


hm2 <- ComplexHeatmap::Heatmap(exprs,
                               row_names_gp = gpar(fontsize = 10),
                               cluster_columns = F,
                               col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                              # top_annotation = col_anno,
                               row_dend_width = unit(0.5, "cm")
)


ht_list = hm1 %v% hm2
draw(ht_list)

pdf(here::here("Figures/Figure6", "Figure6H_Fission_Fusion_heatmap.pdf"), height = 5, width =8, colormodel = "srgb")
print(ht_list)
dev.off()
