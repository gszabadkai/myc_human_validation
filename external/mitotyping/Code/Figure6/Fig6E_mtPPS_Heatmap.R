library(tidyverse)
library(ComplexHeatmap)
mitoPPS <- readRDS(here::here("Data/Guo_et_al", "ProcessedData", "mitoPPS_bins.rds"))
cell_cycle_score <- read_csv(here::here("Data/Guo_et_al", "ProcessedData", "seurat_cell_cycle_scores.csv")) %>%
  full_join(read_csv(here::here("Data/Guo_et_al", "ProcessedData", "cell_to_bin.csv"))) %>%
  group_by(bin) %>%
  mutate(mean_s = mean(S.Score, na.rm = TRUE),
         mean_g2m = mean(G2M.Score, na.rm = TRUE)) %>%
  ungroup() 

dominating_phase <- janitor::tabyl(cell_cycle_score %>%
                                        ungroup() %>%
                                        select(bin, Phase) , Phase, bin) %>%
  pivot_longer(cols = -Phase, names_to = "bin", values_to = "n") %>%
  filter(!n==0) %>%
  mutate(bin = as.numeric(bin)) %>%
  group_by(bin) %>%
  mutate(
    dominant_Phase = Phase[which.max(n)]
  ) %>%
  ungroup() %>%
  select(bin, dominant_Phase) %>%
  unique()

colors <- c(
  "SSC" = "#A98FFF",
  "Diff_spermatogonia" = "#F87169",
  "Early_prim_spermatocytes" = "#5AB711",
  "Late_prim_spermatocytes" = "#FB62D8",
  "Round_spermatids" = "#01B6EB",
  "Elongated_spermatids" = "#1AC49B",
  "Late_spermatids" = "#C49B02"
  
)
meta_binned <- read_csv(here::here("Data/Guo_et_al", "ProcessedData", "meta_binned.csv")) 
df <- mitoPPS %>%
  inner_join(dominating_phase, by = "bin") %>%
  inner_join(meta_binned, by = "bin") %>%
  mutate(mitoPPS = log10(mitoPPS)) %>%
  pivot_wider(names_from = Pathway1, values_from = mitoPPS) %>%
  arrange(mean_pt) 
exprs <- t(df[,-c(1:4)])
library(viridisLite)
library(circlize)

rng <- range(df$mean_pt, na.rm = TRUE)
pal    <- rev(viridis(256)) 
breaks <- seq(rng[1], rng[2], length.out = length(pal))
col_fun <- circlize::colorRamp2(breaks, pal)


col_anno <- columnAnnotation(
  pseudotime = anno_simple(df$mean_pt, col = col_fun, border = FALSE),
  `cell cylce phase`      = df$dominant_Phase,
  `cell state`   = df$dominant_cell_type,
  col        = list(`cell state` = colors,
                    `cell cylce phase` = c("G1" ="#803E75" , "S" ="#FFB300" , "G2M" = "#FF6800")),
  annotation_height = unit.c(unit(4, "mm"), unit(4, "mm"), unit(4, "mm")),
  annotation_name_gp = gpar(fontsize =10)
)

pdf(here::here("Figures", "Figure6", "Figure6E_mitoPPS_heatmap.pdf"), width = 7, height = 6 )
ComplexHeatmap::Heatmap(exprs,
                        name  = "log10(mitoPPS)",
                        row_names_gp = gpar(fontsize = 4),
                        cluster_columns = F,
                        top_annotation = col_anno,
                        cluster_rows = T,
                        row_dend_width = unit(0.5, "cm"),
                        # add range to color scale
                        col = circlize::colorRamp2(c(-0.2, 0, 0.2), c("blue", "white", "red")),
                        km = 4
)
dev.off()

pdf(here::here("Figures", "Figure6", "Suppl_FigS14A_mitoPPS_heatmap.pdf"), width = 10, height = 12)
ComplexHeatmap::Heatmap(exprs,
                        name  = "log10(mitoPPS)",
                        row_names_gp = gpar(fontsize = 6),
                        cluster_columns = F,
                        top_annotation = col_anno,
                        cluster_rows = T,
                        row_dend_width = unit(0.5, "cm"),
                        # add range to color scale
                        col = circlize::colorRamp2(c(-0.2, 0, 0.2), c("blue", "white", "red")),
                        km = 4
)
dev.off()