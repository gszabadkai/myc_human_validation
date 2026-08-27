library(tidyverse)
library(mgcv)
library(viridis)
mtPPS_thrsh <- 1.1
pct_thrsh <- 0.3
pathways_below_threshold <- readRDS(here::here("Data/Guo_et_al", "ProcessedData", "mitoPPS_bins.rds")) %>%
  group_by(Pathway1) %>%
  mutate(above_threshold = ifelse(mitoPPS > 1.1, TRUE, FALSE)) %>%
  # mutate(n = n()) %>%
  group_by(Pathway1, above_threshold) %>%
  mutate(n_bins_enriched = n()) %>%
  filter(above_threshold == TRUE) %>%
  ungroup() %>%
  dplyr::select(Pathway1, n_bins_enriched) %>%
  unique()  %>%
  dplyr::filter(n_bins_enriched <10) %>%
  pull(Pathway1)

delta_filter <- readRDS(here::here("Data/Guo_et_al", "ProcessedData", "mitoPPS_bins.rds")) %>%
  group_by(Pathway1) %>%
  mutate(delta = max(mitoPPS) - min(mitoPPS)) %>%
  filter(delta <0.4) %>%
  select(Pathway1, delta) %>%
  unique() %>%
  pull(Pathway1)

mitoPPS <- readRDS(here::here("Data/Guo_et_al", "ProcessedData", "mitoPPS_bins.rds")) %>%
  filter(!Pathway1 %in% pathways_below_threshold) %>%
  filter(!Pathway1 %in% delta_filter) 
meta_binned <- read_csv(here::here("Data/Guo_et_al", "ProcessedData", "meta_binned.csv")) 


non_oxphos_top <- mitoPPS %>%
  #filter(!Pathway1 %in% oxphos) %>%
  inner_join(meta_binned, by = "bin") %>%
  ungroup() %>%
  mutate(n_up = case_when(mitoPPS > mtPPS_thrsh ~ "up", TRUE ~ NA)) %>%
  
  mutate(Phase = case_when(
    dominant_cell_type %in% c("SSC", "Diff_spermatogonia", "Early_prim_spermatocytes", "Late_prim_spermatocytes") ~ "Early",
    dominant_cell_type %in% c("Round_spermatids", "Elongated_spermatids") ~ "Mid",
    dominant_cell_type %in% "Late_spermatids" ~ "Late",
    TRUE ~ "Mid"
  )) %>%
  group_by(Pathway1, Phase) %>%
  mutate(n_bins = n()) %>%
  filter(n_up %in% "up") %>%
  group_by(Pathway1, Phase, n_up) %>%
  mutate(n_bins_up = n()) %>%
  mutate(n_bins_above_threshold = n_bins_up/n_bins) %>%
  ungroup() %>%
  select(Phase, Pathway1, n_bins_above_threshold) %>%
  unique() %>%
  filter(n_bins_above_threshold > pct_thrsh)

pathways_early <- non_oxphos_top %>%
  filter(Phase %in% "Early") %>%
  pull(Pathway1)

pathways_mid <- non_oxphos_top %>%
  filter(Phase %in% "Mid") %>%
  pull(Pathway1)

pathways_late <- non_oxphos_top %>%
  filter(Phase %in% "Late") %>%
  pull(Pathway1)

data <-   mitoPPS %>%
  filter(Pathway1 %in% pathways_early | Pathway1 %in% pathways_mid | Pathway1 %in% pathways_late) %>%
  mutate(Phase = case_when(Pathway1 %in% pathways_early ~ "Early",
                           Pathway1 %in% pathways_mid ~ "Mid",
                           Pathway1 %in% pathways_late ~ "Late",
                           TRUE ~ "Other")) %>%
  inner_join(meta_binned, by = "bin")
  


colors <- c(
  "SSC" = "#A98FFF",
  "Diff_spermatogonia" = "#F87169",
  "Early_prim_spermatocytes" = "#5AB711",
  "Late_prim_spermatocytes" = "#FB62D8",
  "Round_spermatids" = "#01B6EB",
  "Elongated_spermatids" = "#1AC49B",
  "Late_spermatids" = "#C49B02"
  
)

library(ComplexHeatmap)

df_sub <- data %>%
  filter(Pathway1 %in% pathways_early) %>%
  select(bin, Pathway1, mitoPPS) %>%
  full_join(meta_binned, by = "bin") %>%
  pivot_wider(names_from = Pathway1, values_from = mitoPPS) 
exprs <- t(df_sub[,-c(1:3)])
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

hm1 <- ComplexHeatmap::Heatmap(log10(exprs),
                               row_names_gp = gpar(fontsize = 10),
                               cluster_columns = F,
                               top_annotation = col_anno,
                               clustering_distance_rows = "euclidean",
                               clustering_method_rows = "ward.D2",
                               name = "log10(mitoPPS)",
                                col = circlize::colorRamp2(c(-0.2, 0, 0.2), c("blue", "white", "red")),
                               row_dend_width = unit(0.5, "cm")
)



df_sub <- data %>%
  filter(Pathway1 %in% pathways_mid) %>%
  filter(!Pathway1 %in% pathways_early) %>%
  filter(!Pathway1 %in% pathways_late) %>%
  select(bin, Pathway1, mitoPPS) %>%
  full_join(meta_binned, by = "bin") %>%
  pivot_wider(names_from = Pathway1, values_from = mitoPPS) 
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


hm2 <- ComplexHeatmap::Heatmap(log10(exprs),
                               row_names_gp = gpar(fontsize = 10),
                               cluster_columns = F,
                               clustering_distance_rows = "euclidean",
                               clustering_method_rows = "ward.D2",
                               col = circlize::colorRamp2(c(-0.2, 0, 0.2), c("blue", "white", "red")),
                               show_heatmap_legend = FALSE,
                               # top_annotation = col_anno,
                               row_dend_width = unit(0.5, "cm")
)


df_sub <- data %>%
  filter(Pathway1 %in% pathways_late) %>%
  filter(!Pathway1 %in% pathways_early) %>%
  filter(!Pathway1 %in% pathways_mid) %>%
  select(bin, Pathway1, mitoPPS) %>%
  full_join(meta_binned, by = "bin") %>%
  pivot_wider(names_from = Pathway1, values_from = mitoPPS) 
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


hm3 <- ComplexHeatmap::Heatmap(log10(exprs),
                               row_names_gp = gpar(fontsize = 10),
                               cluster_columns = F,
                               clustering_distance_rows = "euclidean",
                               clustering_method_rows = "ward.D2",
                               col = circlize::colorRamp2(c(-0.2, 0, 0.2), c("blue", "white", "red")),
                               show_heatmap_legend = FALSE,
                               # top_annotation = col_anno,
                               row_dend_width = unit(0.5, "cm")
)


ht_list = hm1 %v% hm2 %v% hm3
pdf(here::here("Figures/Figure6", "Suppl_FigS15_mitoPPS_dynamics_heatmap.pdf"), width = 10, height = 8)
draw(ht_list)
dev.off()




# Figure6F_Top pathways ---------------------------------------------------

theme_dynamics <-   theme(
  axis.title = element_text(size = 6),
  axis.text = element_text(size = 5),
  axis.line = element_line(linewidth = 0.3),
  legend.text = element_text(size = 5),
  legend.title = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position = "none")

p1 <- data %>%
  filter(Pathway1 %in% c("Respirasome assembly", "Heme-containing proteins", "Cytochrome C")) %>%
  ggplot(aes(x = bin, y= mitoPPS, color = Phase,fill = Phase, group = Pathway1)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 15),  linewidth =0.3) +
  geom_point(aes(y = mitoPPS), alpha = 0.3, size = 0.1, shape = 16) +
  scale_color_manual(values = c(
    "Early" = "#D7E05B",
    "Mid" = "#4C9388",
    "Late" = "#3F1559"
  )) +
  scale_fill_manual(values = c(
    "Early" = "#D7E05B",
    "Mid" = "#4C9388",
    "Late" = "#3F1559"
  )) +
  labs(y = "mitoPPS", x = "Bin") +
  theme_classic() +
  theme_dynamics 
p1
ggsave(here::here("Figures/Figure6", "Figure6F_Early_top_pathways_dynamics.png"), width =2, height = 1.5, unit = "in", dpi = 1200)
plotly::ggplotly(p1)


p1 <- data %>%
  filter(Pathway1 %in% c("Propanoate metabolism", "Vitamin B12 metabolism")) %>%
  ggplot(aes(x = bin, y= mitoPPS, color = Phase,fill = Phase, group = Pathway1)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 15),  linewidth =0.3) +
  geom_point(aes(y = mitoPPS), alpha = 0.3, size = 0.1, shape = 16) +
  scale_color_manual(values = c(
    "Early" = "#D7E05B",
    "Mid" = "#4C9388",
    "Late" = "#3F1559"
  )) +
  scale_fill_manual(values = c(
    "Early" = "#D7E05B",
    "Mid" = "#4C9388",
    "Late" = "#3F1559"
  )) +
  labs(y = "mitoPPS", x = "Bin") +
  theme_classic() +
  theme_dynamics 
p1
ggsave(here::here("Figures/Figure6", "Figure6F_Mid_top_pathways_dynamics.png"), width =2, height = 1.5, unit = "in", dpi = 1200)
plotly::ggplotly(p1)


p1 <- data %>%
  filter(Pathway1 %in% c("MICOS complex", "Gluconeogenesis", "Glycolysis")) %>%
  ggplot(aes(x = bin, y= mitoPPS, color = Phase,fill = Phase, group = Pathway1)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 15),  linewidth =0.3) +
  geom_point(aes(y = mitoPPS), alpha = 0.3, size = 0.1, shape = 16) +
  scale_color_manual(values = c(
    "Early" = "#D7E05B",
    "Mid" = "#4C9388",
    "Late" = "#3F1559"
  )) +
  scale_fill_manual(values = c(
    "Early" = "#D7E05B",
    "Mid" = "#4C9388",
    "Late" = "#3F1559"
  )) +
  labs(y = "mitoPPS", x = "Bin") +
  theme_classic() +
  theme_dynamics 
p1
ggsave(here::here("Figures/Figure6", "Figure6F_Late_top_pathways_dynamics.png"), width =2, height = 1.5, unit = "in", dpi = 1200)
plotly::ggplotly(p1)



# Fission/Fusion ----------------------------------------------------------

p1 <- readRDS(here::here("Data/Guo_et_al", "ProcessedData", "mitoPPS_bins.rds"))  %>%
  #filter(!Pathway1 %in% oxphos) %>%
  inner_join(meta_binned, by = "bin") %>%
  filter(Pathway1 %in% c("Fission", "Fusion")) %>%
  ggplot(aes(x = bin, y= mitoPPS, color = Pathway1,fill = Pathway1, group = Pathway1)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 15),  linewidth =0.3) +
  geom_point(aes(y = mitoPPS), alpha = 0.3, size = 0.1, shape = 16) +
  scale_color_manual(values = c(
    "Fusion" = "#D7E05B",
    "Fission" = "#3F1559"
  )) +
  scale_fill_manual(values = c(
    "Fusion" = "#D7E05B",
    "Fission" = "#3F1559"
  )) +
  labs(y = "mitoPPS", x = "Bin") +
  theme_classic() +
  theme_dynamics 
p1
ggsave(here::here("Figures/Figure6", "Figure6I_Fission_Fusion.png"), width =3.2, height = 1.6, unit = "in", dpi = 1200)
plotly::ggplotly(p1)



p1 <- readRDS(here::here("Data/Guo_et_al", "ProcessedData", "mitoPPS_bins.rds"))  %>%
  #filter(!Pathway1 %in% oxphos) %>%
  inner_join(meta_binned, by = "bin") %>%
  filter(Pathway1 %in% c("OXPHOS subunits", "Glycolysis")) %>%
  ggplot(aes(x = bin, y= mitoPPS, color = Pathway1,fill = Pathway1, group = Pathway1)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 15),  linewidth =0.3) +
  geom_point(aes(y = mitoPPS), alpha = 0.3, size = 0.1, shape = 16) +
  scale_color_manual(values = c(
    "OXPHOS subunits" = "#D7E05B",
    "Glycolysis" = "#3F1559"
  )) +
  scale_fill_manual(values = c(
    "OXPHOS subunits" = "#D7E05B",
    "Glycolysis" = "#3F1559"
  )) +
  labs(y = "mitoPPS", x = "Bin") +
  theme_classic() +
  theme_dynamics 
p1
ggsave(here::here("Figures/Figure6", "Suppl_FigS15B_OxGlyc.png"), width =3.2, height = 1.6, unit = "in", dpi = 1200)

p1 <- readRDS(here::here("Data/Guo_et_al", "ProcessedData", "mitoPPS_bins.rds"))  %>%
  #filter(!Pathway1 %in% oxphos) %>%
  inner_join(meta_binned, by = "bin") %>%
  filter(Pathway1 %in% c("Fission", "Cristae formation")) %>%
  ggplot(aes(x = bin, y= mitoPPS, color = Pathway1,fill = Pathway1, group = Pathway1)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 15),  linewidth =0.3) +
  geom_point(aes(y = mitoPPS), alpha = 0.3, size = 0.1, shape = 16) +
  scale_color_manual(values = c(
    "Cristae formation" = "#D7E05B",
    "Fission" = "#3F1559"
  )) +
  scale_fill_manual(values = c(
    "Cristae formation" = "#D7E05B",
    "Fission" = "#3F1559"
  )) +
  labs(y = "mitoPPS", x = "Bin") +
  theme_classic() +
  theme_dynamics 
p1

p1 <- readRDS(here::here("Data/Guo_et_al", "ProcessedData", "mitoPPS_bins.rds"))  %>%
  #filter(!Pathway1 %in% oxphos) %>%
  inner_join(meta_binned, by = "bin") %>%
  filter(Pathway1 %in% c("Gluconeogenesis", "Glycolysis")) %>%
  ggplot(aes(x = bin, y= mitoPPS, color = Pathway1,fill = Pathway1, group = Pathway1)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 15),  linewidth =0.3) +
  geom_point(aes(y = mitoPPS), alpha = 0.3, size = 0.1, shape = 16) +
  scale_color_manual(values = c(
    "Gluconeogenesis" = "#4C9388",
    "Glycolysis" = "#3F1559"
  )) +
  scale_fill_manual(values = c(
    "Gluconeogenesis" = "#4C9388",
    "Glycolysis" = "#3F1559"
  )) +
  labs(y = "mitoPPS", x = "Bin") +
  theme_classic() +
  theme_dynamics 
p1
ggsave(here::here("Figures/Figure6", "Suppl_FigS15B_GlucNeoGlyc.png"), width =3.2, height = 1.6, unit = "in", dpi = 1200)
plotly::ggplotly(p1)

p1 <- readRDS(here::here("Data/Guo_et_al", "ProcessedData", "mitoPPS_bins.rds"))  %>%
  #filter(!Pathway1 %in% oxphos) %>%
  inner_join(meta_binned, by = "bin") %>%
  filter(Pathway1 %in% c("MICOS complex", "Cristae formation")) %>%
  ggplot(aes(x = bin, y= mitoPPS, color = Pathway1,fill = Pathway1, group = Pathway1)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 15),  linewidth =0.3) +
  geom_point(aes(y = mitoPPS), alpha = 0.3, size = 0.1, shape = 16) +
  scale_color_manual(values = c(
    "Cristae formation" = "#D7E05B",
    "MICOS complex" = "#3F1559"
  )) +
  scale_fill_manual(values = c(
    "Cristae formation" = "#D7E05B",
    "MICOS complex" = "#3F1559"
  )) +
  labs(y = "mitoPPS", x = "Bin") +
  theme_classic() +
  theme_dynamics 
p1
ggsave(here::here("Figures/Figure6", "Suppl_FigS15B_MICOS_Crist.png"), width =3.2, height = 1.6, unit = "in", dpi = 1200)
# Example pathways early --------------------------------------------------
theme_dynamics <-   theme(
  axis.title = element_text(size = 7),
  axis.text = element_text(size = 6),
  axis.line = element_line(linewidth = 0.3),
  legend.text = element_text(size = 6),
  legend.title = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position = "none", 
  #legend.box = "horizontal",
  title = element_text(size = 7))
p1 <- data %>%
  filter(Pathway1 %in% c("TOM", "TIM23 presequence pathway", "TIM22 carrier pathway", "SAM")) %>%
  ggplot(aes(x = bin, y= mitoPPS, color = Pathway1,fill = Pathway1, group = Pathway1)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 15),  linewidth =0.3) +
  geom_point(aes(y = mitoPPS), alpha = 0.3, size = 0.1, shape = 16) +
  labs(y = "mitoPPS", x = "Bin", title = "Protein transport and sorting systems") +
  theme_classic() +
  theme_dynamics 
p1
ggsave(here::here("Figures/Figure6/Pathway_dynamics", "Protein_import_sorting_dynamics.png"), width =2.8, height = 2, unit = "in", dpi = 1200)
plotly::ggplotly(p1)

p1 <- data %>%
  filter(Pathway1 %in% c("Cristae formation", "Mitochondrial permeability transition pore", "Calcium cycle")) %>%
  ggplot(aes(x = bin, y= mitoPPS, color = Pathway1,fill = Pathway1, group = Pathway1)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 15),  linewidth =0.3) +
  geom_point(aes(y = mitoPPS), alpha = 0.3, size = 0.1, shape = 16) +
  labs(y = "mitoPPS", x = "Bin", title = "Membrane architecture and organization") +
  theme_classic() +
  theme_dynamics 
p1
ggsave(here::here("Figures/Figure6/Pathway_dynamics", "Membrane_architecture_and_organization_dynamics.png"), width =2.8, height = 2, unit = "in", dpi = 1200)
plotly::ggplotly(p1)

p1 <- data %>%
  filter(Pathway1 %in% c("Transcription", "mtDNA modifications", "Chaperones")) %>%
  ggplot(aes(x = bin, y= mitoPPS, color = Pathway1,fill = Pathway1, group = Pathway1)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 15),  linewidth =0.3) +
  geom_point(aes(y = mitoPPS), alpha = 0.3, size = 0.1, shape = 16) +
  labs(y = "mitoPPS", x = "Bin", title = "Mitochondrial central dogma") +
  theme_classic() +
  theme_dynamics 
p1
ggsave(here::here("Figures/Figure6/Pathway_dynamics", "central_dogma_dynamics.png"), width =2.8, height = 2, unit = "in", dpi = 1200)
plotly::ggplotly(p1)

