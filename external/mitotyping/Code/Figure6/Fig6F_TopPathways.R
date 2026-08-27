library(tidyverse)
library(mgcv)

pathways_below_threshold <- readRDS(here::here("Data/Guo_et_al", "ProcessedData", "mitoPPS_bins.rds")) %>%
  group_by(Pathway1) %>%
  mutate(above_threshold = ifelse(mitoPPS > 1.1, TRUE, FALSE)) %>%
 # mutate(n = n()) %>%
  group_by(Pathway1, above_threshold) %>%
  mutate(n_bins_enriched = n()) %>%
  filter(above_threshold == TRUE) %>%
  filter(n_bins_enriched <3) %>%
  ungroup() %>%
  select(Pathway1, n_bins_enriched) %>%
  unique() %>%
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




# oxphos <- readxl::read_xls(here::here("Data/MitoCarta/OriginalData/HumanMitoCarta3_0.xls"), sheet = 4) %>%
#   select(MitoPathway, `MitoPathways Hierarchy`) %>%
#   separate(`MitoPathways Hierarchy`, into = c("Pathway_Level1", "Pathway_Level2", "Pathway_Level3"), sep = " > ") %>%
#   mutate(Level = case_when(
#     MitoPathway == Pathway_Level1 ~ "Pathway_Level1",
#     MitoPathway == Pathway_Level2 ~ "Pathway_Level2",
#     MitoPathway == Pathway_Level3 ~ "Pathway_Level3"
#   )) %>%
#   unique() %>%
#   dplyr::rename(Pathway = MitoPathway) %>%
#   filter(Pathway_Level1 %in% "OXPHOS") %>%
#   pull(Pathway)

non_oxphos_top <- mitoPPS %>%
 #filter(!Pathway1 %in% oxphos) %>%
  inner_join(meta_binned, by = "bin") %>%
  ungroup() %>%
  mutate(n_up = case_when(mitoPPS > 1.1 ~ "up", TRUE ~ NA)) %>%
  
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
  filter(n_bins_above_threshold > 0.3)

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
                           TRUE ~ "Other"))
all_pathways <- c(pathways_early, pathways_mid, pathways_late)
for (i in 1:length(all_pathways)) {
  phase <- data %>%
    filter(Pathway1 == all_pathways[i]) %>%
    pull(Phase) %>%
    unique()
  data %>%
    filter(Pathway1 == all_pathways[i]) %>%
    ggplot(aes(x = bin, y= mitoPPS, color = Phase, group = Pathway1, shape = Pathway1)) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 15),  linewidth =0.5) +
    geom_point(aes(y = mitoPPS), alpha = 0.4, size = 0.2) +
    scale_color_manual(values = c(
      "Early" = "#D7E05B",
      "Mid" = "#4C9388",
      "Late" = "#3F1559"
    )) +
    labs(y = "mitoPPS", x = "Pseudotime", title = all_pathways[i]) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 5),
      axis.line = element_line(linewidth = 0.4),
      legend.text = element_blank(),
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      title = element_text(size = 7),
      legend.position = "none")
  ggsave(here::here("Figures/Figure6/Pathway_dynamics", paste0(phase, "_", all_pathways[i], "_mitoPPS_trend.png")), width = 3, height = 2)
}

# library(viridisLite)
# library(circlize)
# library(ComplexHeatmap)

  colors <- c(
    "SSC" = "#A98FFF",
    "Diff_spermatogonia" = "#F87169",
    "Early_prim_spermatocytes" = "#5AB711",
    "Late_prim_spermatocytes" = "#FB62D8",
    "Round_spermatids" = "#01B6EB",
    "Elongated_spermatids" = "#1AC49B",
    "Late_spermatids" = "#C49B02"
    
  )

df <- data %>% 
    select(bin, Pathway1, mitoPPS) %>%
    full_join(meta_binned, by = "bin") %>%
    pivot_wider(names_from = Pathway1, values_from = mitoPPS) 
  
exprs <- t(df[,-c(1:3)])

rng <- range(df$mean_pt, na.rm = TRUE)
pal    <- rev(viridis(256)) 
breaks <- seq(rng[1], rng[2], length.out = length(pal))
col_fun <- circlize::colorRamp2(breaks, pal)


col_anno <- columnAnnotation(
  pseudotime = anno_simple(df$mean_pt, col = col_fun, border = FALSE),
  celltype   = df$dominant_cell_type,
  col        = list(celltype = colors)
)

#pdf(here::here("Figures", "Figure6", "Figure6E_mitoPPS_heatmap.pdf"), width = 7, height = 6)
ComplexHeatmap::Heatmap(log10(exprs),
                        row_names_gp = gpar(fontsize = 4),
                        cluster_columns = F,
                        top_annotation = col_anno,
                        cluster_rows = T,
                        row_dend_width = unit(0.5, "cm")
                        # add range to color scale
                    #    col = circlize::colorRamp2(c(-0.2, 0, 0.2), c("blue", "white", "red")),
                       # km = 3
)
#dev.off()
# 
     
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

hm1 <- ComplexHeatmap::Heatmap(exprs,
                               row_names_gp = gpar(fontsize = 10),
                               cluster_columns = F,
                               top_annotation = col_anno,
                              # col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
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


hm2 <- ComplexHeatmap::Heatmap(exprs,
                               row_names_gp = gpar(fontsize = 10),
                               cluster_columns = F,
                             #  col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
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


hm3 <- ComplexHeatmap::Heatmap(exprs,
                               row_names_gp = gpar(fontsize = 10),
                               cluster_columns = F,
                               #  col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                               # top_annotation = col_anno,
                               row_dend_width = unit(0.5, "cm")
)


ht_list = hm1 %v% hm2 %v% hm3
draw(ht_list)











pdf(here::here("Figures/Figure6", "Figure6H_Fission_Fusion_heatmap.pdf"), height = 5, width =8, colormodel = "srgb")
print(ht_list)
dev.off()




oxphos_max <- mitoPPS %>%
  filter(Pathway1 %in% oxphos) %>%
  inner_join(meta_binned, by = "bin") %>%
  group_by(Pathway1) %>%
  mutate(max = max(mitoPPS, na.rm = TRUE)) %>%
  mutate(Group = case_when(mitoPPS == max ~ dominant_cell_type, TRUE ~ NA_character_)) %>%
  filter(!is.na(Group)) %>%
  filter(mitoPPS > 1.2) %>%
  mutate(Phase = case_when(
    dominant_cell_type %in% c("SSC", "Diff_spermatogonia", "Early_prim_spermatocytes", "Late_prim_spermatocytes") ~ "Early",
    dominant_cell_type %in% c("Round_spermatids", "Elongated_spermatids") ~ "Mid",
    dominant_cell_type %in% "Late_spermatids" ~ "Late",
    TRUE ~ "Mid"
  ))


pathways_early_ox <- oxphos_max %>%
  filter(Phase %in% "Early") %>%
  pull(Pathway1)

pathways_mid_ox <- oxphos_max %>%
  filter(Phase %in% "Mid") %>%
  pull(Pathway1)

pathways_late_ox <- oxphos_max %>%
  filter(Phase %in% "Late") %>%
  pull(Pathway1)

nonoxphos_max <- mitoPPS %>%
  filter(!Pathway1 %in% oxphos) %>%
  inner_join(meta_binned, by = "bin") %>%
  group_by(Pathway1) %>%
  mutate(max = max(mitoPPS, na.rm = TRUE)) %>%
  mutate(Group = case_when(mitoPPS == max ~ dominant_cell_type, TRUE ~ NA_character_)) %>%
  filter(!is.na(Group)) %>%
  filter(mitoPPS > 1.2) %>%
  mutate(Phase = case_when(
    dominant_cell_type %in% c("SSC", "Diff_spermatogonia", "Early_prim_spermatocytes", "Late_prim_spermatocytes") ~ "Early",
    dominant_cell_type %in% c("Round_spermatids", "Elongated_spermatids") ~ "Mid",
    dominant_cell_type %in% "Late_spermatids" ~ "Late",
    TRUE ~ "Mid"
  ))



pathways_early <- nonoxphos_max %>%
  filter(Phase %in% "Early") %>%
  pull(Pathway1)

pathways_mid <- nonoxphos_max %>%
  filter(Phase %in% "Mid") %>%
  pull(Pathway1)

pathways_late <- nonoxphos_max %>%
  filter(Phase %in% "Late") %>%
  pull(Pathway1)

data <-   mitoPPS %>%
  filter(Pathway1 %in% pathways_early | Pathway1 %in% pathways_mid | Pathway1 %in% pathways_late) %>%
  mutate(Phase = case_when(Pathway1 %in% pathways_early ~ "Early",
                           Pathway1 %in% pathways_mid ~ "Mid",
                           Pathway1 %in% pathways_late ~ "Late",
                           TRUE ~ "Other")) %>%
  group_by(Pathway1) %>%
  #mutate(smooth = zoo::rollmean(mitoPPS, k = 20, fill = NA, align = "right")) %>%
  #mutate(residuals = smooth - mitoPPS) %>%
  mutate(
    gam_fit = predict(
      gam(mitoPPS ~ s(bin, bs = "cs", k = 15))
    )) %>%
  mutate(gam_resid = gam_fit - mitoPPS) 
ggplot(data, aes(x = gam_fit, y = gam_resid)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  labs(title = "Residuals vs Fitted Values", x = "Fitted", y = "Residuals")

ggplot(data, aes(x = bin, y = gam_resid)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Residuals vs Predictor", y = "Residuals")


test <- data %>%
  mutate(resid = abs(residuals_gam)) %>%
  group_by(Pathway1) %>%
  summarise(mean_resid = mean(resid, na.rm = TRUE), 
            sd_resid = sd(resid, na.rm = TRUE),
            cv_resid = sd_resid/mean_resid,
            max_resid = max(resid, na.rm = T)) %>%
  arrange(desc(cv_resid))
plotly::ggplotly(data %>%
                   unique() %>%
    filter(Pathway1 %in% pathways_early) %>%
    ggplot(aes(x = residuals_gam, color = Pathway1)) +
    geom_density())
data %>%
  filter(Pathway1 %in% "fMet processing") %>%
   ggplot(aes(x = residuals, y = residuals_gam, color = bin)) + 
  geom_point()


data %>%
  filter(Pathway1 %in% "fMet processing") %>%
  ggplot(aes(x = bin, y= mitoPPS, group = Pathway1, shape = Pathway1)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 15),  linewidth =0.5) +
  #geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 15),  linewidth =0.5) +
  geom_point(aes(y = mitoPPS), alpha = 0.4, size = 0.2) + 
  labs(y = "mitoPPS", x = "Pseudotime") + 
  theme_classic() +
  theme(
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 5),
    axis.line = element_line(linewidth = 0.4),
    legend.text = element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    title = element_text(size = 7),
    legend.position = "none")

data %>%
  filter(Pathway1 %in% "fMet processing") %>%
  ggplot(aes(x = smooth, y = gam_fit)) + 
 # ggplot(aes(x = residuals, y = residuals_gam)) + 
  geom_point()
  mutate(
    gam_fit = predict(
      gam(mitoPPS ~ s(bin, bs = "cs", k = 15), data = mitoPPS),
      newdata = mitoPPS
    ))
  inner_join(meta_binned, by = "bin") %>%
  filter(Pathway1 %in% pathways_early | Pathway1 %in% pathways_mid | Pathway1 %in% pathways_late) 
all_pathways <- c(pathways_early, pathways_mid, pathways_late)
for (i in 1:length(all_pathways)) {
  phase <- data %>%
    filter(Pathway1 == all_pathways[i]) %>%
    pull(Phase) %>%
    unique()
  data %>%
    filter(Pathway1 == all_pathways[i]) %>%
    ggplot(aes(x = bin, y= smooth, color = Phase, group = Pathway1, shape = Pathway1)) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 15),  linewidth =0.5) +
    geom_point(aes(y = mitoPPS), alpha = 0.4, size = 0.2) + 
    scale_color_manual(values = c(
      "Early" = "#D7E05B",
      "Mid" = "#4C9388",
      "Late" = "#3F1559"
    )) +
    labs(y = "mitoPPS", x = "Pseudotime", title = all_pathways[i]) + 
    theme_classic() +
    theme(
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 5),
      axis.line = element_line(linewidth = 0.4),
      legend.text = element_blank(),
      legend.title = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      title = element_text(size = 7),
      legend.position = "none")
  ggsave(here::here("Figures/Figure6/Pathway_dynamics", paste0(phase, "_", all_pathways[i], "_mitoPPS_trend.png")), width = 3, height = 2)
}


p <- data %>%
ggplot(aes(x = bin, y= smooth, color = Phase, group = Pathway1, shape = Pathway1)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 15),  linewidth =0.5) +
  scale_color_manual(values = c(
    "Early" = "#D7E05B",
    "Mid" = "#4C9388",
    "Late" = "#3F1559"
  )) +
  labs(y = "mitoPPS", x = "Pseudotime") + 
  theme_classic() +
  theme(
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 5),
    axis.line = element_line(linewidth = 0.4),
    legend.text = element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()) +
  facet_wrap(~fct_relevel(Phase, c("Early", "Mid", "Late")), ncol = 1, scales = "free_y")
p
plotly::ggplotly(p)

p <- data %>%
  filter(Pathway1 %in% c("Cytochromes", "Heme-containing proteins",  "PentosePhosphatePathway")) %>%
  #filter(Pathway1 %in% c( "Propanoate metabolism", "Vitamin B12 metabolism"))
#filter(Pathway1 %in% c("Glycolysis", "MICOS complex", "Selenoproteins")) %>%
  ggplot(aes(x = bin, y= smooth, color = Phase, group = Pathway1, shape = Pathway1)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 15),  linewidth =0.5) +
  scale_color_manual(values = c(
    "Early" = "#D7E05B",
    "Mid" = "#4C9388",
    "Late" = "#3F1559"
  )) +
  labs(y = "mitoPPS", x = "Pseudotime") + 
  theme_classic() +
  theme(
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 5),
    axis.line = element_line(linewidth = 0.4),
    legend.text = element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()) +
  facet_wrap(~fct_relevel(Phase, c("Early", "Mid", "Late")), ncol = 1, scales = "free_y")
p
plotly::ggplotly(p)
