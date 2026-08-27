library(Seurat)
library(tidyverse)
mitogenes <- readxl::read_xls(here::here("Data/MitoCarta/OriginalData/HumanMitoCarta3_0.xls"), sheet = 2) %>%
  pull(Symbol)

seurat <- readRDS(here::here("Data/Guo_et_al", "ProcessedData", "seurat_object.rds"))
mat <- GetAssayData(seurat, layer = "data")
meta <- seurat@meta.data

mitodata <- as.data.frame(mat) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Cell", values_to = "exprs") %>%
  filter(Gene %in% mitogenes) %>%
  full_join(meta %>% rownames_to_column("Cell"), by = "Cell")



# Binned ------------------------------------------------------------------
mitodata %>%
  mutate(Genome = case_when(
    grepl("MT-", Gene) ~ "mtDNA",
    TRUE ~ "nucDNA")) %>%
  na.omit() %>%
  group_by(Genome) %>%
  select(Genome, Gene) %>%
  unique() %>%
  mutate(n = n()) %>%
  select(Genome, n) %>%
  unique()
mt_nuc_binned <- mitodata %>%
  mutate(Genome = case_when(
    grepl("MT-", Gene) ~ "mtDNA",
    TRUE ~ "nucDNA")) %>%
  na.omit() %>%
  group_by(Genome, Cell, Pseudotime) %>%
  mutate(mean = mean(exprs)) %>%
  select(Genome, Cell, Pseudotime, mean) %>%
  unique() %>%
  inner_join(read_csv(here::here("Data/Guo_et_al", "ProcessedData", "cell_to_bin.csv")) %>%
              mutate(Cell = as.character(Cell)), by = "Cell") %>%
  select(Genome, bin, mean_pt, mean) %>%
  unique() %>%
  group_by(Genome) %>%
  arrange(bin) %>%
  mutate(
    smooth = zoo::rollmean(mean, k = 80, fill = NA, align = "right")
  ) %>%
  filter(!is.na(smooth)) %>%
  mutate(smooth = as.numeric(scale(smooth))) %>%
  ungroup() %>%
  unique()  #%>%
  # group_by(Genome) %>%
  # mutate(mean = as.numeric(scale(mean)))
colors <- c(
  "SSC" = "#A98FFF",
  "Diff_spermatogonia" = "#F87169",
  "Early_prim_spermatocytes" = "#5AB711",
  "Late_prim_spermatocytes" = "#FB62D8",
  "Round_spermatids" = "#01B6EB",
  "Elongated_spermatids" = "#1AC49B",
  "Late_spermatids" = "#C49B02"
)

v_lines <- mt_nuc_binned %>%
  inner_join(read_csv(here::here("Data/Guo_et_al", "ProcessedData", "meta_binned.csv"))) %>%
  select(bin, dominant_cell_type) %>%
  unique() %>%
  group_by(dominant_cell_type) %>%
  mutate(max = max(bin)) %>%
  select(dominant_cell_type, max) %>%
  unique() %>%
  arrange(max) %>%
  inner_join(data.frame(colors) %>%
               rownames_to_column("dominant_cell_type"), by = "dominant_cell_type")

mt_nuc_binned %>%
  ggplot(aes(x = bin, y = smooth, color = Genome, group = Genome)) +
  geom_line(alpha =0.9, linewidth = 0.3) +
 # geom_point(alpha = 0.2, size = 0.2) +
  theme_classic() +
  #geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 31), linewidth = 0.5, alpha = 0.3) +
  scale_color_manual(values = c("nucDNA" = "#474B87", "mtDNA" = "#AE525E"))+
  ggnewscale::new_scale_color() +
  geom_vline(data = v_lines, aes(xintercept = max, color = dominant_cell_type),
             inherit.aes = F, linewidth = 0.3, linetype = "dotted") +

  scale_color_manual(values = colors, guide = "none") +
 # geom_vline(xintercept = 27, linetype = "dashed", color = "grey") +
  labs(x = "Binned pseudotime", y = "Mean expression") +
  theme(legend.position = "top") +
  theme(axis.text = element_text(size = 5),
       axis.title = element_text(size = 6),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       legend.position = "none",
       axis.line = element_line(size = 0.3)
  )
ggsave(here::here("Figures", "Figure6",  "Figure6C_mtDNA_nucDNA_pt.png"), width = 4, height = 2, units = "in", dpi = 1200)




# Actual pseudotime -------------------------------------------------------


mt_nuc_pseudotime <- mitodata %>%
  mutate(Genome = case_when(
    grepl("MT-", Gene) ~ "mtDNA",
    TRUE ~ "nucDNA")) %>%
  na.omit() %>%
  group_by(Genome, Cell, Pseudotime) %>%
  mutate(mean = mean(exprs)) %>%
  select(Genome, Cell, Pseudotime, mean) %>%
  unique()  %>%
  group_by(Genome) %>%
  # arrange(Pseudotime) %>%
  # mutate(
  #   smooth = zoo::rollmean(mean, k = 100, fill = NA, align = "right")
  # ) %>%
  # filter(!is.na(smooth)) %>%
  # mutate(smooth = as.numeric(scale(smooth))) %>%
  # ungroup() %>%
  # unique() #%>%
  group_by(Genome) %>%
  mutate(mean = scale(mean))

colors <- c(
  "SSC" = "#A98FFF",
  "Diff_spermatogonia" = "#F87169",
  "Early_prim_spermatocytes" = "#5AB711",
  "Late_prim_spermatocytes" = "#FB62D8",
  "Round_spermatids" = "#01B6EB",
  "Elongated_spermatids" = "#1AC49B",
  "Late_spermatids" = "#C49B02"
  
)


as.data.frame(colors) %>%
  rownames_to_column("CellSubtype")
v_lines <- mt_nuc_pseudotime %>%
  inner_join(readRDS(here::here("Data/Guo_et_al", "ProcessedData", "seurat_meta.rds")) %>%
               rownames_to_column("Cell")) %>%
  select(Cell,Pseudotime, CellSubtype) %>%
  unique() %>%
  group_by(CellSubtype) %>%
  mutate(max = max(Pseudotime)) %>%
  select(CellSubtype, max) %>%
  unique() %>%
  arrange(max) %>%
  mutate(max = case_when(CellSubtype == "SSC" ~ 4.43765148, TRUE ~ max)) %>% ## SSC has 2 outlier cells with high pseudotime
  inner_join(data.frame(colors) %>%
               rownames_to_column("CellSubtype"), by = "CellSubtype")

mt_nuc_pseudotime %>%
  inner_join(readRDS(here::here("Data/Guo_et_al", "ProcessedData", "seurat_meta.rds")) %>%
               rownames_to_column("Cell")) %>%
  select(Cell,Pseudotime, CellSubtype) %>%
  unique() %>%
  filter(CellSubtype %in% "SSC") %>%
  arrange(Pseudotime) %>%
  pull(Pseudotime)


mt_nuc_pseudotime %>%
  ggplot(aes(x = Pseudotime, y = mean, color = Genome, group = Genome)) +
  theme_classic() +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 31), linewidth = 0.5, alpha = 0.3) + 
  scale_color_manual(values = c("nucDNA" = "#474B87", "mtDNA" = "#AE525E"))+
  ggnewscale::new_scale_color() +
  geom_vline(data = v_lines, aes(xintercept = max, color = CellSubtype), 
             inherit.aes = F, linewidth = 0.5, linetype = "dotted") +
  
  scale_color_manual(values = colors, guide = "none") +
  # geom_vline(xintercept = 27, linetype = "dashed", color = "grey") +
  labs(x = "Pseudotime", y = "Mean expression") +
  theme(legend.position = "top") + 
  theme(axis.text = element_text(size = 5),
        axis.title = element_text(size = 6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 0.5)
  )





# bin to pseudotime -------------------------------------------------------


meta <- readRDS(here::here("Data/Guo_et_al", "ProcessedData", "seurat_meta.rds"))  %>%
  rownames_to_column("Cell") %>%
  select(Cell, Pseudotime)
cell_to_bin <- read_csv(here::here("Data/Guo_et_al", "ProcessedData", "cell_to_bin.csv")) %>%
              mutate(Cell = as.character(Cell))
pseudotime_to_bin <- meta %>% inner_join(cell_to_bin, by = "Cell") %>%
  select(Pseudotime, bin) %>%
  arrange(Pseudotime)
library(viridisLite)
rng <- range(pseudotime_to_bin$Pseudotime, na.rm = TRUE)
pal    <- rev(viridis(256))  # "yellow -> green -> teal -> blue -> purple"
breaks <- seq(rng[1], rng[2], length.out = length(pal))
col_fun <- circlize::colorRamp2(breaks, pal)

df_with_colors <- pseudotime_to_bin %>%
  mutate(color = col_fun(Pseudotime))

test <- df_with_colors %>%
  arrange(Pseudotime) %>%
  mutate(round = round(Pseudotime, 0)) %>%
  filter(round %in% c(10, 20, 30, 40, 50)) %>%
  select(round, bin, color)
