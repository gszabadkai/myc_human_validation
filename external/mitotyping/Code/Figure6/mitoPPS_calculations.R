
library(tidyverse)
library(ComplexHeatmap)
library(Seurat)
# library(SeuratDisk)
# library(SeuratData)
library(broom)
library(reshape2)
library(mgcv)
library(abind)
library(data.table)


# Load the data -----------------------------------------------------------

seurat <- readRDS(here::here("Data/Guo_et_al", "ProcessedData", "seurat_object.rds"))
umap_coordinates <- read_csv(here::here("Data/Guo_et_al", "ProcessedData", "seurat_umap.csv")) %>%
  mutate(Cell = as.character(Cell))
meta <- as.data.frame(readRDS(here::here("Data/Guo_et_al", "ProcessedData","seurat_meta.rds"))) %>%
  rownames_to_column("Cell")

gene_sets <- readRDS(here::here("Data/Guo_et_al", "ProcessedData", "gene_sets.rds"))
gene_to_pathway <- stack(gene_sets) %>%
  filter(!grepl("random", ind)) %>%
  dplyr::rename(Gene = values, Pathway = ind)

DefaultAssay(seurat) <- "RNA" 
seurat <- Seurat::AddModuleScore(
  seurat,
  features = gene_sets,
  name     = "PW"               
)

newcols <- paste0("PW", seq_along(gene_sets))
names(seurat@meta.data)[match(newcols, names(seurat@meta.data))] <-
  paste0(names(gene_sets))

# -------------------------------------------------------------------------


data_pathways <- seurat@meta.data %>%
  rownames_to_column("Cell") %>%
  pivot_longer(cols = gene_sets %>% names(), names_to = "Pathway", values_to = "score") %>%
  select(Cell, Pathway, score) 

min <- min(data_pathways$score, na.rm = TRUE)
data_pathways <- seurat@meta.data %>%
  rownames_to_column("Cell") %>%
  pivot_longer(cols = gene_sets %>% names(), names_to = "Pathway", values_to = "score") %>%
  select(Cell, Pathway, score) %>%
  mutate(score_new = score - min + 1e-6) # make all scores positive (note: sd stays, mean changes)

## compare old to new score
data_pathways %>%
  ggplot(aes(x = score, y = score_new) ) +
  geom_point()

data_pathways <- seurat@meta.data %>%
  rownames_to_column("Cell") %>%
  pivot_longer(cols = gene_sets %>% names(), names_to = "Pathway", values_to = "score") %>%
  filter(Pathway %in% gene_to_pathway$Pathway) %>%
  select(Cell, Pathway, score) %>%
  mutate(score = score - min + 1e-6) %>%
  inner_join(read_csv(here::here("Data/Guo_et_al", "ProcessedData", "cell_to_bin.csv")) %>%
               mutate(Cell = as.character(Cell)), by = "Cell") %>%
  group_by(bin, Pathway) %>%
  summarize(score = mean(score, na.rm = TRUE), .groups = "drop")

# Pathway score - based mitoPPS -------------------------------------------

data_ratios <- data_pathways %>%
  group_by(bin) %>%
  nest() %>%
  mutate(pairs = map(data, ~ {
    df <- .x
    expand.grid(
      Pathway1 = df$Pathway,
      Pathway2 = df$Pathway,
      stringsAsFactors = FALSE
    ) %>%
      filter(Pathway1 != Pathway2) %>%
      left_join(df, by = c("Pathway1" = "Pathway")) %>%
      dplyr::rename(value1 = score) %>%
      left_join(df, by = c("Pathway2" = "Pathway")) %>%
      dplyr::rename(value2 = score) %>%
      mutate(ratio = value1 / value2)
  })) %>%
  dplyr::select(-data) %>%
  unnest(pairs)

data_ratios2 <- data_ratios %>%
  group_by(Pathway1, Pathway2) %>%
  mutate(average_ratio = mean(ratio, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(corrected = ratio / average_ratio) 

mitoPPS  <- data_ratios2 %>%
  filter(!is.na(corrected)) %>%
  group_by(bin, Pathway1) %>%
  summarize(mitoPPS = mean(corrected, na.rm = TRUE), .groups = "drop") 

saveRDS(mitoPPS, here::here("Data/Guo_et_al", "ProcessedData", "mitoPPS_bins.rds"))
