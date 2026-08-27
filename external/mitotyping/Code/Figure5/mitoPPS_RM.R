library(tidyverse)

mitogenes <- readxl::read_xls(here::here("Data/MitoCarta/OriginalData/HumanMitoCarta3_0.xls"), sheet = 2) %>%
  pull(Symbol)

data <- read.csv(here::here("Data", "MBM_ROSMAP", "OriginalData", 
                            "unified-cpm.log2cpm.exclude-Broad.rename-gene-and-donor.csv")) 
IDs <- data.frame(group = colnames(data[,3:ncol(data)])) %>%
  separate(group, into = c("Celltype", "ID")) %>%
  select(ID) %>%
  unique()

data <- data %>%
  select(-gene_id) %>%
  pivot_longer(
    cols = !gene_symbol, 
    names_to = "Celltype", 
    values_to = "exprs"
  )%>%
  separate(Celltype, into = c("Celltype", "ID")) %>%
  filter(Celltype %in% c("Inh", "Exc", "OPC", "Oli", "Ast", "Mic", "Peri", "End")) %>%
  filter(gene_symbol %in% mitogenes) %>%
  pivot_wider(names_from = gene_symbol, values_from = exprs) 

write_csv(data, here::here("Data", "MBM_ROSMAP", "ProcessedData", "ROSMAP_mitogenes.csv"))
data <- read_csv(here::here("Data", "MBM_ROSMAP", "ProcessedData", "ROSMAP_mitogenes.csv"))
genes_found <- colnames(data)[-c(1:2)]
sort(mitogenes[which(!mitogenes %in% genes_found)])
## Mitocarta
mitocarta_sheet4 <- readxl::read_xls(here::here("Data/MitoCarta/OriginalData/HumanMitoCarta3_0.xls"), sheet = 4) %>%
  dplyr::select(MitoPathway, Genes) %>%
  na.omit() 
gene_to_pathway <- splitstackshape::cSplit(mitocarta_sheet4, 'Genes', ',') %>%
  column_to_rownames("MitoPathway") %>%
  t() %>%
  as.data.frame() 
gene_to_pathway <- gene_to_pathway %>%
  pivot_longer(cols = colnames(gene_to_pathway), names_to = "Pathway", values_to = "Gene") %>%
  na.omit() 

data_pathways <- read_csv(here::here("Data", "MBM_ROSMAP", "ProcessedData", "ROSMAP_mitogenes.csv")) %>%
#  dplyr::rename(Gene = gene_symbol) %>%
  pivot_longer(cols = -c("ID", "Celltype"), names_to = "Gene") %>%
  full_join(gene_to_pathway) %>%
  filter(!is.na(Pathway)) %>%
  filter(!is.na(Celltype)) %>%
  filter(!is.na(ID)) %>%
  group_by(Celltype, ID, Pathway) %>%
  summarize(score = mean(2^value, na.rm = TRUE), .groups = "drop")

# MitoPPS -----------------------------------------------------------------


all_pathways <- unique(data_pathways$Pathway)

data_ratios <- data_pathways %>%
  unite(Group, Celltype, ID) %>%
  group_by(Group) %>%
  nest() %>%
  mutate(pairs = map(data, ~ {
    df <- .x
    expand.grid(
      Pathway1 = df$Pathway,
      Pathway2 = df$Pathway,
      stringsAsFactors = FALSE
    ) %>%
       dplyr::filter(Pathway1 != Pathway2) %>%
       left_join(df, by = c("Pathway1" = "Pathway")) %>%
       dplyr::rename(value1 = score) %>%
       left_join(df, by = c("Pathway2" = "Pathway")) %>%
       dplyr::rename(value2 = score) %>%
       mutate(ratio = value1 / value2)
  })) %>%
  select(-data) %>%
  unnest(pairs)

# Step 5: Normalize ratio by average across all groups
data_ratios <- data_ratios %>%
  group_by(Pathway1, Pathway2) %>%
  mutate(average_ratio = mean(ratio, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(corrected = ratio / average_ratio)

# Step 6: Compute mitoPPS score per pathway per group
mitoPPS <- data_ratios %>%
  filter(!is.na(corrected)) %>%
  group_by(Group, Pathway1) %>%
  summarize(mitoPPS = mean(corrected, na.rm = TRUE), .groups = "drop") 

write_csv(mitoPPS, here::here("Data", "MBM_ROSMAP", "ProcessedData", "ROSMAP_mitoPPS.csv"))
