library(tidyverse)
library(Seurat)

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

mitogenes <- readxl::read_xls(here::here("Data/MitoCarta/OriginalData/HumanMitoCarta3_0.xls"), sheet = 2) %>%
  pull(Symbol)



mitobrain_data <- readRDS(here::here("Data", "MBM_ROSMAP", "OriginalData",  "MitoBrainMap_SingleCellData_Downsampled_clean.rds"))
mitobrain_data <- subset(mitobrain_data, features = mitogenes)
mitobrain_data <- AverageExpression(mitobrain_data, group.by = c("orig.ident", "CellType"), slot = "data", return.seurat=F)
data_pathways <- mitobrain_data[["RNA"]] %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene)  %>%
  #filter(!grepl("MT-", Gene)) %>%
  full_join(gene_to_pathway) %>%
  filter(!is.na(Pathway)) %>%
  filter(!is.na(name)) %>%
  group_by(name, Pathway) %>%
  mutate(score = mean(value)) %>%
  ungroup() %>%
  select(-c(Gene, value)) %>%
  unique() %>%
  separate(name, into = c("Voxel", "CellType")) %>%
  mutate(Voxel = case_when(Voxel == "voxel1" ~ "Cort",
                           Voxel == "voxel2" ~ "Hipp",
                           Voxel == "voxel3" ~ "Put",
                           Voxel == "voxel4" ~ "WM") )  %>%
  select(CellType, Voxel, Pathway, score) %>%
  filter(CellType %in% c("Inhibitory", "Excitatory", "OPC", "Oligodendrocytes", "Astrocytes", "Microglia",
                         "Pericytes", "Endothelial")) 
all_pathways <- unique(data_pathways$Pathway)

# MitoPPS -----------------------------------------------------------------

mitoPPS_per_voxel <- function(voxel) {
  data_ratios <- data_pathways %>%
    filter(Voxel %in% voxel) %>%
    unite(Group, CellType, Voxel) %>%
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
  
  data_ratios <- data_ratios %>%
    group_by(Pathway1, Pathway2) %>%
    mutate(average_ratio = mean(ratio, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(corrected = ratio / average_ratio)
  
  mitoPPS <- data_ratios %>%
    filter(!is.na(corrected)) %>%
    group_by(Group, Pathway1) %>%
    summarize(mitoPPS = mean(corrected, na.rm = TRUE), .groups = "drop") 
  return(mitoPPS)
}
unique(data_pathways$Voxel)
vox1 <- mitoPPS_per_voxel("Cort")
vox2 <- mitoPPS_per_voxel("Hipp")
vox3 <- mitoPPS_per_voxel("Put")
vox4 <- mitoPPS_per_voxel("WM")
mitoPPS <- rbind(vox1, vox2, vox3, vox4) %>%
  separate(Group, into = c("CellType", "Voxel"))

# 
# mitoPPS_per_voxel <- function() {
#   data_ratios <- data_pathways %>%
#    filter(!Voxel %in% "WM") %>%
#     unite(Group, CellType, Voxel) %>%
#     group_by(Group) %>%
#     nest() %>%
#     mutate(pairs = map(data, ~ {
#       df <- .x
#       expand.grid(
#         Pathway1 = df$Pathway,
#         Pathway2 = df$Pathway,
#         stringsAsFactors = FALSE
#       ) %>%
#         dplyr::filter(Pathway1 != Pathway2) %>%
#         left_join(df, by = c("Pathway1" = "Pathway")) %>%
#         dplyr::rename(value1 = score) %>%
#         left_join(df, by = c("Pathway2" = "Pathway")) %>%
#         dplyr::rename(value2 = score) %>%
#         mutate(ratio = value1 / value2)
#     })) %>%
#     select(-data) %>%
#     unnest(pairs)
#   
#   data_ratios <- data_ratios %>%
#     group_by(Pathway1, Pathway2) %>%
#     mutate(average_ratio = mean(ratio, na.rm = TRUE)) %>%
#     ungroup() %>%
#     mutate(corrected = ratio / average_ratio)
#   
#   mitoPPS <- data_ratios %>%
#     filter(!is.na(corrected)) %>%
#     group_by(Group, Pathway1) %>%
#     summarize(mitoPPS = mean(corrected, na.rm = TRUE), .groups = "drop") 
#   return(mitoPPS)
# }
# # unique(data_pathways$Voxel)
# # vox1 <- mitoPPS_per_voxel("Cort")
# # vox2 <- mitoPPS_per_voxel("Hipp")
# # vox3 <- mitoPPS_per_voxel("Put")
# # vox4 <- mitoPPS_per_voxel("WM")
# mitoPPS <- mitoPPS_per_voxel() %>%
#   separate(Group, into = c("CellType", "Voxel"))
 write_csv(mitoPPS, here::here("Data", "MBM_ROSMAP", "ProcessedData", "MBM_mitoPPS.csv"))
