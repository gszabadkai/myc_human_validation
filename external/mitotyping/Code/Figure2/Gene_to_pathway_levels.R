library(tidyverse)
# Gene to pathway mouse ---------------------------------------------------

## Create gene-to-pathway-list
mitocarta_sheet4 <- readxl::read_xls(here::here("Data","MitoCarta", "OriginalData","MouseMitoCarta3_0.xls"), sheet = 4) %>%
  dplyr::select(MitoPathway, Genes) %>%
  na.omit() 
gene_to_pathway <- splitstackshape::cSplit(mitocarta_sheet4, 'Genes', ',') %>%
  column_to_rownames("MitoPathway") %>%
  t() %>%
  as.data.frame() 
gene_to_pathway <- gene_to_pathway %>%
  pivot_longer(cols = colnames(gene_to_pathway), names_to = "Pathway", values_to = "Gene") %>%
  na.omit() 

## add levels
levels <- readxl::read_xls(here::here("Data", "MitoCarta","OriginalData","MouseMitoCarta3_0.xls"), sheet = 4) %>%
  select(MitoPathway, `MitoPathway Hierarchy`) %>%
  separate(`MitoPathway Hierarchy`, into = c("Pathway_Level1", "Pathway_Level2", "Pathway_Level3"), sep = " > ") %>%
  mutate(Level = case_when(
    MitoPathway == Pathway_Level1 ~ "Pathway_Level1", 
    MitoPathway == Pathway_Level2 ~ "Pathway_Level2", 
    MitoPathway == Pathway_Level3 ~ "Pathway_Level3"
  )) %>%
  unique() %>%
  dplyr::rename(Pathway = MitoPathway) %>%
  select(Pathway, Level) %>%
  unique() 

gene_to_pathway_levels <- gene_to_pathway %>%
  full_join(levels)
write.csv(gene_to_pathway_levels, "Data/MitoCarta/ProcessedData/gene_to_pathway_levels_mouse.csv")

# Gene to pathway human ---------------------------------------------------

## Create gene-to-pathway-list
mitocarta_sheet4 <- readxl::read_xls(here::here("Data","MitoCarta", "OriginalData","HumanMitoCarta3_0.xls"), sheet = 4) %>%
  dplyr::select(MitoPathway, Genes) %>%
  na.omit() 
gene_to_pathway <- splitstackshape::cSplit(mitocarta_sheet4, 'Genes', ',') %>%
  column_to_rownames("MitoPathway") %>%
  t() %>%
  as.data.frame() 
gene_to_pathway <- gene_to_pathway %>%
  pivot_longer(cols = colnames(gene_to_pathway), names_to = "Pathway", values_to = "Gene") %>%
  na.omit() 

## add levels
levels <- readxl::read_xls(here::here("Data", "MitoCarta","OriginalData","HumanMitoCarta3_0.xls"), sheet = 4) %>%
  select(MitoPathway, `MitoPathways Hierarchy`) %>%
  separate(`MitoPathways Hierarchy`, into = c("Pathway_Level1", "Pathway_Level2", "Pathway_Level3"), sep = " > ") %>%
  mutate(Level = case_when(
    MitoPathway == Pathway_Level1 ~ "Pathway_Level1", 
    MitoPathway == Pathway_Level2 ~ "Pathway_Level2", 
    MitoPathway == Pathway_Level3 ~ "Pathway_Level3"
  )) %>%
  unique() %>%
  dplyr::rename(Pathway = MitoPathway) %>%
  select(Pathway, Level) %>%
  unique() 

gene_to_pathway_levels <- gene_to_pathway %>%
  full_join(levels)
write.csv(gene_to_pathway_levels, "Data/MitoCarta/ProcessedData/gene_to_pathway_levels_human.csv")