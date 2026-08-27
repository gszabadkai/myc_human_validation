rm(list = ls())
library(tidyverse)

# Data preparation --------------------------------------------------------------------

samples <- read_csv(here::here("Data", "Fibroblasts", "OriginalData", "RNAseq_meta.csv")) %>%
  mutate(tmp = "Sample_") %>%
  unite(RNAseq_sampleID, tmp, RNAseq_sampleID, sep = "") %>%
  as.data.frame() 

data_mito <- read_csv(here::here("Data","Fibroblasts", "ProcessedData", "RNAseqData_MitoGenes.csv")) %>%
  full_join(samples) %>%
  mutate(Condition = case_when(
    (Clinical_condition %in% "Normal" & Treatments %in% "Control" & Percent_oxygen == 21) ~ "Control",
    (Clinical_condition %in% "Normal" & Treatments %in% "Oligomycin" & Percent_oxygen == 21) ~ "Oligomycin",
    (Clinical_condition %in% "Normal" & Treatments %in% "Oligomycin+DEX" & Percent_oxygen == 21) ~ "Oligomycin_DEX",
    (Clinical_condition %in% "SURF1_Mutation" & Treatments %in% "Control" & Percent_oxygen == 21) ~ "SURF1",
    (Clinical_condition %in% "SURF1_Mutation" & Treatments %in% "Control" & Percent_oxygen == 3) ~ "SURF1_Hypoxia",
    (Clinical_condition %in% "SURF1_Mutation" & Treatments %in% "DEX" & Percent_oxygen == 21) ~ "SURF1_DEX",
    (Clinical_condition %in% "Normal" & Treatments %in% "2DG" & Percent_oxygen == 21) ~ "2DG",
    (Clinical_condition %in% "Normal" & Treatments %in% "Galactose" & Percent_oxygen == 21) ~ "Galactose",
    (Clinical_condition %in% "Normal" & Treatments %in% "Contact_Inhibition" & Percent_oxygen == 21) ~ "Contact_Inhibition",
    (Clinical_condition %in% "Normal" & Treatments %in% "Control" & Percent_oxygen %in% 3) ~ "Hypoxia",
    (Clinical_condition %in% "Normal" & Treatments %in% "DEX" & Percent_oxygen == 21) ~ "DEX",
    (Clinical_condition %in% "Normal" & Treatments %in% "betahydroxybutyrate" & Percent_oxygen == 21) ~ "BHB",
    (Clinical_condition %in% "Normal" & Treatments %in% "Contact_Inhibition" & Percent_oxygen == 3) ~ "Contact_Inhibition_Hypoxia",
    (Clinical_condition %in% "Normal" & Treatments %in% "mitoNUITs" & Percent_oxygen == 21) ~ "MitoNUITs",
    (Clinical_condition %in% "Normal" & Treatments %in% "mitoNUITs+DEX" & Percent_oxygen == 21) ~ "MitoNUITs_DEX",
    TRUE ~ NA
  )) %>%
  filter(!Condition %in% NA)

## MitoCarta ---------------------------------------------------------------

# Create gene-to-pathway-list
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
mitogenes <- readxl::read_xls(here::here("Data","MitoCarta", "OriginalData","HumanMitoCarta3_0.xls"), sheet = 2) %>%
  select(Symbol)
mitogenes <- mitogenes$Symbol

levels <- readxl::read_xls(here::here("Data", "MitoCarta","OriginalData","HumanMitoCarta3_0.xls"), sheet = 4) %>%
  select(MitoPathway, `MitoPathways Hierarchy`) %>%
  separate(`MitoPathways Hierarchy`, into = c("Pathway_Level1", "Pathway_Level2", "Pathway_Level3"), sep = " > ") %>%
  mutate(Level = case_when(
    MitoPathway == Pathway_Level1 ~ "Pathway_Level1", 
    MitoPathway == Pathway_Level2 ~ "Pathway_Level2", 
    MitoPathway == Pathway_Level3 ~ "Pathway_Level3"
  )) %>%
  unique() %>%
  dplyr::rename(Pathway = MitoPathway)


## Pathway data ------------------------------------------------------------

data_pathways <- data_mito %>%
  full_join(gene_to_pathway, by = "Gene") %>%
  filter(!is.na(Pathway)) %>%
  group_by(RNAseq_sampleID, Pathway) %>%
  mutate(nTPM_score = mean(nTPM, na.rm = T)) %>%
  select(-c(Gene, nTPM, Genome)) %>%
  unique() 
write_csv(data_pathways, 
          here::here("Data", "Fibroblasts", "ProcessedData", "MitoPathwayScore_All.csv"))

all_pathways <- unique(data_pathways$Pathway)

data_pathways <- data_pathways %>%
  select(RNAseq_sampleID, Condition, Pathway, nTPM_score) %>%
  pivot_wider(names_from = "Pathway", values_from = "nTPM_score")


## Calculate ratios + mitoPPS --------------------------------------------------------


pw_pairs <- data.frame(pw1 = rep(all_pathways, each = 149)) %>%
  group_by(pw1) %>%
  nest() %>%
  mutate(data = map(data, ~mutate(.x, pw2 = unique(all_pathways)))) %>%
  unnest() %>%
  filter(!pw1==pw2) 

colnames(pw_pairs) <- c("pw1", "pw2")

data_long  <- data_pathways %>%
  pivot_longer(cols = -c(RNAseq_sampleID, Condition), names_to = "Pathway" ) 

janitor::tabyl(data_long, Pathway)
data_added1 <- pw_pairs %>%
  ungroup() %>%
  dplyr::select(pw1) %>%
  dplyr::rename(Pathway = pw1) %>%
  inner_join(data_long, by = "Pathway") %>%
  unique() %>%
  dplyr::rename(value1 = value, Pathway1 = Pathway) %>%
  dplyr::slice(rep(row_number(), 148)) %>%
  arrange(Pathway1) %>%
  group_by(RNAseq_sampleID, Condition) %>%
  nest() %>%
  dplyr::rename(pw1 = data) 

data_added2 <- pw_pairs %>%
  ungroup() %>%
  dplyr::select(pw2) %>%
  dplyr::rename(Pathway = pw2) %>%
  inner_join(data_long, by = "Pathway") %>%
  unique() %>%
  dplyr::rename(value2 = value, Pathway2 = Pathway) %>%
  dplyr::slice(rep(row_number(), 148)) %>%
  arrange(desc(Pathway2)) %>%
  unique() %>%
  group_by(RNAseq_sampleID, Condition) %>%
  nest() %>%
  dplyr::rename(pw2 = data) 

data_ratios =full_join(data_added1, data_added2, by = c("RNAseq_sampleID", "Condition")) %>%
  mutate(combined =  map2(pw1, pw2, ~ cbind(.x,  .y))) %>%
  select(-c(pw1, pw2)) %>%
  unnest(cols = combined) %>%
  unique() %>%
  mutate(ratio = value1 / value2)

ctrl_ratios = data_ratios %>%
  filter(Condition %in% "Control") %>% 
  group_by(Pathway1, Pathway2) %>%
  mutate(ctrl_ratio_av = mean(ratio)) %>%
  select(Pathway1, Pathway2, ctrl_ratio_av) %>%
  unique() 

data_ratios_corrected <- data_ratios %>%
  #filter(!Condition %in% "Control") %>%
  full_join(ctrl_ratios) %>%
  select(-c(value1, value2)) %>%
  dplyr::mutate(corrected_ratio = ratio / ctrl_ratio_av)

meta <- data_mito %>%
  dplyr::rename(Line= Cell_line_inhouse) %>%
  dplyr::select(RNAseq_sampleID, Line, Study_part, Passage, Days_grown_Udays, Treatments, Clinical_condition, Condition) %>%
  unique()

mitoPPS <-data_ratios_corrected %>%
  group_by(RNAseq_sampleID, Pathway1) %>%
  mutate(mitoPPS = mean(corrected_ratio)) %>%
  select(-c(Pathway2, corrected_ratio, ratio, ctrl_ratio_av)) %>%
  unique() %>%
  full_join(meta) %>%
  pivot_wider(names_from = "Pathway1", values_from = "mitoPPS")  %>%
  filter(!RNAseq_sampleID %in% "Sample_152" )  ## This sample has infinity values in Creatine metabolism. Could be an artifact, or induced by Oligomycin treatment
  ## Check if Creatine metabolism is in general higher in oligo treatment
  

write_csv(mitoPPS, 
          here::here("Data", "Fibroblasts", "ProcessedData", "mitoPPS_All.csv"))

