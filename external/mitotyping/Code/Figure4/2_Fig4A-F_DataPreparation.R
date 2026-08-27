rm(list = ls())
library(tidyverse)

samples <- read_csv(here::here("Data", "Fibroblasts", "OriginalData", "RNAseq_meta.csv")) %>%
  mutate(tmp = "Sample_") %>%
  unite(RNAseq_sampleID, tmp, RNAseq_sampleID, sep = "") %>%
  as.data.frame() 

data_mito <- read_csv(here::here("Data","Fibroblasts", "ProcessedData", "RNAseqData_MitoGenes.csv")) %>%
  full_join(samples) %>%
  filter(Clinical_condition %in% "Normal" & Treatments %in% "Control" & Percent_oxygen == 21) %>%
  filter(Study_part %in% 2) %>%
  dplyr::rename(Line = Cell_line_inhouse) %>%
  mutate(Condition = case_when(
    (Clinical_condition %in% "Normal" & Treatments %in% "Control" & Percent_oxygen == 21) ~ "Control",
    TRUE ~ NA
  )) %>%
  filter(!Condition %in% NA)

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


# Calculate MitoPathway scores --------------------------------------------

data_pathways <- data_mito %>%
  full_join(gene_to_pathway, by = "Gene") %>%
  filter(!is.na(Pathway)) %>%
  unique()%>%
  group_by(RNAseq_sampleID, Pathway) %>%
  mutate(score = mean(nTPM, na.rm = T)) %>%
  dplyr::rename(Cell_line_inhouse = Line) %>%
  select(RNAseq_sampleID, Condition, Pathway, score) %>%
  unique()
all_pathways <- unique(data_pathways$Pathway)

#Calculate ratios --------------------------------------------------------


pw_pairs <- data.frame(pw1 = rep(all_pathways, each = 149)) %>%
  group_by(pw1) %>%
  nest() %>%
  mutate(data = map(data, ~mutate(.x, pw2 = unique(all_pathways)))) %>%
  unnest() %>%
  filter(!pw1==pw2) 
colnames(pw_pairs) <- c("pw1", "pw2")

data_long  <- data_pathways

janitor::tabyl(data_long, Pathway)
data_added1 <- pw_pairs %>%
  ungroup() %>%
  dplyr::select(pw1) %>%
  dplyr::rename(Pathway = pw1) %>%
  inner_join(data_long, by = "Pathway") %>%
  unique() %>%
  dplyr::rename(value1 = score, Pathway1 = Pathway) %>%
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
  dplyr::rename(value2 = score, Pathway2 = Pathway) %>%
  dplyr::slice(rep(row_number(), 148)) %>%
  arrange(desc(Pathway2)) %>%
  unique() %>%
  group_by(RNAseq_sampleID, Condition) %>%
  nest() %>%
  dplyr::rename(pw2 = data) 

data_ratios =full_join(data_added1, data_added2, by = c("RNAseq_sampleID", "Condition")) %>%
  dplyr::mutate(combined =  map2(pw1, pw2, ~ cbind(.x,  .y))) %>%
  dplyr::select(-c(pw1, pw2)) %>%
  unnest(cols = combined) %>%
  unique() %>%
  dplyr::mutate(ratio = value1 / value2)
sapply(data_ratios, class)
ctrl_ratios = data_ratios %>%
  group_by(Pathway1, Pathway2) %>%
  mutate(ctrl_ratio_av = mean(ratio)) %>%
  select(Pathway1, Pathway2, ctrl_ratio_av) %>%
  unique() 

data_ratios_corrected <- data_ratios %>%
  full_join(ctrl_ratios) %>%
  select(-c(value1, value2)) %>%
  dplyr::mutate(corrected_ratio = ratio / ctrl_ratio_av)

meta <- data_mito %>%
  select(RNAseq_sampleID, Line, Study_part, Passage, Days_grown_Udays, Treatments, Clinical_condition, Condition) %>%
  unique()

mitoPPS <-data_ratios_corrected %>%
  group_by(RNAseq_sampleID, Pathway1) %>%
  mutate(mitoPPS = mean(corrected_ratio)) %>%
  select(-c(Pathway2, corrected_ratio, ratio, ctrl_ratio_av)) %>%
  unique() %>%
  full_join(meta) %>%
  pivot_wider(names_from = "Pathway1", values_from = "mitoPPS")  

write_csv(mitoPPS, 
          here::here("Data", "Fibroblasts", "ProcessedData", "mitoPPS_Controls_SP2.csv"))

