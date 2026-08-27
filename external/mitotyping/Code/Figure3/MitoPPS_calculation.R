rm(list = ls())


# Setup -------------------------------------------------------------------

library(tidyverse)

# Assing tissue groups
tissue_to_group_hm <- function(x) {
  mutate(x, Group = case_when(
    (Tissue == "bone marrow" | Tissue =="thymus" | Tissue =="tonsil" | Tissue =="lymph node" | Tissue =="spleen")~"Immune",
    (Tissue == "placenta" | Tissue == "ovary" | Tissue == "seminal vesicle" | Tissue == "cervix" | Tissue == "prostate" | 
       Tissue == "endometrium" | Tissue == "fallopian tube" | Tissue == "vagina" | Tissue =="epididymis" | Tissue == "testis") ~ "Reproductive",
    (Tissue == "gallbladder" | Tissue == "duodenum" | Tissue == "small intestine" | 
       Tissue == "stomach" | Tissue == "rectum" | Tissue == "colon"| Tissue == "esophagus"| Tissue == "appendix") ~ "Digestive",
    (Tissue == "cerebellum" | Tissue ==  "white matter" | Tissue ==  "cerebral cortex" | Tissue ==  "choroid plexus" | 
       Tissue ==  "thalamus" | Tissue ==  "hypothalamus" | Tissue ==  "medulla oblongata" | Tissue ==  "basal ganglia" | 
       Tissue ==  "pons" | Tissue ==  "spinal cord" | Tissue ==  "midbrain" | Tissue ==  "amygdala" | Tissue ==  "hippocampal formation" | 
       Tissue == "olfactory bulb") ~"CNS",
    (Tissue == "heart muscle" | Tissue == "skeletal muscle") ~ "Contractile",
    (Tissue == "adrenal gland" | Tissue == "parathyroid gland" | Tissue == "pituitary gland"|
       Tissue == "salivary gland" | Tissue == "thyroid gland" | Tissue == "adipose tissue" ) ~ "Secretory",
    (Tissue == "liver" | Tissue == "kidney") ~ "Anabolic",
    (Tissue == "lung" | Tissue == "pancreas"| Tissue == "skin"| Tissue == "smooth muscle" | 
       Tissue == "tongue" | Tissue == "urinary bladder"| Tissue == "breast" | Tissue == "retina") ~ "Other"), 
    .after = "Tissue")
}

## Get human mitocarta genes
mitogenes <- readxl::read_xls(here::here("Data","MitoCarta","OriginalData","HumanMitoCarta3_0.xls"), sheet = 2) %>%
  dplyr::rename(Gene = Symbol) %>%
  select(Gene)
mitogenes <- mitogenes$Gene
# Load hpa data
hpa_data <- read_tsv(here::here("Data","HumanProteinAtlas","OriginalData","rna_human_tissue_consensus.tsv"),show_col_types = FALSE)
# Processed dataset 
processed_data <-  hpa_data%>%
  dplyr::rename(GeneName = `Gene name`) %>%
  filter(GeneName %in% mitogenes) %>%
  select(GeneName, Tissue, nTPM) %>%
  na.omit() %>%
  pivot_wider(names_from = "GeneName", values_from = "nTPM") %>%
  select(-TSTD3) %>% # has too many NAs
  mutate(Tissue = str_to_lower(Tissue)) %>%
  tissue_to_group_hm() %>%
  mutate(Tissue = str_to_sentence(Tissue))
# Assing tissue groups

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


data_pathways <- processed_data %>%
  pivot_longer(cols = colnames(processed_data[3:ncol(processed_data)]), names_to = "Gene") %>%
  full_join(gene_to_pathway, by = "Gene") %>%
  filter(!is.na(Pathway)) %>%
  unique()%>%
  group_by(Tissue, Pathway) %>%
  mutate(score = mean(value, na.rm = T)) %>%
  select(-c(Gene, value)) %>%
  unique() 

write_csv(data_pathways,here::here("Data", "HumanProteinAtlas","ProcessedData", "Pathway_scores_all_tissues.csv" ))

# Pathway ratios ----------------------------------------------------------

all_pathways <- unique(data_pathways$Pathway)
data_long <- data_pathways %>%
  select(Tissue, Group, Pathway, score) %>%
  unite(Group, Tissue, Group)

source(here::here("Code", "mitoPPS_function.R"))
mitoPPS <- mitoPPS(data = data_long, group_by = Group) %>%
  separate(Group, into = c("Tissue", "Group"), sep = "_") %>%
  pivot_wider(names_from = "Pathway", values_from = "mitoPPS")

write_csv(mitoPPS,here::here("Data", "HumanProteinAtlas","ProcessedData", "mitoPPS_all_tissues.csv" ))

