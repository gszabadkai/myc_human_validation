rm(list = ls())

# Setup -------------------------------------------------------------------

library(limma)
library(Glimma)
library(edgeR)
library(tidyverse)
library(tximport)
library(DESeq2)
library(RColorBrewer)
library(broom)

color_tissue <- c(
  "Adipose tissue"= "#A8A8A8","Adrenal gland"= "#A8A8A8", "Amygdala"= "#1DB100",             
  "Appendix"="#265A8C", "Basal ganglia"="#1DB100","Bone marrow"="#00A89D", 
  "Breast"="khaki", "Cerebellum"= "#1DB100", "Cerebral cortex"= "#1DB100",    "Cortex"= "#1DB100",       
  "Cervix"="#EE220C", "Choroid plexus"= "#1DB100", "Colon"="#265A8C", "Duodenum" ="#265A8C",           
  "Endometrium"="#EE220C", "Epididymis"="#EE220C", "Esophagus"="#265A8C", "Fallopian tube"="#EE220C",       
  "Gallbladder" ="#265A8C", "Heart muscle"="#F8BA00", "Hippocampal formation"= "#1DB100", 
  "Hypothalamus" = "#1DB100", "Kidney"="#EF5FA7", "Liver"="#EF5FA7","Lung"="khaki",                
  "Lymph node"="#00A89D", "Medulla oblongata"= "#1DB100", "Midbrain"= "#1DB100",
  "Olfactory bulb"= "#1DB100", "Ovary" ="#EE220C","Pancreas"="khaki", 
  "Parathyroid gland"="#A8A8A8", "Pituitary gland"="#A8A8A8","Placenta"="#EE220C",              
  "Pons"= "#1DB100","Prostate"="#EE220C","Rectum"= "#265A8C","Retina"="khaki",# "#1DB100",             
  "Salivary gland"="#A8A8A8", "Seminal vesicle" ="#EE220C", "Skeletal muscle"="#F8BA00",       
  "Skin" ="khaki","Small intestine"="#265A8C","Smooth muscle" ="khaki","Spinal cord"= "#1DB100",            
  "Spleen"="#00A89D","Stomach"="#265A8C", "Testis"="#EE220C","Thalamus"= "#1DB100",         
  "Thymus"="#00A89D", "Thyroid gland"="#A8A8A8","Tongue"="khaki","Tonsil"="#00A89D",                
  "Urinary bladder"="khaki", "Vagina"="#EE220C", "White matter"= "#1DB100"         
)

# Data preparation --------------------------------------------------------

## Mitocarta ---------------------------------------------------------------

mitocarta <- readxl::read_xls(here::here("Data", "MitoCarta", "OriginalData",
                                         "HumanMitoCarta3_0.xls"), sheet = 4) %>%
  select(MitoPathway, Genes) %>%
  na.omit() 
gene_to_pathway <- splitstackshape::cSplit(mitocarta, 'Genes', ',') %>%
  column_to_rownames("MitoPathway") %>%
  t() %>%
  as.data.frame() 
gene_to_pathway <- gene_to_pathway %>%
  pivot_longer(cols = colnames(gene_to_pathway), names_to = "Pathway", values_to = "Gene") %>%
  na.omit() 

## Read data ---------------------------------------------------------------

counts_raw <- read.delim("Data/GTEx/OriginalData/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", skip = 2)
counts <- counts_raw %>%
  column_to_rownames("Name") %>%
  select(-Description)
pheno <- read.delim("Data/GTEx/OriginalData/GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.tsv") %>%
  separate(SUBJID, into = c("GTEX", "DonorID")) %>%
  select(-GTEX)
meta <- read.delim("Data/GTEx/OriginalData/GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.tsv") %>%
  filter(SMAFRZE == "RNASEQ") %>%
  separate(SAMPID, into = c("GTEX", "DonorID", "TissueSideID", "SM", "AliquotID"),sep = "-", remove=F ) %>%
  select(-c(GTEX, TissueSideID, SM, AliquotID)) %>%
  mutate(SAMPID = gsub("\\-", ".", SAMPID)) %>%
  inner_join(pheno, by = "DonorID")
geneset <- counts_raw %>%
  select(Name, Description)



mitocarta <- readxl::read_xls(here::here("Data","MitoCarta", "OriginalData", "HumanMitoCarta3_0.xls"), sheet = 2) 
mitogenes <- unique(mitocarta$Symbol)
## Mitocarta symbol, synonyms and ensembl
mitogenes_df <- mitocarta %>%
  dplyr::select(Symbol, Synonyms, EnsemblGeneID_mapping_version_20200130) %>%
  dplyr::rename(Symbol_MC = Symbol, Synonyms_MC = Synonyms, Ensembl_MC = EnsemblGeneID_mapping_version_20200130) %>%
  mutate(Ensembl_MC = gsub("|", " ", Ensembl_MC, fixed=TRUE)) %>%
  separate(Ensembl_MC, into = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j"), sep = " ") %>%
  pivot_longer(cols = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j"), values_to = "Ensembl_MC") %>%
  dplyr::select(-name) %>%
  filter(!is.na(Ensembl_MC)) %>%
  mutate(Synonyms_MC = gsub("|", " ", Synonyms_MC, fixed=TRUE)) %>%
  separate(Synonyms_MC, into = as.character(seq(1:15)), sep = " ") %>%
  pivot_longer(cols = as.character(seq(1:15)), values_to = "Synonyms_MC") %>%
  dplyr::select(-name) %>%
  filter(!is.na(Synonyms_MC)) %>%
  mutate(merge = Ensembl_MC) 


## DGE object -------------------------------------------------------

dge <- DGEList(counts[,meta$SAMPID], group = factor(meta$SMTSD))

## Visualize library size --------------------------------------------------
dge$samples %>%
  rownames_to_column("SAMPID") %>%
  filter(group %in% c("Brain - Cortex", "Liver", "Whole Blood", "Muscle - Skeletal")) %>%
  ggplot(aes(x = reorder(SAMPID, lib.size, decreasing = T ), y = lib.size, fill= lib.size)) +
  geom_col(mapping = NULL) +
  scale_fill_viridis_c(direction = -1) + 
  facet_wrap(~group, scales= "free") +
  theme_bw() +
  theme(#legend.position = "none",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    strip.text = element_text(size = 7),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_rect("white"),
    axis.title = element_text(size = 6))
ggsave(here::here("Figures", "Figure2", "Suppl_FigS3A_Libraries.png"), width = 7.5, height = 3.5, units = "in")

dge$samples %>%
  rownames_to_column("SAMPID") %>%
  filter(group %in% c("Brain - Cortex", "Liver", "Whole Blood", "Muscle - Skeletal")) %>%
  group_by(group) %>%
  mutate(max = max(lib.size), min = min(lib.size), fc = max/min) %>%
  ungroup() %>%
  mutate(global_max = max(lib.size), global_min = min(lib.size)) %>%
  select(-c(SAMPID, lib.size, norm.factors)) %>%
  unique() %>%
  mutate(global_fc = global_max/global_min)
  





## TMM normalization -------------------------------------------------------


keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
TMM <- cpm(dge, log = FALSE)

TMM <- as.data.frame(TMM) %>% 
  rownames_to_column("Name")
fst::write.fst(TMM, here::here("Data", "GTEx", "ProcessedData", "GTExTMM_all_tissues.fst"))
TMM <- fst::read.fst(here::here("Data", "GTEx", "ProcessedData", "GTExTMM_all_tissues.fst"))
TMM$samples$lib.size

# MitoPathway scores per tissue  ------------------------------------------------------------

## TMM ---------------------------------------------------------------------


tissue = unique(meta$SMTS)
janitor::tabyl(meta, SMTS) %>%
  arrange(n)

genes_to_keep <- TMM %>%
  select(Name) %>%
  unique() %>%
  full_join(geneset) %>%
  mutate(name2 = sub("\\..*", "",Name)) %>%
  dplyr::rename(Ensembl_gtex_iso = Name, 
                Ensembl_gtex = name2,
                Symbol_gtex = Description) %>%
  dplyr::mutate(merge = Ensembl_gtex) %>%
  inner_join(mitogenes_df, by = "merge") %>%
  filter(!is.na(Ensembl_MC)) %>%
  filter(!is.na(Ensembl_gtex)) %>%
  unique()  
for (i in 1:length(tissue)) {
  
tissue_this = tissue[i]
samples <- meta %>%
  filter(SMTS %in% tissue_this) %>%
  pull(SAMPID)
MitoTMM <- TMM %>%
  filter(Name %in% genes_to_keep$Ensembl_gtex_iso) %>%
  select(Name, all_of(samples)) %>%
  pivot_longer(cols = -Name, names_to = "SAMPID") %>%
  dplyr::rename(Ensembl_gtex_iso = Name) %>%
  inner_join(genes_to_keep) %>%
  group_by(SAMPID,Ensembl_MC) %>%
  mutate(value = mean(value, na.rm = T)) %>%
  ungroup() %>%
  select(SAMPID, Symbol_MC, value) %>%
  unique() %>%
  dplyr::rename(Gene = Symbol_MC) %>%
  inner_join(gene_to_pathway) %>%
  group_by(SAMPID, Pathway) %>%
  mutate(score = mean(value, na.rm = T)) %>%
  select(SAMPID, Pathway, score) %>%
  unique()

fst::write.fst(MitoTMM, here::here("Data", "GTEx", "ProcessedData","MitoTMM_Tissue", paste0(tissue_this, "_MitoPathway_scores_TMM.fst")))
}


## TPM ---------------------------------------------------------------------

TPM <- read.delim("Data/GTEx/OriginalData/bulk-gex_v8_rna-seq_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", skip = 2)

geneset <- TPM %>%
  select(Name, Description)

genes_to_keep <- TPM %>%
  select(Name) %>%
  unique() %>%
  full_join(geneset) %>%
  mutate(name2 = sub("\\..*", "",Name)) %>%
  dplyr::rename(Ensembl_gtex_iso = Name, 
                Ensembl_gtex = name2,
                Symbol_gtex = Description) %>%
  dplyr::mutate(merge = Ensembl_gtex) %>%
  inner_join(mitogenes_df, by = "merge") %>%
  filter(!is.na(Ensembl_MC)) %>%
  filter(!is.na(Ensembl_gtex)) %>%
  unique()  

for (i in 1:length(tissue)) {
  
  tissue_this = tissue[i]
  samples <- meta %>%
    filter(SMTS %in% tissue_this) %>%
    pull(SAMPID)
  MitoTPM <- TPM %>%
    filter(Name %in% genes_to_keep$Ensembl_gtex_iso) %>%
    select(Name, all_of(samples)) %>%
    pivot_longer(cols = -Name, names_to = "SAMPID") %>%
    dplyr::rename(Ensembl_gtex_iso = Name) %>%
    inner_join(genes_to_keep) %>%
    group_by(SAMPID,Ensembl_MC) %>%
    mutate(value = mean(value, na.rm = T)) %>%
    ungroup() %>%
    select(SAMPID, Symbol_MC, value) %>%
    unique() %>%
    dplyr::rename(Gene = Symbol_MC) %>%
    inner_join(gene_to_pathway) %>%
    group_by(SAMPID, Pathway) %>%
    mutate(score = mean(value, na.rm = T)) %>%
    select(SAMPID, Pathway, score) %>%
    unique()
  
  fst::write.fst(MitoTPM, here::here("Data", "GTEx", "ProcessedData","MitoTMM_Tissue", paste0(tissue_this, "_MitoPathway_scores_TPM.fst")))
}


# Blood, Brain, Liver, Muscle TPM vs TMM ----------------------------------

TPM <- list.files(here::here("Data", "GTEx", "ProcessedData","MitoTMM_Tissue"), pattern = "TPM", full.names = T) %>%
  lapply(fst::read.fst) %>%
  bind_rows() %>%
  inner_join(meta) %>%
  filter(SMTSD %in% c("Whole Blood", "Brain - Cortex", "Liver", "Muscle - Skeletal")) %>%
  select(SAMPID, SMTSD, Pathway, score) %>%
  unique() %>%
  dplyr::rename(TPM = score)

TMM <- list.files(here::here("Data", "GTEx", "ProcessedData","MitoTMM_Tissue"), pattern = "TMM", full.names = T) %>%
  lapply(fst::read.fst) %>%
  bind_rows() %>%
  inner_join(meta) %>%
  filter(SMTSD %in% c("Whole Blood", "Brain - Cortex", "Liver", "Muscle - Skeletal")) %>%
  select(SAMPID, SMTSD, Pathway, score) %>%
  unique() %>%
  dplyr::rename(TMM = score)

combined <- inner_join(TPM, TMM) %>%
  dplyr::rename(Tissue = SMTSD)
## Complex I -------------------------------------------------------------------

combined %>%
  filter(Pathway %in% "Complex I") %>%
  pivot_longer(cols = c("TPM", "TMM"), names_to = "norm_method") %>%
  ggplot(aes(x = log10(value), color = Tissue, fill = Tissue)) +
  geom_density(alpha = 0.4) +
  scale_color_manual(values = c(
    "Brain - Cortex" = "#1DB100",
    "Whole Blood"="#00A89D",
    "Liver"="#EF5FA7",
    "Muscle - Skeletal" = "#F8BA00"
  )) + 
  scale_fill_manual(values = c(
    "Brain - Cortex" = "#1DB100",
    "Whole Blood"="#00A89D",
    "Liver"="#EF5FA7",
    "Muscle - Skeletal" = "#F8BA00"
  )) + 
  labs(x = "log10(Complex I score)") +
  facet_wrap(~fct_relevel(norm_method, c("TPM", "TMM"))) +
  theme_bw() +
  theme(axis.text = element_text(size = 6), 
        axis.title = element_text(size = 6), 
        strip.text = element_text(size = 7),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6))
ggsave(here::here("Figures", "Figure2", "Suppl_FigS3_CI_TPM_TMM.png"), 
       width = 3.64, height = 2, units = "in", dpi = 1200)

## FAO -------------------------------------------------------------------

combined %>%
  filter(Pathway %in% "Fatty acid oxidation") %>%
  pivot_longer(cols = c("TPM", "TMM"), names_to = "norm_method") %>%
  ggplot(aes(x = log10(value), color = Tissue, fill = Tissue)) +
  geom_density(alpha = 0.4) +
  scale_color_manual(values = c(
    "Brain - Cortex" = "#1DB100",
    "Whole Blood"="#00A89D",
    "Liver"="#EF5FA7",
    "Muscle - Skeletal" = "#F8BA00"
  )) + 
  scale_fill_manual(values = c(
    "Brain - Cortex" = "#1DB100",
    "Whole Blood"="#00A89D",
    "Liver"="#EF5FA7",
    "Muscle - Skeletal" = "#F8BA00"
  )) + 
  labs(x = "log10(FAO score)") +
  facet_wrap(~fct_relevel(norm_method, c("TPM", "TMM"))) +
  theme_bw() +
  theme(axis.text = element_text(size = 6), 
        axis.title = element_text(size = 6), 
        strip.text = element_text(size = 7),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6))
ggsave(here::here("Figures", "Figure2", "Suppl_FigS3_FAO_TPM_TMM.png"), 
       width = 3.64, height = 2, units = "in", dpi = 1200)

## Complex II -------------------------------------------------------------------

combined %>%
  filter(Pathway %in% "Complex II") %>%
  pivot_longer(cols = c("TPM", "TMM"), names_to = "norm_method") %>%
  ggplot(aes(x = log10(value), color = Tissue, fill = Tissue)) +
  geom_density(alpha = 0.4) +
  scale_color_manual(values = c(
    "Brain - Cortex" = "#1DB100",
    "Whole Blood"="#00A89D",
    "Liver"="#EF5FA7",
    "Muscle - Skeletal" = "#F8BA00"
  )) + 
  scale_fill_manual(values = c(
    "Brain - Cortex" = "#1DB100",
    "Whole Blood"="#00A89D",
    "Liver"="#EF5FA7",
    "Muscle - Skeletal" = "#F8BA00"
  )) + 
  labs(x = "log10(Complex II score)") +
  facet_wrap(~fct_relevel(norm_method, c("TPM", "TMM"))) +
  theme_bw() +
  theme(axis.text = element_text(size = 6), 
        axis.title = element_text(size = 6), 
        strip.text = element_text(size = 7),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6))
ggsave(here::here("Figures", "Figure2", "Suppl_FigS3_CII_TPM_TMM.png"), 
       width = 3.64, height = 2, units = "in", dpi = 1200)

