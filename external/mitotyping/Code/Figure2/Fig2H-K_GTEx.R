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
  filter(SMTS %in% c("Brain", "Liver")) %>%
  inner_join(pheno, by = "DonorID")

meta_all <- read.delim("Data/GTEx/OriginalData/GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.tsv") %>%
  filter(SMAFRZE == "RNASEQ") %>%
  separate(SAMPID, into = c("GTEX", "DonorID", "TissueSideID", "SM", "AliquotID"),sep = "-", remove=F ) %>%
  select(-c(GTEX, TissueSideID, SM, AliquotID)) %>%
  mutate(SAMPID = gsub("\\-", ".", SAMPID)) %>%
  inner_join(pheno, by = "DonorID")
length(unique(meta_all$DonorID))
length(unique(meta_all$SMTSD))
geneset <- counts_raw %>%
  select(Name, Description) 

# liver <- read.delim("Data/GTEx/OriginalData/GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.tsv") %>%
#   filter(SMAFRZE == "RNASEQ") %>%
#   separate(SAMPID, into = c("GTEX", "DonorID", "TissueSideID", "SM", "AliquotID"),sep = "-", remove=F ) %>%
#   select(-c(GTEX, TissueSideID, SM, AliquotID)) %>%
#   mutate(SAMPID = gsub("\\-", ".", SAMPID)) %>%
#   filter(SMTS %in% c("Liver"))

## TMM normalization -------------------------------------------------------
## !! Only on liver and cns samples
dge <- DGEList(counts[,meta$SAMPID], group = factor(meta$SMTSD))
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
TMM <- cpm(dge, log = FALSE)


fst::write_fst(as.data.frame(TMM) %>% rownames_to_column("Gene"), here::here("Data", "GTEx", "ProcessedData", "TMM.fst"))
TMM <- fst::read_fst(here::here("Data", "GTEx", "ProcessedData", "TMM.fst")) %>%
  column_to_rownames("Gene")



random_samples_for_plotting <- round(runif(10, 0, 2868))
col <- brewer.pal(10, "Paired")
par(mfrow=c(1,2))
cpm <- cpm(counts[,-c(1:2)], log=F)
boxplot(log2(cpm+2/sapply(as.data.frame(cpm), sum))[,random_samples_for_plotting],col=col,  las=2, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")

boxplot(log2(TMM+2/sapply(as.data.frame(TMM), sum))[,random_samples_for_plotting], las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")

#write_csv(TMM, here::here("Data", "GTEx", "ProcessedData", "TMM.csv"))
rm(list = setdiff(ls(), c("counts", "pheno","meta", "TMM", "geneset", "color_tissue")))


# Mitotyping ------------------------------------------------------------

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


gtex_genes <- geneset %>%
  mutate(name2 = sub("\\..*", "",Name)) %>%
  dplyr::rename(Ensembl_gtex_iso = Name, 
                Ensembl_gtex = name2,
                Symbol_gtex = Description) %>%
  dplyr::mutate(merge = Ensembl_gtex) 

match_by_ensembl <- full_join(gtex_genes, mitogenes_df, by = "merge") %>%
  filter(!is.na(Ensembl_MC)) %>%
  filter(!is.na(Ensembl_gtex)) %>%
  unique() 

## Filter the dissimilar ones out
match_by_ensembl_dissimilar <- full_join(gtex_genes, mitogenes_df, by = "merge") %>%
  filter(!is.na(Ensembl_MC)) %>%
  filter(!is.na(Ensembl_gtex)) %>%
  unique() %>%
  mutate(same_ID = case_when(
    Symbol_MC == Symbol_gtex ~ 1,
    TRUE~0
  )) %>%
  filter(same_ID %in% 0) %>%
  group_by(Symbol_gtex) %>%
  mutate(syn_1 = case_when(
    Symbol_gtex == all_of(Synonyms_MC) ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(passed_1 = sum(syn_1))


genes_not_found <- mitogenes[which(!mitogenes %in% match_by_ensembl$Symbol_MC)]
genes_not_found
#to_check<- c("SOC2", "CMC4", "ATP5MF-PTCD1")




## MitoTMM -----------------------------------------------------------------

MitoTMM <- as.data.frame(TMM) %>% 
  rownames_to_column("Name") %>%
  full_join(geneset) %>%
  mutate(name2 = sub("\\..*", "",Name)) %>%
  dplyr::rename(Ensembl_gtex_iso = Name, 
                Ensembl_gtex = name2,
                Symbol_gtex = Description) %>%
  dplyr::mutate(merge = Ensembl_gtex) %>%
  full_join(mitogenes_df, by = "merge") %>%
  filter(!is.na(Ensembl_MC)) %>%
  filter(!is.na(Ensembl_gtex)) %>%
  unique()  %>%
  pivot_longer(cols = -c(Ensembl_gtex, Symbol_gtex, Ensembl_gtex_iso, merge, Symbol_MC, Ensembl_MC, Synonyms_MC)) %>%
  group_by(name,Ensembl_MC) %>%
  mutate(value = mean(value, na.rm = T)) %>%
  select(name, Symbol_MC, value) %>%
  unique()

## n_Donors, n_Tissues
tmp <- MitoTMM %>%
  dplyr::rename(Gene = Symbol_MC, SAMPID = name) %>%
  inner_join(meta) 
length(unique(tmp$DonorID))

## MitoPathways ------------------------------------------------------------

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

levels <- readxl::read_xls(here::here("Data","MitoCarta", "OriginalData","HumanMitoCarta3_0.xls"), sheet = 4) %>%
  select(MitoPathway, `MitoPathways Hierarchy`) %>%
  separate(`MitoPathways Hierarchy`, into = c("Pathway_Level1", "Pathway_Level2", "Pathway_Level3"), sep = " > ") %>%
  mutate(Level = case_when(
    MitoPathway == Pathway_Level1 ~ "Pathway_Level1", 
    MitoPathway == Pathway_Level2 ~ "Pathway_Level2", 
    MitoPathway == Pathway_Level3 ~ "Pathway_Level3"
  )) %>%
  unique() %>%
  dplyr::rename(Pathway = MitoPathway)


## MitoPathway scores
MitoPathwaysTMM<- MitoTMM %>%
  dplyr::rename(Gene = Symbol_MC, SAMPID = name) %>%
  full_join(gene_to_pathway) %>%
  group_by(SAMPID, Pathway) %>%
  mutate(score_simple = mean(value, na.rm = T)) %>%
  select(-c(Gene, value)) %>%
  unique() %>%
  full_join(levels, by = "Pathway") %>%
  inner_join(meta)

fst::write_fst(MitoPathwaysTMM, here::here("Data", "GTEx", "ProcessedData", "MitoPathwaysTMM.fst"))
rm(list = setdiff(ls(), c("counts", "pheno","meta", "TMM", "geneset", "MitoTMM", "MitoPathwaysTMM", "gene_to_pathway", "color_tissue", "levels")))

MitoPathwaysTMM <- fst::read_fst(here::here("Data", "GTEx", "ProcessedData", "MitoPathwaysTMM.fst"))

# Figure 2H Density FAO Cortex vs Liver ------------------------------------

MitoPathwaysTMM %>%
  filter(Pathway %in% "Fatty acid oxidation") %>%
  filter(grepl("Brain - Cortex", SMTSD) | grepl("Liver", SMTSD)) %>%
  mutate(SMTSD = case_when(SMTSD == "Brain - Cortex" ~ "Cerebral cortex",
                           TRUE ~ SMTSD)) %>%
  select(Pathway, SAMPID, SMTSD, score_simple) %>%
  unique() %>%
  ggplot(aes(x = score_simple, group = SMTSD, fill = SMTSD, color = SMTSD)) +
  geom_histogram(alpha = 0.7, bins = 100, linewidth = 0.1) +
  scale_fill_manual(values = unlist(color_tissue)) + 
  scale_color_manual(values = unlist(color_tissue)) + 
  theme_bw() +
  labs(y = "Frequency", x = "Fatty acid oxidation score") + 
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())

ggsave(here::here("Figures", "Figure2", "Fig2H_Cortex_Liver_Density.png"), width= 2.2, height = 2.2, units = "in", dpi = 1200)

## Stats -------------------------------------------------------------------

test <- MitoPathwaysTMM %>%
  filter(Pathway %in% "Fatty acid oxidation") %>%
  filter(grepl("Brain - Cortex", SMTSD) | grepl("Liver", SMTSD)) %>%
  select(SMTSD, score_simple) %>%
  unique() 

brain <- test %>% 
  filter(SMTSD %in% "Brain - Cortex") %>%
  pull(score_simple)

liver <- test %>% 
  filter(SMTSD %in% "Liver") %>%
  pull(score_simple)

wilcox.test(brain, liver)
mean(liver) / mean(brain) 


# Figure 2I Matching samples ----------------------------------------------

test <- MitoPathwaysTMM %>%
  ungroup() %>%
  filter(Pathway %in% "Fatty acid oxidation") %>%
  filter(grepl("Brain - Cortex", SMTSD) | grepl("Liver", SMTSD)) %>%
  mutate(SMTSD = case_when(SMTSD == "Brain - Cortex" ~ "Cortex",
                           TRUE ~ SMTSD)) %>%
  select(SMTSD, DonorID, score_simple) %>%
  unique() %>%
  pivot_wider(names_from = SMTSD, values_from = score_simple) %>%
  na.omit() %>%
  pivot_longer(cols = c("Cortex", "Liver"), names_to = "Tissue") 

length(unique(test$DonorID))
test %>%
  ggplot(aes(x = Tissue, y = value, color = Tissue, group = DonorID)) +
  geom_point(alpha = 0.4) +
  geom_line(color = "grey", alpha = 0.6, linewidth = 0.2) +
  scale_color_manual(values = unlist(color_tissue)) + 
  theme_bw() +
  labs(y = "Fatty acid oxidation score") + 
  theme(
    axis.title = element_text(size = 7),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())

ggsave(here::here("Figures", "Figure2", "Fig2I_Cortex_Liver_matchingSamples_FAO.png"), width= 1.8, height = 2.1, units = "in", dpi = 1200)

## Average fold change -----------------------------------------------------

test <- MitoPathwaysTMM %>%
  ungroup() %>%
  filter(Pathway %in% "Fatty acid oxidation") %>%
  filter(grepl("Brain - Cortex", SMTSD) | grepl("Liver", SMTSD)) %>%
  mutate(SMTSD = case_when(SMTSD == "Brain - Cortex" ~ "Cortex",
                           TRUE ~ SMTSD)) %>%
  select(SMTSD, DonorID, score_simple) %>%
  unique() %>%
  pivot_wider(names_from = SMTSD, values_from = score_simple) %>%
  mutate(foldchange = Liver / Cortex) %>%
  filter(!is.na(foldchange))
min(test$foldchange)
max(test$foldchange)
mean(test$foldchange)


# Figure 2J ETFDH Complex I-------------------------------------------------------------

ETFDH <- as.data.frame(TMM) %>% 
  rownames_to_column("Name") %>%
  full_join(geneset) %>%
  dplyr::rename(Ensembl_gtex_iso = Name, 
                Symbol_gtex = Description) %>%
  filter(Symbol_gtex %in% c("ETFDH")) %>%
  unique()  %>%
  pivot_longer(cols = -c(Symbol_gtex, Ensembl_gtex_iso), names_to = "SAMPID") %>%
  inner_join(meta) %>%
  select(SAMPID, SMTSD,  value, Symbol_gtex) %>%
  filter(grepl("Brain - Cortex", SMTSD) | grepl("Liver", SMTSD))%>%
  unique() %>%
  pivot_wider(names_from = "Symbol_gtex", values_from = value)

CI <- MitoPathwaysTMM %>%
  ungroup() %>%
  filter(Pathway %in% c("Complex I")) %>%
  filter(grepl("Brain - Cortex", SMTSD) | grepl("Liver", SMTSD))%>%
  select(SAMPID, SMTSD, Pathway, score_simple) %>%
  unique()%>%
  pivot_wider(names_from = "Pathway", values_from = score_simple)

tmp <- CI %>%
  full_join(ETFDH) 

tmp %>%
  ggplot(aes(x = `Complex I`, y = `ETFDH`,color = SMTSD, fill = SMTSD)) + #color = as.numeric(TRISCH))) + #
  geom_point(alpha = 0.4) +
  scale_fill_manual(values = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100")) + 
  scale_color_manual(values = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100"))  + 
  theme_classic() +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.line = element_line(linewidth = 0.25),
        legend.position = "none")
ggsave(here::here("Figures", "Figure2", "Fig2J_Cortex_Liver_CI_ETFDH.png"), width=2, height = 2, units = "in", dpi = 1200)


## ETFDH boxplot ------------------------------------------------------------



### Plot in logarithmic scale
tmp %>%
  ggplot(aes(x = SMTSD, y = log10(ETFDH), fill = SMTSD, color = SMTSD)) + 
  geom_boxplot(alpha = 0.6, size = 0.1, outlier.size = 0.2) +
  ylab("log10 ETFDH expression") + 
  scale_fill_manual(values = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100")) + 
  scale_color_manual(values = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100"))  + 
  theme_classic() +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.line = element_line(linewidth = 0.25),
        legend.position = "none")

ggsave(here::here("Figures", "Figure2", "Fig2J_Cortex_Liver_ETFDH_boxplot.png"), width=1.2, height = 1.2, units = "in", dpi = 1200)

## Stats -------------------------------------------------------------------
tmp %>%
  group_by(SMTSD) %>%
  nest() %>%
  mutate(spearman = map_dbl(data, ~cor(.$ETFDH, .$`Complex I`, method = "spearman"))) %>%
  mutate(cortest = map(data, ~tidy(cor.test(.$ETFDH, .$`Complex I`, method = "spearman", exact = FALSE)))) %>%
  mutate(pvalue = map_dbl(cortest, "p.value")) %>%
  arrange(desc(spearman)) %>%
  ungroup() %>%
  mutate(padjust = p.adjust(pvalue, method = "bonferroni"))

tmp %>%
  group_by(SMTSD) %>%
  mutate(av_exprs  = mean(ETFDH, na.rm = T)) %>%
  select(SMTSD, av_exprs) %>%
  unique() %>%
  pivot_wider(names_from = SMTSD, values_from = av_exprs) %>%
  mutate(fc = `Liver`/`Brain - Cortex`)


test_split <- split(tmp$ETFDH, tmp$SMTSD)
lapply(test_split, shapiro.test) # p< 0.05, not normally distributed
fligner.test(ETFDH ~ SMTSD, data = tmp) # p <0.05, variance not equal
brunnermunzel::brunnermunzel.test(test_split[[1]], test_split[[2]], perm = T)
## HedgesG
effsize::cohen.d(tmp$ETFDH, tmp$SMTSD,pooled=TRUE,paired=FALSE,
                 na.rm=FALSE, mu=0, hedges.correction=TRUE,
                 conf.level=0.95,noncentral=FALSE, 
                 within=TRUE)



# Figure 2K ETFDH Complex II-------------------------------------------------------------

ETFDH <- as.data.frame(TMM) %>% 
  rownames_to_column("Name") %>%
  full_join(geneset) %>%
  dplyr::rename(Ensembl_gtex_iso = Name, 
                Symbol_gtex = Description) %>%
  filter(Symbol_gtex %in% c("ETFDH")) %>%
  unique()  %>%
  pivot_longer(cols = -c(Symbol_gtex, Ensembl_gtex_iso), names_to = "SAMPID") %>%
  inner_join(meta) %>%
  select(SAMPID, SMTSD,  value, Symbol_gtex) %>%
  filter(grepl("Brain - Cortex", SMTSD) | grepl("Liver", SMTSD))%>%
  unique() %>%
  pivot_wider(names_from = "Symbol_gtex", values_from = value)

CII <- MitoPathwaysTMM %>%
  ungroup() %>%
  filter(Pathway %in% c("Complex II")) %>%
  filter(grepl("Brain - Cortex", SMTSD) | grepl("Liver", SMTSD))%>%
  select(SAMPID, SMTSD, Pathway, score_simple) %>%
  unique()%>%
  pivot_wider(names_from = "Pathway", values_from = score_simple)

tmp <- CII %>%
  full_join(ETFDH) 

tmp %>%
  ggplot(aes(x = `Complex II`, y = `ETFDH`,color = SMTSD, fill = SMTSD)) + #color = as.numeric(TRISCH))) + #
  geom_point(alpha = 0.4) +
  scale_fill_manual(values = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100")) + 
  scale_color_manual(values = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100"))  + 
  theme_classic() +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.line = element_line(linewidth = 0.25),
        legend.position = "none")
ggsave(here::here("Figures", "Figure2", "Fig2K_Cortex_Liver_CII_ETFDH.png"), width=2, height = 2, units = "in", dpi = 1200)



## Stats -------------------------------------------------------------------
tmp %>%
  group_by(SMTSD) %>%
  mutate(av_exprs  = mean(`Complex II`, na.rm = T)) %>%
  select(SMTSD, av_exprs) %>%
  unique() %>%
  pivot_wider(names_from = SMTSD, values_from = av_exprs) %>%
  mutate(fc = `Liver`/`Brain - Cortex`)
tmp %>%
  group_by(SMTSD) %>%
  nest() %>%
  mutate(spearman = map_dbl(data, ~cor(.$ETFDH, .$`Complex II`, method = "spearman"))) %>%
  mutate(cortest = map(data, ~tidy(cor.test(.$ETFDH, .$`Complex II`, method = "spearman", exact = FALSE)))) %>%
  mutate(pvalue = map_dbl(cortest, "p.value")) %>%
  arrange(desc(spearman)) %>%
  ungroup() %>%
  mutate(padjust = p.adjust(pvalue, method = "bonferroni"))



# # Figure 2I Bivariate plot CI CII ------------------------------------------------
# 
# tmp <- MitoPathwaysTMM %>%
#   ungroup() %>%
#   filter(Pathway %in% c("Complex I", "Complex II")) %>%
#   filter(grepl("Brain - Cortex", SMTSD) | grepl("Liver", SMTSD)) %>%
#   select(DonorID, SMTSD, Pathway, score_simple, TRISCH) %>%
#   unique() %>%
#   mutate(TRISCH = lubridate::hm(TRISCH)) %>%
#   mutate(TRISCH = lubridate::hour(TRISCH)*60 + lubridate::minute(TRISCH)) %>%
#   pivot_wider(names_from = "Pathway", values_from = "score_simple") %>%
#   select(DonorID, SMTSD, `Complex I`, `Complex II`, TRISCH) %>%
#   unique() 
# 
# tmp %>%
#   ggplot(aes(x = `Complex I`, y = `Complex II`,color = SMTSD, fill = SMTSD)) + #color = as.numeric(TRISCH))) + #
#   geom_point(alpha = 0.4) +
#   scale_fill_manual(values = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100")) + 
#   scale_color_manual(values = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100"))  + 
#   theme_classic() +
#   theme(axis.text = element_text(size = 6),
#         axis.title = element_text(size = 7),
#         axis.line = element_line(linewidth = 0.25),
#         legend.position = "none")
# ggsave(here::here("Figures", "Figure2", "Fig2H_Cortex_Liver_CICII.png"), width=2, height = 2, units = "in", dpi = 1200)
# 
# 
# ## CI CII Ratio boxplot ------------------------------------------------------------
# 
# ratio <- tmp %>%
#   dplyr::mutate(CIICIratio = (`Complex II`/`Complex I`)) %>%
#   dplyr::group_by(SMTSD) %>%
#   dplyr::mutate(mean = mean(CIICIratio)) %>%
#   dplyr::mutate(sd = sd(CIICIratio)) %>%
#   dplyr::mutate(count =n()) %>%
#   dplyr::mutate(sem = sd/sqrt(count)) %>%
#   dplyr::mutate(ymin = mean-sem) %>%
#   dplyr::mutate(ymax = mean+sem) 
# 
# 
# ratio %>%
#   ggplot(aes(x = CIICIratio, color = SMTSD)) + 
#   geom_density() +
#   theme_classic()
# 
# ### Plot in logarithmic scale
# ratio %>%
#   ggplot(aes(x = SMTSD, y = log10(CIICIratio), fill = SMTSD, color = SMTSD)) + 
#   geom_boxplot(alpha = 0.6, size = 0.1, outlier.size = 0.2) +
#   ylab("CII / CI [log10]") + 
#   scale_fill_manual(values = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100")) + 
#   scale_color_manual(values = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100"))  + 
#   theme_classic() +
#   theme(axis.text = element_text(size = 6),
#         axis.title = element_text(size = 6),
#         axis.line = element_line(linewidth = 0.25),
#         legend.position = "none")
# 
# ggsave(here::here("Figures", "Figure2", "Fig2H_Cortex_Liver_CICII_ratio_boxplot.png"), width=1.2, height = 1.2, units = "in", dpi = 1200)
# 
# ## Stats -------------------------------------------------------------------
# 
# cns <- ratio %>%
#   filter(SMTSD %in% "Brain - Cortex") %>%
#   pull(CIICIratio)
# liver <- ratio %>%
#   filter(SMTSD %in% "Liver") %>%
#   pull(CIICIratio)
# 
# test_split <- split(ratio$CIICIratio, ratio$SMTSD)
# lapply(test_split, shapiro.test) # p< 0.05, not normally distributed
# fligner.test(CIICIratio ~ SMTSD, data = ratio) # p <0.05, variance not equal
# brunnermunzel::brunnermunzel.test(test_split[[1]], test_split[[2]], perm = T)
# ## HedgesG
# effsize::cohen.d(ratio$CIICIratio, ratio$SMTSD,pooled=TRUE,paired=FALSE,
#                  na.rm=FALSE, mu=0, hedges.correction=TRUE,
#                  conf.level=0.95,noncentral=FALSE, 
#                  within=TRUE)
# 
# ### Stats on raw values
# t.test(cns, liver)
# mean(liver) / mean(cns)
# 
# p <- ratio %>%
#   select(-c(`Complex I`, `Complex II`, CIICIratio) )%>% 
#   ggplot(aes(fill = SMTSD)) +
#   geom_segment(aes(x=SMTSD,xend=SMTSD,y=0, yend=mean,
#                    color = SMTSD), size=12, alpha =0.6) +
#   geom_errorbar(aes(x=SMTSD,  ymin=ymin, ymax=ymax, color = SMTSD), width=0.3,
#                 alpha=0.8, size=0.3) +
#   scale_fill_manual(values = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100")) + 
#   scale_color_manual(values = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100"))  + 
#   labs(y= "CII / CII") +
#   theme_bw() +
#   theme(#legend.position = "none", 
#     axis.title.x = element_blank(),
#     axis.text.x = element_blank(),
#     axis.text.y = element_text(size = 8),
#     axis.title.y = element_text(size = 10),
#     axis.ticks.x = element_blank(),
#     panel.grid.major = element_blank(), 
#     panel.grid.minor = element_blank())
# p

# Supplementary Figure 2 -------------------------------------------------

## Suppl. Fig. S2C Mean difference --------------------------------------------------------------------

MitoPathwaysTMM %>%
  filter(Pathway %in% "Fatty acid oxidation") %>%
  filter(grepl("Brain - Cortex", SMTSD) | grepl("Liver", SMTSD)) %>%
  select(SMTSD, score_simple) %>%
  unique() %>%
  group_by(SMTSD) %>%
  mutate(mean = mean(score_simple),sd = sd(score_simple), count = n()) %>%
  select(SMTSD,mean, sd, count) %>%
  unique() %>%
  mutate(sem = sd/sqrt(count)) %>%
  mutate(ymin = mean-sem) %>%
  mutate(ymax = mean+sem) %>%
  mutate(Group = SMTSD) %>%
  ggplot(aes(color = SMTSD, fill = SMTSD, alpha = 0.6)) +
  geom_segment(aes(x=Group,xend=Group,y=0 , yend=mean), size=4, alpha =0.6) +
  geom_errorbar(aes(x=Group, ymin=ymin, ymax=ymax, color = Group), width=0.4, 
                alpha=0.8, size=0.3) +
  scale_fill_manual(values = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100")) + 
  scale_color_manual(values = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100"))  + 
  theme_classic() +
  ylab("Fatty acid oxidation score") + 
  scale_y_continuous(limits = c(0,245), expand = c(0, 0)) +
  theme(axis.text = element_text(size = 6),
        axis.text.x = element_text(angle =45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 7),
        axis.line = element_line(linewidth = 0.25),
        legend.position = "none")
ggsave(here::here("Figures", "Figure2", "Suppl_FigS2C_Cortex_Liver_Mean.png"), width=1, height = 2.17, units = "in", dpi = 1200)


## Suppl. Fig. S2D Min to max --------------------------------------------------------------------

test <- MitoPathwaysTMM %>%
  filter(Pathway %in% "Fatty acid oxidation") %>%
  filter(grepl("Brain - Cortex", SMTSD) | grepl("Liver", SMTSD)) %>%
  select(SMTSD, score_simple) %>%
  unique() %>%
  group_by(SMTSD) %>%
  mutate(min = min(score_simple), max = max(score_simple),
         median = median(score_simple), mean = mean(score_simple),sd = sd(score_simple), count = n()) %>%
  select(SMTSD, min, max, median, count, mean) %>%
  unique() %>%
  mutate(max_by_min = max/min)

test %>%
  select(-c(median, mean, count, max_by_min)) %>%
  pivot_longer(cols = -SMTSD) %>%
  mutate(Group = case_when(
    (SMTSD == "Liver" & name == "min")~ "3 Liver min",
    (SMTSD == "Liver" & name == "max")~ "4 Liver max",
    (SMTSD == "Brain - Cortex" & name == "min")~ "1 Cortex min",
    (SMTSD == "Brain - Cortex" & name == "max")~ "2 Cortex max"#,
  )) %>%
  ggplot(aes(color = SMTSD, fill = SMTSD, alpha = 0.6)) +
  geom_segment(aes(x=Group,xend=Group,y=0 , yend=value), size=8, alpha =0.6) +
  scale_fill_manual(values = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100")) + 
  scale_color_manual(values = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100"))  + 
  theme_classic() +
  ylab("Fatty acid oxidation score") + 
  scale_y_continuous(limits = c(0,470), expand = c(0, 0)) +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.line = element_line(linewidth = 0.25),
        legend.position = "none")
ggsave(here::here("Figures", "Figure2", "Suppl_FigS2D_Cortex_Liver_MinMax.png"), width=1.8, height = 2, units = "in", dpi = 1200)


## Suppl. Fig. S2E Heatmap FAD - dependent enzymes --------------------------------------------------
FAD <- as.data.frame(TMM) %>% 
  rownames_to_column("Name") %>%
  full_join(geneset) %>%
  dplyr::rename(Ensembl_gtex_iso = Name, 
                Symbol_gtex = Description) %>%
  filter(Symbol_gtex %in% c("ETFDH", "DHODH", "CHDH", "PRODH", "PRODH2", "GPDH", "SQOR", 
                            "SQR", "SQRDL", "CGI-44", "PRO1975",
                            "SDHA", "SDHB", "SDHC", "SDHD")) %>%
  unique()  %>%
  pivot_longer(cols = -c(Symbol_gtex, Ensembl_gtex_iso), names_to = "SAMPID") %>%
  inner_join(meta) %>%
  select(SAMPID, SMTSD,  value, Symbol_gtex) %>%
  unique() 
saveRDS(FAD,here::here("Data", "GTEx", "ProcessedData","FAD_TMM_GTEx.rds") )
test <- FAD %>%
  filter(grepl("Brain - Cortex", SMTSD) | grepl("Liver", SMTSD)) %>%
  pivot_wider(names_from = "Symbol_gtex", values_from = "value") 

exprs <- test[,-c(1:2)] 

library(ComplexHeatmap)
## Row annotation
col_anno = columnAnnotation(
  Tissue=test$SMTSD,
  col=list(Tissue = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100")),
  show_annotation_name = F,
  show_legend =  T,
  simple_anno_size = unit(0.2, "cm")
  # annotation_legend_param = list(nrow=8)
)
range(exprs)
exprs <- t(scale(exprs))
distance    = dist(t(as.matrix(exprs)), method="euclidean" ) 
coldend     = hclust(distance, method="ward.D2")  
row_dist    = dist(as.matrix(exprs), method="euclidean" )
rowdend     = hclust(row_dist, method="ward.D2")

#col_fun = circlize::colorRamp2(c(-4, 0, 4), c("#fff899", "#ffa76b", "#ff4242"))
HM2 <- Heatmap(exprs, 
               name = "Rel. exprs", 
               cluster_columns = coldend,
               #col = col_fun,
               show_column_names =T,
               show_column_dend = T,
               row_title="",
               row_dend_side = "right",
               row_dend_width = unit(0.1, "in"),
               column_dend_height =  unit(0.1, "in"),
               row_names_side = "left",
               show_row_dend=T,
               show_row_names=T,
               top_annotation=col_anno,
               row_names_gp = grid::gpar(fontsize = 5),
               column_names_gp = grid::gpar(fontsize = 6),
               heatmap_legend_param = list(
                 #title = "mtPPS (log10)",
                 title_position = "lefttop-rot",
                 legend_label_gp = element_text(size = 6),
                 legend_height = unit(0.2, "in")
               )
)
png(here::here("Figures", "Figure2", "Suppl_FigS2E_Heatmap_FAD.png"),height=2,width=4,units="in",
    res=1200, 
    pointsize = 6)
draw(HM2, heatmap_legend_side = "right", annotation_legend_side = "right")#, height = unit(0.05, "mm")*1135)
dev.off()


## Suppl. Fig. S2F FAD-score / CI subunits  -----------------------------------------------------------------




FAD_score <- FAD %>%
  filter(grepl("Brain - Cortex", SMTSD) | grepl("Liver", SMTSD)) %>%
  group_by(SMTSD, SAMPID) %>%
  mutate(score = mean(value)) %>%
  select(c(SMTSD, SAMPID, score)) %>%
  unique() 


pathway_subset <- MitoPathwaysTMM %>% 
  filter(grepl("Brain - Cortex", SMTSD) | grepl("Liver", SMTSD)) %>%
  filter(Pathway %in% c("Fatty acid oxidation", 
                        "CI subunits", "CII subunits")) %>%
  select(SAMPID, SMTSD,Pathway, score_simple) %>%
  unique() 

test <- FAD_score %>%
  inner_join(pathway_subset) 

test %>%
  ggplot(aes(y = score, x = score_simple, color = SMTSD)) +
  geom_point(alpha = 0.4) +
  scale_fill_manual(values = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100")) + 
  scale_color_manual(values = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100"))  + 
  theme_classic() +
  ylab("FAD-dependent enzymes score") + 
  xlab("Pathway score") + 
  #  scale_y_continuous(limits = c(0,245), expand = c(0, 0)) +
  theme(axis.text = element_text(size = 6),
        axis.text.x = element_text(angle =45, hjust = 1),
        axis.title.x = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.line = element_line(linewidth = 0.25),
        strip.text = element_text(size = 8),
        strip.background = element_rect(fill = "white"),
        panel.background=element_rect(color = "black", fill = "white"),
        legend.position = "none") + 
  facet_grid(~fct_relevel(Pathway, c("CI subunits",
                                     "CII subunits", "Fatty acid oxidation")), scales = "free")

ggsave(here::here("Figures", "Figure2", "Suppl_FigS2F_Cortex_Liver_FAD.png"), width=5.1, height = 2.29, units = "in", dpi = 1200)

### Stats -------------------------------------------------------------------

test %>%
  select(-SAMPID) %>%
  group_by(SMTSD, Pathway) %>%
  nest() %>%
  mutate(spearman = map_dbl(data, ~cor(.$score, .$score_simple, method = "spearman"))) %>%
  mutate(cortest = map(data, ~tidy(cor.test(.$score, .$score_simple, method = "spearman", exact = FALSE)))) %>%
  mutate(pvalue = map_dbl(cortest, "p.value")) %>%
  arrange(desc(spearman)) %>%
  ungroup() %>%
  mutate(padjust = p.adjust(pvalue, method = "bonferroni"))
## Suppl. Fig. S2G Bivariate CII FAO ------------------------------------------------------------

tmp <- MitoPathwaysTMM %>%
  ungroup() %>%
  filter(Pathway %in% c("Fatty acid oxidation", "Complex II")) %>%
  filter(grepl("Brain - Cortex", SMTSD) | grepl("Liver", SMTSD)) %>%
  select(DonorID, SMTSD, Pathway, score_simple, TRISCH) %>%
  unique() %>%
  mutate(TRISCH = lubridate::hm(TRISCH)) %>%
  mutate(TRISCH = lubridate::hour(TRISCH)*60 + lubridate::minute(TRISCH)) %>%
  pivot_wider(names_from = "Pathway", values_from = "score_simple") %>%
  select(DonorID, SMTSD, `Fatty acid oxidation`, `Complex II`, TRISCH) %>%
  unique() 
tmp %>%
  ggplot(aes(x = `Fatty acid oxidation`, y = `Complex II`,color = SMTSD, fill = SMTSD)) + #color = as.numeric(TRISCH))) + #
  geom_point(alpha = 0.4) +
  scale_fill_manual(values = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100")) + 
  scale_color_manual(values = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100"))  + 
  theme_classic() +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.line = element_line(linewidth = 0.25),
        legend.position = "none")
ggsave(here::here("Figures", "Figure2", "Suppl_Fig2G_Cortex_Liver_FAOCII.png"), width=2, height = 2, units = "in", dpi = 1200)

tmp %>%
  dplyr::rename(Tissue = SMTSD) %>%
  ggplot(aes(x = `Fatty acid oxidation`, y = `Complex II`,color = Tissue, fill = Tissue)) + #color = as.numeric(TRISCH))) + #
  geom_point(alpha = 0.4) +
  scale_fill_manual(values = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100")) + 
  scale_color_manual(values = c("Liver"="#EF5FA7", "Brain - Cortex"= "#1DB100"))  + 
  theme_classic() +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.line = element_line(linewidth = 0.25),
        legend.title = element_text(size =7),
        legend.text = element_text(size = 7))
ggsave(here::here("Figures", "Figure2", "Suppl_FigS2G_Cortex_Liver_FAOCII_legend.png"), width=3, height = 2, units = "in", dpi = 1200)


### Stats -------------------------------------------------------------------

tmp %>%
  select(-DonorID) %>%
  group_by(SMTSD) %>%
  nest() %>%
  mutate(spearman = map_dbl(data, ~cor(.$`Fatty acid oxidation`, .$`Complex II`, method = "spearman"))) %>%
  mutate(cortest = map(data, ~tidy(cor.test(.$`Fatty acid oxidation`, .$`Complex II`, method = "spearman", exact = FALSE)))) %>%
  mutate(pvalue = map_dbl(cortest, "p.value")) %>%
  arrange(desc(spearman))









