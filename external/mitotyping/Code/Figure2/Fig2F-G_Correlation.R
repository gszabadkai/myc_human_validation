rm(list = ls())


# Setup -------------------------------------------------------------------

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggfortify)

color_groups <- c( "CNS"= "#1DB100",
                   "Contractile"="#F8BA00",
                   "Reproductive"="#EE220C",
                   "Digestive"="#265A8C",
                   "Anabolic"="#EF5FA7",
                   "Secretory"="#A8A8A8",
                   "Other"="khaki",
                   "Immune"="#00A89D"
)

color_tissue <- c(
  "Adipose tissue"= "#A8A8A8","Adrenal gland"= "#A8A8A8", "Amygdala"= "#1DB100",             
  "Appendix"="#265A8C", "Basal ganglia"="#1DB100","Bone marrow"="#00A89D", 
  "Breast"="khaki", "Cerebellum"= "#1DB100", "Cerebral cortex"= "#1DB100",        
  "Cervix"="#EE220C", "Choroid plexus"= "#1DB100", "Colon"="#265A8C", "Duodenum" ="#265A8C",           
  "Endometrium"="#EE220C", "Epididymis"="#EE220C", "Esophagus"="#265A8C", "Fallopian tube"="#EE220C",       
  "Gallbladder" ="#265A8C", "Heart muscle"="#F8BA00", "Hippocampal formation"= "#1DB100", 
  "Hypothalamus" = "#1DB100", "Kidney"="#EF5FA7", "Liver"="#EF5FA7","Lung"="khaki",                
  "Lymph node"="#00A89D", "Medulla oblongata"= "#1DB100", "Midbrain"= "#1DB100",
  "Olfactory bulb"= "#1DB100", "Ovary" ="#EE220C","Pancreas"="khaki", 
  "Parathyroid gland"="#A8A8A8", "Pituitary gland"="#A8A8A8","Placenta"="#EE220C",              
  "Pons"= "#1DB100","Prostate"="#EE220C","Rectum"= "#265A8C","Retina"= "khaki",             
  "Salivary gland"="#A8A8A8", "Seminal vesicle" ="#EE220C", "Skeletal muscle"="#F8BA00",       
  "Skin" ="khaki","Small intestine"="#265A8C","Smooth muscle" ="khaki","Spinal cord"= "#1DB100",            
  "Spleen"="#00A89D","Stomach"="#265A8C", "Testis"="#EE220C","Thalamus"= "#1DB100",         
  "Thymus"="#00A89D", "Thyroid gland"="#A8A8A8","Tongue"="khaki","Tonsil"="#00A89D",                
  "Urinary bladder"="khaki", "Vagina"="#EE220C", "White matter"= "#1DB100"         
)

theme_pca <-   theme(
  axis.title = element_text(size = 7),
  axis.text = element_text(size = 6),
  axis.text.x = element_text(angle = 45, hjust = 1),
  legend.text = element_blank(),
  legend.title = element_blank(),
  legend.position = "none",
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank())


# Data preparation --------------------------------------------------------

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
  select(-TSTD3) # has too many NAs
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
processed_data <- processed_data %>%
  tissue_to_group_hm() %>%
  mutate(Tissue = str_to_sentence(Tissue)) #%>%
 #filter(!Tissue %in% "Olfactory bulb") ## test if reproducible without OB (version 21.0 accessed, but 21.1 lacks OB, see also https://www.proteinatlas.org/about/releases#21.0)


# MitoPathways ------------------------------------------------------------

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


# MitoPathway score -------------------------------------------------------


data_pathways <- processed_data %>%
  pivot_longer(cols = -c(Tissue, Group), names_to = "Gene") %>%
  full_join(gene_to_pathway) %>%
  filter(!is.na(Pathway)) %>%
  group_by(Pathway, Tissue) %>%
  mutate(score = mean(value)) %>%
  select(-c(Gene, value)) %>%
  unique() %>%
  pivot_wider(names_from = Pathway, values_from = score)



# Figure 2F Correlation heatmap -------------------------------------------


cormat <- cor(data_pathways[,-c(1:2)], method = "spearman")

distance    = dist(t(cormat), method="euclidean" )
coldend     = hclust(distance, method="ward.D2")
row_dist    = dist(cormat, method="euclidean" )
rowdend     = hclust(row_dist, method="ward.D2")

col_fun = colorRamp2(c(-0.5, 0, 0.5, 1), c("blue", "white", "orange", "red"))
hm <- ComplexHeatmap::Heatmap(cormat, name = "mat",
                              rect_gp = gpar(type = "none"),
                              col = col_fun,
                              cluster_rows = rowdend,
                              cluster_columns = coldend,
                              cell_fun = function(j, i, x, y, w, h, fill) {
                                if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
                                  grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                                }
                              },
                              row_dend_width = unit(0.15, "in"),
                              show_column_dend = F,
                              show_row_names = F,
                              show_column_names = F
                              )

png(here::here("Figures", "Figure2", "Fig2F_CorrHeatmap_hm.png"), width = 4, height = 3.6, units = "in", res = 1200)
draw(hm)

dev.off()


## Figure 2F Example bivariates ------------------------------------------------------


p <-data_pathways %>%
  mutate(ratio = `Complex III` / `Complex V`) %>%
  ggplot(aes(x = `Complex V`, y = `Complex III` , label = ratio)) + 
  geom_point(shape = 21, fill = "#797979", size = 1,stroke = 0.2) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 4))
print(p)
plotly::ggplotly(p)
cor.test(data_pathways$`Complex III`, data_pathways$`Complex IV`, "two.sided","spearman")

ggsave(here::here("Figures", "Figure2", "Fig2F_Bivariate1.png"), width = 1.2, height = 1.15, units = "in", dpi = 1200)
plotly::ggplotly(p)

p <-data_pathways %>%
  mutate(ratio = `mt-tRNA modifications` / `Complex I`) %>%
  ggplot(aes(x = `Complex I`, y = `mt-tRNA modifications` , label = ratio)) + 
  geom_point(shape = 21, fill = "#797979", size = 1,stroke = 0.2) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 4))
print(p)
ggsave(here::here("Figures", "Figure2", "Fig2F_Bivariate2.png"), width = 1.17, height = 1.15, units = "in", dpi = 1200)
plotly::ggplotly(p)
cor.test(data_pathways$`mt-tRNA modifications`, data_pathways$`Complex I`, "two.sided","spearman")


p <-data_pathways %>%
  mutate(ratio = `mtDNA modifications` / `SAM`) %>%
  ggplot(aes(x = `SAM`, y = `mtDNA modifications` , label = ratio)) + 
  geom_point(shape = 21, fill = "#797979", size = 1,stroke = 0.2) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 4))
print(p)
ggsave(here::here("Figures", "Figure2", "Fig2F_Bivariate3.png"), width = 1.2, height = 1.15, units = "in", dpi = 1200)
plotly::ggplotly(p)
cor.test(data_pathways$`mtDNA modifications`, data_pathways$`SAM`, "two.sided","spearman")





# Figure 2G network graph -------------------------------------------------

data_tissue_cor <- processed_data %>%
  pivot_longer(cols = -c(Tissue, Group), names_to = "Gene") %>%
  full_join(gene_to_pathway) %>%
  filter(!is.na(Pathway)) %>%
  group_by(Pathway, Tissue) %>%
  mutate(score = mean(value)) %>%
  select(-c(Gene, value,Group)) %>%
  unique() %>%
  mutate(Tissue = case_when(
    Tissue == "Adipose tissue" ~ "Adipose",
    Tissue == "Adrenal gland" ~ "Adrenal",
    Tissue == "Bone marrow" ~ "Bone mar",
    Tissue == "Cerebral cortex" ~ "Cortex",
    Tissue == "Choroid plexus" ~ "Choroid pl",
    Tissue == "Fallopian tube" ~ "Fallopian",
    Tissue == "Heart muscle" ~ "Heart",
    Tissue == "Hippocampal formation" ~ "Hippoc",
    Tissue == "Medulla oblongata" ~ "Medulla obl",
    Tissue == "Olfactory bulb" ~ "Olf bulb",
    Tissue == "Parathyroid gland" ~ "Parathyroid",
    Tissue == "Pituitary gland" ~ "Pituitary",
    Tissue == "Salivary gland" ~ "Salivary",
    Tissue == "Seminal vesicle" ~ "Seminal v",
    Tissue == "Skeletal muscle" ~ "Sk muscle",
    Tissue == "Small intestine" ~ "Small int",
    Tissue == "Smooth muscle" ~ "Sm muscle",
    Tissue == "Urinary bladder" ~ "Ur bladder",
    Tissue == "White matter" ~ "WM",
    TRUE ~ Tissue 
  )) %>%
  pivot_wider(names_from = Tissue, values_from = score) %>%
  column_to_rownames("Pathway")

library(igraph)
cormat <- as.matrix(as.dist(cor(data_tissue_cor, method = "spearman")))

col_fun = colorRamp2(c(0.6, 0.8, 1), c("yellow", "orange", "red"))
hm <- ComplexHeatmap::Heatmap(cormat, name = "mat",
                              rect_gp = gpar(type = "none"),
                              col = col_fun,
                              cluster_rows = T,#rowdend,
                              cluster_columns = T,#coldend,
                              cell_fun = function(j, i, x, y, w, h, fill) {
                                if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
                                  grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                                }
                              },
                              row_dend_width = unit(0.15, "in"),
                              show_column_dend = F,
                              show_row_names = F,
                              show_column_names = F
)
draw(hm)
png(here::here("Figures", "Figure2", "Fig2G_CorrHeatmap_hm_forNetwork.png"), width = 4, height = 3.6, units = "in", res = 1200)
draw(hm)
dev.off()

g <-graph_from_adjacency_matrix(
  cormat,
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)
g <- delete_edges(g, E(g)[which(E(g)$weight<0.85)])
color_tissue <- as.data.frame(as.factor(V(g))) %>%
  rownames_to_column("Tissue") %>%
  mutate(Group = case_when(
    (Tissue == "Bone mar" | Tissue =="Thymus" | Tissue =="Tonsil" | Tissue =="Lymph node" | Tissue =="Spleen")~"Immune",
    (Tissue == "Placenta" | Tissue == "Ovary" | Tissue == "Seminal v" | Tissue == "Cervix" | Tissue == "Prostate" | 
       Tissue == "Endometrium" | Tissue == "Fallopian" | Tissue == "Vagina" | Tissue =="Epididymis" | Tissue == "Testis") ~ "Reproductive",
    (Tissue == "Gallbladder" | Tissue == "Duodenum" | Tissue == "Small int" | 
       Tissue == "Stomach" | Tissue == "Rectum" | Tissue == "Colon"| Tissue == "Esophagus"| Tissue == "Appendix") ~ "Digestive",
    (Tissue == "Cerebellum" | Tissue ==  "WM" | Tissue ==  "Cortex" | Tissue ==  "Choroid pl" | 
       Tissue ==  "Thalamus" | Tissue ==  "Hypothalamus" | Tissue ==  "Medulla obl" | Tissue ==  "Basal ganglia" | 
       Tissue ==  "Pons" | Tissue ==  "Spinal cord" | Tissue ==  "Midbrain" | Tissue ==  "Amygdala" | Tissue ==  "Hippoc" | 
       Tissue == "Olf bulb") ~"CNS",
    (Tissue == "Heart" | Tissue == "Sk muscle") ~ "Contractile",
    (Tissue == "Adrenal" | Tissue == "Parathyroid" | Tissue == "Pituitary"|
       Tissue == "Salivary" | Tissue == "Thyroid gland" | Tissue == "Adipose" ) ~ "Secretory",
    (Tissue == "Liver" | Tissue == "Kidney") ~ "Anabolic",
    (Tissue == "Lung" | Tissue == "Pancreas"| Tissue == "Skin"| Tissue == "Sm muscle" | 
       Tissue == "Tongue" | Tissue == "Ur bladder"| Tissue == "Breast" | Tissue == "Retina") ~ "Other"), 
    .after = "Tissue") %>%
  mutate(color = case_when(
    Group ==  "CNS"~ "#1DB100",
    Group ==  "Contractile"~"#F8BA00",
    Group ==  "Reproductive"~"#EE220C",
    Group ==  "Digestive"~"#265A8C",
    Group ==  "Anabolic"~"#EF5FA7",
    Group ==  "Secretory"~"#A8A8A8",
    Group ==  "Other"~"khaki",
    Group ==  "Immune"~"#00A89D")) %>%
  select(Tissue, color) 

png(here::here("Figures", "Figure2", "Fig2G_Network_graph.png"), width = 6.6, height = 6.6, units = "in", res = 1200)
plot(g, layout = layout.kamada.kawai, vertex.color = color_tissue$color, vertex.label.cex = 0.4, 
     vertex.size=12, vertex.lwd = 1, edge.width=log10(edge.betweenness(g)), 
     vertex.label.color = "black", vertex.label.family = "Arial")
# print(plot(g, layout = layout.kamada.kawai, vertex.color = color_tissue$color, vertex.label.cex = 0.4, 
#      vertex.size=12, vertex.lwd = 1, edge.width=log10(edge.betweenness(g)), 
#      vertex.label.color = "black", vertex.label.family = "Arial"))
dev.off()


# Bivariate example -------------------------------------------------------

data_tissue_cor %>%
  ggplot(aes(x = log10(Testis), y = log10(Liver))) + 
  geom_point()

p <- data_tissue_cor %>%
  rownames_to_column("Pathway") %>%
  ggplot(aes(x = Cortex, y = Liver, label = Pathway))+ 
  geom_point(shape = 21, fill = "#797979", size = 1,stroke = 0.2) +
  theme_bw() +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 7))
print(p)
ggsave(here::here("Figures", "Figure2", "Fig2G_Bivariate.png"), width = 1.2, height = 1, units = "in" , dpi = 1200)
plotly::ggplotly(p)

cor.test(data_tissue_cor$Liver, data_tissue_cor$Cortex, method = "spearman")
