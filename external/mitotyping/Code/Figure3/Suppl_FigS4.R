rm(list = ls())

# Setup -------------------------------------------------------------------

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggfortify)

color_groups <- c( "Brain"= "#1DB100",
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
       Tissue == "olfactory bulb") ~"Brain",
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
  mutate(Tissue = str_to_sentence(Tissue)) #%>%
#filter(!Tissue %in% "Olfactory bulb") ## test if reproducible without OB (version 21.0 accessed, but 21.1 lacks OB, see also https://www.proteinatlas.org/about/releases#21.0)
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


# Suppl Figure S4A Fold change ---------------------------------------------------


## Combine data
data_pathways <- processed_data %>%
  pivot_longer(cols = colnames(processed_data[3:ncol(processed_data)]), names_to = "Gene") %>%
  full_join(gene_to_pathway, by = "Gene") %>%
  full_join(levels) %>%
  filter(!is.na(Pathway)) %>%
  unique() %>%
  
  filter(Group %in% c("Brain", "Anabolic")) %>%
  filter(!Tissue %in% "Kidney") %>%
  select(-c(Level)) %>%
  mutate(Group = case_when(Group == "Anabolic" ~ "Liver", 
                           TRUE ~ "Brain")) %>%
  group_by(Pathway,Group) %>%
  mutate(mean = mean(value)) %>%
  select(-c(Gene, value, Tissue)) %>%
  unique() %>%
  pivot_wider(names_from = "Group", values_from = "mean") %>%
  select(Pathway_Level1, Pathway_Level2, Pathway_Level3, Brain, Liver)

fc <- data_pathways  %>%
  mutate(log2_fc = log2(Brain/Liver)) %>%
  mutate(test = Liver / Brain)
### Calculate log2 fc
fc <- data_pathways  %>%
  mutate(log2_fc = log2(Brain/Liver)) %>%
  arrange(log2_fc)%>% #arrange data ascending
  column_to_rownames("Pathway") %>%
  mutate(xaxis = seq(1:149)) %>% #add x-axis to sort ascending from left to right
  rownames_to_column("Pathway") %>%
  select(Pathway, Pathway_Level1, Pathway_Level3, log2_fc, xaxis)
p <- fc %>% 
  ggplot(aes(x =xaxis, y = log2_fc, color = Pathway_Level1, label = Pathway_Level3 )) +
  geom_point(alpha = 0.7, size = 2, shape = 18) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(y = "Log2 fold change Brain to Liver", x = "Ranked Mitopathway scores") +
  theme_bw() +
  geom_vline(aes(xintercept =123), linetype = "dotted") +
  scale_color_brewer(palette = "Dark2", name = "Pathway group")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(size=6),
        legend.position = "none",
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size=6),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

png(here::here("Figures", "Figure3", "Suppl_FigS4A_Brain_Liver.png"),width=3.4,height=1.7,
    units="in",res=1200, pointsize = 10) 

print(p)
dev.off()

plotly::ggplotly(p)


# Suppl Figure S4B-C Bivariates --------------------------------------------------


data_pathways <- processed_data %>%
  pivot_longer(cols = colnames(processed_data[3:ncol(processed_data)]), names_to = "Gene") %>%
  full_join(gene_to_pathway, by = "Gene") %>%
  filter(!is.na(Pathway)) %>%
  unique()%>%
  group_by(Tissue, Pathway) %>%
  mutate(score = mean(value, na.rm = T)) %>%
  select(-c(Gene, value)) %>%
  unique() %>%
  pivot_wider(names_from = "Pathway", values_from = "score")

theme_bivariate <-   theme(
  legend.title = element_blank(),
  legend.position = "none",
  axis.title = element_text(size = 6), 
  axis.text = element_text(size =6)
)


## B 
p <- data_pathways %>%
  ggplot(aes(x = log10(`Glycerol phosphate shuttle`), y = log10(`Urea cycle`))) +
  geom_point(aes(color = Group, fill = Group, label = Tissue), alpha = 0.6, size = 2) +
  theme_bw() +
  labs(x = "G3P shuttle (log10)", y = "Urea cycle (log10)") + 
  scale_color_manual(values = color_groups) +
  scale_fill_manual(values = color_groups) +
  theme_bivariate

print(p)
ggsave(here::here("Figures", "Figure3", "Suppl_FigS4B_G3P_Urea.png"), width = 1.77, height =1.75, 
       units = "in", dpi = 1200, pointsize = 10)

data_pathways %>%
  select(Tissue, `Urea cycle`, `Glycerol phosphate shuttle`,
         `Serine metabolism`, `Tetrahydrobiopterin synthesis`) %>%
  mutate(Group2 = case_when(Tissue == "Liver" ~ "Liver",
                            TRUE ~ "All_others")) %>%
  pivot_longer(cols = -c(Group2, Tissue), names_to = "Pathway") %>%
  group_by(Group2, Pathway) %>%
  mutate(average = mean(value),
         min = min(value),
         max = max(value)) %>%
  select(Pathway, Group2, average, min, max) %>%
  unique() %>%
  pivot_longer(cols = -c(Pathway, Group2)) %>% 
  pivot_wider(names_from = "Group2", values_from = "value") %>%
  mutate(FC_liver_all_others = Liver / All_others) %>%
  mutate(FC_all_others_liver = All_others / Liver)
## C 
p <- data_pathways %>%
  ggplot(aes(x = log10(`Serine metabolism`), y = log10(`Tetrahydrobiopterin synthesis`))) +
  geom_point(aes(color = Group, fill = Group, label = Tissue), alpha = 0.6, size = 2) +
  theme_bw() +
  labs(x = "Serine metabolism (log10)", y = "BH4 synthesis (log10)") + 
  scale_color_manual(values = color_groups) +
  scale_fill_manual(values = color_groups) +
  theme_bivariate
print(p)
ggsave(here::here("Figures", "Figure3", "Suppl_FigS4C_Serine_BH4.png"), width = 1.77, height =1.75, 
       units = "in", dpi = 1200, pointsize = 10)


# Suppl. Figure S4D Radar chart ------------------------------------------------------

pathways <- data_pathways %>%
  pivot_longer(cols = -c(Tissue, Group), names_to = "Pathway", values_to = "value") %>%
  group_by(Pathway) %>%
  mutate(zscore = (value - mean(value)) / sd(value)) %>%
  select(-value) %>%
  arrange(Tissue, desc(zscore)) %>%
  group_by(Tissue) %>%
  dplyr::slice(1:3) %>%
  filter(Tissue %in% c("Liver", "Skeletal muscle", "Cerebral cortex", "Small intestine")) %>%
  pull(Pathway)

data <- data_pathways %>%
  pivot_longer(cols = -c(Tissue, Group), names_to = "Pathway", values_to = "value") %>%
  group_by(Pathway) %>%
  mutate(zscore = (value - mean(value)) / sd(value)) %>%
  select(-value) %>%
  filter(Pathway %in% pathways) %>%
  pivot_wider(names_from = "Pathway", values_from = "zscore") %>%
  select(-Group) %>%
  filter(Tissue %in% c("Liver", "Skeletal muscle", "Cerebral cortex", "Small intestine")) 

colMax <- function(data) sapply(data, max, na.rm = TRUE)
colMin <- function(data) sapply(data, min, na.rm = TRUE)
upper <- max(colMax(data[,2:ncol(data)]))
upper <- upper+upper/100
lower <- min(colMin(data[,2:ncol(data)]))
lower <- lower+lower/100
data <- rbind(rep(upper,3) , rep(lower,3),data) %>%
  column_to_rownames("Tissue")
feature_order <- c( "mt-mRNA modifications","Complex I", "CI subunits",
                    "Serine metabolism", "Urea cycle", "Xenobiotic metabolism", 
                    "Heme synthesis and processing",
                                  "Proline metabolism", "Vitamin A metabolism",     
                                  "mtDNA stability and decay", "Coenzyme Q metabolism", 
                                  "Intramitochondrial membrane interactions" 
)
data  <- data [, feature_order]
library(fmsb)

colors_border=c(
                 rgb(118/255, 206/255, 102/255,0.9),
                 rgb(245/255, 142/255, 192/255,0.9),
                 rgb(250/255, 209/255, 107/255,0.9) ,
                 rgb(52/255, 89/255, 136/255,0.9))
colors_in=c(
             rgb(118/255, 206/255, 102/255,0.4),
             rgb(245/255, 142/255, 192/255,0.4),
             rgb(250/255,209/255, 107/255,0.4),
             rgb(52/255, 89/255, 136/255,0.4))

png("Figures/Figure3/Suppl_FigS4D_radar.png", width = 6, height = 6, units = "in", res = 1200)
radarchart(data  , axistype=0, pcol = colors_border, pfcol = colors_in,
           plwd=c(0.5,1.5),  pty=16,plty=c(1,3),
           cglcol="grey",cglwd = 1, cglty=3, axislabcol="grey",
           vlcex=0.8, title="xxx", calcex=0.8, seg = 3
) 
legend(x=1.1, y=0.75, legend = rownames(data[-c(1,2),]), bty = "n", 
       pch=16 , col=colors_border , text.col = "black", cex=1, pt.cex=1.5)
  dev.off()

  png("Figures/Figure3/Suppl_FigS4D_radar_noLabels.png", width = 6, height = 6, units = "in", res = 1200)
  radarchart(data  , axistype=0, pcol = colors_border, pfcol = colors_in,
           plwd=c(0.5,1.5),  pty=16,plty=c(1,3),
           cglcol="grey",cglwd = 1, cglty=3, axislabcol="grey",
           vlcex=0.001, title="xxx", calcex=0.8, seg = 3
) 
legend(x=1.1, y=0.75, legend = rownames(data[-c(1,2),]), bty = "n", 
       pch=16 , col=colors_border , text.col = "black", cex=1, pt.cex=1.5)
 dev.off()



# Suppl Figure S4E barplots ------------------------------------------------------

 png(here::here("Figures", "Figure3", "Suppl_FigS4E_barplots.png"), width = 3, 
     height = 2.3,  unit = "in", res = 1200)
 
 p <- processed_data %>%
   pivot_longer(cols = colnames(processed_data[3:ncol(processed_data)]), names_to = "Gene") %>%
   full_join(gene_to_pathway, by = "Gene") %>%
   na.omit() %>%
   unique() %>%
   group_by(Pathway, Tissue) %>%
   mutate(score = mean(value)) %>%
   group_by(Pathway) %>%
   mutate(mean_score = mean(value)) %>%
   select(-c(value, Gene, Group)) %>%
   unique() %>%
   filter(Pathway %in% c("Complex I", "Serine metabolism", "Vitamin A metabolism", "Coenzyme Q metabolism")) %>%
   filter(Tissue %in% c("Cerebral cortex", "Liver", "Skeletal muscle", "Small intestine")) %>%
   pivot_wider(names_from = "Tissue", values_from = "score") %>%
   pivot_longer(cols = c("mean_score", "Cerebral cortex", "Liver", "Skeletal muscle", "Small intestine"), names_to = "Tissue") %>%
   mutate(Pathway = case_when(
     Pathway == "Complex I" ~ "Complex I",
     Pathway == "Serine metabolism" ~ "Serine metab.",
     Pathway == "Vitamin A metabolism" ~ "VitA metab.",
     Pathway == "Coenzyme Q metabolism" ~ "CoQ metab."
   )) %>%
   ggplot(aes(x = fct_relevel(Tissue, c("mean_score", "Cerebral cortex", "Skeletal muscle", "Liver", "Small intestine")), 
              y = value, Group = Pathway)) +
   geom_segment(aes(y=1, yend=value, xend = Tissue,
                    color = Tissue), size=2, alpha =0.6) +
   scale_color_manual(values = c("mean_score" = "grey", 
                                 "Cerebral cortex" = "#1DB100", 
                                 "Liver" = "#EF5FA7", 
                                 "Skeletal muscle" = "#F8BA00", 
                                 "Small intestine"="#265A8C")) + 
   labs(y = "Pathway score", x = "Tissue") +
     facet_wrap(~Pathway, scales= "free") + 
   theme_bw() +
   theme(axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.title = element_text(size = 7),
         strip.text = element_text(size = 6),
         axis.text.y = element_text(size = 6),
         legend.position = "none")
 print(p)
 dev.off()
 

