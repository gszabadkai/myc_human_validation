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


# Figure 3A OxPhos -------------------------------------------------------------------
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
 
 sub <- data_pathways %>%
   select(Tissue, Group, `Complex I`, `Complex II`, `Complex III`, `Complex IV`, `Complex V`, `Fatty acid oxidation`)  %>%
   mutate_if(is.numeric, log10)
 
 p0 <- sub %>%
   ggplot(aes(y = `Complex I`, x = `Fatty acid oxidation`)) +
   geom_point(aes(color = Group, fill = Group),
              alpha = 0.6, size = 2)+
   ylim(min(sub$`Complex I`),max(sub$`Complex I`)) +
   scale_color_manual(values = color_groups) +
   scale_fill_manual(values = color_groups) +
   labs(y = "Complex I (log10)",x = "Fatty acid oxidation (log10)") +
   theme_bw()
 
 
 p1 <- sub %>%
   ggplot(aes(y = `Complex I`, x = `Complex II`)) +
   geom_point(aes(color = Group, fill = Group),
              alpha = 0.6, size = 2)+
   ylim(min(sub$`Complex I`),max(sub$`Complex I`)) +
   scale_color_manual(values = color_groups) +
   scale_fill_manual(values = color_groups) +
   labs(y = "Complex I (log10)",x = "Complex II (log10)") +
   theme_bw() 
 p2 <- sub %>%
   mutate(Group2 = case_when(
     Group == "Brain" ~ "A",
     TRUE ~ "B"
   )) %>%
   ggplot(aes(y = `Complex I`, x = `Complex III`)) +
   geom_point(aes(color = Group, fill = Group),
              alpha = 0.6, size = 2)+
   geom_smooth(aes(fill = Group2), method = "lm", linewidth = 0.3, 
               color = "black", linetype ="dotted", se = F) +
   ylim(min(sub$`Complex I`),max(sub$`Complex I`)) +
   scale_color_manual(values = color_groups) +
   scale_fill_manual(values = color_groups) +
   labs(y = "Complex I (log10)",x = "Complex III (log10)") +
   theme_bw() 

 p3 <- sub %>%
   mutate(Group2 = case_when(
     Group == "Brain" ~ "A",
     TRUE ~ "B"
   )) %>%
   ggplot(aes(y = `Complex I`, x = `Complex IV`)) +
   geom_point(aes(color = Group, fill = Group),
              alpha = 0.6, size = 2)+
   geom_smooth(aes(fill = Group2), method = "lm", linewidth = 0.3,  
               color = "black", linetype ="dotted", se = F) +
   ylim(min(sub$`Complex I`),max(sub$`Complex I`)) +
   scale_color_manual(values = color_groups) +
   scale_fill_manual(values = color_groups) +
   labs(y = "Complex I (log10)",x = "Complex IV (log10)") +
   theme_bw() 
 
 p <- egg::ggarrange(p0 +   
                       theme(
                         legend.position = "none" ,
                         axis.text = element_text(size = 6), 
                         axis.title = element_text(size = 6), 
                         strip.background = element_blank(),
                         strip.text.x = element_blank(),
                         plot.margin = margin(r = 0.5, l = 0.5) ),
                     p1 +   
                     theme(axis.text.y = element_blank(),
                           axis.text = element_text(size = 6), 
                           axis.title = element_text(size = 6), 
                           axis.ticks.y = element_blank(),
                           axis.title.y = element_blank(),
                           strip.background = element_blank(),
                           strip.text.x = element_blank(),
                           legend.position = "none" ,
                           plot.margin = margin(r = 0.5, l = 0.5) ),
                     p2 +  
                       theme(axis.text.y = element_blank(),
                             axis.text = element_text(size = 6), 
                             axis.title = element_text(size = 6), 
                             axis.ticks.y = element_blank(),
                             axis.title.y = element_blank(),
                             strip.background = element_blank(),
                             strip.text.x = element_blank(),
                             legend.position = "none" ,
                             plot.margin = margin(r = 0.5, l = 0.5) ),
                     p3 +   
                       theme(axis.text.y = element_blank(),
                             axis.text = element_text(size = 6), 
                             axis.title = element_text(size = 6), 
                             axis.ticks.y = element_blank(),
                             axis.title.y = element_blank(),
                             strip.background = element_blank(),
                             strip.text.x = element_blank(),
                             legend.position = "none" ,
                             plot.margin = margin(r = 0.5, l = 0.5) ),#+
                      ncol = 4)

 
 png(here::here("Figures", "Figure3",  "Fig3A_CI_FAO_II_III_IV_Grid.png"),width=6,
     height=1.6,
     units="in",res=1200, pointsize = 10)
 print(p)
 dev.off()
 
 
 
 
 # p <- data_pathways %>%
 #   ggplot(aes(x = log10(`Fatty acid oxidation`), y = log10(`Complex I`))) +
 #   geom_point(aes(color = Group, fill = Group, label = Tissue), alpha = 0.6, size = 2) +
 #   theme_bw() +
 #   labs(x = "Fatty acid oxidation (log10)", y = "Complex I (log10)") + 
 #   scale_color_manual(values = color_groups) +
 #   scale_fill_manual(values = color_groups) +
 #   theme_bivariate
 # plotly::ggplotly(p)
 # p
 # ggsave(here::here("Figures", "Figure4", "Figure4B_FAO_CI.png"), width = 1.77, height =1.75, 
 #        units = "in", dpi = 1200, pointsize = 10)
 

# Figure 3B CI CIII Ratio ------------------------------------------------------------

 data_pathways %>%
   select(Tissue, Group, `Complex I`, `Complex II`, `Complex III`, `Complex IV`, `Complex V`)  %>%
   mutate(Group2 = case_when(
     Group == "Brain" ~ "Brain",
     TRUE ~ "all_others"
   )) %>%
   mutate(CI_CIII_ratio = `Complex I` / `Complex III`) %>%
   group_by(Group2) %>%
   mutate(av_ratio = mean(CI_CIII_ratio)) %>%
   select(Group2, av_ratio) %>%
   unique()
 
test <-  data_pathways %>%
   select(Tissue, Group, `Complex I`, `Complex III`)  %>%
   mutate(CICIIIratio = `Complex I` / `Complex III`) %>%
   group_by(Group) %>%
   unique() %>%
   mutate(mean = mean(CICIIIratio)) %>%
   mutate(sd = sd(CICIIIratio)) %>%
   mutate(count = n()) %>%
   mutate(sem = sd/sqrt(count)) %>%
   mutate(ymin = mean-sem) %>%
   mutate(ymax = mean+sem) %>%
   select(-c(Tissue, `Complex I`, `Complex III`, CICIIIratio) )%>% #Group, mean, sem, ymin, ymax) %>% 
   unique() 

p <- test %>%
   ggplot(aes(fill = Group)) +
   geom_segment(aes(x=dplyr::desc(reorder(Group, mean)),
                    xend=dplyr::desc(reorder(Group, mean)),y=1, yend=mean,
                    color = Group), size=3.5, alpha =0.6) +
   geom_errorbar(aes(x=dplyr::desc(reorder(Group, mean)), ymin=ymin, ymax=ymax, color = Group), width=0.4, 
                 alpha=0.8, size=0.3) +
   scale_fill_manual(values = color_groups) +
   scale_color_manual(values = color_groups) +
   labs(y= "CI / CIII") +
   theme_bw() +
   theme(legend.position = "none", 
         axis.title.x = element_blank(),
         axis.text.x = element_blank(),
         axis.text.y = element_text(size = 6),
         axis.title.y = element_text(size = 6),
         axis.ticks.x = element_blank(),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())

 png(here::here("Figures", "Figure3",  "Fig3B_CI_III_ratio.png"),width=1.7,height=1.54,
     units="in",res=1200, pointsize = 10)
 print(p)
 dev.off()

 

## Stats -------------------------------------------------------------------


 
test <-  data_pathways %>%
   select(Tissue, Group, `Complex I`, `Complex II`, `Complex III`, `Complex IV`, `Complex V`)  %>%
   mutate(Group2 = case_when(
     Group == "Brain" ~ "Brain",
     TRUE ~ "all_others"
   )) %>%
   mutate(CI_CIII_ratio = `Complex I` / `Complex III`) %>%
   ungroup() %>%
   select(Group2, CI_CIII_ratio) 
 
 test %>%
   group_by(Group2) %>%
   mutate(mean = mean(CI_CIII_ratio)) %>%
   select(Group2, mean) %>%
   unique() %>%
   pivot_wider(names_from = "Group2", values_from = "mean") %>%
   mutate(FC =Brain / all_others)

test_split <- split(test$CI_CIII_ratio, test$Group2)
lapply(test_split, shapiro.test) # p> 0.05, normally distributed
fligner.test(CI_CIII_ratio ~ Group2, data = test) # p <0.05, variance not equal
t.test(test$CI_CIII_ratio ~ test$Group2, var.equal = F)
## HedgesG
effsize::cohen.d(test$CI_CIII_ratio, test$Group2,pooled=TRUE,paired=FALSE,
        na.rm=FALSE, mu=0, hedges.correction=TRUE,
        conf.level=0.95,noncentral=FALSE, 
        within=TRUE, subject=NA)



# Supplemental Figure S5A CI CII Ratio ------------------------------------------------------------

data_pathways %>%
  select(Tissue, Group, `Complex I`, `Complex II`, `Complex III`, `Complex IV`, `Complex V`)  %>%
  mutate(Group2 = case_when(
    Group == "Brain" ~ "Brain",
    TRUE ~ "all_others"
  )) %>%
  mutate(CI_CII_ratio = `Complex I` / `Complex II`) %>%
  group_by(Group2) %>%
  mutate(av_ratio = mean(CI_CII_ratio)) %>%
  select(Group2, av_ratio) %>%
  unique()

test <-  data_pathways %>%
  select(Tissue, Group, `Complex I`, `Complex II`)  %>%
  mutate(CICIIratio = `Complex I` / `Complex II`) %>%
  group_by(Group) %>%
  unique() %>%
  mutate(mean = mean(CICIIratio)) %>%
  mutate(sd = sd(CICIIratio)) %>%
  mutate(count = n()) %>%
  mutate(sem = sd/sqrt(count)) %>%
  mutate(ymin = mean-sem) %>%
  mutate(ymax = mean+sem) %>%
  select(-c(Tissue, `Complex I`, `Complex II`, CICIIratio) )%>% #Group, mean, sem, ymin, ymax) %>% 
  unique() 

p <- test %>%
  ggplot(aes(fill = Group)) +
  geom_segment(aes(x=dplyr::desc(reorder(Group, mean)),
                   xend=dplyr::desc(reorder(Group, mean)),y=0, yend=mean,
                   color = Group), size=3.5, alpha =0.6) +
  geom_errorbar(aes(x=dplyr::desc(reorder(Group, mean)), ymin=ymin, ymax=ymax, color = Group), width=0.4, 
                alpha=0.8, size=0.3) +
  scale_fill_manual(values = color_groups) +
  scale_color_manual(values = color_groups) +
  labs(y= "CI / CII") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

png(here::here("Figures", "Figure3",  "Suppl_FigS5A_CI_II_ratio.png"),width=1.7,height=1.54,
    units="in",res=1200, pointsize = 10)
print(p)
dev.off()



## Stats -------------------------------------------------------------------



test <-  data_pathways %>%
  select(Tissue, Group, `Complex I`, `Complex II`, `Complex III`, `Complex IV`, `Complex V`)  %>%
  mutate(Group2 = case_when(
    Group == "Brain" ~ "Brain",
    TRUE ~ "all_others"
  )) %>%
  mutate(CI_CII_ratio = `Complex I` / `Complex II`) %>%
  ungroup() %>%
  select(Group2, CI_CII_ratio) 

test_split <- split(test$CI_CII_ratio, test$Group2)
lapply(test_split, shapiro.test) # p> 0.05, normally distributed
fligner.test(CI_CII_ratio ~ Group2, data = test) # p >0.05, variance equal
t.test(test$CI_CII_ratio ~ test$Group2, var.equal = T)
## HedgesG
effsize::cohen.d(test$CI_CII_ratio, test$Group2,pooled=TRUE,paired=FALSE,
                 na.rm=FALSE, mu=0, hedges.correction=TRUE,
                 conf.level=0.95,noncentral=FALSE, 
                 within=TRUE, subject=NA)



# Supplemental Figure S5B CI CIV Ratio ------------------------------------------------------------

data_pathways %>%
  select(Tissue, Group, `Complex I`, `Complex II`, `Complex III`, `Complex IV`, `Complex V`)  %>%
  mutate(Group2 = case_when(
    Group == "Brain" ~ "Brain",
    TRUE ~ "all_others"
  )) %>%
  mutate(CI_CIV_ratio = `Complex I` / `Complex IV`) %>%
  group_by(Group2) %>%
  mutate(av_ratio = mean(CI_CIV_ratio)) %>%
  select(Group2, av_ratio) %>%
  unique()

test <-  data_pathways %>%
  select(Tissue, Group, `Complex I`, `Complex IV`)  %>%
  mutate(CICIVratio = `Complex I` / `Complex IV`) %>%
  group_by(Group) %>%
  unique() %>%
  mutate(mean = mean(CICIVratio)) %>%
  mutate(sd = sd(CICIVratio)) %>%
  mutate(count = n()) %>%
  mutate(sem = sd/sqrt(count)) %>%
  mutate(ymin = mean-sem) %>%
  mutate(ymax = mean+sem) %>%
  select(-c(Tissue, `Complex I`, `Complex IV`, CICIVratio) )%>% #Group, mean, sem, ymin, ymax) %>% 
  unique() 

p <- test %>%
  ggplot(aes(fill = Group)) +
  geom_segment(aes(x=fct_relevel(Group, c("Brain", "Other", "Reproductive",
                                          "Secretory", "Contractile", "Anabolic", 
                                          "Digestive", "Immune")),
                   xend=fct_relevel(Group, c("Brain", "Other", "Reproductive",
                                             "Secretory", "Contractile", "Anabolic", 
                                             "Digestive", "Immune")),y=0, yend=mean,
                   color = Group), size=3.5, alpha =0.6) +
  geom_errorbar(aes(x=fct_relevel(Group, c("Brain", "Other", "Reproductive",
                                           "Secretory", "Contractile", "Anabolic", 
                                           "Digestive", "Immune")), ymin=ymin, ymax=ymax, color = Group), width=0.4, 
                alpha=0.8, size=0.3) +
  scale_fill_manual(values = color_groups) +
  scale_color_manual(values = color_groups) +
  labs(y= "CI / CIV") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

png(here::here("Figures", "Figure3",  "Suppl_FigS5B_CI_IV_ratio.png"),width=1.7,height=1.54,
    units="in",res=1200, pointsize = 10)
print(p)
dev.off()



## Stats -------------------------------------------------------------------



test <-  data_pathways %>%
  select(Tissue, Group, `Complex I`, `Complex II`, `Complex III`, `Complex IV`, `Complex V`)  %>%
  mutate(Group2 = case_when(
    Group == "Brain" ~ "Brain",
    TRUE ~ "all_others"
  )) %>%
  mutate(CI_CIV_ratio = `Complex I` / `Complex IV`) %>%
  ungroup() %>%
  select(Group2, CI_CIV_ratio) 

test_split <- split(test$CI_CIV_ratio, test$Group2)
lapply(test_split, shapiro.test) # p> 0.05, normally distributed
fligner.test(CI_CIV_ratio ~ Group2, data = test) # p <0.05, variance not equal
t.test(test$CI_CIV_ratio ~ test$Group2, var.equal = F)

## HedgesG
effsize::cohen.d(test$CI_CIV_ratio, test$Group2,pooled=TRUE,paired=FALSE,
                 na.rm=FALSE, mu=0, hedges.correction=TRUE,
                 conf.level=0.95,noncentral=FALSE, 
                 within=TRUE, subject=NA)

# Supplemental Figure 4B CI CV Ratio ------------------------------------------------------------

data_pathways %>%
  #pivot_wider(names_from = "Pathway", values_from = "score") %>%
  select(Tissue, Group, `Complex I`, `Complex II`, `Complex III`, `Complex IV`, `Complex V`)  %>%
  mutate(Group2 = case_when(
    Group == "Brain" ~ "Brain",
    TRUE ~ "all_others"
  )) %>%
  mutate(CI_CV_ratio = `Complex I` / `Complex V`) %>%
  group_by(Group2) %>%
  mutate(av_ratio = mean(CI_CV_ratio)) %>%
  select(Group2, av_ratio) %>%
  unique()

test <-  data_pathways %>%
#  pivot_wider(names_from = "Pathway", values_from = "score") %>%
  select(Tissue, Group, `Complex I`, `Complex V`)  %>%
  mutate(CICVratio = `Complex I` / `Complex V`) %>%
  group_by(Group) %>%
  unique() %>%
  mutate(mean = mean(CICVratio)) %>%
  mutate(sd = sd(CICVratio)) %>%
  mutate(count = n()) %>%
  mutate(sem = sd/sqrt(count)) %>%
  mutate(ymin = mean-sem) %>%
  mutate(ymax = mean+sem) %>%
  select(-c(Tissue, `Complex I`, `Complex V`, CICVratio) )%>% #Group, mean, sem, ymin, ymax) %>% 
  unique() 

p <- test %>%
  ggplot(aes(fill = Group)) +
  geom_segment(aes(x=dplyr::desc(reorder(Group, mean)),
                   xend=dplyr::desc(reorder(Group, mean)),y=0, yend=mean,
                   color = Group), size=3.5, alpha =0.6) +
  geom_errorbar(aes(x=dplyr::desc(reorder(Group, mean)), ymin=ymin, ymax=ymax, color = Group), width=0.4, 
                alpha=0.8, size=0.3) +
  scale_fill_manual(values = color_groups) +
  scale_color_manual(values = color_groups) +
  labs(y= "CI / CV") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

png(here::here("Figures", "Figure3",  "Suppl_FigS5C_CI_V_ratio.png"),width=1.7,height=1.54,
    units="in",res=1200, pointsize = 10)
print(p)
dev.off()



## Stats -------------------------------------------------------------------



test <-  data_pathways %>%
  select(Tissue, Group, `Complex I`, `Complex II`, `Complex III`, `Complex IV`, `Complex V`)  %>%
  mutate(Group2 = case_when(
    Group == "Brain" ~ "Brain",
    TRUE ~ "all_others"
  )) %>%
  mutate(CI_CV_ratio = `Complex I` / `Complex V`) %>%
  ungroup() %>%
  select(Group2, CI_CV_ratio) 

test_split <- split(test$CI_CV_ratio, test$Group2)
lapply(test_split, shapiro.test) # p< 0.05, not normally distributed
fligner.test(CI_CV_ratio ~ Group2, data = test) # p <0.05, variance not equal
brunnermunzel::brunnermunzel.test(test_split[[1]], test_split[[2]], perm = T)
## HedgesG
effsize::cohen.d(test$CI_CV_ratio, test$Group2,pooled=TRUE,paired=FALSE,
                 na.rm=FALSE, mu=0, hedges.correction=TRUE,
                 conf.level=0.95,noncentral=FALSE, 
                 within=TRUE, subject=NA)



# Supplemental Figure S5D CI ratios across all pathways --------------------

## Pathway ratios ----------------------------------------------------------
data_pathways <- processed_data %>%
  pivot_longer(cols = colnames(processed_data[3:ncol(processed_data)]), names_to = "Gene") %>%
  full_join(gene_to_pathway, by = "Gene") %>%
  filter(!is.na(Pathway)) %>%
  unique()%>%
  group_by(Tissue, Pathway) %>%
  mutate(score = mean(value, na.rm = T)) %>%
  select(-c(Gene, value)) %>%
  unique() 
all_pathways <- unique(data_pathways$Pathway)

pw_pairs <- data.frame(pw1 = rep(all_pathways, each = 149)) %>%
  group_by(pw1) %>%
  nest() %>%
  mutate(data = map(data, ~mutate(.x, pw2 = unique(all_pathways)))) %>%
  unnest() %>%
  filter(!pw1==pw2) 
colnames(pw_pairs) <- c("pw1", "pw2")

janitor::tabyl(data_pathways, Pathway)
data_added1 <- pw_pairs %>%
  ungroup() %>%
  dplyr::select(pw1) %>%
  dplyr::rename(Pathway = pw1) %>%
  inner_join(data_pathways, by = "Pathway") %>%
  unique() %>%
  dplyr::rename(value1 = score, Pathway1 = Pathway) %>%
  dplyr::slice(rep(row_number(), 148)) %>%
  arrange(Pathway1) %>%
  group_by(Tissue, Group) %>%
  nest() %>%
  dplyr::rename(pw1 = data) 

data_added2 <- pw_pairs %>%
  ungroup() %>%
  dplyr::select(pw2) %>%
  dplyr::rename(Pathway = pw2) %>%
  inner_join(data_pathways, by = "Pathway") %>%
  unique() %>%
  dplyr::rename(value2 = score, Pathway2 = Pathway) %>%
  dplyr::slice(rep(row_number(), 148)) %>%
  arrange(desc(Pathway2)) %>%
  unique() %>%
  group_by(Tissue, Group) %>%
  nest() %>%
  dplyr::rename(pw2 = data) 

data_ratios =full_join(data_added1, data_added2, by = c("Tissue", "Group")) %>%
  mutate(combined =  map2(pw1, pw2, ~ cbind(.x,  .y))) %>%
  select(-c(pw1, pw2)) %>%
  unnest(cols = combined) %>%
  unique() %>%
  mutate(ratio = value1 / value2)


## Plot CI ratios ----------------------------------------------------------

p <- data_ratios %>%
  filter(Pathway1 %in% "Complex I") %>%
  ungroup() %>%
  mutate(mean_ratio = mean(ratio)) %>%
  ggplot(aes(x = dplyr::desc(reorder(Tissue, ratio)), y= ratio, color = Group, label = Pathway2)) +
  geom_point(alpha = 0.7) +
  geom_hline(aes(yintercept = mean_ratio), linetype = "dotted") + 
  labs(y = "CI ratio to all other pathways") + 
  scale_color_manual(values = color_groups) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 6),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 7))

png(here::here("Figures", "Figure3",  "Suppl_FigS5D_CI_vs_all_other.png"),width=8.5,height=2,
    units="in",res=1200, pointsize = 8)
print(p)
dev.off()




# Supplemental Figure 4E Average CI ratio per tissue group ----------------


test <-  data_ratios %>%
  filter(Pathway1 %in% "Complex I") %>%
  group_by(Tissue) %>%
  mutate(mean = mean(ratio)) %>%
  mutate(sd = sd(ratio)) %>%
  mutate(count = n()) %>%
  mutate(sem = sd/sqrt(count)) %>%
  mutate(ymin = mean-sem) %>%
  mutate(ymax = mean+sem) %>%
  select(-c(ratio, value1, value2, Pathway1, Pathway2) )%>% #Group, mean, sem, ymin, ymax) %>% 
  unique() 

p <- test %>%
  ggplot(aes(fill = Group)) +
  geom_segment(aes(x=dplyr::desc(reorder(Tissue, mean)),
                   xend=dplyr::desc(reorder(Tissue, mean)),y=0, yend=mean,
                   color = Group), size=3, alpha =0.6) +
  geom_errorbar(aes(x=dplyr::desc(reorder(Tissue, mean)), ymin=ymin, ymax=ymax, color = Group), width=0.4, 
                alpha=0.8, size=0.3) +
  geom_hline(aes(yintercept = mean(mean)), linetype = "dotted" ) +
  scale_fill_manual(values = color_groups) +
  scale_color_manual(values = color_groups) +
  labs(y= "Average CI ratios") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

png(here::here("Figures", "Figure3",  "Suppl_FigS5E_Average_CIratios.png"),width=7.23,height=2,
    units="in",res=1200, pointsize = 8)
print(p)
dev.off()


## stats -------------------------------------------------------------------

test %>%
  mutate(Group2 = case_when(Group == "Brain" ~ "Above_average",
                            Tissue == "Heart" ~ "Above_average",
                            TRUE ~ "All_others")) %>%
  ungroup() %>%
  select(Group2, mean) %>%
  group_by(Group2) %>%
  mutate(Groupmean = mean(mean)) %>%
  select(-mean) %>%
  unique() %>%
  pivot_wider(names_from = "Group2", values_from = "Groupmean") %>%
  mutate(fc = Above_average / All_others)

         