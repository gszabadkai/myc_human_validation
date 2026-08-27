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

mitoPPS <- read_csv(here::here("Data", "HumanProteinAtlas","ProcessedData", "mitoPPS_all_tissues.csv" ))

# Run PCA
pca <- prcomp(mitoPPS[, -c(1:2)], scale. = T)
summary(pca)


# Figure 3D PC1 PC2 ---------------------------------------------------------------

p <- autoplot(pca, data =mitoPPS, elipse = TRUE, colour = 'Tissue',loadings = F, loadings.label=F,
              loadings.color = "black", #scale = F,
              fill = 'Group', size = 2.8, alpha = 0.6, shape = 21)+
  scale_color_manual(values = color_tissue) +
  scale_fill_manual(values = color_groups) +
  theme_bw() +
  theme_pca 

png(here::here("Figures", "Figure3", "Fig3D_PC1_PC2.png"),width=2,height=2,units="in",res=1200)
print(p)
dev.off()

df <- as.data.frame(p$data)%>% 
  select(PC1, PC2, Group, Tissue) %>%
  group_by(Group) %>%
  mutate(PC1_average = mean(PC1), PC2_average = mean(PC2))
# Figure 3D PC1 PC3 ---------------------------------------------------------------

p <- autoplot(pca, x = 3, y =1, data =mitoPPS, elipse = TRUE, colour = 'Tissue',
              loadings = F, loadings.label=F, #scale = F,
              loadings.color = "black", 
              fill = 'Group', size = 2.8, alpha = 0.6, shape = 21)+
  scale_color_manual(values = color_tissue) +
  scale_fill_manual(values = color_groups) +
  theme_bw() +
  theme_pca + scale_x_reverse()
p
png(here::here("Figures", "Figure3", "Fig3D_PC1_PC3.png"),width=2,height=2,units="in",res=1200)
print(p)
dev.off()

df <- as.data.frame(p$data)%>% 
  select(PC1, PC3, Group, Tissue) %>%
  group_by(Group) %>%
  mutate(PC1_average = mean(PC1), PC3_average = mean(PC3))


# Figure 3D 3D plot -----------------------------------------------------------------
## Commented out because it requires developer tools to be installed
# summary_pca <- summary(pca)
# var_1 <- round(summary_pca$importance[2,1]*100,2)
# var_2 <- round(summary_pca$importance[2,2]*100,2)
# var_3 <- round(summary_pca$importance[2,3]*100,2)
# group_color_df <- data.frame(Color = color_groups) %>%
#   rownames_to_column("Group")
# df <- pca$x
# df <- data.frame(PC1=df[,1], PC2=df[,2], PC3=df[,3], Tissue = mitoPPS$Tissue, Group=as.factor(mitoPPS$Group)) %>%
#   full_join(group_color_df, by = "Group") %>%
#   na.omit() 
# 
# rgl::par3d(cex=0.001)
# rgl::bg3d(color = "white")
# with(df, rgl::plot3d(PC1,PC2,PC3, col= Color, alpha = 0.6, size = 15, type = "p", 
#                      xlab = "", ylab = "", zlab = "",
#                      ann = FALSE, axes = FALSE
# ))
# rgl::box3d()


# Figure 3E loadings ------------------------------------------------------


loadings <- pca$rotation
write_csv(as.data.frame(loadings) %>% rownames_to_column("Pathway"), 
          here::here("Supplemental/Figure3E_loadings.csv"))
PC1_top <- as.data.frame(loadings) %>%
  arrange(dplyr::desc(PC1)) %>%
  dplyr::slice(1:5) 
PC1_bottom <- as.data.frame(loadings) %>%
  arrange(PC1) %>%
  dplyr::slice(1:5) 

PC2_top <- as.data.frame(loadings) %>%
  arrange(dplyr::desc(PC2)) %>%
  dplyr::slice(1:5) 
PC2_bottom <- as.data.frame(loadings) %>%
  arrange(PC2) %>%
  dplyr::slice(1:5) 

PC3_top <- as.data.frame(loadings) %>%
  arrange(dplyr::desc(PC3)) %>%
  dplyr::slice(1:5) 
PC3_bottom <- as.data.frame(loadings) %>%
  arrange(PC3) %>%
  dplyr:: slice(1:5) 



## Color-code by level1

mitocarta <- readxl::read_xls(here::here("Data","MitoCarta", "OriginalData",
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

levels <- readxl::read_xls(here::here("Data","MitoCarta","OriginalData","HumanMitoCarta3_0.xls"), sheet = 4) %>%
  select(MitoPathway, `MitoPathways Hierarchy`) %>%
  separate(`MitoPathways Hierarchy`, into = c("Pathway_Level1", "Pathway_Level2", "Pathway_Level3"), sep = " > ") %>%
  mutate(Level = case_when(
    MitoPathway == Pathway_Level1 ~ "Pathway_Level1", 
    MitoPathway == Pathway_Level2 ~ "Pathway_Level2", 
    MitoPathway == Pathway_Level3 ~ "Pathway_Level3"
  )) %>%
  unique() %>%
  dplyr::rename(Pathway = MitoPathway) 



p <- as.data.frame(loadings) %>%
  rownames_to_column("Pathway") %>%
  select(Pathway, PC1) %>%
  full_join(levels) %>%
  ggplot(aes(x = dplyr::desc(reorder(Pathway, PC1)), y = PC1, fill = Pathway_Level1, color = Pathway_Level1)) +
  geom_bar(stat = "identity", alpha =0.8, size =0.3, width = 0.4)+ coord_flip() +
  scale_fill_brewer(palette = "Dark2", name = "Pathway group")+
  scale_color_brewer(palette = "Dark2", name = "Pathway group")+
  theme_classic() +
  geom_hline(yintercept = 0) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none",
        axis.text =  element_text(size = 6),
        axis.title =  element_text(size = 7))
png(here::here("Figures", "Figure3", "Fig3E_PC1_loadings.png"), width = 1.5, height = 2, units = "in", res = 1200)
print(p)
dev.off()


p <- as.data.frame(loadings) %>%
  rownames_to_column("Pathway") %>%
  select(Pathway, PC2) %>%
  full_join(levels) %>%
  ggplot(aes(x = dplyr::desc(reorder(Pathway, PC2)), y = PC2, fill = Pathway_Level1, color = Pathway_Level1)) +
  geom_bar(stat = "identity", alpha =0.8, size =0.3, width = 0.4)+ coord_flip() +
  scale_fill_brewer(palette = "Dark2", name = "Pathway group")+
  scale_color_brewer(palette = "Dark2", name = "Pathway group")+
  theme_classic() +
  geom_hline(yintercept = 0) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none",
        axis.text =  element_text(size = 6),
        axis.title =  element_text(size = 7))
png(here::here("Figures", "Figure3", "Fig3E_PC2_loadings.png"), width = 1.5, height = 2, units = "in", res = 1200)
print(p)
dev.off()

p <- as.data.frame(loadings) %>%
  rownames_to_column("Pathway") %>%
  select(Pathway, PC3) %>%
  full_join(levels) %>%
  ggplot(aes(x = dplyr::desc(reorder(Pathway, PC3)), y = PC3, fill = Pathway_Level1, color = Pathway_Level1)) +
  geom_bar(stat = "identity", alpha =0.8, size =0.3, width = 0.4)+ coord_flip() +
  scale_fill_brewer(palette = "Dark2", name = "Pathway group")+
  scale_color_brewer(palette = "Dark2", name = "Pathway group")+
  theme_classic() +
  geom_hline(yintercept = 0) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none",
        axis.text =  element_text(size = 6),
        axis.title =  element_text(size = 7))
png(here::here("Figures", "Figure3", "Fig3E_PC3_loadings.png"), width = 1.5, height = 2, units = "in", res = 1200)
print(p)
dev.off()

