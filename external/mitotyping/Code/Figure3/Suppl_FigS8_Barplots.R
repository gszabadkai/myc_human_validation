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

mitoPPS <- read_csv(here::here("Data", "HumanProteinAtlas","ProcessedData", "mitoPPS_all_tissues.csv" ))




subset <- mitoPPS %>%
  pivot_longer(cols = -c(Tissue, Group), names_to = "Pathway") %>%
  filter(Pathway %in% c(
    "Apoptosis", 
    "Cholesterol, bile acid, steroid synthesis",
    "Complex I",
    "Complex II",
    "Fatty acid oxidation",
    "Calcium homeostasis",
    "mtDNA repair",
    "Protein import and sorting",
    "Glycerol phosphate shuttle",
    "Malate-aspartate shuttle"
  )) %>%
  mutate(value = log10(value)) 

pathway_selection <- unique(subset$Pathway)
for (i in 1:length(pathway_selection)) {
  pathway = pathway_selection[i]
  
  subset %>%
    filter(Pathway %in% pathway) %>%
    ggplot(aes(x = desc(reorder(Tissue,value)), y = value, fill = Group, color = Group, label = Tissue)) +
    geom_segment(aes(x=dplyr::desc(reorder(Tissue, value)),xend=dplyr::desc(reorder(Tissue, value)),y=0, yend=value,
                     color = Group), size=0.8, alpha =0.7) +
    scale_fill_manual(values = color_groups) +
    scale_color_manual(values = color_groups) +
    labs(title= pathway) +
    theme_bw() +
    theme(legend.position = "none",  
          plot.title = element_text(size =8),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 6),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) 
  ggsave(here::here("Figures","Figure3", paste0("Suppl_FigS8_", pathway, "_PPS_barplot.png")), height = 2.8, width =2, dpi = 1200)
}

subset_wide <- subset %>%
  pivot_wider(names_from = "Pathway", values_from = "value")

