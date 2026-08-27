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

MitoPathwayScores <- read_csv(here::here("Data", "HumanProteinAtlas","ProcessedData", "Pathway_scores_all_tissues.csv" ))
mitoPPS <- read_csv(here::here("Data", "HumanProteinAtlas","ProcessedData", "mitoPPS_all_tissues.csv" )) 

range(mitoPPS[,-c(1:2)])
range(MitoPathwayScores$score)

# Figure 3C Heatmap mitoPPS -----------------------------------------------------------

## Meta and expression data
meta  <- mitoPPS %>% select(Tissue, Group)
exprs <- mitoPPS %>%
  column_to_rownames("Tissue") %>%
  select(-Group)
exprs = log10(exprs)


## kmeans rows
set.seed(123)
k.max <- 15
wss <- sapply(1:k.max,
              function(k){kmeans(exprs, k, nstart=50,iter.max = 15 )$tot.withinss})
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE,
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")


## Row annotation
col_anno = columnAnnotation(
  `Tissue Group`=meta$Group,
  col=list(`Tissue Group`  = color_groups),
  show_annotation_name = F,
  show_legend =  T,
  annotation_legend_param = list(nrow=8)
)
range(exprs)

## Transpose 
exprs <- t(exprs)
## Row and column dendrogram
distance    = dist(t(as.matrix(exprs)), method="euclidean" )
coldend     = hclust(distance, method="ward.D2")
row_dist    = dist(as.matrix(exprs), method="euclidean" )
rowdend     = hclust(row_dist, method="ward.D2")

## Heatmap
set.seed(456)
HM2 <- Heatmap(exprs, 
               name = "Pathway weight", 
               row_km =5,
               cluster_columns = coldend,
               show_column_names =T,
               show_column_dend = T,
               row_title="",
               row_dend_side = "right",
               row_dend_width = unit(0.25, "in"),
               column_dend_height =  unit(0.25, "in"),
               row_names_side = "left",
               show_row_dend=T,
               show_row_names=T,
               top_annotation=col_anno,
               row_names_gp = grid::gpar(fontsize = 5),
               column_names_gp = grid::gpar(fontsize = 6),
               heatmap_legend_param = list(
                 title = "mitoPPS (log10)",
                 title_position = "lefttop-rot",
                 legend_height = unit(1, "in")
               )
)

png(here::here("Figures", "Figure3", "Fig3C_mitoPPS_Heatmap.png"),
    #height =12, width = 8, units = 'in', res = 1200)
    height =6.5, width = 8, units = 'in', res = 1200)
set.seed(123)
ht <- draw(HM2)
ht
dev.off()
row_order <- ht@ht_list[["Pathway weight"]]@row_order


# Suppl Figure S7A Heatmap mitoPPS detailed---------------------------------------------------------

exprs_ordered <- exprs[row_order,]
HM2 <- Heatmap(exprs_ordered, 
               name = "Pathway weight", 
               cluster_columns = coldend,
               cluster_rows = F,
               show_column_names =T,
               show_column_dend = T,
               row_title="",
               column_dend_height =  unit(0.25, "in"),
               row_names_side = "left",
               show_row_dend=F,
               show_row_names=T,
               top_annotation=col_anno,
               row_names_gp = grid::gpar(fontsize = 7),
               column_names_gp = grid::gpar(fontsize = 7),
               heatmap_legend_param = list(
                 title = "mitoPPS (log10)",
                 title_position = "lefttop-rot",
                 legend_height = unit(1, "in")
               )
)



png(here::here("Figures", "Figure3", "Suppl_FigS7A_mitoPPS_Heatmap_detailed.png"),
    height =14, width = 9, units = 'in', res = 1200)
set.seed(123)
ht <- draw(HM2)
ht
dev.off()


# Suppl Figure S7B Heatmap MitoPathwayScore  ------------------------------------------------

## Meta and expression data
meta  <- MitoPathwayScores %>% 
  pivot_wider(names_from = "Pathway", values_from = "score") %>%
  select(Tissue, Group)
exprs <- MitoPathwayScores %>% 
  group_by(Pathway) %>%
  mutate(zscore = (score - mean(score)) / sd(score)) %>%
  select(-score) %>%
  arrange(Pathway) %>%
  pivot_wider(names_from = "Pathway", values_from = "zscore") %>%
  column_to_rownames("Tissue") %>%
  select(-Group) 
exprs <- t(exprs)
range(exprs)
exprs_ordered <- exprs[row_order,]

## Row annotation
col_anno = columnAnnotation(
  `Tissue Group`=meta$Group,
  col=list(`Tissue Group`  = color_groups),
  show_annotation_name = F,
  show_legend =  T,
  annotation_legend_param = list(nrow=8)
)

## Heatmap
HM2 <- Heatmap(exprs_ordered, 
               name = "Pathway weight", 
               cluster_columns = coldend,
               cluster_rows = F,
               show_column_names =T,
               show_column_dend =T,
               row_title="",
               column_dend_height =  unit(0.25, "in"),
               row_names_side = "left",
               show_row_dend=F,
               show_row_names=T,
               top_annotation=col_anno,
               row_names_gp = grid::gpar(fontsize = 7),
               column_names_gp = grid::gpar(fontsize = 7),
               heatmap_legend_param = list(
                 title = "Pathway zscore",
                 title_position = "lefttop-rot",
                 legend_height = unit(1, "in")
               )
)




ht <- draw(HM2)
png(here::here("Figures", "Figure3", "Suppl_FigS7B_PathwayZscore_Heatmap_detailed.png"),
    height =14, width = 9, units = 'in', res = 1200)
draw(ht)
dev.off()


# Suppl Figure S7C Average mito expression ---------------------------------
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
  group_by(Tissue) %>%
  mutate(Avg_mito_exprs = mean(nTPM)) %>%
  select(Tissue, Avg_mito_exprs) %>%
  unique() %>%
  tissue_to_group_hm()%>%
  mutate(Tissue = str_to_sentence(Tissue))

processed_data %>%
  ggplot(aes(color = Group)) +
  geom_segment(aes(x=dplyr::desc(reorder(Tissue, Avg_mito_exprs)),
                   xend=dplyr::desc(reorder(Tissue,  Avg_mito_exprs)),
                   y=0, yend= Avg_mito_exprs), 
               size=1.3, alpha =0.7) +
  labs(title= "Average mitochondrial gene expression") +
  labs(y= "Average mito gene expression (nTPM)") +
  scale_color_manual(values = color_groups) + 
  geom_hline(yintercept = 0) +

  theme_classic() +
  theme(legend.position = "none",  
        plot.title = element_text(size =7),
        axis.title.y  = element_text(size =6),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        #axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave(here::here("Figures", "Figure3", "Suppl_FigS7C_Avg_mito_exprs_allTissues.png"), width = 4, height = 2, units = "in", dpi = 1200)


# Suppl Figure S7D Score vs mitoPPS bivariate --------------------------------

test <- mitoPPS %>%
  pivot_longer(cols = -c("Tissue", "Group"), names_to = "Pathway", values_to = "mitoPPS") %>%
  inner_join(MitoPathwayScores) 

p <- test %>%
  ggplot(aes(x = log10(mitoPPS), y = log10(score), fill =Group, 
             color = Group, label = Pathway)) +
  geom_point(size = 1.2, alpha = 0.6, shape =21) +
  scale_fill_manual(values = color_groups) +
  scale_color_manual(values = color_groups) +
  theme_classic() + 
  theme(legend.position = "none",  
        plot.title = element_text(size =7),
        axis.title  = element_text(size =6),
        axis.text = element_text(size = 6),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p
ggsave(here::here("Figures", "Figure3", "Suppl_FigS7D_Score_vs_mitoPPS_bivariate.png"), width = 2.3, height = 2.3, units = "in", dpi = 1200)


##  Correlation overall mitoPPS - MitoPathwayScore ------------------------------------
cor.test(test$mitoPPS, test$score, method = "spearman")

##  Correlation detail mitoPPS - MitoPathwayScore ------------------------------------
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
mtgene_pathways <- gene_to_pathway %>%
  filter(grepl("MT-", Gene)) %>%
  select(Pathway) %>%
  unique() %>%
  pull(Pathway)


### mtDNA pathways -------------------------------------------------------------


mtgene_data <- test %>%
  filter(Pathway %in% mtgene_pathways) %>%
  filter(!is.na(Pathway)) %>%
  select(Tissue, Group, Pathway, mitoPPS, score) %>%
  unique() 
cor.test(mtgene_data$mitoPPS, mtgene_data$score, method = "spearman")

### Non-mtDNA pathways -------------------------------------------------------------

non_mtgene_data <- test %>%
  filter(!Pathway %in% mtgene_pathways) %>%
  filter(!is.na(Pathway)) %>%
  select(Tissue, Group, Pathway, mitoPPS, score) %>%
  unique() 
cor.test(non_mtgene_data$mitoPPS, non_mtgene_data$score, method = "spearman")


# In text: Fold difference for some pathways ---------------------------------------

## OxPhos brain vs immune

test <- mitoPPS %>%
  filter(Group %in% c("CNS", "Immune")) %>%
  select(Tissue, Group, OXPHOS) %>%
  group_by(Group) %>%
  mutate(mean = mean(OXPHOS))
test_split <- split(test$OXPHOS, test$Group)
lapply(test_split, shapiro.test) # p close to 0.05 --> assume not normal distribution
fligner.test(OXPHOS ~ Group, data = test) # equal variance
wilcox.test(test$OXPHOS~ test$Group)
mean(test_split$CNS) / mean(test_split$Immune)

## Heme synthesis and processing
test <- mitoPPS %>%
  mutate(Group2 = case_when(
    (Tissue == "Liver" | Tissue == "Colon" |
      Tissue == "Rectum" | Tissue == "Duodenum" |
     Tissue == "Small intestine" 
    ) ~ "A",
    TRUE ~ "B"
  )) %>%
  select(Tissue, Group2, `Heme synthesis and processing`)  %>%
  group_by(Group2) %>%
  mutate(mean = mean(`Heme synthesis and processing`))

test_split <- split(test$`Heme synthesis and processing`, test$Group2)
lapply(test_split, shapiro.test) # not normally distributed
fligner.test(`Heme synthesis and processing` ~ Group2, data = test) # variance not equal
brunnermunzel::brunnermunzel.test(test_split[[1]], test_split[[2]], perm = T)
mean(test_split$A) / mean(test_split$B)

## Creatine metabolism
test <- mitoPPS %>%
  mutate(Group2 = case_when(
    (Group == "Contractile" | Group == "Anabolic" |
       Tissue == "Tongue" | Tissue == "Pancreas") ~ "A",
    TRUE ~ "B"
  )) %>%
  select(Tissue, Group2, `Creatine metabolism`) %>%
  group_by(Group2) %>%
  mutate(mean = mean(`Creatine metabolism`))
test_split <- split(test$`Creatine metabolism`, test$Group2)
lapply(test_split, shapiro.test) # not normally distributed
fligner.test(`Creatine metabolism` ~ Group2, data = test) # variance not equal
brunnermunzel::brunnermunzel.test(test_split[[1]], test_split[[2]], perm = T)
mean(test_split$A) / mean(test_split$B)

## Folate 1C metabolism
test <- mitoPPS %>%
  mutate(Group2 = case_when(
    (Group == "CNS" | 
       Tissue == "Liver" ) ~ "A",
    TRUE ~ "B"
  )) %>%
  select(Tissue, Group2, `Folate and 1-C metabolism`) %>%
  group_by(Group2) %>%
  mutate(mean = mean(`Folate and 1-C metabolism`))
test_split <- split(test$`Folate and 1-C metabolism`, test$Group2)
lapply(test_split, shapiro.test) # normal distribution
test %>%
  ggplot(aes(x = Group2, y= `Folate and 1-C metabolism`)) +
  geom_boxplot()
fligner.test(`Folate and 1-C metabolism` ~ Group2, data = test) # equal variance
t.test(test$`Folate and 1-C metabolism`~ test$Group2, var.equal = T)
mean(test_split$A) / mean(test_split$B)

## Cholesterol in adrenal
test <- mitoPPS %>%
  mutate(Group2 = case_when(
    (Tissue == "Adrenal gland" ) ~ "A",
    TRUE ~ "B"
  )) %>%
  select(Tissue, Group2, `Cholesterol, bile acid, steroid synthesis`)  %>%
  group_by(Group2) %>%
  mutate(mean = mean(`Cholesterol, bile acid, steroid synthesis`))
test_split <- split(test$`Cholesterol, bile acid, steroid synthesis`, test$Group2)
mean(test_split$A) / mean(test_split$B)

## VitaminD in adrenal
test <- mitoPPS %>%
  mutate(Group2 = case_when(
    (Tissue == "Adrenal gland" ) ~ "A",
    TRUE ~ "B"
  )) %>%
  select(Tissue, Group2, `Vitamin D metabolism`)  %>%
  group_by(Group2) %>%
  mutate(mean = mean(`Vitamin D metabolism`))
test_split <- split(test$`Vitamin D metabolism`, test$Group2)
mean(test_split$A) / mean(test_split$B)

## FeS in adrenal
test <- mitoPPS %>%
  mutate(Group2 = case_when(
    (Tissue == "Adrenal gland" ) ~ "A",
    TRUE ~ "B"
  )) %>%
  select(Tissue, Group2, `Fe-S cluster biosynthesis`)  %>%
  group_by(Group2) %>%
  mutate(mean = mean(`Fe-S cluster biosynthesis`))
test_split <- split(test$`Fe-S cluster biosynthesis`, test$Group2)
mean(test_split$A) / mean(test_split$B)

## CII, G3P, FAO in liver
test <- mitoPPS %>%
  filter(Tissue %in% "Liver") %>%
  pivot_longer(cols = -c(Tissue, Group, `Glycerol phosphate shuttle`), names_to = "Pathway", values_to = "mitoPPS") %>%
  filter(Pathway %in% c("Complex II","Fatty acid oxidation")) %>%
  mutate(FC = mitoPPS / `Glycerol phosphate shuttle`)

