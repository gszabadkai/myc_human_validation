rm(list = ls())
library(tidyverse)
library(ggfortify)
library(Rtsne)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(broom)


# Figure 4AB Proliferation, Senescence ----------------------------------------



data <- readRDS(here::here("Data", "Fibroblasts", "ProcessedData", "nTPM_Lifespan.rds")) %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "RNAseq_sampleID", values_to = "nTPM") %>%
  filter(Gene %in% c("CDKN2A",
                     "CDKN1A",
                     "TP53",
                     "MKI67", 
                     "TOP2A", 
                     "RRM2"
  ))
meta <- read.csv(here::here("Data","Fibroblasts","OriginalData", "RNAseq_meta.csv")) %>%
  mutate(RNAseq_sampleID = as.character(RNAseq_sampleID)) %>%
  mutate(Sample = "Sample") %>%
  unite(RNAseq_sampleID, Sample, RNAseq_sampleID) %>%
  filter(Study_part %in% 2) %>%
  filter(Clinical_condition %in% "Normal") %>%
  filter(Percent_oxygen %in% 21) %>%
  filter(Treatments %in% "Control")

data_senesc_prolif <- inner_join(data, meta) %>%
  mutate(Pathway = case_when((Gene == "CDKN2A" | Gene == "CDKN1A" | Gene == "TP53") ~"Senescence",
                             TRUE ~ "Proliferation"))  %>%
  group_by(RNAseq_sampleID, Pathway) %>%
  mutate(score = mean(nTPM)) %>%
  dplyr::select(RNAseq_sampleID, Cell_line_inhouse, Passage, Days_grown_Udays, Pathway, score) %>% 
  unique() %>%
  pivot_wider(names_from = "Pathway", values_from = "score")



## Plot proliferation passage ----------------------------------------------


data_senesc_prolif %>%
  select(Passage, Proliferation) %>%
  unique() %>%
  ggplot(aes(x = Passage, y = Proliferation, color = Passage)) + 
  geom_point(alpha = 0.7) +
  theme_bw() +
  scale_color_viridis_c(direction = -1) +
  theme(legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))
ggsave(here::here("Figures", "Figure4",  "Fig4A_Proliferation_passage.png"), width = 1.7, height = 1.1, units = "in", dpi = 1200)

## Plot senescence passage ----------------------------------------------

data_senesc_prolif %>%
  select(Passage, Senescence) %>%
  unique() %>%
  ggplot(aes(x = Passage, y = Senescence, color = Passage)) + 
  geom_point(alpha = 0.7) +
  theme_bw() +
  scale_color_viridis_c(direction = -1) +
  theme(legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))
ggsave(here::here("Figures", "Figure4",  "Fig4B_Senescence_passage.png"), width = 1.7, height = 1.1, units = "in", dpi = 1200)


rm(list = ls())

theme_pca <- theme(axis.text = element_text(size = 6),
                   axis.title = element_text(size = 7),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.position = "none"
)

mitoPPS <- read_csv(here::here("Data", "Fibroblasts", "ProcessedData", "mitoPPS_Controls_SP2.csv"))


# Figure 4C PCAs ---------------------------------------------------------------

## hFB12 -------------------------------------------------------------------

data_pca <- mitoPPS %>%
  filter(Line %in% "hFB12") %>%
  filter(Study_part %in% 2)
pca <- prcomp(data_pca[,-c(1:8)] , scale. = T)
loadings <- pca$rotation
write_csv(as.data.frame(loadings) %>% rownames_to_column("Pathway"), 
          here::here("Supplemental/Figure4C_loadings_HC1.csv"))

autoplot(pca, data = data_pca, colour = 'Passage', alpha = 0.6,
         fill = 'Passage', size = 2.3, shape = 21, lwd = 0.5)+
  viridis::scale_color_viridis(direction = -1) +
  viridis::scale_fill_viridis(direction = -1) +
  theme_minimal() +
  theme_pca

ggsave(here::here("Figures", "Figure4", "Fig4C_PCA_mitoPPS_hFB12.png"), height = 1.35, width = 1.5,unit = "in", dpi = 1200)


## hFB13 -------------------------------------------------------------------

data_pca <- mitoPPS %>%
  filter(Line %in% "hFB13") %>%
  filter(Study_part %in% 2)
pca <- prcomp(data_pca[,-c(1:8)] , scale. = T)
loadings <- pca$rotation
write_csv(as.data.frame(loadings) %>% rownames_to_column("Pathway"), 
          here::here("Supplemental/Figure4C_loadings_HC2.csv"))
autoplot(pca, data = data_pca, colour = 'Passage', alpha = 0.6,
         fill = 'Passage', size = 2.3, shape = 21, lwd = 0.5)+
  viridis::scale_color_viridis(direction = -1) +
  viridis::scale_fill_viridis(direction = -1) +
  theme_minimal() +
  theme_pca

ggsave(here::here("Figures", "Figure4", "Fig4C_PCA_mitoPPS_hFB13.png"), height = 1.35, width = 1.5,unit = "in", dpi = 1200)

## hFB14 -------------------------------------------------------------------

data_pca <- mitoPPS %>%
  filter(Line %in% "hFB14") %>%
  filter(Study_part %in% 2)
pca <- prcomp(data_pca[,-c(1:8)] , scale. = T)
loadings <- pca$rotation
write_csv(as.data.frame(loadings) %>% rownames_to_column("Pathway"), 
          here::here("Supplemental/Figure4C_loadings_HC3.csv"))
autoplot(pca, data = data_pca, colour = 'Passage', alpha = 0.6,
         fill = 'Passage', size = 2.3, shape = 21, lwd = 0.5)+
  viridis::scale_color_viridis(direction = -1) +
  viridis::scale_fill_viridis(direction = -1) +
  theme_minimal() +
  theme_pca

ggsave(here::here("Figures", "Figure4", "Fig4C_PCA_mitoPPS_hFB14.png"), height = 1.35, width = 1.5,unit = "in", dpi = 1200)


# Figure 4D Heatmap -------------------------------------------------------
all_pathways <- colnames(mitoPPS[,9:ncol(mitoPPS)])
cor_res <- mitoPPS %>%
  pivot_longer(cols = all_of(all_pathways), names_to = "Pathway") %>%
  group_by(Pathway) %>%
  nest(-Pathway) %>%
  mutate(cor=map(data,~cor.test(.x$value, .x$Passage, method = "sp"))) %>%
  mutate(tidied = map(cor, tidy)) %>%
  unnest(tidied, .drop = T) %>%
  ungroup() %>%
  mutate(padj = p.adjust(p.value, method = "BH"))

## Top and bottom 10 pathways
top_10 <- cor_res %>%
  dplyr::filter( (padj < 0.05) | (padj = 0.05) ) %>%
  ungroup() %>%
  arrange(desc(estimate)) %>%
  dplyr::slice(1:10) %>%
  pull(Pathway)

bottom_10 <- cor_res %>%
  dplyr::filter( (padj < 0.05) | (padj = 0.05) ) %>%
  ungroup() %>%
  arrange(estimate) %>%
  dplyr::slice(1:10) %>%
  pull(Pathway)

test <- mitoPPS  %>%
  ungroup() %>%
  dplyr::select(Passage, Days_grown_Udays, Line, all_of(c(top_10, bottom_10))) %>%
  arrange(Line, Passage) %>%
  dplyr::rename(`Mt dynamics and surveillance` = `Mitochondrial dynamics and surveillance`,
         `NAD biosynthesis and metab.` = `NAD biosynthesis and metabolism`)

# Meta and expression data
meta  <- test[,1:3] #%>% column_to_rownames("Unique_variable_name")

exprs <- test[4:ncol(test)]
exprs = log10(exprs)
exprs = t(exprs)
rownames(exprs)
colnames(exprs)

# Col annotation
col_anno = columnAnnotation(
  `Cell_line`=meta$Line,
  Passage = anno_barplot(meta$Passage, axis_param=list(gp=gpar(fontsize = 4))),
  col=list(`Cell_line`  = c("hFB12" = "#636363", "hFB13" = "#969696", "hFB14" = "#D9D9D9")),
  show_annotation_name = F,
  show_legend =  T
)



col_anno <- re_size(col_anno, annotation_height = unit(c(0.1,0.3), "in"))
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue","white", "#ff4242")) #"#ffa76b",
row_dist    = dist(as.matrix(exprs), method="euclidean" )
rowdend     = hclust(row_dist, method="ward.D2")
HM2 <- Heatmap(exprs,
               name = "Pathway weight",
               cluster_rows = rowdend,
               show_column_names =F,
               cluster_columns =F,
               row_title="",
               row_dend_side = "right",
               row_names_side = "left",
               show_row_dend=F,
               show_row_names=T,
               top_annotation=col_anno,
               row_names_gp = grid::gpar(fontsize = 5),
               column_names_gp = grid::gpar(fontsize = 6),
               heatmap_legend_param = list(
                 title = "mitoPPS (log10)",
                 title_position = "lefttop-rot",
                 legend_height = unit(0.9, "in"),
                 labels_gp = gpar(fontsize = 8),
                 title_gp = gpar(fontsize = 7, fontface= "bold")
               )
)


png(here::here("Figures","Figure4", "Fig4D_mitoPPS_Age_CorrHeatmap.png"), width = 3.8, height = 2.7, units = "in", res = 600)
draw(HM2)
dev.off()


# Figure 4E Scatter plots ---------------------------------------------------------------


mitoPPS %>%
  dplyr::select(Line, Passage, `Fission`, `Mitophagy`, `mtDNA repair`, `mtDNA maintenance`) %>%
  pivot_longer(cols = -c(Line, Passage), names_to = "Pathway", values_to = "mitoPPS") %>%
  ggplot(aes(x = Passage, y = mitoPPS, color = Pathway)) +
  geom_point(size =0.8, alpha = 0.7) +
  scale_color_manual(values = c("Fission" ="#CE9332", "Mitophagy" = "#55BC82",
                                "mtDNA repair" ="#56BCC2", "mtDNA maintenance" = "#4FADF0")) +
  theme_bw() +
  facet_wrap(~Pathway, scales = "free", ncol = 4) +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.position =  "none",
        strip.text = element_text(size = 6))
ggsave(here::here("Figures", "Figure4", "Fig4E_Scatterplots.png"), width = 4.5, height = 1.2, units = "in", dpi = 1200)

cor_res %>%
  filter(Pathway %in% c("Fission", "Mitophagy", "mtDNA repair", "mtDNA maintenance"))


# Figure 4F PCA Tissues + Fibroblasts -------------------------------------

tissue_mitoPPS <- read_csv(here::here("Data", "HumanProteinAtlas", "ProcessedData", "mitoPPS_all_tissues.csv"))
fibroblast_mitoPPS <- read_csv(here::here("Data", "Fibroblasts", "ProcessedData", "mitoPPS_Controls_SP2.csv"))
meta <- read.csv(here::here("Data","Fibroblasts","OriginalData", "RNAseq_meta.csv")) %>%
  mutate(RNAseq_sampleID = as.character(RNAseq_sampleID)) %>%
  mutate(Sample = "Sample") %>%
  unite(RNAseq_sampleID, Sample, RNAseq_sampleID)
lifespan_tissues <-  fibroblast_mitoPPS %>%
  mutate(Group = "Fibroblasts") %>%
  mutate(Tissue = "Control") %>%
  inner_join(meta) %>%
  select(RNAseq_sampleID, Cell_line_inhouse, Study_part, Passage, Days_grown_Udays, Tissue, Group, all_of(all_pathways)) %>%
  pivot_longer(cols = -c(RNAseq_sampleID, Cell_line_inhouse, Study_part, Passage, Days_grown_Udays, Tissue, Group), names_to = "Pathway") %>%
  full_join(tissue_mitoPPS %>%pivot_longer(cols = -c(Tissue, Group) , names_to = "Pathway")%>%
              mutate(Passage =NA)) %>%
  pivot_wider(names_from = "Pathway", values_from = "value") #%>%

lifespan_tissues <-  fibroblast_mitoPPS %>%
  mutate(Group = "Fibroblasts") %>%
  mutate(Tissue = "Control") %>%
  inner_join(meta) %>%
  select(RNAseq_sampleID, Cell_line_inhouse, Study_part, Passage, Days_grown_Udays, Tissue, Group, all_of(all_pathways)) %>%
  pivot_longer(cols = -c(RNAseq_sampleID, Cell_line_inhouse, Study_part, Passage, Days_grown_Udays, Tissue, Group), names_to = "Pathway") %>%
  full_join(tissue_mitoPPS %>%pivot_longer(cols = -c(Tissue, Group) , names_to = "Pathway")%>%
              mutate(Passage =NA)) %>%
  pivot_wider(names_from = "Pathway", values_from = "value") %>%
  arrange(desc(Group))



## Tissue + Fibroblast PCA -------------------------------------------------


pca <- prcomp(lifespan_tissues[, -c(1:7)], scale = T)
summary(pca)
loadings <- pca$rotation
autoplot(pca, data = lifespan_tissues , colour = 'Passage', alpha = 0.6, scale = 1,
         fill = 'Passage', size = 2, shape = 21, lwd = 0.5, )+
  viridis::scale_color_viridis(direction = -1, na.value = "grey") +
  viridis::scale_fill_viridis(direction = -1, na.value = "grey") +
  theme_minimal() +
  theme_pca + 
  scale_y_reverse() ## Note the switch in sign in PCA, axis reversed for visualization purposes. 
  ## Loadings comparable to Figure 5B, what was positive before is now negative
#ggsave(here::here("Figures", "Figure6", "Fig6E_PCA_Tissue_Fibroblasts.png"), , height = 1.35, width = 1.5,unit = "in", dpi = 1200)
ggsave(here::here("Figures", "Figure4", "Fig4F_PCA_Tissue_Fibroblasts.png"), height = 1.56, width = 1.68,unit = "in", dpi = 1200)



## Fibroblast subset -------------------------------------------------------------


pca$x <- pca$x[which(lifespan_tissues$Group == "Fibroblasts"), ]
autoplot(pca, data = lifespan_tissues %>%
           filter(Group %in% "Fibroblasts"), colour = 'Passage', alpha = 0.6, scale = 1,
         fill = 'Passage', size = 2.3, shape = 21, lwd = 0.5)+
  viridis::scale_color_viridis(direction = -1) +
  viridis::scale_fill_viridis(direction = -1) +
  theme_minimal() +
  theme_pca + 
  scale_y_reverse()
ggsave(here::here("Figures", "Figure4", "Fig4F_PCA_Tissue_Fibroblasts_subset.png"), 
       height = 1.56, width = 1.68,unit = "in", dpi = 1200)




## get R session info
sessionInfo()
