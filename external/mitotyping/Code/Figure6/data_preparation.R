library(SingleCellExperiment)
library(slingshot)
library(tidyverse)
library(ComplexHeatmap)
library(Seurat)
# library(SeuratDisk)
# library(SeuratData)


# Data preparation ---------------------------------------------------------------

######______ Read cluster data for annotation (from HPA cluster)
clusters <- read.delim(here::here("Data","Guo_et_al", "OriginalData", "rna_single_cell_cluster.tsv")) %>%
  filter(Tissue %in% "testis")
 
genes <- clusters %>%
  select(Gene, Gene.name) %>%
  unique() %>%
  dplyr::rename(Ensembl = Gene, Gene = Gene.name)

######______ Read raw per cell count data  
counts <- read.delim(here::here("Data","Guo_et_al", "OriginalData","read_count.tsv"), row.names = 1, header = T,
                     check.names = F) %>%
  as.data.frame() %>%
  rownames_to_column("Ensembl") %>%
  inner_join(genes) %>%
  select(-Ensembl) %>%
  group_by(Gene) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  column_to_rownames("Gene") %>%
  as.matrix()
#
######______ Meta data 
cluster_to_celltype <- clusters %>%
  mutate(cluster = as.numeric(gsub("c-", "", Cluster))) %>%
  select(cluster, Cell.type) %>%
  dplyr::rename(CellType = Cell.type) %>%
  unique()

meta <- read.delim("/Users/annamonzel/Desktop/R/R_Testing/HPA_CellType/Data/Original/rna_single_cell_read_count/testis/cell_data.tsv") %>%
  column_to_rownames("cell_id") %>%
  full_join(cluster_to_celltype, by = "cluster") %>%
  unite(CellType_cluster, CellType, cluster, remove = F) 



seurat <- Seurat::CreateSeuratObject(counts = counts, assay = "RNA", project = "testis", meta.data = meta)


# Find MitoGenes ----------------------------------------------------------


mitocarta_sheet4 <- readxl::read_xls(here::here("Data/MitoCarta/OriginalData/HumanMitoCarta3_0.xls"), sheet = 4) %>%
  dplyr::select(MitoPathway, Genes) %>%
  na.omit() 
gene_to_pathway <- splitstackshape::cSplit(mitocarta_sheet4, 'Genes', ',') %>%
  column_to_rownames("MitoPathway") %>%
  t() %>%
  as.data.frame() 
gene_to_pathway <- gene_to_pathway %>%
  pivot_longer(cols = colnames(gene_to_pathway), names_to = "Pathway", values_to = "Gene") %>%
  na.omit() 
mitogenes <- unique(gene_to_pathway$Gene)
mitogenes_found <- mitogenes[mitogenes %in% rownames(seurat)]
mitogenes_missing <- unique(gene_to_pathway$Gene)[!unique(gene_to_pathway$Gene) %in% rownames(seurat)]
rm(mitocarta_sheet4)



## Find via MC synonyms
synonyms <- readxl::read_xls(here::here("Data/MitoCarta/OriginalData/HumanMitoCarta3_0.xls"), sheet = 2) %>%
  dplyr::select(Symbol, Synonyms) %>%
  mutate(Synonyms = strsplit(Synonyms, "\\|")) %>%  
  unnest(Synonyms) %>%                            
  mutate(Synonyms = trimws(Synonyms))  %>%
  full_join(gene_to_pathway, by = c("Symbol" = "Gene"))

gene_synonym_match <- synonyms %>%
  filter(Symbol %in% mitogenes_missing) 
synonyms_found <- gene_synonym_match$Synonyms[gene_synonym_match$Synonyms %in% rownames(seurat)] 
found_via_MC <- synonyms %>%
  filter(Synonyms %in% synonyms_found) %>%
  select(Symbol, Synonyms) %>%
  unique() %>%
  pull(Symbol)
mitogenes_still_missing <- mitogenes_missing[!mitogenes_missing %in% found_via_MC]


## Find via HGNC helper
hgnc_check <- HGNChelper::checkGeneSymbols(mitogenes_still_missing,
                                           species = "human",
                                           unmapped.as.na = FALSE)
hgnc_synonyms_found <- hgnc_check$Suggested.Symbol[hgnc_check$Suggested.Symbol %in% rownames(seurat)] 
found_via_HGNC <- hgnc_check %>%
  filter(Suggested.Symbol %in% hgnc_synonyms_found) %>%
  select(x, Suggested.Symbol) %>%
  unique() %>%
  pull(x)
mitogenes_still_missing <- mitogenes_missing[!mitogenes_missing %in% c(found_via_MC, found_via_HGNC)]


gene_to_pathway <- gene_to_pathway %>%
  mutate(Gene = case_when(
    Gene == "PHB" ~ "PHB1",
     Gene == "C12orf65" ~ "MTRFR",
     Gene == "ZADH2" ~ "PTGR3",
     Gene == "ATP5MD" ~ "ATP5MK", 
     Gene == "ATP5MPL" ~ "ATP5MJ",
    TRUE ~ Gene
  ))
mitogenes <- gene_to_pathway$Gene %>% unique()
rm(list = setdiff(ls(), c("counts", "seurat", "gene_to_pathway", "mitogenes", "meta")))


# Filter by expression ----------------------------------------------------
## Filter genes out that are detected in less than 1% of cells
library(Matrix)
cts <- seurat@assays[["RNA"]]@layers[["counts"]]
n_cells <- ncol(cts) #total number of cells in the dataset
min_prop <- 0.001   #minimum fraction of cells that must express a gene (1%)
min_cells <- max(10, ceiling(0.001 * n_cells))  # ~0.3% or at least 10 cells
min_total <- 10 #minimum total UMI count per gene (sum across all cells)
detected_in <- Matrix::rowSums(cts > 0) #number of cells where each gene is expressed (nonzero counts)
total_umis  <- Matrix::rowSums(cts) #total UMI count per gene across all cells
keep <- detected_in >= min_cells & total_umis >= min_total
table(keep)   

seurat_filtered <- seurat[keep, ]    # subset genes

mitogenes_before_filter <- mitogenes[mitogenes %in% rownames(seurat)]
mitogenes_after_filter <- mitogenes[mitogenes %in% rownames(seurat_filtered)]

mitogenes_lost <- mitogenes_before_filter[!mitogenes_before_filter %in% mitogenes_after_filter]
test <- gene_to_pathway %>%
  group_by(Pathway) %>%
  mutate(n_genes_pathway = n()) %>%
  filter(Gene %in% mitogenes_lost) %>%
  group_by(Pathway) %>%
  mutate(n_genes_lost = n())

seurat <- seurat_filtered
rm(list = setdiff(ls(), c("counts", "gene_to_pathway",  "seurat", "meta")))


theme_umap <-   theme(legend.text = element_text(size = 14),
                      axis.text = element_text(size = 12),
                      axis.title = element_text(size = 14),
                      axis.line = element_line(linewidth = 0.5),
                      legend.position = "none")

umap <- meta[,c("umap_x", "umap_y")]
colnames(umap) <- c("UMAP_1", "UMAP_2")
rownames(umap) <- rownames(meta)
umap_dr <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_",assay = DefaultAssay(seurat))
seurat[["umap"]] <- umap_dr
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
Idents(seurat) <- "CellType"
DimPlot(seurat, label = F, raster.dpi = c(1024, 1024))+
  theme_umap

ggsave(here::here("Figures", "Figure6", "Suppl_FigS13A_seurat_all_cells.png"), width =4.6,height =4.6,units = "in")

######______ Assign CellSubtypes based on Germ cell markers of Guo et al 

seurat <- subset(seurat, subset = CellType %in% c("spermatogonia", "spermatocytes", "early spermatids", "late spermatids"))
Idents(seurat) <- "cluster"
DimPlot(seurat, label = T, raster.dpi = c(1024, 1024))+
  theme_umap
ggsave(here::here("Figures", "Figure6", "Suppl_FigS13B_seurat_subset.png"), width =4.6,height =4.6,units = "in")
FeaturePlot(seurat, 
            features = c("UTF1",  "FGFR3","ID4",
                         "KIT", "DMRT1", "DMRTB1", "SYCP3",
                         "SPO11", "MLH3", "ZPBP", 
                         "TNP1", "PRM2"))
ggsave(here::here("Figures", "Figure6", "Suppl_FigS13D_FeaturePlot_GermMarkers.png"), width =14, height = 10)



seurat$CellSubtype <- NA
seurat$CellSubtype[seurat$cluster == 14] <- "SSC"
seurat$CellSubtype[seurat$cluster == 4] <- "Diff_spermatogonia"
seurat$CellSubtype[seurat$cluster == 1] <- "Early_prim_spermatocytes"
seurat$CellSubtype[seurat$cluster == 17] <- "Late_prim_spermatocytes"
seurat$CellSubtype[seurat$cluster == 16] <- "Late_prim_spermatocytes"
seurat$CellSubtype[seurat$cluster == 8] <- "Round_spermatids"
seurat$CellSubtype[seurat$cluster == 3] <- "Elongated_spermatids"
seurat$CellSubtype[seurat$cluster == 15] <- "Late_spermatids"
seurat$CellSubtype[seurat$cluster == 9] <- "Late_spermatids"
seurat$CellSubtype[seurat$cluster == 18] <- "Late_spermatids"
seurat$CellSubtype[seurat$cluster == 5] <- "Late_spermatids"
seurat$CellSubtype[seurat$cluster == 13] <- "Late_spermatids"
seurat$CellSubtype[seurat$cluster == 11] <- "Late_spermatids"
seurat$CellSubtype[seurat$cluster == 2] <- "Late_spermatids"
seurat$CellSubtype[seurat$cluster == 6] <- "Late_spermatids"
seurat$CellSubtype_cluster <- paste0(seurat$CellSubtype, '-',seurat$cluster)
Idents(seurat) <- "CellSubtype"
DimPlot(seurat, label = T)+
  theme_umap
ggsave(here::here("Figures", "Figure6", "Suppl_FigS13C_seurat_subset_assigned.png"), width =5, height = 5)



######______  Re-normalize and cluster for downstream analyses 
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
ElbowPlot(seurat, ndims = 50)
seurat <- RunUMAP(seurat, dims = 1:30, reduction = "pca", min.dist = 0.3, n.neighbors = 30)
Idents(seurat) <- "CellSubtype"
DimPlot(seurat, label = F)+
  theme_umap
ggsave(here::here("Figures", "Figure6", "Figure6A_seurat_subset_final.png"), width =5, height = 5)
DimPlot(seurat, label = T)+
  theme_umap
ggsave(here::here("Figures", "Figure6", "Suppl_FigS13E_seurat_subset_final.png"), width =5, height = 5)



# Pseudotime --------------------------------------------------------------

sce <- as.SingleCellExperiment(seurat)
sce <- slingshot(sce, reducedDim = "UMAP")
slingshot(sce, reducedDim = "UMAP", )
pt <- slingPseudotime(sce)[, 1]  # NA for cells not on the lineage

umap <- reducedDims(sce)$UMAP
df <- data.frame(UMAP_1 = umap[,1],
                 UMAP_2 = umap[,2],
                 Pseudotime = pt)

ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = Pseudotime)) +
  geom_point(size =0.5) +
  scale_color_viridis_c(na.value = "gray80", direction = -1) +
  theme_classic() +
  labs(title = "Pseudotime by Slingshot",
       color = "Pseudotime")
ggsave(here::here("Figures", "Figure6", "Figure6B_seurat_pseudotime.png"), width =5.5, height = 5)

pt <- as.data.frame(pt) %>% 
  rownames_to_column("Cell") %>%
  dplyr::rename(Pseudotime = pt)
write_csv(pt, here::here("Data/Guo_et_al", "ProcessedData", "seurat_pseudotime.csv"))


# Cell cycle scoring ------------------------------------------------------

seurat <- CellCycleScoring(seurat,g2m.features = cc.genes$g2m.genes,
                           s.features = cc.genes$s.genes)
cell_cycle_scores <- seurat@meta.data %>%
  rownames_to_column("Cell") %>%
  select(Cell, G2M.Score, S.Score, Phase)

write_csv(cell_cycle_scores, here::here("Data/Guo_et_al", "ProcessedData", "seurat_cell_cycle_scores.csv"))

# Export data -------------------------------------------------------------

meta <- seurat@meta.data %>%
  rownames_to_column("Cell") %>%
  full_join(pt, by = "Cell")%>%
  column_to_rownames("Cell")
seurat <-AddMetaData(seurat, metadata = meta)
saveRDS(seurat, here::here("Data/Guo_et_al", "ProcessedData", "seurat_object.rds"))

meta <- seurat@meta.data
mat <- GetAssayData(seurat, layer = "count")
saveRDS(mat, here::here("Data/Guo_et_al", "ProcessedData", "seurat_counts.rds"))
saveRDS(meta, here::here("Data/Guo_et_al", "ProcessedData", "seurat_meta.rds"))

umapCoord <- as.data.frame(Embeddings(object = seurat[["umap"]])) %>%
  rownames_to_column("Cell")
write_csv(umapCoord, here::here("Data/Guo_et_al", "ProcessedData", "seurat_umap.csv"))


# Cell to bin -------------------------------------------------------------

# genes <- unique(gene_to_pathway$Gene)
# exprs <- as.matrix(GetAssayData(seurat, assay = "RNA", slot = "data"))
# exprs <- exprs[rownames(exprs) %in% genes, ]
cell_to_bin <- meta %>%
  rownames_to_column("Cell") %>%
  filter(!Pseudotime ==0) %>%
  arrange(Pseudotime) %>%
  mutate(
    bin = (row_number() - 1) %/% 10 + 1  
  ) %>%
  group_by(bin) %>%
  mutate(mean_pt = mean(Pseudotime)) %>%
  select(Cell, bin, mean_pt, CellSubtype) %>%
  unique() 
write_csv(cell_to_bin, here::here("Data/Guo_et_al", "ProcessedData", "cell_to_bin.csv"))


dominating_celltype <- janitor::tabyl(cell_to_bin %>%
                                        ungroup() %>%
                                        select(bin, CellSubtype) , CellSubtype,bin) %>%
  pivot_longer(cols = -CellSubtype, names_to = "bin", values_to = "n") %>%
  filter(!n==0) %>%
  mutate(bin = as.numeric(bin)) %>%
  group_by(bin) %>%
  mutate(
    dominant_cell_type = CellSubtype[which.max(n)]
  ) %>%
  ungroup() %>%
  select(bin, dominant_cell_type)

meta_binned <- cell_to_bin %>%
  select(bin, mean_pt) %>%
  full_join(dominating_celltype, by = "bin") %>%
  filter(!is.na(mean_pt)) %>%
  unique()

write_csv(meta_binned, here::here("Data/Guo_et_al", "ProcessedData", "meta_binned.csv"))



# Gene sets ---------------------------------------------------------------
mat <- readRDS(here::here("Data/Guo_et_al", "ProcessedData", "seurat_counts.rds"))

mat <- GetAssayData(seurat, layer = "count")

## Manually defined pathways
glyc <- tibble(Pathway = "Glycolysis", Gene = c(
  "HK1", "HK2", "HK3", "GCK", "GPI", "PFKL", "PFKM", "PFKP",
  "ALDOA", "ALDOB", "ALDOC", "TPI1", "GAPDH", "PGK1",
  "PGAM1", "PGAM2", "ENO1", "ENO2", "ENO3", "PKM", "PKLR"))

ppp <- tibble(Pathway = "PentosePhosphatePathway", Gene = c(
  "G6PD", "PGLS", "PGD", "RPIA", "RPE", "TKT", "TALDO1"))


extra_pathways <- bind_rows(glyc, ppp
)
gene_to_pathway <- bind_rows(gene_to_pathway, extra_pathways)
gene_sets <- split(gene_to_pathway$Gene, gene_to_pathway$Pathway)
saveRDS(gene_sets, here::here("Data/Guo_et_al", "ProcessedData", "gene_sets.rds"))
saveRDS(gene_sets, here::here("Code", "Figure6/Docker_scGSEA", "gene_sets.rds"))
saveRDS(gene_sets, here::here("/Users/annamonzel/Desktop/Docker/scGSEA/scGSEA_Guo", "gene_sets.rds"))


