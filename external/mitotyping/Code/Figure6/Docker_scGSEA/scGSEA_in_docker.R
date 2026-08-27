options(timeout = 1200)
require(msigdbr)
library(gficf)
library(ggplot2)
counts <- readRDS("project/Final/seurat_counts.rds")
meta <- readRDS("project/Final/seurat_meta.rds")

# Data normalization and gene filtering
data <- gficf( M = counts)

# Create PCA-subspace using overdispersed genes
data <- runPCA(data = data,dim = 10,use.odgenes = T)

# Create t-UMAP space
data <-runReduction(data = data,reduction = "umap",nt = 2,verbose = T,n_neighbors=150)

# Cell meta-data can stored in the data$embedded data.frame
# Let' add the info about the cell-line, stripping this information
# from the name of the cell and storing it into ccl column.
data$embedded$ccl = sapply(
  strsplit(x = rownames(data$embedded),
           split = "_",fixed = T)
  ,function(x) x[1]
)

gene_set <- readRDS("project/Final/gene_sets.rds")
#data <- resetScGSEA(data)
#sample(1:1000, 10)
seeds <- c(430, 803, 374, 362, 489,  23, 712, 826, 772, 361)
list <- list()
for (i in 1:length(seeds)){
#set.seed(seeds[i])
data = runScGSEA(data = data,
                 geneID = "symbol",
                 species = "human",
                 #category = "H",
                 pathway.list = gene_set,
                 minSize = 1, 
                 nmf.k = 100,
                 fdr.th = .05,
                 seed = seeds[i],
                 rescale = "byGS",
                 verbose = T)
list[[i]] <- data
names(list)[i] <- paste0("Seed_", seeds[i])
}

#final <- unlist(list)

saveRDS(list, here::here("project/Final/scGSEA_testis_multiseeds10_list_5pct.rds"))




