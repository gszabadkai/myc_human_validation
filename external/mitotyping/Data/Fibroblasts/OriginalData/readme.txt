The human lifespan study dataset: 
- Download from GEO GSE179848 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179848 
- Place the gene_count_lengthScaledTPM.csv in the "Data/Fibroblasts/OriginalData" folder. 
- The meta data can be extracted from "SraRunTable.csv" or from the shiny app https://columbia-picard.shinyapps.io/MitotypeExplorer/. I'll leave a compressed copy of the metadata in the "Data/Fibroblasts/OriginalData" folder.
- Note that the data on GEO was prepared by Sturm et al. 2022, with different R software and package versions.
It might not 100% match the data prepared using the markdown code above (kallisto, tximport). 
This does not affect the figures overall, but some statistics might differ slightly (e.g. the PCA variances, the top pathways in passage correlation). The raw data that was used to create the figures can be shared upon request. ALso happy to rerun the analysis in the future just with the GEO dataset. 
