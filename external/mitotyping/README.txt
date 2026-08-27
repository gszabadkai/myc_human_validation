###__________ Data availability __________###

In order to run the code, the data needs to be downloaded from their respective public repositories:

Mouse and human MitoCarta3.0 datasets: 
- Download from https://www.broadinstitute.org/mitocarta
- Place "HumanMitoCarta3_0.xls" and "MouseMitoCarta3_0.xls" in the "Data/MitoCarta/OriginalData" folder.

Human Protein Atlas dataset version 21.0: 
- Download from https://v21.proteinatlas.org/about/download. 
- Important note: the publicly available version 21.1 no longer includes the olfactory bulb.  For robustness analysis, we re-ran the code without the olfactory bulb, and found no major impact on the results, but the dataset only has 54 tissues. Some stats might be slightly different. 
- Place "rn_human_tissue_consensus.tsv" in the "Data/HumanProteinAtlas/OriginalData" folder.

GTEx dataset version 8: 
- Download from https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression
- No protected data needed
- Place all files (gene_tpm, sample attributes, subject phenotypes, gene read, gene_median_tpm) in the folder "Data/GTEx/OriginalData" and extract the .gz files in the same folder.

The human lifespan study dataset: 
- Download from GEO GSE179848 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179848 
- Place the gene_count_lengthScaledTPM.csv in the "Data/Fibroblasts/OriginalData" folder. 
- The meta data can be extracted from "SraRunTable.csv" or from the shiny app https://columbia-picard.shinyapps.io/MitotypeExplorer/. I'll leave a compressed copy of the metadata in the "Data/Fibroblasts/OriginalData" folder.
- Note that the data on GEO was prepared by Sturm et al. 2022, with different R software and package versions.
It might not 100% match the data prepared using the markdown code above (kallisto, tximport). 
This does not affect the figures overall, but some statistics might differ slightly (e.g. the PCA variances, the top pathways in passage correlation). The raw data that was used to create the figures can be shared upon request. ALso happy to rerun the analysis in the future just with the GEO dataset. 

###__________ General __________###

The code can be executed from the main.R file. The code is divided into different sections, each corresponding to a different figure or analysis. The main file will also generate all subfolders necessary. 

###__________ Processed Data __________###

When running the code, data will be taken from the "OriginalData" folders and saved in the "ProcessedData" folders. This is done to avoid running the same data processing steps multiple times.

###__________ Figures __________###

Figures will be automatically saved in the "Figures" folders. Some figures created with pdf() and dev.off() might not appear in the figures folder when sourced in the "main.R" file. In this case, run script individually.
The code for the 3D plots is also added, however, users might experience issues with the execution. For macOS users that have XQuartz installed, the code should work without any issues. To avoid any errors, the code for the 3D plots has been commented out. This affects the code for figure 5 and supplemental figure s8.
