generate_heatmaps_for_gene_set <- function(
    set_name,
    gene_symbols,
    dds_obj,
    lfc_variant = c("shrunk", "raw"),
    output_dir = "outputs/heatmaps_int/",
    show_plot = TRUE
) {
  require(ComplexHeatmap)
  require(circlize)
  require(dplyr)
  require(tibble)
  require(grid)
  select <- dplyr::select
  
  lfc_variant <- match.arg(lfc_variant)
  
  # === Load appropriate LFC data ===
  lfc_df <- switch(
    lfc_variant,
    shrunk = readRDS("results/combined_df_annotated.rds"),
    raw = readRDS("results/combined_df_annotated_raw.rds")
  )
  
  # === Select LFC columns ===
  lfc_columns <- if (lfc_variant == "shrunk") {
    c("myc_6W_log2FC", "myc_12W_log2FC", "timepoint_neg_log2FC")
  } else {
    c("myc_6W_log2FC_raw", "myc_12W_log2FC_raw", "timepoint_neg_log2FC_raw")
  }
  
  # === Filter LFC matrix ===
  
  lfc_mat <- lfc_df %>%
    filter(mgi_symbol %in% gene_symbols) %>%
    distinct(mgi_symbol, .keep_all = TRUE) %>%   # REMOVE duplicates
    select(mgi_symbol, all_of(lfc_columns)) %>%
    column_to_rownames("mgi_symbol") %>%
    as.matrix()
  
  
  if (nrow(lfc_mat) == 0) {
    warning(paste("No genes found for set:", set_name))
    return(NULL)
  }
  
  # === Expression matrix ===
  norm_counts <- counts(dds_obj, normalized = TRUE)
  log_expr <- log2(norm_counts + 1)
  
  gene_map <- lfc_df %>% select(gene, mgi_symbol) %>% distinct()
  log_expr <- log_expr[rownames(log_expr) %in% gene_map$gene, ]
  rownames(log_expr) <- gene_map$mgi_symbol[match(rownames(log_expr), gene_map$gene)]
  log_expr <- log_expr[!is.na(rownames(log_expr)), ]
  
  coldata <- as.data.frame(colData(dds_obj))
  group_levels <- levels(coldata$group)
  group_means <- sapply(group_levels, function(g) {
    samples <- rownames(coldata)[coldata$group == g]
    rowMeans(log_expr[, samples, drop = FALSE])
  })
  rownames(group_means) <- rownames(log_expr)
  
  # === Filter genes with complete data
  common_genes <- intersect(rownames(lfc_mat), rownames(group_means))
  lfc_mat <- lfc_mat[common_genes, , drop = FALSE]
  expr_mat <- group_means[common_genes, , drop = FALSE]
  expr_scaled <- t(scale(t(expr_mat)))
  
  # === Row annotations
  row_meta <- lfc_df %>%
    filter(mgi_symbol %in% rownames(lfc_mat)) %>%
    select(mgi_symbol, starts_with("time_effect_category"), group_sig_status) %>%
    distinct() %>%
    column_to_rownames("mgi_symbol") %>%
    .[rownames(lfc_mat), , drop = FALSE]
  
  # Define colors
  time_colors <- c(
    "direct_Myc_reduction" = "firebrick",
    "baseline_driven_reduction" = "purple",
    "baseline_driven_increase" = "steelblue",
    "direct_Myc_increase" = "darkgreen",
    "no_change" = "lightgrey"
  )
  
  sig_colors <- c(
    up = "blue",
    down = "orange",
    ns = "grey90"
  )
  
  # Build row annotation
  ra_list <- list()
  time_col <- grep("time_effect_category", colnames(row_meta), value = TRUE)
  if (length(time_col) > 0) {
    ra_list$Time_Effect <- row_meta[[time_col]]
  }
  if ("group_sig_status" %in% colnames(row_meta)) {
    ra_list$Group_Change <- row_meta$group_sig_status
  }
  
  row_ha <- rowAnnotation(
    df = as.data.frame(ra_list),
    col = list(
      Time_Effect = time_colors,
      Group_Change = sig_colors
    ),
    show_annotation_name = TRUE
  )
  
  # Output path
  safe_name <- gsub("[/\\|:*?\"<>]", "_", set_name)
  pdf_path <- file.path(output_dir, paste0(safe_name, "_", lfc_variant, ".pdf"))
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # === Heatmaps ===
  # === Dynamic height based on number of genes ===
  n_genes <- nrow(lfc_mat)
  
  row_font_max <- 11
  row_font_min <- 6
  fontsize_row <- max(row_font_min, min(row_font_max, 12 - 0.05 * n_genes))
  
  ht1 <- Heatmap(
    lfc_mat,
    name = "log2FC",
    cluster_rows = TRUE,
    clustering_distance_rows = "pearson",     
    clustering_method_rows = "ward.D2", 
    cluster_columns = FALSE,
    column_order = lfc_columns,
    show_column_names = TRUE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = fontsize_row),
    column_names_gp = gpar(fontsize = 6),
    cell_fun = function(j, i, x, y, ...) {
      grid.text(sprintf("%.2f", lfc_mat[i, j]), x, y, gp = gpar(fontsize = 8))
    }
  )
  
  ht2 <- Heatmap(
    expr_scaled,
    name = "log2norm_z",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_order = c("6W_neg", "6W_pos", "12W_neg", "12W_pos"),
    show_column_names = TRUE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = fontsize_row),
    column_names_gp = gpar(fontsize = 6),
    cell_fun = function(j, i, x, y, ...) {
      grid.text(sprintf("%.2f", expr_mat[i, j]), x, y, gp = gpar(fontsize = 8))
    }
  )
  
  ht_list <- ht1 + ht2
  
  # Save and display
  height_per_gene <- 0.25  # adjust for font size + cell padding
  min_height <- 4
  max_height <- 20
  pdf_height <- max(min_height, min(n_genes * height_per_gene, max_height))
  
  pdf(pdf_path, width = 7, height = pdf_height)
  grid.newpage()
  draw(row_ha + ht1 + ht2, merge_legend = TRUE)
  dev.off()
  
  if (show_plot) {
    grid.newpage()
    draw(row_ha + ht1 + ht2, merge_legend = TRUE)
  }
  
  return(list(
    pdf_file = pdf_path,
    lfc_matrix = lfc_mat
  ))
}
