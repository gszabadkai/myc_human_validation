# 03_build_covariate_table.R
# =============================================================================
# One patient x covariate table for TCGA-BRCA, assembled from scripts 01 and 02
# plus the GISTIC and clinical snapshots, with a coverage report.
#
# This is the join where things go wrong quietly, so every merge is asserted and
# every source's coverage is reported rather than assumed. Nothing is imputed
# and nothing is silently dropped.
#
# SCALE DISCIPLINE: not applicable. No expression values are read - only the
# sample list from script 01's VST object.
#
# SPECIES: human. TCGA-BRCA. See CLAUDE.md.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
suppressPackageStartupMessages(library(data.table))

message("\n03: TCGA-BRCA covariate table\n", strrep("=", 78))

PATH_ALIQUOTS <- here::here("data", "tcga_clinical",
                            "gdc_brca_rnaseq_aliquots_2026-08-28.tsv")

# -----------------------------------------------------------------------------
# 1. Sample universe, from script 01
# -----------------------------------------------------------------------------
vst <- readRDS(file.path(DIR_RESULTS, "tcga_brca_vst.rds"))
if (is.null(vst$sample_map)) {
  stop("script 01's output has no sample_map. Re-source script 01 - it is ",
       "cached, so this is quick - then re-run 03.", call. = FALSE)
}
sample_map <- vst$sample_map
patients   <- colnames(vst$mat)
rm(vst); invisible(gc())
stopifnot(identical(sort(patients), sort(sample_map$patient)))
message("1. patients from script 01: ", length(patients))

cov <- tibble::tibble(patient = patients) %>%
  dplyr::left_join(sample_map, by = "patient")

# -----------------------------------------------------------------------------
# 2. Sequencing plate
# -----------------------------------------------------------------------------
# Plate is a pre-specified Block C covariate (plan section 8, "Plate / batch,
# TCGA barcode"). It lives in field 6 of the ALIQUOT barcode, which the Xena
# matrix does not carry - its columns are sample-level. So it is joined from a
# dated GDC metadata snapshot, keyed on the sample id script 01 actually kept.
ali <- data.table::fread(PATH_ALIQUOTS, sep = "\t", header = TRUE,
                         showProgress = FALSE)
ali <- ali[sample_id %in% cov$sample_id]

# Five sample ids carry two RNA-seq aliquots each, one of them re-sequenced on
# plate A277. These are exactly the five Xena averaged (confirmed independently:
# the GDC aliquot count and the half-integer signature identify the same five).
# Script 01 already excluded the averaged columns, so anything still duplicated
# here is a genuine two-aliquot sample; take the lowest plate code for a
# deterministic answer and report it.
dup_ali <- ali$sample_id[duplicated(ali$sample_id)]
if (length(dup_ali)) {
  message("   ", length(unique(dup_ali)),
          " sample id(s) with multiple aliquots; taking lowest plate code")
  data.table::setorder(ali, sample_id, plate)
  ali <- ali[!duplicated(sample_id)]
}

cov <- cov %>%
  dplyr::left_join(tibble::as_tibble(ali[, .(sample_id, plate, aliquot_barcode)]),
                   by = "sample_id")
message("2. plate assigned for ", sum(!is.na(cov$plate)), " of ", nrow(cov),
        " patients across ", dplyr::n_distinct(cov$plate, na.rm = TRUE), " plates")

# -----------------------------------------------------------------------------
# 3. Genomics, from script 02
# -----------------------------------------------------------------------------
gpath <- file.path(DIR_RESULTS, "tcga_brca_genomics.rds")
if (!file.exists(gpath)) {
  stop("results/tcga_brca_genomics.rds not found - run script 02 first.",
       call. = FALSE)
}
gen <- readRDS(gpath)

cov <- cov %>%
  dplyr::left_join(gen$mutations,  by = "patient") %>%
  dplyr::left_join(gen$purity,     by = "patient") %>%
  dplyr::left_join(gen$leukocyte,  by = "patient") %>%
  dplyr::left_join(gen$aneuploidy, by = "patient") %>%
  dplyr::left_join(gen$mc3_covered, by = "patient")
message("3. genomics joined")

# -----------------------------------------------------------------------------
# 4. GISTIC copy number
# -----------------------------------------------------------------------------
# ISAR is PRIMARY, matching G2 (see data/gistic_tcga_brca/README.md). Firehose
# is carried alongside as the pre-specified sensitivity source. Both are read
# with the SAME sample rule G2 used, so the CNV columns here are the same calls
# G2 tested - not a re-derivation that might drift.
G2_GENES <- c("MYC", "MCL1", "BCL2L1", "BBC3", "BAX", "PTEN", "TP53", "PIK3CA")

.read_gistic_subset <- function(path, cmd = NULL, label) {
  dt <- if (is.null(cmd)) {
    data.table::fread(path, sep = "\t", header = TRUE, showProgress = FALSE)
  } else {
    data.table::fread(cmd = cmd, sep = "\t", header = TRUE, showProgress = FALSE)
  }
  data.table::setnames(dt, 1:3, c("gene_symbol", "locus_id", "cytoband"))
  bc  <- setdiff(names(dt), c("gene_symbol", "locus_id", "cytoband"))
  pat <- vapply(strsplit(bc, "-"), function(p) paste(p[1:3], collapse = "-"),
                character(1))
  sty <- substr(vapply(strsplit(bc, "-"), `[`, character(1), 4), 1, 2)
  keep <- bc[sty == "01" & pat %in% patients]
  kpat <- vapply(strsplit(keep, "-"), function(p) paste(p[1:3], collapse = "-"),
                 character(1))
  stopifnot(!any(duplicated(kpat)))
  sub <- dt[gene_symbol %in% G2_GENES, c("gene_symbol", keep), with = FALSE]
  m <- as.matrix(sub[, keep, with = FALSE]); storage.mode(m) <- "integer"
  rownames(m) <- sub$gene_symbol; colnames(m) <- kpat
  message("   ", label, ": ", ncol(m), " of ", length(patients), " patients")
  m
}

message("4. GISTIC")
g_isar <- .read_gistic_subset(
  PATH_GISTIC_ISAR,
  cmd = paste("gzip -dc", shQuote(PATH_GISTIC_ISAR)), label = "ISAR (primary)")
g_fh <- .read_gistic_subset(PATH_GISTIC_FH_TXT, label = "Firehose (sensitivity)")

.cnv_col <- function(m, gene) {
  v <- rep(NA_integer_, length(patients)); names(v) <- patients
  present <- intersect(patients, colnames(m))
  v[present] <- m[gene, present]
  unname(v)
}
for (g in G2_GENES) {
  cov[[paste0("cnv_", g)]]    <- .cnv_col(g_isar, g)
  cov[[paste0("cnvfh_", g)]]  <- .cnv_col(g_fh,   g)
}

# -----------------------------------------------------------------------------
# 5. Derived variables - the pre-registered definitions
# -----------------------------------------------------------------------------
# Thresholds from the G2 design note section 2.3 (gene- AND direction-specific).
# Alteration rules from the script 02 header, agreed 2026-08-28.
cov <- cov %>%
  dplyr::mutate(
    # H2 / H3 strata
    TP53_status    = dplyr::if_else(is.na(TP53_mut), NA_character_,
                                    dplyr::if_else(TP53_mut, "mutant", "wildtype")),
    PIK3CA_altered = PIK3CA_mut,
    # PTEN: mutation OR GISTIC deep deletion. Either alone misclassifies.
    PTEN_altered   = dplyr::case_when(
      is.na(PTEN_mut) & is.na(cnv_PTEN) ~ NA,
      (!is.na(PTEN_mut) & PTEN_mut) | (!is.na(cnv_PTEN) & cnv_PTEN == -2L) ~ TRUE,
      TRUE ~ FALSE
    ),
    PI3K_pathway_intact = dplyr::case_when(
      is.na(PIK3CA_altered) | is.na(PTEN_altered) ~ NA,
      !PIK3CA_altered & !PTEN_altered ~ TRUE,
      TRUE ~ FALSE
    ),
    # BUFFER, at the G2 primary thresholds
    buffer_MCL1   = cnv_MCL1   ==  2L,
    buffer_BCL2L1 = cnv_BCL2L1 >=  1L,
    BUFFER = dplyr::case_when(
      is.na(buffer_MCL1) & is.na(buffer_BCL2L1) ~ NA,
      (!is.na(buffer_MCL1) & buffer_MCL1) |
        (!is.na(buffer_BCL2L1) & buffer_BCL2L1) ~ TRUE,
      TRUE ~ FALSE
    ),
    MYC_amp   = cnv_MYC  == 2L,
    BBC3_loss = cnv_BBC3 <= -1L
  )

# -----------------------------------------------------------------------------
# 6. Clinical: PAM50, receptor status, TNBC
# -----------------------------------------------------------------------------
clin <- readr::read_tsv(PATH_TCGA_CLINICAL, show_col_types = FALSE, na = "NA")
.clean_rec <- function(x) ifelse(x %in% c("Positive", "Negative"), x, NA_character_)

clin <- clin %>%
  dplyr::transmute(
    patient  = patient_barcode,
    PAM50    = dplyr::na_if(PAM50_SUBTYPE, "NA"),
    er_call  = .clean_rec(ER_STATUS_BY_IHC),
    pr_call  = .clean_rec(PR_STATUS_BY_IHC),
    her2_call = dplyr::coalesce(.clean_rec(HER2_FISH_STATUS), .clean_rec(IHC_HER2)),
    aneuploidy_cbio = as.numeric(ANEUPLOIDY_SCORE),
    FGA             = as.numeric(FRACTION_GENOME_ALTERED)
  ) %>%
  dplyr::mutate(
    TNBC = dplyr::case_when(
      is.na(er_call) | is.na(pr_call) | is.na(her2_call) ~ NA,
      er_call == "Negative" & pr_call == "Negative" & her2_call == "Negative" ~ TRUE,
      TRUE ~ FALSE
    )
  )
cov <- dplyr::left_join(cov, clin, by = "patient")
message("5. clinical joined: PAM50 for ", sum(!is.na(cov$PAM50)),
        ", TNBC callable for ", sum(!is.na(cov$TNBC)), " patients")

# -----------------------------------------------------------------------------
# 7. Coverage report
# -----------------------------------------------------------------------------
# The number that matters is not how many patients exist, but how many survive
# the primary model's complete-case requirement. Reported here so the loss is a
# known quantity before any model is fitted, not a surprise afterwards.
BLOCK_C_VARS <- c("purity", "leukocyte_fraction", "PAM50", "TP53_status", "plate")

coverage <- tibble::tibble(
  variable = names(cov),
  n_present = vapply(cov, function(x) sum(!is.na(x)), integer(1)),
  pct = round(100 * vapply(cov, function(x) mean(!is.na(x)), numeric(1)), 1)
) %>% dplyr::arrange(n_present)

message("\n6. coverage, least complete first:")
coverage %>%
  dplyr::filter(n_present < nrow(cov)) %>%
  as.data.frame() %>% head(15) %>% print()

complete_C <- stats::complete.cases(cov[, BLOCK_C_VARS])
message("\n   Block C complete cases (", paste(BLOCK_C_VARS, collapse = ", "), "):")
message("   ", sum(complete_C), " of ", nrow(cov),
        sprintf("  (%.1f%%)", 100 * mean(complete_C)))
cov$complete_block_c <- complete_C

# >>> PRE-REGISTRATION DECISION, NEEDS SIGN-OFF BEFORE ANY MODEL IS FITTED.
# Missing-data policy. Proposed: COMPLETE CASES for the primary Block C model,
# with the n reported alongside every estimate, and no imputation. Reasons:
# the missingness is in covariates rather than the exposure or endpoint; it is
# driven by assay coverage rather than by biology; and multiple imputation of
# PAM50 or plate would invent structure the design depends on.
# Pre-specified sensitivity: refit on all patients with the covariates that are
# fully observed, to show the complete-case restriction is not doing the work.

message("\n7. sample-set reconciliation")
message("   expression (script 01)      : ", length(patients))
message("   with ISAR GISTIC            : ", sum(!is.na(cov$cnv_MYC)))
message("   with Firehose GISTIC        : ", sum(!is.na(cov$cnvfh_MYC)))
message("   with MC3 coverage           : ", sum(cov$mc3_covered, na.rm = TRUE))
message("   with ABSOLUTE purity        : ", sum(!is.na(cov$purity)))
message("   Block C complete            : ", sum(complete_C))

# -----------------------------------------------------------------------------
# 8. Save
# -----------------------------------------------------------------------------
saveRDS(list(covariates = cov, coverage = coverage,
             block_c_vars = BLOCK_C_VARS,
             gistic_isar = g_isar, gistic_firehose = g_fh,
             built = Sys.time()),
        file.path(DIR_RESULTS, "tcga_brca_covariates.rds"))
readr::write_csv(cov, file.path(DIR_TABLES, "tcga_brca_covariate_table.csv"))
readr::write_csv(coverage, file.path(DIR_TABLES, "tcga_brca_covariate_coverage.csv"))

message("\n03: done. results/tcga_brca_covariates.rds + 2 tables")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  cv <- readRDS(file.path(DIR_RESULTS, "tcga_brca_covariates.rds"))$covariates

  # --- do the strata match published BRCA frequencies? ---------------------
  # TP53 ~30-35%, PIK3CA ~35%, PTEN-altered ~10-15%. A large deviation means a
  # rule is wrong, not that TCGA is unusual.
  round(100 * c(
    TP53   = mean(cv$TP53_status == "mutant", na.rm = TRUE),
    PIK3CA = mean(cv$PIK3CA_altered, na.rm = TRUE),
    PTEN   = mean(cv$PTEN_altered, na.rm = TRUE),
    BUFFER = mean(cv$BUFFER, na.rm = TRUE),
    MYCamp = mean(cv$MYC_amp, na.rm = TRUE)
  ), 1)

  # --- does BUFFER reproduce G2 in this joined table? ----------------------
  # It should: same GISTIC calls, same thresholds. A mismatch means the join
  # lost or shifted samples.
  with(cv, table(MYC_amp, buffer_BCL2L1, useNA = "ifany"))
  with(cv, table(MYC_amp, buffer_MCL1,   useNA = "ifany"))

  # --- what is complete-case actually costing, and is it random? -----------
  # If the dropped patients differ by subtype, complete-case is not innocuous
  # and the sensitivity refit matters.
  with(cv, table(PAM50, complete_block_c, useNA = "ifany"))
  with(cv, table(TNBC,  complete_block_c, useNA = "ifany"))

  # --- plate: enough levels, and none tiny enough to be unfittable? --------
  sort(table(cv$plate))

  # --- purity vs leukocyte fraction, the two stromal covariates ------------
  cor(cv$purity, cv$leukocyte_fraction, use = "complete.obs", method = "spearman")

  # --- ISAR vs Firehose disagreement on the calls that matter --------------
  with(cv, table(ISAR = cnv_MYC, Firehose = cnvfh_MYC, useNA = "ifany"))

}
