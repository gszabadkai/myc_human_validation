# 01_fetch_tcga_expression.R
# =============================================================================
# TCGA-BRCA expression: fetch GDC harmonised STAR counts, build the two
# normalised matrices the rest of the pipeline needs, and discharge the
# outstanding half of gate G1.
#
# =============================================================================
# SCALE DISCIPLINE - the most likely silent error in this project.
# =============================================================================
# This script produces TWO matrices with OPPOSITE scale requirements. They are
# saved to separate files and MUST NOT be merged, coerced into one another, or
# passed to the wrong consumer.
#
#   results/tcga_brca_vst.rds       LOG scale (VST).  For GSVA / ssGSEA, with
#                                   kcdf = "Gaussian".  Consumers: 06, 07, 08.
#
#   results/tcga_brca_linear.rds    LINEAR scale (DESeq2 size-factor normalised
#                                   counts, NOT logged).  For mitoPPS only.
#                                   Consumer: 07, mitoPPS block only.
#
# If you find yourself writing log2() on the linear object or expm1() on the VST
# object, stop: you are about to reintroduce the error this split exists to
# prevent.  See CLAUDE.md, "Scale discipline".
#
# SPECIES: human.  TCGA-BRCA, human ENSG/HGNC.  See CLAUDE.md.
#
# COHORT-RELATIVITY: these matrices are TCGA-only.  GSVA scores derived from
# them are cohort-relative and are never comparable to scores from METABRIC,
# SCAN-B or the neoadjuvant cohorts.  Meta-analyse effect estimates, not scores.
#
# DOWNLOAD: several GB into data/raw/, which is gitignored and NOT on origin.
# The raw download is cached; re-sourcing this script does not re-download.
# data/raw/ needs its own backup outside git.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

suppressPackageStartupMessages({
  library(SummarizedExperiment)   # for assay() on the DESeqTransform
  library(DESeq2)
  library(data.table)
})

message("\n01: TCGA-BRCA expression\n", strrep("=", 78))

DIR_TCGA_RAW <- file.path(DIR_DATA, "raw", "tcga_brca_expression")
.ensure_dir(DIR_TCGA_RAW)

PATH_SE      <- file.path(DIR_TCGA_RAW, "tcga_brca_star_counts.rds")
PATH_VST     <- file.path(DIR_RESULTS, "tcga_brca_vst.rds")
PATH_LINEAR  <- file.path(DIR_RESULTS, "tcga_brca_linear.rds")
PATH_G1CORR  <- file.path(DIR_RESULTS, "g1_correlation_criterion.rds")

# -----------------------------------------------------------------------------
# 1. Fetch: UCSC Xena STAR counts matrix
# -----------------------------------------------------------------------------
# SOURCE DECISION, 2026-08-28. The plan says "GDC harmonised STAR counts via
# TCGAbiolinks (or recount3)". The GDC per-file route was tried first and is
# impractical here: 5.2 GB across 1231 files, measured at 0.05 files/s, i.e.
# about 7 hours, over a transport that failed twice on the way.
#
# Xena re-hosts the SAME GDC STAR-Counts quantification as a single 138 MB
# matrix. Same pipeline, same GENCODE v36, same versioned ENSG ids.
#
# What this costs, stated rather than waved away:
#   1. Provenance is one hop longer - Xena's snapshot of a GDC release, not the
#      GDC API directly. Mitigated by snapshotting with a SHA-256 and a README,
#      exactly as done for CollecTRI and GISTIC.
#   2. Values arrive as log2(count + 1) and must be inverted. Verified lossless:
#      the stored values are exact log2 of integers (3.321928... = log2(10),
#      3.0 = log2(8), 0.0 = log2(1)), so round(2^x - 1) recovers raw counts
#      exactly. The inversion is asserted below, not assumed.
#   3. Xena columns are sample-level with vial (TCGA-B6-A1KC-01B), so some
#      aliquot collapsing is already baked in upstream and is not visible to us.
#      Section 2 therefore uses a slightly different rule from the GDC path;
#      it is documented there.
#
# Gene annotation does NOT come from Xena. Xena ships only ENSG ids, and its
# published probeMap is GENCODE v22 against v36 data. Instead the annotation is
# read from one of the GDC per-file downloads already on disk, which carries
# gene_id / gene_name / gene_type at the matching GENCODE v36. Same build, local,
# already provenanced.

PATH_XENA <- file.path(DIR_DATA, "raw", "tcga_expression_xena",
                       "TCGA-BRCA.star_counts.tsv.gz")
XENA_SHA256 <- "058d79121460c535b73312247a55d3108d18ea6ddd0ccc1070b5dd93ceeedaa4"

if (!file.exists(PATH_XENA)) {
  stop("Xena matrix not found at ", PATH_XENA,
       "\nSee data/tcga_expression_xena/README.md for the URL and checksum.",
       call. = FALSE)
}

if (file.exists(PATH_SE)) {
  message("1. cached counts matrix found, skipping rebuild")
  cached <- readRDS(PATH_SE)
  counts_raw   <- cached$counts
  gene_ann     <- cached$annotation
  averaged_col <- cached$averaged_col
} else {
  message("1. reading Xena STAR counts matrix (138 MB, log2 scale on disk)")
  xz <- data.table::fread(cmd = paste("gzip -dc", shQuote(PATH_XENA)),
                          sep = "\t", header = TRUE, showProgress = FALSE)
  gid <- xz[[1]]
  xm  <- as.matrix(xz[, -1, with = FALSE])
  rownames(xm) <- gid
  message("   ", nrow(xm), " genes x ", ncol(xm), " samples")

  # --- invert log2(count + 1) -----------------------------------------------
  # Verified against 4 GDC per-file downloads: after inversion and rounding the
  # values match GDC "unstranded" raw counts EXACTLY across 4,000 genes, on a
  # blind match (each GDC sample identified exactly one Xena column).
  #
  # BUT the inversion is not integral everywhere, and the exception matters.
  # Deviations are either ~0 or EXACTLY 0.5 - never in between - and they are
  # confined to 5 of 1,226 columns, in which ~49% of genes are half-integers.
  # That is the signature of averaging two integers, and every affected column
  # is a "-01A" for a patient who also has a "-01B"/"-01C" column.
  #
  # Xena has AVERAGED replicate aliquots into those five columns. They are not
  # raw counts, DESeq2 must not be given them, and section 2 must not select
  # them. They are detected here rather than assumed, so a future Xena release
  # that averages different samples is caught rather than silently accepted.
  message("   inverting log2(count + 1) -> raw counts")
  lin <- 2^xm - 1
  frac <- abs(lin - round(lin))

  avg_flag <- apply(frac, 2L, function(v) any(v > 0.25, na.rm = TRUE))
  if (any(avg_flag)) {
    message("   ", sum(avg_flag), " column(s) contain half-integer values, ",
            "i.e. Xena aliquot averages, and are flagged:")
    for (nm in colnames(lin)[avg_flag]) message("     ", nm)
  }
  # Any residual non-integrality OUTSIDE those columns would mean the values are
  # not log2(raw count + 1) at all, and that is fatal.
  resid <- max(frac[, !avg_flag, drop = FALSE], na.rm = TRUE)
  message(sprintf("   max deviation outside those columns: %.3e", resid))
  if (resid > 1e-3) {
    stop("log2 inversion is not returning integers even outside the averaged ",
         "columns (max deviation ", resid, "). The Xena values are not ",
         "log2(raw count + 1) as assumed; stop and re-check the source.",
         call. = FALSE)
  }

  counts_raw <- round(lin)
  storage.mode(counts_raw) <- "integer"
  # Held as a PLAIN VARIABLE, not an attribute on the matrix. Attributes are
  # silently dropped by `m[rows, ]`, and the row-alignment step below does
  # exactly that - which on the first run left this empty and let all five
  # averaged columns through into the output.
  averaged_col <- colnames(lin)[avg_flag]
  rm(xz, xm, lin, frac); invisible(gc())

  # --- annotation from a local GDC file, GENCODE v36 ------------------------
  gdc_any <- list.files(DIR_TCGA_RAW, pattern = "augmented_star_gene_counts\\.tsv$",
                        recursive = TRUE, full.names = TRUE)
  if (!length(gdc_any)) {
    stop("no GDC per-file download found to take gene annotation from. ",
         "Keep at least one *.augmented_star_gene_counts.tsv under ",
         DIR_TCGA_RAW, call. = FALSE)
  }
  ann <- data.table::fread(gdc_any[1], sep = "\t", skip = 1L,
                           showProgress = FALSE)
  gene_ann <- as.data.frame(ann[, c("gene_id", "gene_name", "gene_type")])
  gene_ann <- gene_ann[!startsWith(gene_ann$gene_id, "N_"), ]
  message("   annotation from ", basename(gdc_any[1]), ": ",
          nrow(gene_ann), " genes (GENCODE v36)")

  saveRDS(list(counts = counts_raw, annotation = gene_ann,
               averaged_col = averaged_col,
               source = "UCSC Xena GDC hub, TCGA-BRCA.star_counts.tsv.gz",
               sha256 = XENA_SHA256, built = Sys.time()), PATH_SE)
  message("   cached to ", PATH_SE)
}

# Align matrix rows to the annotation. Xena and GDC share GENCODE v36 versioned
# ids, so this is an exact join; anything unmatched is a real inconsistency and
# is reported rather than silently dropped.
common <- intersect(rownames(counts_raw), gene_ann$gene_id)
message("   gene ids matched to annotation: ", length(common),
        " of ", nrow(counts_raw))
if (length(common) < 0.95 * nrow(counts_raw)) {
  stop("fewer than 95% of Xena gene ids matched the GDC annotation. ",
       "The two are probably on different GENCODE builds.", call. = FALSE)
}
counts_raw <- counts_raw[common, , drop = FALSE]
gene_ann   <- gene_ann[match(common, gene_ann$gene_id), ]
stopifnot(identical(rownames(counts_raw), gene_ann$gene_id))

# -----------------------------------------------------------------------------
# 2. Samples: primary tumours only, one per patient
# -----------------------------------------------------------------------------
# Same INTENT as G2 (see data/gistic_tcga_brca/README.md): primary solid tumour
# only, one column per patient, so the expression and copy-number analyses stay
# aligned.
#
# The mechanics differ from the GDC per-file route, and the difference matters.
# Xena columns are SAMPLE-level with the vial letter (TCGA-B6-A1KC-01B) rather
# than full aliquot barcodes.
#
# TIE-BREAK, and why it is not simply "vial A". Section 1 finds that Xena has
# AVERAGED replicate aliquots into 5 columns, all of them "-01A" for patients
# that also carry a clean "-01B"/"-01C". A naive "keep the lowest vial" rule
# would therefore select the averaged, non-integer column for every one of those
# patients and discard the clean one - the exact wrong choice, silently, in 5
# cases.
#
# So the rule is, in order:
#   1. drop averaged columns outright where the patient has a clean alternative
#   2. among what remains, keep the LOWEST vial letter (A before B before C)
#
# A patient whose ONLY column is averaged would be kept and reported, since
# dropping the patient entirely is worse than carrying one averaged sample. On
# the current release no patient is in that position.

bc  <- colnames(counts_raw)
pat <- vapply(strsplit(bc, "-"), function(p) paste(p[1:3], collapse = "-"), character(1))
sty <- substr(vapply(strsplit(bc, "-"), `[`, character(1), 4), 1, 2)
vial <- substr(vapply(strsplit(bc, "-"), `[`, character(1), 4), 3, 3)

message("\n2. sample types present: ",
        paste(sprintf("%s=%d", names(table(sty)), table(sty)), collapse = " "))

stopifnot(exists("averaged_col"))
averaged <- colnames(counts_raw) %in% averaged_col
# Guard against the failure this rule exists to prevent: if the flag list is
# non-empty but nothing matches, the list has been lost somewhere upstream and
# the averaged columns would sail through unnoticed.
if (length(averaged_col) > 0L && !any(averaged)) {
  stop(length(averaged_col), " averaged column(s) were flagged in section 1 but ",
       "none match the matrix columns. The flag list has been lost - do not ",
       "proceed, the averaged aliquots would be selected silently.", call. = FALSE)
}

keep       <- which(sty == "01")
counts_raw <- counts_raw[, keep, drop = FALSE]
pat        <- pat[keep]; vial <- vial[keep]; averaged <- averaged[keep]

message("   Xena-averaged columns among primary tumours: ", sum(averaged))

# Drop an averaged column only when that patient has a clean alternative.
has_clean <- pat %in% pat[!averaged]
drop_avg  <- averaged & has_clean
if (any(drop_avg)) {
  message("   dropping ", sum(drop_avg),
          " averaged column(s) in favour of a clean aliquot for the same patient")
  counts_raw <- counts_raw[, !drop_avg, drop = FALSE]
  pat <- pat[!drop_avg]; vial <- vial[!drop_avg]; averaged <- averaged[!drop_avg]
}
if (any(averaged)) {
  message("   NOTE: ", sum(averaged), " patient(s) retained with an averaged ",
          "column because no clean aliquot exists: ",
          paste(pat[averaged], collapse = ", "))
}

if (any(duplicated(pat))) {
  ord        <- order(pat, vial)          # A before B before C
  counts_raw <- counts_raw[, ord, drop = FALSE]
  pat        <- pat[ord]; vial <- vial[ord]
  first      <- !duplicated(pat)
  message("   ", sum(!first), " further duplicate column(s) dropped (kept lowest vial)")
  counts_raw <- counts_raw[, first, drop = FALSE]
  pat        <- pat[first]
}
# Record WHICH sample id survived for each patient, before the columns are
# renamed to patient barcodes. Script 03 needs it to join the sequencing plate,
# which lives in the aliquot barcode and cannot be recovered from a patient id.
# Re-deriving this in 03 by re-applying the rules would put the same logic in
# two files, which is how they drift apart.
sample_map <- tibble::tibble(patient = pat, sample_id = colnames(counts_raw))

colnames(counts_raw) <- pat
stopifnot(!any(duplicated(colnames(counts_raw))))
message("   primary tumours, one per patient: ", ncol(counts_raw))

# -----------------------------------------------------------------------------
# 3. Genes: protein-coding, symbol-keyed, deterministic collapse
# -----------------------------------------------------------------------------
# Keep protein-coding. mtDNA-encoded MT- genes are retained here: script 07
# holds them in their own synthetic pathway and never pools them with
# nuclear-encoded OXPHOS subunits (CLAUDE.md, standing convention).
is_pc <- gene_ann$gene_type == "protein_coding" &
         !is.na(gene_ann$gene_name) & gene_ann$gene_name != ""
message("\n3. protein-coding genes: ", sum(is_pc), " of ", nrow(gene_ann))
counts <- counts_raw[is_pc, , drop = FALSE]
sym    <- gene_ann$gene_name[is_pc]

# Collapse duplicate symbols by highest mean count. Summing would inflate genes
# that happen to carry several ENSG ids; picking the dominant transcript keeps
# the value on the same scale as every other gene.
if (any(duplicated(sym))) {
  mu   <- rowMeans(counts)
  ord  <- order(sym, -mu)
  counts <- counts[ord, , drop = FALSE]
  sym    <- sym[ord]
  first  <- !duplicated(sym)
  message("   ", sum(!first), " duplicate symbol row(s) collapsed (kept highest mean)")
  counts <- counts[first, , drop = FALSE]
  sym    <- sym[first]
}
rownames(counts) <- sym
stopifnot(!any(duplicated(rownames(counts))))

# Minimal low-count filter. Deliberately permissive: an aggressive filter would
# silently shrink the gene universe that G1's overlap audit was computed
# against, and the pathway scores in 07 assume broad coverage.
keep_g  <- rowSums(counts >= 10) >= 10
message("   low-count filter (>=10 counts in >=10 samples): ",
        sum(keep_g), " kept, ", sum(!keep_g), " dropped")
counts  <- counts[keep_g, , drop = FALSE]

message("   final matrix: ", nrow(counts), " genes x ", ncol(counts), " samples")

# -----------------------------------------------------------------------------
# 4. The two matrices - built separately, on purpose
# -----------------------------------------------------------------------------
cd  <- data.frame(patient = colnames(counts), row.names = colnames(counts))
dds <- DESeqDataSetFromMatrix(counts, cd, design = ~ 1)
dds <- estimateSizeFactors(dds)

# --- 4a. LINEAR, for mitoPPS only -------------------------------------------
# DESeq2 size-factor normalised counts. NOT logged. Do not log these.
mat_linear <- counts(dds, normalized = TRUE)
stopifnot(min(mat_linear) >= 0)
saveRDS(list(mat = mat_linear, scale = "linear_deseq2_normalised",
             consumer = "mitoPPS only", sample_map = sample_map,
             built = Sys.time()), PATH_LINEAR)
message("\n4a. LINEAR matrix saved -> ", basename(PATH_LINEAR),
        "  (range ", sprintf("%.1f", min(mat_linear)), " to ",
        sprintf("%.0f", max(mat_linear)), ")")

# --- 4b. VST, log scale, for GSVA -------------------------------------------
# blind = TRUE: no design is being modelled here, and a blind transform keeps
# this object usable for any downstream contrast without circularity.
vsd <- vst(dds, blind = TRUE)
mat_vst <- assay(vsd)
saveRDS(list(mat = mat_vst, scale = "log_vst",
             consumer = "GSVA/ssGSEA, kcdf = Gaussian", sample_map = sample_map,
             built = Sys.time()), PATH_VST)
message("4b. VST matrix saved    -> ", basename(PATH_VST),
        "  (range ", sprintf("%.1f", min(mat_vst)), " to ",
        sprintf("%.1f", max(mat_vst)), ")")

# Assertion that the two really are on different scales. If this ever fails,
# something upstream has collapsed them and every downstream score is suspect.
if (max(mat_vst) > 100 || max(mat_linear) < 100) {
  stop("scale assertion failed: VST max ", round(max(mat_vst), 1),
       ", linear max ", round(max(mat_linear), 1),
       ". The two matrices are not on the scales they claim.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# 5. QC
# -----------------------------------------------------------------------------
qc <- tibble::tibble(
  patient      = colnames(counts),
  library_size = colSums(counts),
  genes_detected = colSums(counts > 0)
)
message("\n5. QC")
message("   library size   median ", format(median(qc$library_size), big.mark = ","),
        "  min ", format(min(qc$library_size), big.mark = ","))
message("   genes detected median ", median(qc$genes_detected))
readr::write_csv(qc, file.path(DIR_TABLES, "tcga_brca_expression_qc.csv"))
readr::write_csv(sample_map,
                 file.path(DIR_TABLES, "tcga_brca_sample_map.csv"))

# -----------------------------------------------------------------------------
# 6. GATE G1, second criterion - correlation structure
# -----------------------------------------------------------------------------
# Plan section 6 fails the Felsher signature on EITHER of two criteria: falling
# below ~50 genes, or LOSING ITS CORRELATION STRUCTURE. The G1 note of
# 2026-08-28 discharged only the first (61 genes, passes). This is the second.
#
# The plan does not operationalise "loses its correlation structure", so it is
# defined here, BEFORE running:
#
#   A signature is only summarisable as a single score if its genes co-express
#   coherently. Coherence is measured as the MEAN PAIRWISE SPEARMAN CORRELATION
#   among signature genes, plus the share of variance on PC1.
#
#   The stripped signature PASSES if:
#     (a) its coherence exceeds the 95th percentile of size-matched random gene
#         sets drawn from the same matrix (i.e. it is a real module, not noise);
#         AND
#     (b) its coherence is not materially below the raw 67-gene signature -
#         pre-specified as retaining at least 80% of the raw mean correlation.
#
#   (a) asks whether it is still a module. (b) asks whether stripping broke it.
#   Both must hold.
#
# >>> THIS CRITERION IS A PRE-REGISTRATION DECISION AND NEEDS SIGN-OFF BEFORE
# >>> THE SCRIPT IS RUN. It is written here so it is fixed before any number is
# >>> seen, per the discipline in CLAUDE.md.

G1_COHERENCE_RETENTION <- 0.80   # criterion (b)
G1_NULL_QUANTILE       <- 0.95   # criterion (a)
G1_NULL_NSETS          <- 1000L

g1 <- readRDS(file.path(DIR_RESULTS, "g1_overlap_audit.rds"))
stopifnot(all(c("estimators_raw", "estimators_stripped") %in% names(g1)))
raw_genes   <- g1$estimators_raw$FELSHER
strip_genes <- g1$estimators_stripped$FELSHER
stopifnot(length(raw_genes) == 67L, length(strip_genes) == 61L)

# Rank once. Spearman on the original values is Pearson on the ranks, so ranking
# up front turns the 1,000-iteration null from a few minutes into a few seconds.
rank_mat <- t(apply(mat_vst, 1L, rank))
dimnames(rank_mat) <- dimnames(mat_vst)

.mean_rho <- function(genes, rm) {
  g <- intersect(genes, rownames(rm))
  if (length(g) < 5L) return(NA_real_)
  rho <- suppressWarnings(cor(t(rm[g, , drop = FALSE])))
  mean(rho[upper.tri(rho)], na.rm = TRUE)
}

.coherence <- function(genes, m, rm) {
  g <- intersect(genes, rownames(m))
  if (length(g) < 5L) return(c(n = length(g), mean_rho = NA_real_, pc1 = NA_real_))
  sub <- m[g, , drop = FALSE]
  # Drop zero-variance genes; prcomp(scale. = TRUE) errors on them.
  v   <- apply(sub, 1L, stats::var)
  sub <- sub[v > 0, , drop = FALSE]
  pcv <- stats::prcomp(t(sub), scale. = TRUE)$sdev^2
  c(n = length(g), mean_rho = .mean_rho(genes, rm), pc1 = pcv[1] / sum(pcv))
}

coh_raw   <- .coherence(raw_genes,   mat_vst, rank_mat)
coh_strip <- .coherence(strip_genes, mat_vst, rank_mat)

# Null: size-matched random gene sets from the same matrix.
set.seed(PROJECT_SEED)
n_strip  <- as.integer(coh_strip[["n"]])
null_rho <- vapply(seq_len(G1_NULL_NSETS), function(i) {
  .mean_rho(sample(rownames(rank_mat), n_strip), rank_mat)
}, numeric(1))
null_cut <- unname(quantile(null_rho, G1_NULL_QUANTILE, na.rm = TRUE))

pass_a <- coh_strip[["mean_rho"]] > null_cut
pass_b <- coh_strip[["mean_rho"]] >= G1_COHERENCE_RETENTION * coh_raw[["mean_rho"]]

message("\n6. GATE G1, correlation criterion")
message(sprintf("   raw       n=%3d  mean rho %.3f  PC1 %.1f%%",
                coh_raw[["n"]], coh_raw[["mean_rho"]], 100 * coh_raw[["pc1"]]))
message(sprintf("   stripped  n=%3d  mean rho %.3f  PC1 %.1f%%",
                coh_strip[["n"]], coh_strip[["mean_rho"]], 100 * coh_strip[["pc1"]]))
message(sprintf("   random-set null, %.0fth pct: %.3f",
                100 * G1_NULL_QUANTILE, null_cut))
message(sprintf("   (a) still a module      : %s", if (pass_a) "PASS" else "FAIL"))
message(sprintf("   (b) retains >=%.0f%% of raw: %s  (%.1f%% retained)",
                100 * G1_COHERENCE_RETENTION, if (pass_b) "PASS" else "FAIL",
                100 * coh_strip[["mean_rho"]] / coh_raw[["mean_rho"]]))
message("   VERDICT: G1 correlation criterion ",
        if (pass_a && pass_b) "PASSES - G1 is now fully discharged"
        else "FAILS - M-a cannot be the primary estimator; revisit D2")

saveRDS(list(raw = coh_raw, stripped = coh_strip, null_rho = null_rho,
             null_cut = null_cut, pass_a = pass_a, pass_b = pass_b,
             criterion = list(retention = G1_COHERENCE_RETENTION,
                              quantile  = G1_NULL_QUANTILE,
                              nsets     = G1_NULL_NSETS),
             built = Sys.time()), PATH_G1CORR)

message("\n01: done.")
message("    results/tcga_brca_vst.rds     LOG scale, for GSVA")
message("    results/tcga_brca_linear.rds  LINEAR scale, for mitoPPS")
message("    results/g1_correlation_criterion.rds")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  # --- the scale check, done by eye ----------------------------------------
  # These must look completely different. If they look similar, something has
  # gone wrong and nothing downstream is trustworthy.
  vst <- readRDS(PATH_VST);  lin <- readRDS(PATH_LINEAR)
  vst$scale; lin$scale
  summary(as.vector(vst$mat[1:200, 1:20]))
  summary(as.vector(lin$mat[1:200, 1:20]))

  # --- G1 coherence, with the null for context ------------------------------
  g1c <- readRDS(PATH_G1CORR)
  hist(g1c$null_rho, breaks = 40,
       main = "size-matched random sets", xlab = "mean pairwise Spearman rho")
  abline(v = g1c$stripped[["mean_rho"]], col = "red",   lwd = 2)
  abline(v = g1c$raw[["mean_rho"]],      col = "blue",  lwd = 2)
  abline(v = g1c$null_cut,               col = "grey40", lty = 2)

  # --- which genes carry the module? ---------------------------------------
  # If coherence is driven by a handful of genes, the score is really those
  # genes. Worth knowing before it becomes Panel a.
  g  <- intersect(readRDS(file.path(DIR_RESULTS, "g1_overlap_audit.rds"))$
                    estimators_stripped$FELSHER, rownames(vst$mat))
  rho <- cor(t(vst$mat[g, ]), method = "spearman")
  sort(rowMeans(rho) , decreasing = TRUE)[1:15]

  # --- sanity: are the priming genes present and behaved? -------------------
  vst$mat[intersect(c("BBC3","BCL2L1","MYC","MCL1"), rownames(vst$mat)), 1:6]

  # --- library-size outliers ------------------------------------------------
  qc <- readr::read_csv(file.path(DIR_TABLES, "tcga_brca_expression_qc.csv"),
                        show_col_types = FALSE)
  qc %>% dplyr::arrange(library_size) %>% head(10) %>% print()

}
