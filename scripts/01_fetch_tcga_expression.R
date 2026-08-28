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

.check_packages("TCGAbiolinks", "TCGA download")
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
})

message("\n01: TCGA-BRCA expression\n", strrep("=", 78))

DIR_TCGA_RAW <- file.path(DIR_DATA, "raw", "tcga_brca_expression")
.ensure_dir(DIR_TCGA_RAW)

PATH_SE      <- file.path(DIR_TCGA_RAW, "tcga_brca_star_counts_se.rds")
PATH_VST     <- file.path(DIR_RESULTS, "tcga_brca_vst.rds")
PATH_LINEAR  <- file.path(DIR_RESULTS, "tcga_brca_linear.rds")
PATH_G1CORR  <- file.path(DIR_RESULTS, "g1_correlation_criterion.rds")

# -----------------------------------------------------------------------------
# 1. Fetch (cached)
# -----------------------------------------------------------------------------
# GDCdownload writes into a GDCdata/ directory under the working directory.
# It is pointed at data/raw/ explicitly so nothing lands in the repo root.

# TCGAbiolinks' `directory` argument controls only where the FINAL per-file
# folders are placed. It writes its chunk tarballs, the extracted UUID
# directories and MANIFEST.txt into getwd() regardless. Run from the project
# root, that fills the repo root with hundreds of UUID folders and a multi-
# hundred-MB tarball. So the download is fenced inside its own directory and the
# working directory is restored on exit, including on error.
#
# The download is resumable: GDCdownload skips files already present, so an
# interrupted run continues rather than starting over.
#
# If the API route keeps failing on truncated chunks, the robust alternative is
# the GDC Data Transfer Tool: install gdc-client, then add method = "client".
.fetch_gdc <- function(query, dest, per_chunk, timeout_s) {
  old_wd <- setwd(dest)                    # setwd() returns the PREVIOUS wd
  on.exit(setwd(old_wd), add = TRUE)
  old_to <- options(timeout = timeout_s)
  on.exit(options(old_to), add = TRUE)
  TCGAbiolinks::GDCdownload(query, directory = "GDCdata",
                            files.per.chunk = per_chunk)
  TCGAbiolinks::GDCprepare(query, directory = "GDCdata")
}

# Small chunks on purpose. The default packs ~1 GB into one tarball, and a
# single truncated stream throws the whole chunk away and re-downloads it.
# 50 files is roughly 200 MB, so a failure costs a fifth as much.
GDC_FILES_PER_CHUNK <- 50L

# THE ONE THAT ACTUALLY MATTERS. R's default download timeout is 60 SECONDS,
# which no multi-GB download can meet. It presents as
#   "truncated gzip input: Unknown error: -1"
# followed, after TCGAbiolinks' retry also times out, by
#   "Error in if (ret == 1) ... : argument is of length zero"
# The giveaway that it is a clock and not a size limit: the first run truncated
# at 63 MB and the retry at 140 MB. A size cap fails at the same point every
# time; a wall-clock timeout fails wherever the transfer has reached.
# Set generously - this is a ceiling, not a delay.
GDC_TIMEOUT_SECONDS <- 14400L   # 4 hours

if (file.exists(PATH_SE)) {
  message("1. cached SummarizedExperiment found, skipping download")
  se <- readRDS(PATH_SE)
} else {
  message("1. querying GDC (~5.2 GB over ~1231 files; resumable, runs once)")
  q <- TCGAbiolinks::GDCquery(
    project       = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  se <- .fetch_gdc(q, DIR_TCGA_RAW, GDC_FILES_PER_CHUNK, GDC_TIMEOUT_SECONDS)
  saveRDS(se, PATH_SE)
  message("   cached to ", PATH_SE)
}

stopifnot(inherits(se, "SummarizedExperiment"))
message("   raw object: ", nrow(se), " genes x ", ncol(se), " samples")
message("   assays: ", paste(assayNames(se), collapse = ", "))

# -----------------------------------------------------------------------------
# 2. Samples: primary tumours only, one per patient
# -----------------------------------------------------------------------------
# The SAME rule as G2 (see data/gistic_tcga_brca/README.md): sample type -01,
# one column per patient.  Using a different rule here would silently
# de-align the expression analysis from the copy-number analysis.
#
# Unlike GISTIC, TCGA expression has genuine within-patient replicate aliquots.
# Those need a deterministic tie-break, so the rule is stated rather than left
# to whichever column happens to come first.

bc  <- colnames(se)
pat <- vapply(strsplit(bc, "-"), function(p) paste(p[1:3], collapse = "-"), character(1))
sty <- substr(vapply(strsplit(bc, "-"), `[`, character(1), 4), 1, 2)

message("\n2. sample types present: ",
        paste(sprintf("%s=%d", names(table(sty)), table(sty)), collapse = " "))

keep <- which(sty == "01")
se   <- se[, keep]
pat  <- pat[keep]

# Deterministic tie-break: for a patient with several -01 aliquots, keep the one
# with the highest library size.  Highest depth is the defensible choice and,
# being a property of the data rather than of barcode ordering, it is stable
# across re-runs and across GDC re-releases.
if (any(duplicated(pat))) {
  libsize <- colSums(assay(se, "unstranded"))
  ord     <- order(pat, -libsize)
  se      <- se[, ord]
  pat     <- pat[ord]
  first   <- !duplicated(pat)
  message("   ", sum(!first), " replicate aliquot(s) dropped (kept highest depth)")
  se  <- se[, first]
  pat <- pat[first]
}
colnames(se) <- pat
stopifnot(!any(duplicated(colnames(se))))
message("   primary tumours, one per patient: ", ncol(se))

# -----------------------------------------------------------------------------
# 3. Genes: protein-coding, symbol-keyed, deterministic collapse
# -----------------------------------------------------------------------------
rd <- as.data.frame(rowData(se))
stopifnot(all(c("gene_name", "gene_type") %in% names(rd)))

# "unstranded" is the raw count assay of the STAR-Counts workflow. GDCprepare
# also returns tpm_/fpkm_ assays; those are already normalised and must not be
# fed to DESeq2, so the assay is named explicitly and asserted rather than
# taken positionally.
if (!"unstranded" %in% assayNames(se)) {
  stop("expected an 'unstranded' raw-count assay, found: ",
       paste(assayNames(se), collapse = ", "), call. = FALSE)
}
counts <- assay(se, "unstranded")
mode(counts) <- "integer"

# Keep protein-coding. mtDNA-encoded MT- genes are retained here: script 07
# holds them in their own synthetic pathway and never pools them with
# nuclear-encoded OXPHOS subunits (CLAUDE.md, standing convention).
is_pc <- rd$gene_type == "protein_coding" & !is.na(rd$gene_name) & rd$gene_name != ""
message("\n3. protein-coding genes: ", sum(is_pc), " of ", nrow(rd))
counts <- counts[is_pc, , drop = FALSE]
sym    <- rd$gene_name[is_pc]

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
             consumer = "mitoPPS only", built = Sys.time()), PATH_LINEAR)
message("\n4a. LINEAR matrix saved -> ", basename(PATH_LINEAR),
        "  (range ", sprintf("%.1f", min(mat_linear)), " to ",
        sprintf("%.0f", max(mat_linear)), ")")

# --- 4b. VST, log scale, for GSVA -------------------------------------------
# blind = TRUE: no design is being modelled here, and a blind transform keeps
# this object usable for any downstream contrast without circularity.
vsd <- vst(dds, blind = TRUE)
mat_vst <- assay(vsd)
saveRDS(list(mat = mat_vst, scale = "log_vst",
             consumer = "GSVA/ssGSEA, kcdf = Gaussian", built = Sys.time()),
        PATH_VST)
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
