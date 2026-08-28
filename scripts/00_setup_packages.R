# 00_setup_packages.R
# =============================================================================
# Package loader, path constants and snapshot-shape expectations.
# Sourced by every numbered script. No side effects beyond attaching core
# packages and writing one session_info file per session.
#
# POLICY: this script CHECKS packages, it does not install them. Auto-installing
# behind a source() is how an analysis environment silently drifts between runs.
# A missing package produces one stop() naming everything that is absent.
#
# No renv; packages are installed system-wide. See CLAUDE.md.
# =============================================================================

# -----------------------------------------------------------------------------
# Package tiers
# -----------------------------------------------------------------------------
# core     : attached here, used everywhere
# analysis : checked for presence only, attached by the scripts that need them
#            (GSVA and DESeq2 are slow to attach and most scripts do not use them)
.pkg_core <- c(
  "here", "dplyr", "tibble", "tidyr", "readr", "readxl",
  "stringr", "purrr", "ggplot2"
)

.pkg_analysis <- c(
  "GSVA", "DESeq2", "limma", "decoupleR", "msigdbr", "survival", "data.table"
)

# Needed by scripts 01 and 02 only. Recorded, not fatal here: gates G1 and G2
# do not touch TCGA expression or genomics.
.pkg_deferred <- c("TCGAbiolinks")

.check_packages <- function(pkgs, tier) {
  have <- vapply(pkgs, function(p) requireNamespace(p, quietly = TRUE), logical(1))
  if (any(!have)) {
    stop(
      "Missing ", tier, " packages: ",
      paste(pkgs[!have], collapse = ", "),
      "\nInstall them, then re-source this script. Nothing is auto-installed.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.check_packages(.pkg_core, "core")
.check_packages(.pkg_analysis, "analysis")

suppressPackageStartupMessages(
  invisible(lapply(.pkg_core, library, character.only = TRUE))
)

# Deferred packages: warn once, never stop.
local({
  have <- vapply(.pkg_deferred, function(p) requireNamespace(p, quietly = TRUE),
                 logical(1))
  if (any(!have)) {
    message("00: deferred package(s) not installed, needed by scripts 01/02 only: ",
            paste(.pkg_deferred[!have], collapse = ", "))
  }
})

# -----------------------------------------------------------------------------
# Options
# -----------------------------------------------------------------------------
options(stringsAsFactors = FALSE)

# Seed constant. Deliberately NOT applied globally here: a seed set at source()
# time gives false reproducibility, because the result then depends on how many
# random draws happened earlier in the session. Scripts that use randomness call
# set.seed(PROJECT_SEED) immediately before the draw.
PROJECT_SEED <- 20260827L

# NOTE on OmnipathR: attaching it writes an omnipathr-log/ directory into the
# working directory. This project does not call it for retrieval - CollecTRI is
# consumed from the dated snapshot in data/collectri_human/ - so it is not
# attached here. omnipathr-log/ is gitignored in case a session attaches it.
# See data/collectri_human/README.md for why the snapshot exists.

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
# All relative to the project root via here::here(). No absolute paths.
DIR_DATA    <- here::here("data")
DIR_RESULTS <- here::here("results")
DIR_OUTPUTS <- here::here("outputs")
DIR_TABLES  <- here::here("outputs", "tables")

PATH_MITOCARTA <- here::here("data", "mitocarta_human", "Human.MitoCarta3.0.xls")
PATH_COLLECTRI <- here::here("data", "collectri_human",
                             "collectri_human_omnipath.tsv.gz")
PATH_FELSHER   <- here::here("data", "genesets_from_library_human",
                             "felsher_integrative_signature.csv")

# G2 copy-number inputs. The two GISTIC files are large and live under
# data/raw/, which is gitignored and NOT on origin. Provenance, URLs and
# SHA-256 sums are in data/gistic_tcga_brca/README.md; re-download from there
# rather than hunting for a backup.
DIR_GISTIC_RAW <- here::here("data", "raw", "gistic_tcga_brca")

PATH_GISTIC_ISAR <- file.path(DIR_GISTIC_RAW,
                              "ISAR_GISTIC.all_thresholded.by_genes.txt.gz")
PATH_GISTIC_FH_TAR <- file.path(
  DIR_GISTIC_RAW,
  "gdac.broadinstitute.org_BRCA-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0.tar.gz"
)
# Extracted from the tarball on first use by script 05.
PATH_GISTIC_FH_TXT <- file.path(DIR_GISTIC_RAW, "firehose",
                                "all_thresholded.by_genes.txt")

PATH_GDC_BRCA_CASES <- here::here("data", "tcga_clinical",
                                  "gdc_brca_cases_2026-08-28.tsv")
PATH_TCGA_CLINICAL  <- here::here("data", "tcga_clinical",
                                  "tcga_brca_clinical_snapshot_2026-08-28.tsv")

.ensure_dir <- function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  invisible(p)
}

invisible(lapply(c(DIR_RESULTS, DIR_OUTPUTS, DIR_TABLES), .ensure_dir))

# -----------------------------------------------------------------------------
# Snapshot shape expectations
# -----------------------------------------------------------------------------
# Transcribed from the provenance READMEs. Scripts assert against these so that
# a file swapped underneath the pipeline fails loudly instead of quietly
# changing a result. If one of these legitimately changes, re-snapshot and
# update BOTH the README and this block.
EXPECT_MITOCARTA_GENES        <- 1136L   # Sheet 2, "A Human MitoCarta3.0"
EXPECT_MITOCARTA_BACKGRD      <- 19247L  # Sheet 3, "B Human All Genes"
# Sheet 4 ends with 5 entirely blank padding rows. Both numbers are asserted:
# the raw row count catches a swapped file, the pathway count is the real one.
EXPECT_MITOCARTA_PATHWAY_ROWS <- 154L    # Sheet 4 raw, incl. 5 blank rows
EXPECT_MITOCARTA_PATHWAYS     <- 149L    # Sheet 4 after dropping the blanks
EXPECT_MITOCARTA_MTDNA    <- 13L     # MT- prefixed protein-coding
EXPECT_FELSHER_GENES      <- 67L
EXPECT_COLLECTRI_ROWS     <- 64723L
EXPECT_COLLECTRI_TFS      <- 1201L

# G1 decision threshold, plan section 6: below this the Felsher signature cannot
# serve as the primary MYC estimator (D2) and CollecTRI + the 8q24 CNV
# instrument carry the MYC axis instead.
G1_MIN_SIGNATURE_SIZE <- 50L

# --- G2 shape expectations ---------------------------------------------------
EXPECT_GISTIC_ISAR_GENES   <- 24203L
EXPECT_GISTIC_ISAR_SAMPLES <- 9991L   # pan-cancer, all 33 TCGA types
EXPECT_GISTIC_FH_GENES     <- 24776L
EXPECT_GISTIC_FH_SAMPLES   <- 1080L   # BRCA primary tumours only
EXPECT_GDC_BRCA_CASES      <- 1098L
# After the sample rule (GDC case list, then sample type -01, one per patient).
EXPECT_BRCA_N_ISAR         <- 1043L
EXPECT_BRCA_N_FH           <- 1080L
# TNBC rule, see data/tcga_clinical/README.md.
EXPECT_TNBC_CALLABLE       <- 951L
EXPECT_TNBC_N              <- 161L

# --- G2 pre-registered thresholds -------------------------------------------
# Fixed 2026-08-28, see docs/2026-08-28_D4_scoping_and_G2_design.md section 2.3.
# Gene-specific AND direction-specific. BBC3 is the pro-apoptotic PRIME
# numerator, so its PRIME-lowering event is LOSS, not gain; MYC/MCL1/BCL2L1 are
# tested for gain. These are not the same test and must never be described
# collectively as "co-amplification".
#
# BCL2L1 sits at >= +1 rather than +2 because at +2 it is n=28 and
# uninterpretable, and it is this arm's a priori focus (the PRIME denominator,
# and the term Lee et al. 2017 never touch).
# BBC3 sits at <= -1 because homozygous deletion is n=5 and its joint count with
# MYC amplification is ZERO - that test is empty, not merely underpowered.
# BAX is the 19q13 REGIONAL CONTROL for BBC3, not a hypothesis.
G2_THRESHOLDS <- tibble::tribble(
  ~gene,     ~direction, ~primary_rule, ~role,
  "MYC",     "gain",     "eq2",         "exposure",
  "MCL1",    "gain",     "eq2",         "partner_replication",
  "BCL2L1",  "gain",     "ge1",         "partner_novel",
  "BBC3",    "loss",     "le_neg1",     "partner_novel",
  "BAX",     "loss",     "le_neg1",     "regional_control"
)

# Secondary grid: every gene is additionally reported at the OTHER threshold in
# its direction's pair, so the sensitivity of the result to that choice is
# visible. The pairing is keyed on the gene's PRIMARY rule, not on direction --
# keying it on direction alone silently maps BCL2L1 (primary "ge1") back onto
# "ge1" and the grid then reports the primary twice.
G2_ALT_RULE <- c(
  eq2     = "ge1",       # gain: high-level <-> gain-or-amplification
  ge1     = "eq2",
  le_neg1 = "le_neg2",   # loss: any loss <-> homozygous deletion
  le_neg2 = "le_neg1"
)

# GUARD: odds ratios computed at different thresholds are NOT comparable across
# genes. Compare within a threshold only. See design note section 2.3.

# G2 pass criterion, fixed before any statistic was computed (design note 2.8):
# a partner passes if the ANEUPLOIDY-ADJUSTED odds ratio exceeds 1 with a 95% CI
# excluding 1, at its primary threshold and direction, in source A (ISAR).
# Source B must agree in direction; a B disagreement is reported, never
# overridden.
G2_ALPHA <- 0.05

# -----------------------------------------------------------------------------
# Session record
# -----------------------------------------------------------------------------
# One file per R session, not per source(). Cheap provenance for a pipeline that
# is hand-run interactively rather than executed end to end.
local({
  stamp <- format(Sys.time(), "%Y%m%d")
  f <- file.path(DIR_RESULTS, paste0("session_info_", stamp, ".rds"))
  if (!file.exists(f)) {
    saveRDS(utils::sessionInfo(), f)
    message("00: session info written to results/", basename(f))
  }
})

message("00: setup complete (", length(.pkg_core), " core packages attached)")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  # Confirm every snapshot the gates need is present and readable.
  file.exists(PATH_MITOCARTA)
  file.exists(PATH_COLLECTRI)
  file.exists(PATH_FELSHER)

  # Package versions actually in play.
  vapply(c(.pkg_core, .pkg_analysis),
         function(p) as.character(utils::packageVersion(p)), character(1))

  # Paths resolve where you expect.
  DIR_RESULTS
  DIR_TABLES

}
