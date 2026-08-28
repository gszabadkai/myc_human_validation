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
  "GSVA", "DESeq2", "limma", "decoupleR", "msigdbr", "survival"
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
