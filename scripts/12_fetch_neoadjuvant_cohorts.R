# 12_fetch_neoadjuvant_cohorts.R
# =============================================================================
# Fetch, parse and harmonise the three H4 neoadjuvant cohorts. NOTHING IS
# SCORED HERE AND NO ASSOCIATION IS COMPUTED.
#
# Built to:
#   docs/2026-08-28_D5_cohort_selection.md            - which cohorts, and why
#   docs/2026-08-30_STATE_frozen_and_H4_buffer_declaration.md - what H4 needs
#   docs/2026-08-27_human_validation_plan.md sections 3, 11 (script 12)
#
# =============================================================================
# THE DISCIPLINE LINE, WHICH THIS SCRIPT SITS EXACTLY ON
# =============================================================================
# This script reads the pCR column. That is unavoidable - it has to be parsed
# and harmonised before it can be modelled. What it must NOT do, and does not
# do, is relate any score to it. D5 section 2 held that line during cohort
# selection ("no score-versus-pCR association was computed in any candidate
# cohort") and it carries here for the same reason: once an association has been
# seen, no declaration can honestly be amended.
#
# So this script computes marginals only - n, endpoint coding, subtype
# composition, gene coverage, scale. Every pCR count printed below is already
# published in the D5 note (319 of 988 for I-SPY2), so nothing new about the
# outcome is learned by running it.
#
# NO GSVA. NO mitoPPS. NO MYC score. NO STATE. Those belong in script 13, where
# each cohort is scored in its own single cohort-relative run.
#
# =============================================================================
# SCALE - READ THIS BEFORE WRITING SCRIPT 13
# =============================================================================
# ALL THREE COHORTS ARE DEPOSITED ON THE LOG2 SCALE. Verified from the files:
#
#   GSE194040  Agilent gene-level, values ~4-12, ComBat-adjusted (D5 section 4)
#   GSE164458  "RNAseq_log2_Processed", values ~1-15
#   GSE25066   GPL96 series matrix, values ~5-13
#
# CONSEQUENCE, AND IT IS NOT SMALL. GSVA wants log scale and is satisfied.
# mitoPPS wants LINEAR DESeq2-normalised counts (CLAUDE.md), and no such matrix
# exists for any of these three - two are microarrays and the third is deposited
# already logged. Exponentiating recovers a linear scale arithmetically but not
# the quantity mitoPPS was defined on: mitoPPS takes ratios of pathway MEANS
# across genes, so it assumes cross-gene comparability that probe affinity does
# not provide, and I-SPY2's matrix has additionally been ComBat-adjusted on the
# log scale, so 2^x is a batch-corrected intensity rather than an abundance.
#
# This script MEASURES the facts and does NOT make the decision. Section 6
# reports them. The two-instrument rule declared in
# docs/2026-08-30_STATE_frozen_and_H4_buffer_declaration.md section 6.5 has to
# be resolved by a recorded amendment BEFORE script 13, and the reason must be
# stated as FEASIBILITY - the instrument does not exist in these cohorts - not
# as anything about a result. That amendment is legitimate only while it is
# written from this section's output and nothing else.
#
# =============================================================================
# COHORT-RELATIVITY: these three are scored SEPARATELY in script 13, never
# pooled, and their effect estimates are meta-analysed (D5 section 6.2,
# CLAUDE.md). This script therefore keeps them in three separate objects and
# provides no combined matrix, deliberately.
#
# SPECIES: human. Human gene symbols throughout. See CLAUDE.md.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

message("\n12: fetch the three H4 neoadjuvant cohorts\n", strrep("=", 78))

# -----------------------------------------------------------------------------
# 0. Constants
# -----------------------------------------------------------------------------
DIR_NEO <- here::here("data", "raw", "neoadjuvant")

GEO_BASE <- "https://ftp.ncbi.nlm.nih.gov/geo"

# name = c(url-relative-path, minimum bytes). The GEO gateway serves neither
# HEAD nor Content-Range, so an exact byte guard is not obtainable the way
# script 02's is. A FLOOR plus the parsed-shape assertions in sections 2-4 is
# the substitute: the floor catches a truncated or error-page download, and the
# dimension checks catch a wrong file. Floors are ~90% of the sizes shown in the
# GEO directory listing on 2026-08-30.
INPUTS <- list(
  ispy2_expr = c(
    "series/GSE194nnn/GSE194040/suppl/GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt.gz",
    42e6),
  ispy2_sm_a = c(
    "series/GSE194nnn/GSE194040/matrix/GSE194040-GPL20078_series_matrix.txt.gz",
    40e3),
  ispy2_sm_b = c(
    "series/GSE194nnn/GSE194040/matrix/GSE194040-GPL30493_series_matrix.txt.gz",
    25e3),
  brightness_expr = c(
    "series/GSE164nnn/GSE164458/suppl/GSE164458_BrighTNess_RNAseq_log2_Processed_ASTOR.txt.gz",
    65e6),
  brightness_sm = c(
    "series/GSE164nnn/GSE164458/matrix/GSE164458_series_matrix.txt.gz",
    12e3),
  hatzis_sm = c(
    "series/GSE25nnn/GSE25066/matrix/GSE25066_series_matrix.txt.gz",
    53e6),
  gpl96_annot = c(
    "platforms/GPLnnn/GPL96/annot/GPL96.annot.gz",
    3.8e6)
)

# Expected shapes, from D5 and from the GEO metadata read on 2026-08-30.
EXPECT_ISPY2_N       <- 988L
EXPECT_ISPY2_PCR     <- 319L    # D5 section 4; published, not new information
EXPECT_BRIGHTNESS_N  <- 482L
EXPECT_HATZIS_N      <- 508L
EXPECT_GPL96_PROBES  <- 22283L

# Minimum overlap between an expression matrix's column names and a phenotype
# identifier column before that column is accepted as the join key.
MIN_ID_MATCH <- 0.90

# Genes H4 needs to find in every cohort. Not the full sets - those are checked
# by coverage in section 7 - but the ones whose absence would silently remove a
# term from the declared model.
CRITICAL_GENES <- c("MYC", "MCL1", "BCL2L1", "BBC3", "BCL2L11")

# GSE25066 carries four EXPRESSION-DERIVED PUBLISHED PREDICTORS. They are
# recorded because they are in the file, and they are FORBIDDEN as covariates or
# comparators: every one is a function of the same matrix H4's exposure is built
# from, so adjusting for them or benchmarking against them is circular.
HATZIS_FORBIDDEN <- c("set_class", "chemosensitivity_prediction",
                      "dlda30_prediction", "rcb_0_i_prediction")

# =============================================================================
# 1. Inputs - the author downloads these; this script does not fetch
# =============================================================================
# Same contract as script 02: assert and stop with the exact command. Nothing is
# downloaded behind a source().
message("\n1. inputs")

.ensure_dir(DIR_NEO)

.local_path <- function(key) file.path(DIR_NEO, basename(INPUTS[[key]][1]))

missing <- character(0)
short   <- character(0)
for (k in names(INPUTS)) {
  p <- .local_path(k)
  if (!file.exists(p)) {
    missing <- c(missing, k)
  } else if (file.size(p) < as.numeric(INPUTS[[k]][2])) {
    short <- c(short, sprintf("%s (%s bytes, expected at least %s)",
                              basename(p), format(file.size(p), big.mark = ","),
                              format(as.numeric(INPUTS[[k]][2]), big.mark = ",")))
  }
}
if (length(missing) || length(short)) {
  cmds <- vapply(missing, function(k)
    sprintf("  curl -L -o %s \\\n    %s/%s",
            .local_path(k), GEO_BASE, INPUTS[[k]][1]), character(1))
  stop("neoadjuvant inputs are not ready.\n",
       if (length(missing)) paste0("MISSING (", length(missing), "):\n",
                                   paste(cmds, collapse = "\n"), "\n") else "",
       if (length(short)) paste0("TRUNCATED:\n  ",
                                 paste(short, collapse = "\n  "), "\n",
                                 "Delete and re-download these.\n") else "",
       "Provenance and the full list are in data/neoadjuvant/README.md.",
       call. = FALSE)
}
message("   all ", length(INPUTS), " inputs present and above the size floor")

# -----------------------------------------------------------------------------
# Series-matrix helpers
# -----------------------------------------------------------------------------
# A GEO series matrix carries enormous metadata LINES - one field per sample,
# all on one line. Measured in these files:
#
#   GSE194040-GPL20078  !Sample_data_processing etc.   896,657 bytes
#   GSE25066            !Sample_data_processing        185,443 bytes
#   GSE164458           longest line                   124,861 bytes
#
# vroom's connection buffer defaults to 131,072, so GSE25066 stops readr with
# "The size of the connection buffer was not large enough to fit a complete
# line" - which reads like a corrupt download and is not one. GSE164458 passes
# only because it lands 6,211 bytes under the default; that is luck, not
# margin, so the buffer is raised for every read here rather than for the one
# file that happened to fail.
#
# Raised per call and RESTORED, so sourcing this script leaves no global
# environment change behind for whatever runs next in the session.
VROOM_BUF <- 8388608L    # 8 MB, about 9x the largest line seen

.vroom_big <- function(expr) {
  old <- Sys.getenv("VROOM_CONNECTION_SIZE", unset = NA_character_)
  Sys.setenv(VROOM_CONNECTION_SIZE = VROOM_BUF)
  on.exit({
    if (is.na(old)) Sys.unsetenv("VROOM_CONNECTION_SIZE")
    else Sys.setenv(VROOM_CONNECTION_SIZE = old)
  }, add = TRUE)
  force(expr)
}
# A GEO series matrix is a "!"-prefixed metadata block, then
# !series_matrix_table_begin, the table, then !series_matrix_table_end. The
# metadata block is read line by line so a 59 MB expression table is never
# pulled into memory as character.
.sm_header <- function(path, max_lines = 800L) {
  con <- gzfile(path, "rt"); on.exit(close(con))
  out <- character(0)
  repeat {
    l <- readLines(con, n = 1L, warn = FALSE)
    if (!length(l)) break
    out <- c(out, l)
    if (startsWith(l, "!series_matrix_table_begin")) break
    if (length(out) >= max_lines) {
      stop(basename(path), ": no !series_matrix_table_begin in the first ",
           max_lines, " lines. Not a series matrix, or truncated.",
           call. = FALSE)
    }
  }
  out
}

.sm_field <- function(hdr, tag) {
  l <- hdr[startsWith(hdr, tag)]
  if (!length(l)) return(NULL)
  lapply(l, function(x) {
    v <- strsplit(x, "\t", fixed = TRUE)[[1]][-1]
    gsub('^"|"$', "", v)
  })
}

# Characteristic lines are "key: value", one value per sample. THEY ARE NOT
# ALIGNED ACROSS SAMPLES, and assuming they are corrupts the endpoint.
#
# Measured in GSE25066: the series is ragged from the 7th characteristic line
# onward, because its 310- and 198-sample blocks were deposited with different
# field orders. Line 12 carries `pathologic_response_pcr_rd` for 310 samples and
# `pathologic_response_rcb_class` for the other 198. A parser that names the
# line by its majority key - which this one used to do - hands back RCB classes
# labelled as pCR for 198 of 508 patients. Nothing about that is visible
# downstream; it simply changes who is a responder.
#
# So the key is read PER SAMPLE and the table is pivoted on it. A sample that
# never carried a given key gets NA for it, which is the truth.
.sm_characteristics <- function(hdr, n_samples) {
  ch <- .sm_field(hdr, "!Sample_characteristics_ch1")
  if (is.null(ch)) return(tibble::tibble(.rows = n_samples))

  long <- dplyr::bind_rows(lapply(ch, function(v) {
    kv <- grepl(":", v, fixed = TRUE)
    tibble::tibble(
      sample = seq_along(v),
      key    = ifelse(kv, trimws(sub(":.*$", "", v)), NA_character_),
      value  = ifelse(kv, trimws(sub("^[^:]*:[ ]?", "", v)), NA_character_))
  }))
  long <- long[!is.na(long$key) & nzchar(long$key), , drop = FALSE]
  long$value[long$value %in% c("", "NA", "na", "N/A", "n/a", "null")] <-
    NA_character_

  # A sample may legitimately carry the same key twice; number the repeats
  # rather than letting pivot_wider silently make a list column.
  long <- long %>%
    dplyr::group_by(.data$sample, .data$key) %>%
    dplyr::mutate(key = if (dplyr::n() > 1L)
                          paste0(.data$key, "_", dplyr::row_number())
                        else .data$key) %>%
    dplyr::ungroup()

  wide <- tidyr::pivot_wider(long, names_from = "key", values_from = "value")
  out <- dplyr::left_join(tibble::tibble(sample = seq_len(n_samples)), wide,
                          by = "sample")
  dplyr::select(out, -"sample")
}

.sm_pheno <- function(path) {
  hdr <- .sm_header(path)
  gsm <- .sm_field(hdr, "!Sample_geo_accession")
  ttl <- .sm_field(hdr, "!Sample_title")
  if (is.null(gsm) || is.null(ttl)) {
    stop(basename(path), ": missing !Sample_geo_accession or !Sample_title.",
         call. = FALSE)
  }
  gsm <- gsm[[1]]; ttl <- ttl[[1]]
  stopifnot(length(gsm) == length(ttl))
  ch <- .sm_characteristics(hdr, length(gsm))
  if (nrow(ch) && nrow(ch) != length(gsm)) {
    stop(basename(path), ": ", nrow(ch), " characteristic rows against ",
         length(gsm), " samples.", call. = FALSE)
  }
  dplyr::bind_cols(tibble::tibble(gsm = gsm, title = ttl), ch)
}

.sm_table <- function(path) {
  hdr <- .sm_header(path)
  .vroom_big(readr::read_tsv(path, skip = length(hdr), comment = "!",
                             show_col_types = FALSE, progress = FALSE))
}

# The join key is DISCOVERED, not assumed, because none of these three cohorts
# joins on the obvious column. Measured 2026-08-30 against the deposited files:
#
#   GSE194040  expression columns are I-SPY2 research IDs. `patient id` matches
#              653 of GPL20078's 654 - ONE sample disagrees with its own title -
#              while `!Sample_title` minus the "ISPY2_" prefix matches 654/654.
#   GSE164458  `!Sample_title` is "102001_RNAseq" against an expression column
#              "102001". Raw overlap is ZERO; suffix-stripped it is 482/482.
#
# Assuming a key here is how 988 samples silently become 40. Every candidate
# column is scored under each declared transformation and the best is used, but
# only if it clears MIN_ID_MATCH - otherwise the script stops and shows its
# working rather than joining badly.
.ID_TRANSFORMS <- list(
  raw          = function(x) x,
  strip_prefix = function(x) sub("^[A-Za-z][A-Za-z0-9]*_", "", x),
  strip_suffix = function(x) sub("_[A-Za-z][A-Za-z0-9]*$", "", x))

.match_ids <- function(expr_ids, pheno, label) {
  grid <- expand.grid(col = names(pheno), tf = names(.ID_TRANSFORMS),
                      stringsAsFactors = FALSE)
  grid$score <- vapply(seq_len(nrow(grid)), function(i) {
    v <- .ID_TRANSFORMS[[grid$tf[i]]](trimws(as.character(pheno[[grid$col[i]]])))
    mean(expr_ids %in% v)
  }, numeric(1))
  grid <- grid[order(-grid$score), ]
  top <- utils::head(grid, 4)
  msg <- paste(sprintf("%s/%s %.3f", top$col, top$tf, top$score),
               collapse = "  |  ")
  if (grid$score[1] < MIN_ID_MATCH) {
    stop(label, ": no phenotype column matches the expression column names ",
         "at >= ", MIN_ID_MATCH, ", under any declared transformation.",
         "\nBest candidates (column/transform score): ", msg,
         "\nThe join key has changed upstream; inspect before proceeding.",
         call. = FALSE)
  }
  message("     join key: '", grid$col[1], "' under '", grid$tf[1], "' matches ",
          sprintf("%.1f%%", 100 * grid$score[1]), " of expression columns")
  list(col = grid$col[1], tf = grid$tf[1], score = grid$score[1])
}

# Applies what .match_ids found, and asserts the result lines up exactly.
.align_pheno <- function(pheno, expr_ids, m, label) {
  pheno$sample_id <- .ID_TRANSFORMS[[m$tf]](trimws(as.character(pheno[[m$col]])))
  if (anyDuplicated(pheno$sample_id)) {
    stop(label, ": join key '", m$col, "' is not unique after '", m$tf, "'.",
         call. = FALSE)
  }
  out <- pheno[match(expr_ids, pheno$sample_id), ]
  if (!identical(out$sample_id, expr_ids)) {
    stop(label, ": ", sum(is.na(out$sample_id)), " expression column(s) have no",
         " phenotype row after joining on '", m$col, "'/'", m$tf, "'.",
         call. = FALSE)
  }
  out
}

# A matrix is accepted as log2 only if it looks like one. This is a guard, not a
# transformation: nothing here rescales anything.
# A file written by R with row.names = TRUE has ONE FEWER header field than it
# has data fields, and readr does not treat that as row names: it reads 988
# names against 989 fields and MERGES THE OVERFLOW INTO THE LAST COLUMN, so the
# final sample silently becomes the character string "8.427\t8.8056" - two
# patients' values glued together - while every other column parses cleanly.
# read_tsv warns, but the warning says "parsing issues", not "your last sample
# is now two samples".
#
# So the header is read separately and the body is given explicit names. The
# ragged case is detected rather than guessed at: if the header is one short,
# the first data field is row names.
.read_matrix_tsv <- function(path, label) {
  con <- gzfile(path, "rt")
  hdr <- strsplit(readLines(con, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1]]
  close(con)
  hdr <- gsub('^"|"$', "", hdr)

  con <- gzfile(path, "rt")
  n_field <- length(strsplit(readLines(con, n = 2L, warn = FALSE)[2], "\t",
                             fixed = TRUE)[[1]])
  close(con)

  if (n_field == length(hdr) + 1L) {
    message("     row-named layout: ", length(hdr), " header fields, ",
            n_field, " data fields")
    nm <- c("gene", hdr)
  } else if (n_field == length(hdr)) {
    nm <- c("gene", hdr[-1])
  } else {
    stop(label, ": header has ", length(hdr), " fields and the first data row ",
         n_field, ". Neither a row-named nor a plain layout.", call. = FALSE)
  }
  out <- .vroom_big(readr::read_tsv(path, skip = 1L, col_names = nm,
                                    show_col_types = FALSE, progress = FALSE))
  if (ncol(out) != length(nm)) {
    stop(label, ": read ", ncol(out), " columns, expected ", length(nm), ".",
         call. = FALSE)
  }
  out
}

# as.matrix() on a data frame containing one stray character column silently
# returns a CHARACTER matrix, and storage.mode(x) <- "double" then fills it with
# NAs behind a warning. Checked rather than hoped for.
.numeric_matrix <- function(df, idcol, label) {
  ok <- vapply(df[, -1, drop = FALSE], is.numeric, logical(1))
  if (!all(ok)) {
    stop(label, ": ", sum(!ok), " non-numeric expression column(s), first is '",
         names(ok)[!ok][1], "'. The file layout has changed.", call. = FALSE)
  }
  M <- as.matrix(df[, -1, drop = FALSE])
  rownames(M) <- df[[idcol]]
  M
}

.assert_log2 <- function(M, label) {
  rg <- range(M, na.rm = TRUE)
  if (rg[2] > 40) {
    stop(label, " has max ", signif(rg[2], 4), ", which is not a log2 scale. ",
         "The deposited file has changed; every scale statement in this ",
         "script's header must be re-checked before script 13 runs.",
         call. = FALSE)
  }
  tibble::tibble(cohort = label, min = rg[1], max = rg[2],
                 median = stats::median(M, na.rm = TRUE),
                 frac_negative = mean(M < 0, na.rm = TRUE),
                 frac_na = mean(is.na(M)),
                 n_genes = nrow(M), n_samples = ncol(M))
}

RES_SCALE <- list()

# =============================================================================
# 2. GSE194040 - I-SPY2-990. PRIMARY.
# =============================================================================
# Two platforms (GPL20078 n=654, GPL30493 n=334) with ONE combined gene-level
# expression file at n=988. The phenotype therefore comes from both series
# matrices stacked, and the count is asserted.
#
# The matrix is ComBat-adjusted as deposited (D5 section 4). That is a property
# of the source, not a choice made here, and it must be stated in Methods.
message("\n2. GSE194040 - I-SPY2-990 (primary)")

ispy2_e <- .read_matrix_tsv(.local_path("ispy2_expr"), "GSE194040")
ispy2_m <- .numeric_matrix(ispy2_e, "gene", "GSE194040")
message("   expression: ", nrow(ispy2_m), " genes x ", ncol(ispy2_m), " samples")
if (ncol(ispy2_m) != EXPECT_ISPY2_N) {
  stop("I-SPY2 expression has ", ncol(ispy2_m), " samples, expected ",
       EXPECT_ISPY2_N, ".", call. = FALSE)
}

ispy2_p <- dplyr::bind_rows(.sm_pheno(.local_path("ispy2_sm_a")),
                            .sm_pheno(.local_path("ispy2_sm_b")))
message("   phenotype:  ", nrow(ispy2_p), " samples from 2 series matrices")
stopifnot(nrow(ispy2_p) == EXPECT_ISPY2_N)

ispy2_key <- .match_ids(colnames(ispy2_m), ispy2_p, "GSE194040")
ispy2_p <- .align_pheno(ispy2_p, colnames(ispy2_m), ispy2_key, "GSE194040")

# I-SPY2's `arm` field carries ONE typo, and it matters more than a typo should:
# "Paclitaxel + AMG 386" (n = 114) and "Paclitaxel + AMG-386" (n = 1) are the
# same arm spelled two ways. Left alone that is a singleton factor level, which
# either drops a patient or makes any arm-adjusted design rank-deficient - and
# it would be found, if at all, as a confusing model error rather than as a
# spelling difference.
#
# Normalised explicitly rather than by a silent gsub, and the merge is printed,
# because collapsing two treatment labels is a data decision and not a format
# fix. Only separator variants are merged; nothing else is touched.
# A blanket hyphen-to-space rule would be WRONG here: it turns "T-DM1 +
# Pertuzumab" into "T DM1 + Pertuzumab", and that hyphen is part of the drug's
# name. So only labels that genuinely COLLIDE once punctuation and case are
# ignored get merged, and the more common spelling wins. A label with no
# collision partner is never touched.
.norm_arm <- function(x) {
  key <- gsub("[^a-z0-9]+", "", tolower(x))
  tab <- table(x)
  out <- x
  for (k in unique(key)) {
    lv <- unique(x[key == k])
    if (length(lv) < 2L) next
    keep <- lv[which.max(tab[lv])]
    for (l in setdiff(lv, keep)) {
      message("   arm label merged: '", l, "' (n = ", tab[[l]], ") -> '",
              keep, "' (n = ", tab[[keep]], ")")
      out[x == l] <- keep
    }
  }
  out
}

# Endpoint: `pcr` is already 0/1 in this series.
ispy2_p <- ispy2_p %>%
  dplyr::mutate(
    pcr       = suppressWarnings(as.integer(.data$pcr)),
    hr        = suppressWarnings(as.integer(.data$hr)),
    her2      = suppressWarnings(as.integer(.data$her2)),
    subtype   = dplyr::case_when(
      is.na(hr) | is.na(her2) ~ NA_character_,
      hr == 1 & her2 == 0     ~ "HRpos_HER2neg",
      hr == 0 & her2 == 0     ~ "TNBC",
      hr == 1 & her2 == 1     ~ "HRpos_HER2pos",
      hr == 0 & her2 == 1     ~ "HRneg_HER2pos"),
    treatment = .norm_arm(trimws(as.character(.data$arm))))

message("   pCR: ", sum(ispy2_p$pcr, na.rm = TRUE), " of ",
        sum(!is.na(ispy2_p$pcr)),
        sprintf("  (%.1f%%)", 100 * mean(ispy2_p$pcr, na.rm = TRUE)),
        "   [D5 published ", EXPECT_ISPY2_PCR, "]")
if (sum(ispy2_p$pcr, na.rm = TRUE) != EXPECT_ISPY2_PCR) {
  stop("I-SPY2 pCR count is ", sum(ispy2_p$pcr, na.rm = TRUE), ", not the ",
       EXPECT_ISPY2_PCR, " D5 recorded. The deposited file has changed.",
       call. = FALSE)
}
message("   arms: ", dplyr::n_distinct(ispy2_p$treatment),
        "   |  subtypes: ",
        paste(names(table(ispy2_p$subtype)), collapse = ", "))

RES_SCALE[["GSE194040"]] <- .assert_log2(ispy2_m, "GSE194040")

# =============================================================================
# 3. GSE164458 - BrighTNess. REPLICATION.
# =============================================================================
message("\n3. GSE164458 - BrighTNess (replication)")

bt_e <- .read_matrix_tsv(.local_path("brightness_expr"), "GSE164458")
bt_m <- .numeric_matrix(bt_e, "gene", "GSE164458")
message("   expression: ", nrow(bt_m), " genes x ", ncol(bt_m), " samples")

bt_p <- .sm_pheno(.local_path("brightness_sm"))
message("   phenotype:  ", nrow(bt_p), " samples")

bt_key <- .match_ids(colnames(bt_m), bt_p, "GSE164458")
bt_p <- .align_pheno(bt_p, colnames(bt_m), bt_key, "GSE164458")

# Endpoint: "pCR" / "RD".
.pcr01 <- function(x) {
  v <- toupper(trimws(as.character(x)))
  out <- rep(NA_integer_, length(v))
  out[v %in% c("PCR", "1", "YES")] <- 1L
  out[v %in% c("RD", "0", "NO")]   <- 0L
  bad <- setdiff(unique(v[!is.na(v) & nzchar(v)]),
                 c("PCR", "RD", "0", "1", "YES", "NO"))
  if (length(bad)) {
    stop("unrecognised pCR value(s): ", paste(bad, collapse = ", "),
         call. = FALSE)
  }
  out
}
bt_p$pcr <- .pcr01(bt_p$pathologic_complete_response)
bt_p$treatment <- trimws(as.character(bt_p$planned_arm_code))

# BrighTNess is TNBC-only by design, so there is no receptor-subtype variable to
# adjust for and none is invented. The declaration's PAM50 adjustment is met
# here by the trial's own restriction; section 5 records that explicitly.
bt_p$subtype <- "TNBC"

message("   pCR: ", sum(bt_p$pcr, na.rm = TRUE), " of ", sum(!is.na(bt_p$pcr)),
        sprintf("  (%.1f%%)", 100 * mean(bt_p$pcr, na.rm = TRUE)))
message("   arms: ", paste(sort(unique(bt_p$treatment)), collapse = ", "))

RES_SCALE[["GSE164458"]] <- .assert_log2(bt_m, "GSE164458")

# =============================================================================
# 4. GSE25066 - Hatzis. THIRD, and the only one needing probe collapse.
# =============================================================================
# GPL96 (HG-U133A). The deposited supplementary is a 1.1 GB CEL tarball, so the
# SERIES MATRIX is the input: it carries the normalised table and the phenotype
# in one file.
#
# PROBE COLLAPSE RULE, fixed here before any model:
#   1. drop probes whose GPL96 `Gene symbol` is empty;
#   2. drop MULTI-MAPPING probes - symbol containing "///". This is what removes
#      D5's `211692_s_at`, the single BBC3 probe that also maps to MIR3190 and
#      MIR3191. Dropping it is the honest outcome: BBC3 is then absent from this
#      cohort rather than measured by a probe that is not specific to it;
#   3. where several probes survive for one gene, keep the one with the HIGHEST
#      MEAN across samples (maxMean).
message("\n4. GSE25066 - Hatzis (third)")

hz_p <- .sm_pheno(.local_path("hatzis_sm"))
message("   phenotype:  ", nrow(hz_p), " samples")
stopifnot(nrow(hz_p) == EXPECT_HATZIS_N)

hz_t <- .sm_table(.local_path("hatzis_sm"))
names(hz_t)[1] <- "probe"
hz_probe <- .numeric_matrix(hz_t, "probe", "GSE25066")
message("   expression: ", nrow(hz_probe), " probes x ", ncol(hz_probe),
        " samples")
stopifnot(ncol(hz_probe) == EXPECT_HATZIS_N)

# --- GPL96, and the collapse -------------------------------------------------
# NOT read with .sm_header: a platform annot file has no
# !series_matrix_table_begin, so that helper would exhaust its line budget and
# stop. The marker here is !platform_table_begin.
local({
  con <- gzfile(.local_path("gpl96_annot"), "rt"); on.exit(close(con))
  gpl_hdr <<- readLines(con, n = 200L, warn = FALSE)
})
n_skip <- which(startsWith(gpl_hdr, "!platform_table_begin"))
if (length(n_skip) != 1L) {
  stop("GPL96.annot.gz: found ", length(n_skip), " !platform_table_begin ",
       "markers in the first 200 lines, expected 1. Wrong or truncated file.",
       call. = FALSE)
}
gpl <- .vroom_big(readr::read_tsv(.local_path("gpl96_annot"), skip = n_skip,
                                  comment = "!", show_col_types = FALSE,
                                  progress = FALSE))
stopifnot(all(c("ID", "Gene symbol") %in% names(gpl)))
message("   GPL96: ", nrow(gpl), " probes annotated")
if (nrow(gpl) != EXPECT_GPL96_PROBES) {
  message("   NOTE: GPL96 has ", nrow(gpl), " probes, expected ",
          EXPECT_GPL96_PROBES, " - annotation version differs, recorded")
}

sym <- gpl[["Gene symbol"]]
names(sym) <- gpl$ID
sym <- sym[!is.na(sym) & nzchar(trimws(sym))]
n_multi <- sum(grepl("///", sym, fixed = TRUE))
sym <- sym[!grepl("///", sym, fixed = TRUE)]
message("   dropped ", n_multi, " multi-mapping probes and kept ", length(sym),
        " single-gene probes")

keep <- intersect(rownames(hz_probe), names(sym))
hz_sub <- hz_probe[keep, , drop = FALSE]
g <- unname(sym[keep])
rm_ <- rowMeans(hz_sub, na.rm = TRUE)
ord <- order(g, -rm_)                       # per gene, highest mean first
pick <- ord[!duplicated(g[ord])]
hz_m <- hz_sub[pick, , drop = FALSE]
rownames(hz_m) <- g[pick]
hz_m <- hz_m[order(rownames(hz_m)), , drop = FALSE]
message("   collapsed to ", nrow(hz_m), " genes (maxMean probe per gene)")

hz_key <- .match_ids(colnames(hz_m), hz_p, "GSE25066")
hz_p <- .align_pheno(hz_p, colnames(hz_m), hz_key, "GSE25066")

hz_p$pcr <- .pcr01(hz_p$pathologic_response_pcr_rd)
hz_p$treatment <- "taxane-anthracycline"    # single-arm (D5 section 4)
hz_p <- hz_p %>%
  dplyr::mutate(subtype = dplyr::case_when(
    is.na(er_status_ihc) | is.na(her2_status)      ~ NA_character_,
    er_status_ihc == "P" & her2_status == "N"      ~ "HRpos_HER2neg",
    er_status_ihc == "N" & her2_status == "N"      ~ "TNBC",
    er_status_ihc == "P" & her2_status == "P"      ~ "HRpos_HER2pos",
    er_status_ihc == "N" & her2_status == "P"      ~ "HRneg_HER2pos"))

message("   pCR: ", sum(hz_p$pcr, na.rm = TRUE), " of ", sum(!is.na(hz_p$pcr)),
        sprintf("  (%.1f%%)", 100 * mean(hz_p$pcr, na.rm = TRUE)))
message("   FORBIDDEN as covariates (expression-derived published predictors): ",
        paste(intersect(HATZIS_FORBIDDEN, names(hz_p)), collapse = ", "))

RES_SCALE[["GSE25066"]] <- .assert_log2(hz_m, "GSE25066")

# =============================================================================
# 5. Harmonised covariates - and what is NOT available where
# =============================================================================
# The declaration (2026-08-30 note, section 5) makes PAM50 adjustment mandatory,
# on the strength of the TCGA composition. PAM50 exists only in GSE25066. The
# harmonisable subtype variable across all three is RECEPTOR SUBTYPE, and that
# is what is built above.
#
# Recorded rather than silently substituted, because it is a departure from the
# letter of the declaration and must be visible:
#   GSE194040  hr / her2 flags -> 4 receptor subtypes. No PAM50.
#   GSE164458  TNBC only. Subtype is CONSTANT, so it cannot be adjusted for and
#              does not need to be - the trial's restriction does the work.
#   GSE25066   pam50_class present AND receptor subtype. Both carried; the
#              harmonised variable is receptor subtype so the three cohorts
#              adjust for the same thing.
message("\n5. harmonised covariates")

.cohort_pheno <- function(p, cohort) {
  tibble::tibble(
    cohort    = cohort,
    sample_id = p$sample_id,
    gsm       = p$gsm,
    pcr       = p$pcr,
    subtype   = p$subtype,
    treatment = p$treatment,
    pam50     = if ("pam50_class" %in% names(p)) as.character(p$pam50_class)
                else NA_character_)
}
pheno <- dplyr::bind_rows(
  .cohort_pheno(ispy2_p, "GSE194040"),
  .cohort_pheno(bt_p,    "GSE164458"),
  .cohort_pheno(hz_p,    "GSE25066"))

pheno %>% dplyr::count(cohort, subtype) %>%
  tidyr::pivot_wider(names_from = subtype, values_from = n, values_fill = 0L) %>%
  as.data.frame() %>% print(row.names = FALSE)

pheno %>% dplyr::count(cohort, pcr) %>%
  tidyr::pivot_wider(names_from = pcr, values_from = n, values_fill = 0L,
                     names_prefix = "pcr_") %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 6. THE SCALE FACTS - the input to the mitoPPS amendment, and nothing else
# =============================================================================
message("\n6. scale - all three cohorts")

scale_tab <- dplyr::bind_rows(RES_SCALE)
scale_tab %>% as.data.frame() %>% print(row.names = FALSE)

message("\n   GSVA:    usable in all three (log scale, kcdf = 'Gaussian').")
message("   mitoPPS: NO linear DESeq2-normalised matrix exists for any of the")
message("            three. See this script's header. The two-instrument rule")
message("            in the 2026-08-30 declaration section 6.5 needs a RECORDED")
message("            AMENDMENT on FEASIBILITY grounds before script 13 runs.")
message("            Do not resolve it by writing script 13 a different way.")

# =============================================================================
# 7. Gene coverage - would the declared model lose a term?
# =============================================================================
message("\n7. gene coverage")

MATS <- list(GSE194040 = ispy2_m, GSE164458 = bt_m, GSE25066 = hz_m)

crit <- dplyr::bind_rows(lapply(names(MATS), function(cn) {
  tibble::tibble(cohort = cn, gene = CRITICAL_GENES,
                 present = CRITICAL_GENES %in% rownames(MATS[[cn]]))
})) %>%
  tidyr::pivot_wider(names_from = gene, values_from = present)
crit %>% as.data.frame() %>% print(row.names = FALSE)

if (!all(unlist(crit[, -1]))) {
  message("\n   AT LEAST ONE CRITICAL GENE IS ABSENT. That is a fact about the")
  message("   platform, not an error. Record which term it removes from which")
  message("   cohort BEFORE script 13, and report it - do not silently drop it.")
}

# Set-level coverage against the arms H4 actually uses. Read from the snapshot,
# not rebuilt.
gs <- readRDS(file.path(DIR_RESULTS, "tcga_brca_mito_scores.rds"))
arm_sets <- gs$arm_sets

# M-a's genes are taken from the SAME object script 06 scored from - the G1
# audit's STRIPPED set (61 genes), not the raw 67-gene CSV. Reading the CSV here
# would measure coverage of a set this arm never uses, and its first column is
# Ensembl ids rather than symbols.
g1 <- readRDS(file.path(DIR_RESULTS, "g1_overlap_audit.rds"))
felsher_ma <- g1$estimators_stripped$FELSHER
stopifnot(length(felsher_ma) == 61L)

sets <- c(arm_sets[intersect(c("OXPHOS subunits", "OXPHOS assembly factors"),
                             names(arm_sets))],
          list(`Felsher M-a (stripped 61)` = felsher_ma))

cov_tab <- dplyr::bind_rows(lapply(names(MATS), function(cn) {
  dplyr::bind_rows(lapply(names(sets), function(sn) {
    g <- unique(sets[[sn]])
    tibble::tibble(cohort = cn, set = sn, n_set = length(g),
                   n_present = sum(g %in% rownames(MATS[[cn]])),
                   frac = sum(g %in% rownames(MATS[[cn]])) / length(g))
  }))
}))
cov_tab %>% as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 8. Save - three separate objects, deliberately
# =============================================================================
# There is no combined matrix and there will not be one. GSVA is cohort-relative
# (CLAUDE.md) and D5 section 6.2 requires the ESTIMATES to be meta-analysed, not
# the scores pooled. A combined matrix on disk is an invitation to do the wrong
# thing later.
message("\n8. save")

out <- list(
  cohorts = list(
    GSE194040 = list(expr = ispy2_m, pheno = ispy2_p, scale = "log2",
                     platform = "Agilent GPL20078 + GPL30493, gene level",
                     note = paste("ComBat-adjusted as deposited (D5 section 4);",
                                  "state in Methods"),
                     role = "PRIMARY"),
    GSE164458 = list(expr = bt_m, pheno = bt_p, scale = "log2",
                     platform = "RNA-seq, log2 processed as deposited",
                     note = "TNBC only; subtype is constant",
                     role = "REPLICATION"),
    GSE25066  = list(expr = hz_m, pheno = hz_p, scale = "log2",
                     platform = "Affymetrix GPL96 HG-U133A",
                     note = paste("probe-collapsed: multi-mapping dropped",
                                  "(incl. the only BBC3 probe), then maxMean"),
                     role = "THIRD")),
  pheno    = pheno,
  scale    = scale_tab,
  coverage = list(critical = crit, sets = cov_tab),
  join_keys = tibble::tibble(
    cohort    = c("GSE194040", "GSE164458", "GSE25066"),
    column    = c(ispy2_key$col, bt_key$col, hz_key$col),
    transform = c(ispy2_key$tf,  bt_key$tf,  hz_key$tf),
    matched   = c(ispy2_key$score, bt_key$score, hz_key$score)),
  rules = list(
    probe_collapse = paste("GPL96: drop empty symbols, drop '///' multi-mapping",
                           "probes, then highest mean across samples per gene"),
    endpoint       = "pCR = 1, RD = 0; unrecognised values stop the script",
    subtype        = paste("receptor subtype (HR/HER2) harmonised across all",
                           "three; PAM50 exists only in GSE25066 and is carried",
                           "alongside, not used as the harmonised variable"),
    forbidden      = HATZIS_FORBIDDEN,
    no_association = paste("no score-versus-pCR association is computed in this",
                           "script; marginals only (D5 section 2)"),
    arm_labels     = paste("I-SPY2 'AMG-386' merged into 'AMG 386' (n = 1 vs",
                           "114); separator variants only"),
    control_arm    = paste("I-SPY2's 179-patient Paclitaxel control arm is 94",
                           "HRpos_HER2neg + 85 TNBC and contains ZERO HER2+",
                           "patients, so a score x treatment contrast against",
                           "control is only estimable in those two strata"),
    mitopps        = paste("NOT available: no linear DESeq2-normalised matrix",
                           "exists in any of the three. Needs a recorded",
                           "amendment on feasibility grounds before script 13.")),
  built = Sys.time())

saveRDS(out, file.path(DIR_RESULTS, "neoadjuvant_cohorts.rds"))
readr::write_csv(pheno,     file.path(DIR_TABLES, "neoadjuvant_pheno.csv"))
readr::write_csv(scale_tab, file.path(DIR_TABLES, "neoadjuvant_scale.csv"))
readr::write_csv(cov_tab,   file.path(DIR_TABLES, "neoadjuvant_coverage.csv"))

message("\n12: done.")
message("    results/neoadjuvant_cohorts.rds")
message("    outputs/tables/  3 tables")
message("    NEXT: the mitoPPS feasibility amendment, BEFORE script 13.")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  n <- readRDS(file.path(DIR_RESULTS, "neoadjuvant_cohorts.rds"))

  # --- the two things that decide script 13's shape ------------------------
  n$scale %>% as.data.frame() %>% print(row.names = FALSE)
  n$coverage$critical %>% as.data.frame() %>% print(row.names = FALSE)

  # --- is anything actually on a linear scale after all? -------------------
  # If a matrix were linear its histogram would be right-skewed with a long
  # tail; a log2 matrix is roughly symmetric. Look before believing the header.
  op <- par(mfrow = c(1, 3))
  for (cn in names(n$cohorts)) {
    hist(sample(n$cohorts[[cn]]$expr, 2e5), breaks = 80, main = cn, xlab = "")
  }
  par(op)

  # --- the endpoint, per cohort and per subtype ----------------------------
  # MARGINALS ONLY. Do not add a score to this table.
  n$pheno %>% dplyr::count(cohort, subtype, pcr) %>%
    tidyr::pivot_wider(names_from = pcr, values_from = n, values_fill = 0L,
                       names_prefix = "pcr_") %>%
    as.data.frame() %>% print(row.names = FALSE)

  # --- I-SPY2 arms: several are HER2-targeted and confounded with subtype --
  # D5 says arm must be adjusted for or stratified. This is why.
  ip <- n$cohorts$GSE194040$pheno
  table(ip$treatment, ip$subtype)

  # --- did the GSE25066 collapse behave? -----------------------------------
  hz <- n$cohorts$GSE25066$expr
  c(genes = nrow(hz), samples = ncol(hz))
  intersect(c("BBC3", "BCL2L1", "MCL1", "BCL2L11", "MYC"), rownames(hz))
  # BBC3 is EXPECTED to be absent here - its only GPL96 probe multi-maps.

  # --- gene-name overlap across the three ----------------------------------
  # Not a merge. Just how much of the axis is measured in common.
  gg <- lapply(n$cohorts, function(x) rownames(x$expr))
  c(all_three = length(Reduce(intersect, gg)),
    vapply(gg, length, integer(1)))

  # --- what the join keys were ---------------------------------------------
  lapply(n$cohorts, function(x) utils::head(x$pheno$sample_id, 3))

}
