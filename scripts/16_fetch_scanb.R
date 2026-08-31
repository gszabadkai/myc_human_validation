# 16_fetch_scanb.R
# =============================================================================
# Fetch, parse and normalise SCAN-B (GEO GSE202203) for the BIM replication.
# NOTHING IS SCORED HERE AND NO ASSOCIATION IS COMPUTED.
#
# Built to:
#   docs/2026-08-31_scanb_bim_replication_declaration.md   <- read this first
#   docs/2026-08-29_block_c_result_H1_not_supported.md section 9 (the declaration)
#   docs/2026-08-27_human_validation_plan.md section 3 (amended 2026-08-31)
#
# =============================================================================
# WHY THIS ACCESSION AND NOT THE ONE THE PLAN NAMED
# =============================================================================
# The plan said "SCAN-B (Brueffer)". Checked 2026-08-31 against the GEO FTP
# listings and the series-matrix processing fields:
#
#   GSE96058  n=3273  gene_expression_..._transformed.csv.gz  cufflinks FPKM,
#                     deposited as log2(FPKM + 0.1).              NO COUNTS
#   GSE81538  n=405   gene_expression_405_transformed.csv.gz     NO COUNTS
#   GSE202203 n=3207  RawCounts_gene_3207.tsv.gz                 COUNTS
#
# A deposit that arrives already logged cannot carry mitoPPS - that is exactly
# what amendment A1 established for the three neoadjuvant cohorts. So GSE96058
# would have been blocked for the same reason, and GSE202203 is the only SCAN-B
# deposit on which the declaration's co-primary clause can be honoured.
#
# =============================================================================
# SCALE - THE THING THIS SCRIPT EXISTS TO GET RIGHT
# =============================================================================
# The two instruments have OPPOSITE requirements and must not share an object
# (CLAUDE.md). This script is the only place in the SCAN-B arm where the split
# happens, and it happens once:
#
#   RawCounts (NON-INTEGER, StringTie -e estimated counts)
#     -> round()                              integer counts
#     -> DESeq2 size factors
#          |-- counts(dds, normalized = TRUE)   LINEAR -> mitoPPS ONLY
#          |-- vst(dds, blind = TRUE)           LOG    -> GSVA ONLY
#
# The counts really are non-integer: the deposited file carries values like
# "A1BG 412.35820192". That is the prepDE-style output of StringTie run with
# -e, not a corrupted file. DESeq2 requires integers, and gene-level counts are
# all that is deposited, so tximport with transcript lengths is not available.
# Rounding is a DOCUMENTED DEVIATION from TCGA, where GDC counts are already
# integers - declaration section 7. Section 4 below reports how many values it
# changed so the deviation is quantified rather than asserted.
#
# =============================================================================
# COHORT-RELATIVITY
# =============================================================================
# GSVA is cohort-relative and mitoPPS is composition-dependent, so SCAN-B is
# scored in ONE run over ALL of its samples in script 17 and its scores are
# never compared numerically with TCGA's. Only the sign and the CI of the
# interaction transfer. That is why the matrices here are built on every sample
# that passes QC, not on the model's complete cases - exactly as script 01
# builds TCGA on 1,095 while Block C fits on 938.
#
# =============================================================================
# THE FENCE - declaration section 10, and it is enforced below, not trusted
# =============================================================================
# GSE202203 carries overall survival, relapse-free interval, endocrine/chemo
# treatment flags, and two expression-derived ESR1/ESR2 values. SCAN-B was
# fetched for the BIM replication and for nothing else. Those columns are
# NAMED in SCANB_FORBIDDEN and DROPPED, and section 3 asserts that none of them
# reaches the saved object.
#
# They are dropped rather than hidden behind a flag on purpose. A toggle gets
# flipped; a named list has to be edited in a commit with a dated note beside
# it. If F2 (treatment-stratified survival) is ever decided, that decision
# amends SCANB_FORBIDDEN and says so in its own note.
#
# The ESR1/ESR2 values are forbidden for a second and independent reason: they
# are functions of the same matrix the exposure is built from, so using them as
# covariates would be circular. Same rule as HATZIS_FORBIDDEN in script 12.
#
# NO GSVA. NO mitoPPS. NO MYC score. NO MODELS. Those are script 17.
#
# SPECIES: human. Human gene symbols throughout. See CLAUDE.md.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

suppressPackageStartupMessages({
  library(SummarizedExperiment)   # assay() on the DESeqTransform
  library(DESeq2)
})

message("\n16: fetch SCAN-B (GSE202203) for the BIM replication\n", strrep("=", 78))

# =============================================================================
# 0. Constants
# =============================================================================
DIR_SCANB <- here::here("data", "raw", "scanb")

GEO_BASE <- "https://ftp.ncbi.nlm.nih.gov/geo"

PATH_VST    <- file.path(DIR_RESULTS, "scanb_vst.rds")
PATH_LINEAR <- file.path(DIR_RESULTS, "scanb_linear.rds")
PATH_PHENO  <- file.path(DIR_RESULTS, "scanb_pheno.rds")

# name = c(url-relative-path, minimum bytes). Floors are ~90% of the sizes in
# the GEO directory listing read on 2026-08-31. The GEO gateway serves neither
# HEAD nor Content-Range, so a floor plus the shape assertions below is the
# substitute for an exact byte guard - same contract as scripts 02 and 12.
INPUTS <- list(
  counts = c(
    "series/GSE202nnn/GSE202203/suppl/GSE202203_RawCounts_gene_3207.tsv.gz",
    183e6),
  sm_hiseq = c(
    "series/GSE202nnn/GSE202203/matrix/GSE202203-GPL11154_series_matrix.txt.gz",
    230e3),
  sm_nextseq = c(
    "series/GSE202nnn/GSE202203/matrix/GSE202203-GPL18573_series_matrix.txt.gz",
    40e3)
)

# Expected shapes, from the GEO metadata read on 2026-08-31. The samples are
# split across two sequencers and BOTH series matrices are needed; asserting
# only the total means a missing platform file cannot pass unnoticed.
EXPECT_N_TOTAL    <- 3207L
EXPECT_N_NEXTSEQ  <- 294L      # GPL18573; GPL11154 carries the remaining 2,913

# Minimum overlap between expression column names and a phenotype identifier
# column before that column is accepted as the join key. Script 12's value.
MIN_ID_MATCH <- 0.90

# Genes the declared model needs. Absence of any of these is a stop: BCL2L11 is
# the endpoint, BCL2L1 corroborates direction, and BBC3 is the control that
# makes the other two mean anything (declaration section 11).
CRITICAL_GENES <- c("BCL2L11", "BCL2L1", "BBC3", "MYC")

# Coverage floor for the sets the declared model needs. REPORTED here and
# ENFORCED in script 17: this script builds no model, so a low fraction is not
# a reason to throw away a 200 MB download. 0.80 is script 14's MIN_SET_FRAC,
# adopted under amendment A2.
MIN_SET_FRAC <- 0.80

# The fence. Dropped, and asserted absent from the saved object. See the header.
SCANB_FORBIDDEN <- c(
  "overall survival days", "overall survival years", "overall survival event",
  "relapse free interval days", "relapse free interval years",
  "relapse free interval event",
  "endocrine treated", "chemo treated",
  "esr1 expression", "esr2 expression")

# =============================================================================
# 1. Inputs - the author downloads these; this script does not fetch
# =============================================================================
# Same contract as scripts 02 and 12: assert, and stop with the exact command.
# Nothing is downloaded behind a source().
message("\n1. inputs")

.ensure_dir(DIR_SCANB)

.local_path <- function(key) file.path(DIR_SCANB, basename(INPUTS[[key]][1]))

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
  stop("SCAN-B inputs are not ready.\n",
       if (length(missing)) paste0("MISSING (", length(missing), "):\n",
                                   paste(cmds, collapse = "\n"), "\n") else "",
       if (length(short)) paste0("TRUNCATED:\n  ",
                                 paste(short, collapse = "\n  "), "\n",
                                 "Delete and re-download these.\n") else "",
       "Provenance and the full list are in data/scanb/README.md.",
       call. = FALSE)
}
message("   all ", length(INPUTS), " inputs present and above the size floor")

# Checksums are recorded so data/scanb/README.md can carry them and so a silent
# re-download of a changed file is visible. tools::md5sum is base; sha256 needs
# openssl, which may not be installed, so it is attempted and not required.
.sha256 <- function(p) {
  if (!requireNamespace("openssl", quietly = TRUE)) return(NA_character_)
  con <- file(p, "rb")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  paste(openssl::sha256(con))
}
checksums <- tibble::tibble(
  file   = vapply(names(INPUTS), function(k) basename(.local_path(k)), character(1)),
  bytes  = vapply(names(INPUTS), function(k) file.size(.local_path(k)), numeric(1)),
  md5    = vapply(names(INPUTS), function(k) unname(tools::md5sum(.local_path(k))),
                  character(1)),
  sha256 = vapply(names(INPUTS), function(k) .sha256(.local_path(k)), character(1)))
checksums %>% as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 2. Series-matrix helpers
# =============================================================================
# Carried over from script 12 unchanged except where noted. Every one of them
# exists because it caught a real failure in the neoadjuvant cohorts that was
# invisible in the output. See data/neoadjuvant/README.md before editing any of
# them.
#
# The connection buffer: a GEO series matrix puts one field per sample on a
# single line, so !Sample_data_processing runs to hundreds of kilobytes.
# GSE25066's was 185,443 bytes against vroom's 131,072 default, which stops
# readr with a message that reads like a corrupt download and is not one.
# Raised per call and RESTORED, so sourcing leaves no global change behind.
VROOM_BUF <- 8388608L

.vroom_big <- function(expr) {
  old <- Sys.getenv("VROOM_CONNECTION_SIZE", unset = NA_character_)
  Sys.setenv(VROOM_CONNECTION_SIZE = VROOM_BUF)
  on.exit({
    if (is.na(old)) Sys.unsetenv("VROOM_CONNECTION_SIZE")
    else Sys.setenv(VROOM_CONNECTION_SIZE = old)
  }, add = TRUE)
  force(expr)
}

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

# Characteristic lines are "key: value", one value per sample, and THEY ARE NOT
# ALIGNED ACROSS SAMPLES. In GSE25066 that mislabelled 198 of 508 endpoints,
# because two deposit blocks used different field orders and a majority-key
# parser handed back RCB classes named as pCR. So the key is read PER SAMPLE
# and the table pivoted on it; a sample that never carried a key gets NA.
#
# SCAN-B's ESR1/ESR2 lines contain a SECOND colon inside the value
# ("esr1 expression: log2(tpm+0.1): 5.97"), so the key is "esr1 expression" and
# the value keeps the "log2(tpm+0.1): " prefix. Both are in SCANB_FORBIDDEN and
# never used, so the ragged value is not cleaned up - noted here so the next
# reader does not think it is a bug.
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

# The join key is DISCOVERED under a DECLARED transformation set, never assumed.
# In BrighTNess the obvious key had ZERO raw overlap ("102001_RNAseq" against a
# column "102001"), and assuming it would have turned 482 samples into none.
#
# SCAN-B adds a fourth transformation: the external id is
# "Q009012.C009079.S000008.l.r.m.c.lib.g.k2.a" and the expression columns are
# "S000008", so the S-token has to be extractable. All four are declared HERE,
# before the match rate is seen, so the set is not tuned to the answer.
.ID_TRANSFORMS <- list(
  raw            = function(x) x,
  strip_prefix   = function(x) sub("^[A-Za-z][A-Za-z0-9]*_", "", x),
  strip_suffix   = function(x) sub("_[A-Za-z][A-Za-z0-9]*$", "", x),
  scanb_s_token  = function(x) ifelse(grepl("(^|\\.)S[0-9]+(\\.|$)", x),
                                      sub("^.*(^|\\.)(S[0-9]+)(\\.|$).*$", "\\2", x),
                                      NA_character_))

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
  message("   join key: '", grid$col[1], "' under '", grid$tf[1], "' matches ",
          sprintf("%.1f%%", 100 * grid$score[1]), " of expression columns")
  list(col = grid$col[1], tf = grid$tf[1], score = grid$score[1])
}

# A file written with row.names = TRUE has ONE FEWER header field than data
# field, and readr does not treat that as row names: it MERGES THE OVERFLOW
# INTO THE LAST COLUMN, so the final sample silently becomes "8.427\t8.8056" -
# two patients glued together - while everything else parses cleanly. That
# happened in I-SPY2. The ragged case is detected, not guessed at.
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
    message("   row-named layout: ", length(hdr), " header fields, ",
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

# as.matrix() on a frame with one stray character column returns a CHARACTER
# matrix, and storage.mode(x) <- "double" then fills it with NAs behind a
# warning. Checked rather than hoped for.
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

# =============================================================================
# 3. Phenotype - both platforms, then the fence
# =============================================================================
message("\n3. phenotype")

ph_hi <- .sm_pheno(.local_path("sm_hiseq"))
ph_nx <- .sm_pheno(.local_path("sm_nextseq"))
ph_hi$platform <- "GPL11154 (HiSeq 2000)"
ph_nx$platform <- "GPL18573 (NextSeq 500)"
message("   GPL11154: ", nrow(ph_hi), " samples")
message("   GPL18573: ", nrow(ph_nx), " samples")

if (nrow(ph_nx) != EXPECT_N_NEXTSEQ) {
  stop("GPL18573 has ", nrow(ph_nx), " samples, expected ", EXPECT_N_NEXTSEQ,
       ". The deposit has changed; re-read the metadata before proceeding.",
       call. = FALSE)
}

pheno_all <- dplyr::bind_rows(ph_hi, ph_nx)
if (nrow(pheno_all) != EXPECT_N_TOTAL) {
  stop("series matrices give ", nrow(pheno_all), " samples in total, expected ",
       EXPECT_N_TOTAL, ".", call. = FALSE)
}
message("   combined: ", nrow(pheno_all), " samples, ", ncol(pheno_all),
        " metadata columns")

# --- 3.1 the fence, enforced -------------------------------------------------
# Declaration section 10. Dropped here, and asserted absent from the saved
# object in section 8, so the fence is a property of the artefact rather than a
# promise in a comment.
present_forbidden <- intersect(SCANB_FORBIDDEN, names(pheno_all))
message("   FENCE: dropping ", length(present_forbidden), " of ",
        length(SCANB_FORBIDDEN), " forbidden columns -> ",
        paste(present_forbidden, collapse = ", "))
if (!length(present_forbidden)) {
  warning("none of SCANB_FORBIDDEN was found in the deposit. Either the field ",
          "names changed or the wrong file was read - check before trusting ",
          "the fence.", call. = FALSE)
}
pheno_all <- dplyr::select(pheno_all, -dplyr::any_of(SCANB_FORBIDDEN))

# --- 3.2 the analysis variables ----------------------------------------------
# Only what the declared model needs, plus what a Methods table needs. Anything
# not named here stays out; adding a variable is a decision, not a convenience.
#
# PAM50 "Unclassified" is coded NA, not kept as a level. It is not a subtype -
# it records that the nearest-centroid call failed - and the declaration takes
# complete cases on PAM50. "Normal" IS kept as a level, mirroring TCGA, where
# BRCA_Normal (n = 29) is a covariate level that is simply never fitted as a
# stratum (script 09 section on strata).
.num <- function(x) suppressWarnings(as.numeric(x))

# Named here so a renamed GEO field fails with a sentence rather than with
# "Column `pam50 subtype` doesn't exist" from inside a transmute.
NEEDED_CH <- c("scanb external id", "pam50 subtype", "age at diagnosis", "nhg",
               "tumor size", "lymph node group", "er status", "pgr status",
               "her2 status", "ki67 status", "clinical groups",
               "histopathological type")
absent_ch <- setdiff(NEEDED_CH, names(pheno_all))
if (length(absent_ch)) {
  stop("GSE202203: characteristic field(s) missing from the series matrices -> ",
       paste(absent_ch, collapse = ", "),
       "\nThe deposit's field names have changed. Re-read the metadata and ",
       "update NEEDED_CH and section 3.2 before proceeding; do NOT guess a ",
       "substitute column.", call. = FALSE)
}

pheno <- pheno_all %>%
  dplyr::transmute(
    gsm            = .data$gsm,
    title          = .data$title,
    platform       = .data$platform,
    external_id    = .data$`scanb external id`,
    PAM50          = dplyr::na_if(.data$`pam50 subtype`, "Unclassified"),
    age            = .num(.data$`age at diagnosis`),
    NHG            = .data$nhg,
    tumor_size     = .num(.data$`tumor size`),
    lymph_node     = .data$`lymph node group`,
    er_call        = .data$`er status`,
    pgr_call       = .data$`pgr status`,
    her2_call      = .data$`her2 status`,
    ki67_call      = .data$`ki67 status`,
    clinical_group = .data$`clinical groups`,
    histology      = .data$`histopathological type`)

n_unclass <- sum(pheno_all$`pam50 subtype` == "Unclassified", na.rm = TRUE)
message("   PAM50: ", sum(!is.na(pheno$PAM50)), " called, ", n_unclass,
        " 'Unclassified' -> NA, ", sum(is.na(pheno$PAM50)) - n_unclass,
        " already missing")
pheno %>% dplyr::count(PAM50) %>% as.data.frame() %>% print(row.names = FALSE)

# --- 3.3 duplicate subjects --------------------------------------------------
# GSE96058 shipped 136 technical replicates. This file claims none, which is a
# claim and not a guarantee. The external id is Q<qid>.C<case>.S<sample>.<...>,
# so a repeated C-token is the same case sequenced twice.
pheno$case_token <- unname(vapply(pheno$external_id, function(x) {
  m <- regmatches(x, regexpr("(^|\\.)C[0-9]+(\\.|$)", x))
  if (!length(m)) NA_character_ else gsub("[.]", "", m)
}, character(1)))

dup_case <- sum(duplicated(pheno$case_token[!is.na(pheno$case_token)]))
if (dup_case > 0L) {
  message("   ", dup_case, " duplicate case token(s); keeping the first by ",
          "external id order (declaration section 11)")
  ord   <- order(pheno$external_id)
  pheno <- pheno[ord, ]
  pheno <- pheno[!duplicated(pheno$case_token) | is.na(pheno$case_token), ]
} else {
  message("   no duplicate case tokens: ",
          dplyr::n_distinct(pheno$case_token), " distinct cases in ",
          nrow(pheno), " samples")
}

# =============================================================================
# 4. Expression - read, symbols, filter, round
# =============================================================================
message("\n4. expression matrix")

# ~20k rows x 3,208 columns is 64 million cells; expect a couple of minutes
# and roughly 1.5 GB while the frame and the matrix both exist.
message("   reading 204 MB / 3,208 columns - this takes a few minutes")
raw_df <- .read_matrix_tsv(.local_path("counts"), "GSE202203 counts")
message("   read ", nrow(raw_df), " rows x ", ncol(raw_df) - 1L, " samples")

M <- .numeric_matrix(raw_df, "gene", "GSE202203 counts")
rm(raw_df); invisible(gc())

if (ncol(M) != EXPECT_N_TOTAL) {
  stop("counts matrix has ", ncol(M), " samples, expected ", EXPECT_N_TOTAL,
       ".", call. = FALSE)
}

# --- 4.1 the scale guard -----------------------------------------------------
# The mirror image of script 12's .assert_log2. There, a matrix claiming to be
# log2 with a max above 40 meant the deposit had changed. Here a matrix
# claiming to be COUNTS with a small maximum would mean the same thing, and it
# would be far more damaging: DESeq2 would run happily on FPKM and every
# mitoPPS value downstream would be wrong with nothing to show for it.
rg <- range(M, na.rm = TRUE)
message(sprintf("   range %.3f to %s, median %.2f, %.2f%% zeros",
                rg[1], format(round(rg[2]), big.mark = ","),
                stats::median(M, na.rm = TRUE), 100 * mean(M == 0, na.rm = TRUE)))
if (rg[1] < 0) {
  stop("negative values in a counts matrix (min ", signif(rg[1], 4),
       "). This is not RawCounts.", call. = FALSE)
}
if (rg[2] < 1e4) {
  stop("counts matrix has max ", signif(rg[2], 4), ", which is far too small ",
       "for counts - this looks like FPKM or TPM. The wrong supplementary ",
       "file has been downloaded: it must be RawCounts, not TPM_Raw.",
       call. = FALSE)
}
if (anyNA(M)) {
  stop(sum(is.na(M)), " NA value(s) in the counts matrix. Inspect before ",
       "proceeding; DESeq2 will not accept them.", call. = FALSE)
}

# --- 4.2 symbols -------------------------------------------------------------
# The deposit is already collapsed on gene symbols and restricted to GENCODE 27
# protein-coding transcripts (series !Sample_data_processing), so there is no
# gene_type column to filter on and none is needed. What IS needed is the same
# duplicate handling as script 01, in case the collapse left any: keep the row
# with the highest mean rather than summing, so the value stays on the same
# scale as every other gene.
sym <- rownames(M)
bad <- is.na(sym) | !nzchar(sym) | sym %in% c("NA", "-", ".")
if (any(bad)) {
  message("   dropping ", sum(bad), " row(s) with no usable symbol")
  M <- M[!bad, , drop = FALSE]; sym <- sym[!bad]
}
if (any(duplicated(sym))) {
  mu  <- rowMeans(M)
  ord <- order(sym, -mu)
  M   <- M[ord, , drop = FALSE]; sym <- sym[ord]
  first <- !duplicated(sym)
  message("   ", sum(!first), " duplicate symbol row(s) collapsed (kept highest mean)")
  M <- M[first, , drop = FALSE]; sym <- sym[first]
}
rownames(M) <- sym
stopifnot(!any(duplicated(rownames(M))))

# --- 4.3 rounding, and how much it changed -----------------------------------
# The documented deviation, quantified. Reported BEFORE rounding so the number
# is a property of the deposit rather than of this script.
counts          <- round(M)
frac_noninteger <- mean(M != counts)
max_shift       <- max(abs(M - counts))
message(sprintf("   %.1f%% of values are non-integer (StringTie estimated counts)",
                100 * frac_noninteger),
        sprintf("; max rounding shift %.3f", max_shift))
if (max_shift > 0.5 + 1e-9) {
  stop("a rounding shift above 0.5 is arithmetically impossible; the matrix ",
       "is not what it claims to be.", call. = FALSE)
}
rm(M); invisible(gc())

# Integer storage halves the matrix and is what DESeq2 coerces to anyway - but
# it silently yields NA above 2^31-1, so the range is checked first rather than
# discovered as a wall of NAs inside estimateSizeFactors().
if (max(counts) > .Machine$integer.max) {
  stop("a count exceeds .Machine$integer.max (", max(counts), "). DESeq2 ",
       "coerces to integer and would produce NAs. Inspect before proceeding.",
       call. = FALSE)
}
storage.mode(counts) <- "integer"

# --- 4.4 low-count filter ----------------------------------------------------
# Script 01's rule verbatim: >= 10 counts in >= 10 samples. Deliberately
# permissive there, and even more permissive here in proportional terms -
# 10 of 3,207 samples is 0.3% against 0.9% in TCGA. Kept LITERAL rather than
# rescaled, because the point of the filter is to remove genes that are never
# measured, and because a filter that differs between the discovery and the
# replication cohort is one more thing that could explain a null.
keep_g <- rowSums(counts >= 10) >= 10
message("   low-count filter (>=10 counts in >=10 samples): ", sum(keep_g),
        " kept, ", sum(!keep_g), " dropped")
counts <- counts[keep_g, , drop = FALSE]
message("   final matrix: ", nrow(counts), " genes x ", ncol(counts), " samples")

# =============================================================================
# 5. Join expression to phenotype
# =============================================================================
message("\n5. join")

key <- .match_ids(colnames(counts), pheno, "GSE202203")
pheno$sample_id <- .ID_TRANSFORMS[[key$tf]](trimws(as.character(pheno[[key$col]])))
if (anyDuplicated(pheno$sample_id)) {
  stop("join key '", key$col, "' is not unique after '", key$tf, "'.",
       call. = FALSE)
}

# Samples dropped as duplicate cases in 3.3 have no phenotype row, so the
# matrix is restricted to what survived rather than the join being forced.
keep_s <- colnames(counts) %in% pheno$sample_id
if (!all(keep_s)) {
  message("   dropping ", sum(!keep_s), " expression column(s) with no ",
          "phenotype row (duplicate cases removed in 3.3)")
  counts <- counts[, keep_s, drop = FALSE]
}
pheno <- pheno[match(colnames(counts), pheno$sample_id), ]
if (!identical(pheno$sample_id, colnames(counts))) {
  stop("phenotype and expression are not aligned after the join.", call. = FALSE)
}
message("   aligned: ", ncol(counts), " samples")

# =============================================================================
# 6. The two matrices - built separately, on purpose
# =============================================================================
# Script 01 section 4, applied to SCAN-B. blind = TRUE for the same reason: no
# design is being modelled here, and a blind transform keeps the object usable
# downstream without circularity.
#
# MEMORY: ~19k genes x 3,207 samples is about 500 MB per double matrix, and
# there are two of them plus the dds. Expect this section to want ~3-4 GB and a
# few minutes.
message("\n6. normalisation")

cd  <- data.frame(sample = colnames(counts), row.names = colnames(counts))
dds <- DESeqDataSetFromMatrix(counts, cd, design = ~ 1)
dds <- estimateSizeFactors(dds)
sf  <- sizeFactors(dds)
message(sprintf("   size factors: median %.3f, range %.3f to %.3f",
                stats::median(sf), min(sf), max(sf)))

# --- 6a. LINEAR, for mitoPPS only -------------------------------------------
mat_linear <- counts(dds, normalized = TRUE)
stopifnot(min(mat_linear) >= 0)
message(sprintf("6a. LINEAR  range %.1f to %s", min(mat_linear),
                format(round(max(mat_linear)), big.mark = ",")))

# --- 6b. VST, log scale, for GSVA -------------------------------------------
vsd     <- vst(dds, blind = TRUE)
mat_vst <- assay(vsd)
message(sprintf("6b. VST     range %.1f to %.1f", min(mat_vst), max(mat_vst)))

# The assertion that the two really are on different scales. If this ever
# fails, something has collapsed them and every score downstream is suspect.
if (max(mat_vst) > 100 || max(mat_linear) < 100) {
  stop("scale assertion failed: VST max ", round(max(mat_vst), 1),
       ", linear max ", round(max(mat_linear), 1),
       ". The two matrices are not on the scales they claim.", call. = FALSE)
}

# =============================================================================
# 7. Coverage guards
# =============================================================================
message("\n7. coverage")

# --- 7.1 the genes the declared model needs ---------------------------------
crit <- tibble::tibble(
  gene    = CRITICAL_GENES,
  present = CRITICAL_GENES %in% rownames(mat_vst))
crit %>% as.data.frame() %>% print(row.names = FALSE)
if (!all(crit$present)) {
  stop("critical gene(s) absent from SCAN-B: ",
       paste(crit$gene[!crit$present], collapse = ", "),
       ".\nBCL2L11 is the endpoint, BCL2L1 corroborates direction, and BBC3 is ",
       "the control that makes the other two mean anything. The replication ",
       "cannot be interpreted without all of them (declaration section 11).",
       call. = FALSE)
}

# --- 7.2 mtDNA-encoded genes -------------------------------------------------
# GENCODE 27 protein-coding carries the 13 MT- genes, but the deposit is
# symbol-collapsed and that is a claim worth checking rather than assuming.
# They matter because mitoPPS's universe splits them into their own synthetic
# pathway (CLAUDE.md standing convention) and because the matched null skips
# that pathway for abundance reasons.
mt13 <- c("MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
          "MT-CO1", "MT-CO2", "MT-CO3", "MT-ATP6", "MT-ATP8", "MT-CYB")
mt_found <- intersect(mt13, rownames(mat_vst))
message("   mtDNA-encoded protein-coding genes present: ", length(mt_found),
        " of 13")
if (length(mt_found) < 13L) {
  message("   MISSING: ", paste(setdiff(mt13, mt_found), collapse = ", "))
  message("   Not a stop. Script 17 must decide whether the mtDNA pathway is ",
          "constructible; the declared model does not use it, and the null ",
          "skips it in TCGA too.")
}

# --- 7.3 symbol harmonisation to the SCAN-B vocabulary ----------------------
# THIS SECTION IS NOT OPTIONAL AND IT IS NOT A CONVENIENCE.
#
# SCAN-B is annotated against UCSC knownGenes downloaded 22 SEPTEMBER 2014.
# Human MitoCarta 3.0 (2020) carries CURRENT HGNC symbols. The ATP synthase
# subunits were renamed wholesale in 2018, so the two vocabularies disagree on
# nineteen of the eighty-nine genes in `OXPHOS subunits` - which is the
# EXPOSURE of the declared model:
#
#   ATP5F1A <- ATP5A1     ATP5MC1 <- ATP5G1     ATP5PB  <- ATP5F1
#   ATP5F1B <- ATP5B      ATP5MC2 <- ATP5G2     ATP5PD  <- ATP5H
#   ATP5F1C <- ATP5C1     ATP5MC3 <- ATP5G3     ATP5PF  <- ATP5J
#   ATP5F1D <- ATP5D      ATP5MD  <- USMG5      ATP5PO  <- ATP5O
#   ATP5F1E <- ATP5E      ATP5ME  <- ATP5I      ATP5IF1 <- ATPIF1
#   ATP5MPL <- C14orf2    ATP5MF  <- ATP5J2     ATP5MG  <- ATP5L
#   DMAC2L  <- ATP5S
#
# Unharmonised, `OXPHOS subunits` covers 0.787 of its genes here against 1.000
# in TCGA, and script 17's coverage floor would stop the replication. It would
# stop for the WRONG REASON: the genes are all present, under the names they
# had in 2014. Silently accepting 70 of 89 would be worse still - a Complex V
# with no F1 head and no c-ring is a different exposure wearing the same label.
#
# The map is script 07 section 2's, applied to this matrix instead of TCGA's,
# and reproduced here rather than sourced because script 07 is a pipeline script
# with side effects, not a function library. Its properties are what make it
# safe, and they are worth restating:
#
#   - the source is MitoCarta's own curated `Synonyms` column, not a guess and
#     not a new annotation package introduced mid-arm;
#   - it runs FORWARD only, current symbol -> its listed alias. The reverse
#     direction is dangerous and script 07 documents why: COX1/COX2/COX3
#     resolve to the PROSTAGLANDIN SYNTHASES, injecting two abundant
#     inflammatory genes into an OXPHOS set, invisibly;
#   - a candidate that is itself a current MitoCarta symbol for a different
#     gene is refused, and so is any symbol with more than one surviving
#     candidate. Ambiguity is left UNRESOLVED and reported by name.
#
# Verified on this matrix 2026-08-31: all 19 resolve, none ambiguous, no two
# map to the same target, and coverage goes 0.787 -> 1.000.
#
# The Felsher and Hallmark sets are NOT MitoCarta genes, so MitoCarta's
# synonyms cannot resolve their handful of renames (H2AX, POLR1G, VARS1 and a
# few more). They sit at 0.95-0.98, comfortably above the floor, and they are
# REPORTED BY NAME rather than rescued with a general alias source. Introducing
# one for this cohort only would be a new instrument, chosen after seeing which
# genes were missing, and that is exactly the move the coverage floor exists to
# make unnecessary.
message("\n7.3 symbol harmonisation to the SCAN-B vocabulary")

gs       <- readRDS(file.path(DIR_RESULTS, "tcga_brca_mito_scores.rds"))
arm_sets <- gs$arm_sets
cov_sets <- gs$covariate_sets
g1       <- readRDS(file.path(DIR_RESULTS, "g1_overlap_audit.rds"))
felsher_ma <- g1$estimators_stripped$FELSHER
stopifnot(length(felsher_ma) == 61L)

# intersect() would quietly return an empty list for a renamed set, and the
# coverage floor below would then pass because there was nothing to fail.
need_sets <- c("OXPHOS subunits", "PROLIF_DISJOINT", "PROLIF_STD")
have_sets <- c(names(arm_sets), names(cov_sets))
if (!all(need_sets %in% have_sets)) {
  stop("set(s) the declared model needs are absent from the script 07 object: ",
       paste(setdiff(need_sets, have_sets), collapse = ", "), call. = FALSE)
}

sets <- c(arm_sets["OXPHOS subunits"],
          cov_sets[c("PROLIF_DISJOINT", "PROLIF_STD")],
          list(`Felsher M-a (stripped 61)` = felsher_ma))

mitocarta_background <- suppressWarnings(readxl::read_xls(PATH_MITOCARTA, sheet = 3))
stopifnot(all(c("Symbol", "Synonyms") %in% colnames(mitocarta_background)))

MATRIX_SYMBOLS <- rownames(mat_vst)
MC_SYMBOLS     <- unique(mitocarta_background$Symbol)

.syn_map <- local({
  syn <- strsplit(ifelse(is.na(mitocarta_background$Synonyms), "",
                         mitocarta_background$Synonyms), "|", fixed = TRUE)
  tibble::tibble(symbol = rep(mitocarta_background$Symbol, lengths(syn)),
                 alias  = unlist(syn)) %>%
    dplyr::filter(!is.na(alias), alias != "", alias != symbol) %>%
    dplyr::distinct()
})

# ONE map over the union of every symbol the declared model uses, applied to
# each set separately - script 07's rule. Building it per set would let the
# setdiff() guard below see a different "already in the set" universe each
# time.
.build_symbol_map <- function(genes) {
  genes   <- sort(unique(genes[!is.na(genes) & genes != ""]))
  present <- genes %in% MATRIX_SYMBOLS
  out     <- stats::setNames(genes, genes)
  status  <- ifelse(present, "matched", "unresolved")

  for (i in which(!present)) {
    g    <- genes[[i]]
    cand <- .syn_map$alias[.syn_map$symbol == g]
    cand <- cand[cand %in% MATRIX_SYMBOLS]
    cand <- cand[!(cand %in% MC_SYMBOLS & cand != g)]
    cand <- setdiff(cand, genes)
    if (length(cand) == 1L) {
      out[[g]]    <- cand
      status[[i]] <- "resolved"
    } else if (length(cand) > 1L) {
      status[[i]] <- "ambiguous"
    }
  }
  list(map = out,
       report = tibble::tibble(input_symbol = genes, status = status,
                               resolved_to = unname(out[genes])))
}

bm           <- .build_symbol_map(unlist(sets, use.names = FALSE))
symbol_map   <- bm$map
symbol_report <- bm$report
print(table(symbol_report$status))

resolved <- symbol_report[symbol_report$status == "resolved", ]
if (nrow(resolved)) {
  message("   ", nrow(resolved), " symbol(s) resolved through MitoCarta synonyms:")
  resolved %>% as.data.frame() %>% print(row.names = FALSE)
}
ambig <- symbol_report$input_symbol[symbol_report$status == "ambiguous"]
if (length(ambig)) {
  message("   AMBIGUOUS, left unresolved: ", paste(ambig, collapse = ", "))
}

# Two symbols must never resolve onto the same row - that would double-weight a
# gene inside a pathway mean and there would be nothing to see downstream.
mapped_all <- unname(symbol_map)
mapped_in  <- mapped_all[mapped_all %in% MATRIX_SYMBOLS]
if (anyDuplicated(mapped_in)) {
  stop("the symbol map sends two inputs to the same matrix row: ",
       paste(unique(mapped_in[duplicated(mapped_in)]), collapse = ", "),
       call. = FALSE)
}

.remap <- function(g) unname(symbol_map[g])

# --- 7.4 the sets script 17 will score, before and after harmonisation ------
# REPORTED here, ENFORCED in 17. This script builds no model, so a low fraction
# is not a reason to discard a 200 MB download - but it is something the author
# must see before spending time on script 17.
cov_tab <- dplyr::bind_rows(lapply(names(sets), function(sn) {
  g <- unique(sets[[sn]])
  h <- .remap(g)
  tibble::tibble(set = sn, n_set = length(g),
                 n_raw       = sum(g %in% MATRIX_SYMBOLS),
                 frac_raw    = sum(g %in% MATRIX_SYMBOLS) / length(g),
                 n_present   = sum(h %in% MATRIX_SYMBOLS),
                 frac        = sum(h %in% MATRIX_SYMBOLS) / length(g),
                 missing     = paste(g[!h %in% MATRIX_SYMBOLS], collapse = ", "))
}))
cov_tab %>% dplyr::select(-"missing") %>% as.data.frame() %>%
  print(row.names = FALSE)
for (i in seq_len(nrow(cov_tab))) {
  if (nzchar(cov_tab$missing[i])) {
    message("   ", cov_tab$set[i], " still missing: ", cov_tab$missing[i])
  }
}

low <- cov_tab$set[cov_tab$frac < MIN_SET_FRAC]
if (length(low)) {
  warning("set(s) below the ", MIN_SET_FRAC, " coverage floor in SCAN-B ",
          "EVEN AFTER harmonisation: ", paste(low, collapse = ", "),
          ". Script 17 enforces this; decide before running it.",
          call. = FALSE)
} else {
  message("   all sets at or above the ", MIN_SET_FRAC,
          " coverage floor after harmonisation")
}

# --- 7.5 what the declared model will actually have --------------------------
# Complete cases on what the primary specification needs. Reported, not applied
# - the matrices stay whole, because GSVA scores the cohort in one run and only
# the MODEL takes complete cases (declaration section 11).
model_ok <- !is.na(pheno$PAM50)
message("   complete cases on PAM50 (the only covariate that can be missing): ",
        sum(model_ok), " of ", nrow(pheno))
message("   NOTE: n for the fits is settled in script 17, after scoring. ",
        "Nothing is pre-set here.")

# =============================================================================
# 8. Save
# =============================================================================
message("\n8. save")

prov <- list(
  accession   = "GSE202203",
  cohort      = "SCAN-B (Sweden Cancerome Analysis Network - Breast)",
  purpose     = paste("the BIM replication declared in the Block C note",
                      "section 9; see docs/2026-08-31_scanb_bim_replication",
                      "_declaration.md"),
  platform    = "Illumina HiSeq 2000 (GPL11154) + NextSeq 500 (GPL18573)",
  pipeline    = paste("HISAT2 -> StringTie -e, GENCODE 27 protein-coding,",
                      "collapsed on gene symbols by the depositors"),
  deviation   = paste("deposited counts are NON-INTEGER StringTie estimates",
                      "and are rounded for DESeq2;", sprintf("%.1f%%", 100 * frac_noninteger),
                      "of values were non-integer, max shift",
                      sprintf("%.3f", max_shift)),
  filter      = "protein-coding as deposited; >=10 counts in >=10 samples",
  join        = sprintf("'%s' under '%s', %.1f%% matched", key$col, key$tf,
                        100 * key$score),
  pam50       = "'Unclassified' coded NA; 'Normal' kept as a level",
  symbols     = paste("SCAN-B is annotated against UCSC knownGenes 2014 and",
                      "carries pre-2018 HGNC symbols. The declared sets are",
                      "harmonised to it through MitoCarta's own Synonyms",
                      "column, forward direction only (script 07 section 2).",
                      "symbol_map and symbol_report carry the result; script 17",
                      "MUST apply symbol_map before scoring."),
  duplicates  = sprintf("%d duplicate case token(s) removed", dup_case),
  forbidden   = SCANB_FORBIDDEN,
  fence       = paste("survival, treatment and ESR1/ESR2 columns are dropped;",
                      "SCAN-B was fetched for the BIM replication and nothing",
                      "else (declaration section 10)"),
  checksums   = checksums,
  built       = Sys.time())

# The fence, asserted on the artefact rather than promised in a comment.
leaked <- intersect(SCANB_FORBIDDEN, names(pheno))
if (length(leaked)) {
  stop("FENCE BREACH: forbidden column(s) reached the saved phenotype -> ",
       paste(leaked, collapse = ", "), call. = FALSE)
}
if (any(grepl("surviv|relapse|chemo|endocrine|esr[12]", names(pheno),
              ignore.case = TRUE))) {
  stop("FENCE BREACH: a column whose name looks like an outcome or a treatment ",
       "reached the saved phenotype -> ",
       paste(grep("surviv|relapse|chemo|endocrine|esr[12]", names(pheno),
                  ignore.case = TRUE, value = TRUE), collapse = ", "),
       call. = FALSE)
}

saveRDS(list(mat = mat_linear, scale = "linear_deseq2_normalised",
             consumer = "mitoPPS only", pheno = pheno, prov = prov),
        PATH_LINEAR)
saveRDS(list(mat = mat_vst, scale = "log_vst",
             consumer = "GSVA, kcdf = Gaussian", pheno = pheno, prov = prov),
        PATH_VST)
saveRDS(list(pheno = pheno, prov = prov, coverage = cov_tab, critical = crit,
             size_factors = sf, mt_found = mt_found,
             symbol_map = symbol_map, symbol_report = symbol_report), PATH_PHENO)

readr::write_csv(pheno,   file.path(DIR_TABLES, "scanb_pheno.csv"))
readr::write_csv(cov_tab, file.path(DIR_TABLES, "scanb_coverage.csv"))

message("\n16: done.")
message("    results/scanb_linear.rds   LINEAR, for mitoPPS ONLY")
message("    results/scanb_vst.rds      LOG (VST), for GSVA ONLY")
message("    results/scanb_pheno.rds    phenotype + provenance + guards")
message("    outputs/tables/  2 tables")
message("    NEXT: script 17. Its section 2 is the TCGA calibration and it is ",
        "a HARD GATE - if it fails, nothing in SCAN-B is fitted.")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  v <- readRDS(file.path(DIR_RESULTS, "scanb_vst.rds"))
  l <- readRDS(file.path(DIR_RESULTS, "scanb_linear.rds"))
  p <- readRDS(file.path(DIR_RESULTS, "scanb_pheno.rds"))

  # --- the scale split, which is the whole point of this script -------------
  c(vst = v$scale, linear = l$scale)
  summary(as.vector(v$mat[1:200, 1:20]))     # roughly 4-16
  summary(as.vector(l$mat[1:200, 1:20]))     # long right tail, not logged

  # A log matrix is roughly symmetric; a linear one is not. Look, do not trust
  # the label.
  op <- par(mfrow = c(1, 2))
  hist(sample(v$mat, 2e5), breaks = 80, main = "VST (log)", xlab = "")
  hist(log10(sample(l$mat, 2e5) + 1), breaks = 80, main = "linear, log10 axis",
       xlab = "")
  par(op)

  # --- the three endpoint genes and the exposure ---------------------------
  v$mat[intersect(c("BCL2L11", "BCL2L1", "BBC3", "MYC"), rownames(v$mat)), 1:6]

  # --- covariate composition, for Methods ----------------------------------
  p$pheno %>% dplyr::count(PAM50) %>% as.data.frame() %>% print(row.names = FALSE)
  p$pheno %>% dplyr::count(NHG)   %>% as.data.frame() %>% print(row.names = FALSE)
  summary(p$pheno$age)
  table(p$pheno$platform)

  # --- is the platform split confounded with subtype? ----------------------
  # Not a covariate in the declared model. If this table is badly unbalanced,
  # say so in Methods rather than adding a term after the fact.
  table(p$pheno$platform, p$pheno$PAM50)

  # --- coverage, and the join ----------------------------------------------
  p$coverage %>% dplyr::select(-missing) %>% as.data.frame() %>%
    print(row.names = FALSE)
  p$prov$join
  p$prov$deviation

  # --- the symbol map, which is the thing most likely to be wrong ----------
  # Every "resolved" row should be a historical HGNC rename you can recognise.
  # If any of them is not, stop and look before script 17 scores anything.
  table(p$symbol_report$status)
  p$symbol_report %>% dplyr::filter(status == "resolved") %>%
    as.data.frame() %>% print(row.names = FALSE)

  # The exposure, gene by gene, as script 17 will actually see it.
  ox <- unique(readRDS(file.path(DIR_RESULTS,
        "tcga_brca_mito_scores.rds"))$arm_sets[["OXPHOS subunits"]])
  data.frame(mitocarta = ox, in_scanb_as = unname(p$symbol_map[ox]),
             found = unname(p$symbol_map[ox]) %in% rownames(v$mat))

  # --- the fence -----------------------------------------------------------
  # Should be character(0). If it is not, something re-added an outcome column.
  grep("surviv|relapse|chemo|endocrine|esr", names(p$pheno),
       ignore.case = TRUE, value = TRUE)

  # --- size factors: a bimodal distribution would mean a batch effect the
  #     blind VST cannot see. Descriptive only.
  hist(p$size_factors, breaks = 60, main = "DESeq2 size factors", xlab = "")

}
