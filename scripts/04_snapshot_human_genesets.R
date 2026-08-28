# 04_snapshot_human_genesets.R
# =============================================================================
# GATE G1 - MYC signature overlap audit (plan section 6).
#
# SCALE: not applicable. This script handles GENE IDENTITY ONLY. No expression
# data is read, so neither the GSVA log-scale rule nor the mitoPPS linear-scale
# rule applies here. The first script where scale matters is 07.
#
# The snapshotting this script was originally specified to do is already done
# and committed (f63f14d). The three inputs live in data/ with provenance
# READMEs. This script therefore LOADS and VALIDATES them, harmonises gene
# symbols, and runs the audit.
#
# Question G1 answers: how much of each MYC estimator is definitionally
# mitochondrial, or definitionally proliferative? Without this, any
# MYC-to-OXPHOS result is partly true by construction.
#
# DECISION THIS FEEDS (plan section 6, D2): if the MitoCarta-stripped Felsher
# signature falls below ~50 genes it cannot serve as the primary MYC estimator,
# and the CollecTRI regulon plus the 8q24 CNV instrument carry the MYC axis.
#
# WHAT THIS SCRIPT DOES NOT ANSWER. Plan section 6 makes the signature fail on
# either of two criteria: falling below ~50 genes, OR losing its correlation
# structure. Correlation structure needs expression data, which arrives in
# script 01. This script settles the COUNT criterion only. G1 is not fully
# discharged until the correlation half runs. Do not report G1 as closed on the
# strength of this script alone.
#
# Strip set (decided 2026-08-28): the FULL 1,136-gene MitoCarta inventory, not
# just OXPHOS subunits. Conservative - any mitochondrial gene is removed. The
# OXPHOS breakdown is reported as the diagnostic, since OXPHOS-subunit overlap
# is what would create circularity with the section 7.2 OXPHOS score.
#
# Outputs:
#   results/g1_overlap_audit.rds        - audit tibble, stripped gene vectors
#                                         (script 06 consumes the vectors),
#                                         harmonisation report, metadata
#   outputs/tables/g1_overlap_audit.csv
#   outputs/tables/g1_symbol_harmonisation.csv
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

# =============================================================================
# 1. Load and validate the snapshots
# =============================================================================
# Every load asserts the shape recorded in the provenance README. A file swapped
# underneath the pipeline fails here rather than quietly changing a result.

message("04: loading snapshots")

mitocarta_inventory <- suppressWarnings(
  readxl::read_xls(PATH_MITOCARTA, sheet = 2)
)
mitocarta_background <- suppressWarnings(
  readxl::read_xls(PATH_MITOCARTA, sheet = 3)
)
# Sheet 4 carries a leading unnamed index column that readxl names "2".
# Address columns by name, never by position. See the MitoCarta README.
mitocarta_pathways <- suppressWarnings(
  readxl::read_xls(PATH_MITOCARTA, sheet = 4)
)

stopifnot(
  nrow(mitocarta_inventory)  == EXPECT_MITOCARTA_GENES,
  nrow(mitocarta_background) == EXPECT_MITOCARTA_BACKGRD,
  nrow(mitocarta_pathways)   == EXPECT_MITOCARTA_PATHWAYS,
  all(c("Symbol", "Synonyms") %in% colnames(mitocarta_background)),
  all(c("MitoPathway", "Genes") %in% colnames(mitocarta_pathways))
)

felsher_raw_table <- utils::read.csv(PATH_FELSHER, stringsAsFactors = FALSE)
stopifnot(
  nrow(felsher_raw_table) == EXPECT_FELSHER_GENES,
  "Symbol" %in% colnames(felsher_raw_table)
)

collectri <- utils::read.delim(gzfile(PATH_COLLECTRI), stringsAsFactors = FALSE)
stopifnot(
  nrow(collectri) == EXPECT_COLLECTRI_ROWS,
  dplyr::n_distinct(collectri$source_genesymbol) == EXPECT_COLLECTRI_TFS,
  all(c("source_genesymbol", "target_genesymbol",
        "is_stimulation", "is_inhibition") %in% colnames(collectri))
)

message("04: snapshots validated against expected shapes")

# =============================================================================
# 2. Symbol harmonisation
# =============================================================================
# NOT COSMETIC. Five of the 67 Felsher symbols are deprecated and match nothing
# in a 19,247-gene background: CD3EAP, CIRH1A, VARS, KARS, METTL13. One of them
# matters - KARS resolves to KARS1, which IS in MitoCarta. A naive symbol join
# therefore reports the Felsher-MitoCarta overlap as 5 when it is 6.
#
# Two hazards in the alias table constrain the rule:
#   - ~2,700 aliases are ambiguous (map to more than one current symbol)
#   - ~1,000 aliases are ALSO the current symbol of a DIFFERENT gene
#
# So the rule is: only attempt to resolve a symbol that failed to match a
# current symbol, and only accept the resolution when the alias is unambiguous.
# Anything ambiguous or unresolvable is KEPT under its original symbol - it is
# still a member of the signature - and reported by name. Dropping it would
# shrink the estimator artificially and bias the G1 threshold test.
#
# Vocabulary: everything is harmonised to MitoCarta Sheet-3 symbols, so all
# sets in the audit share one vocabulary. MSigDB and CollecTRI ship current
# HGNC; MitoCarta 3.0 is 2020-era. Harmonising only one side would undercount.

.build_alias_map <- function(background) {
  syn <- strsplit(
    ifelse(is.na(background$Synonyms), "", background$Synonyms),
    "|", fixed = TRUE
  )
  tibble::tibble(
    alias  = unlist(syn),
    symbol = rep(background$Symbol, lengths(syn))
  ) %>%
    dplyr::filter(!is.na(alias), alias != "") %>%
    dplyr::distinct()
}

alias_map        <- .build_alias_map(mitocarta_background)
current_symbols  <- unique(mitocarta_background$Symbol)

# Returns the harmonised vector plus a per-gene status report.
# status: matched | resolved | ambiguous | unresolved
.harmonise <- function(symbols, set_name) {
  symbols <- unique(symbols[!is.na(symbols) & symbols != ""])

  matched   <- symbols[symbols %in% current_symbols]
  unmatched <- setdiff(symbols, current_symbols)

  hits <- alias_map %>% dplyr::filter(alias %in% unmatched)

  per_alias <- hits %>%
    dplyr::group_by(alias) %>%
    dplyr::summarise(n_symbol = dplyr::n_distinct(symbol), .groups = "drop")

  unambiguous <- per_alias$alias[per_alias$n_symbol == 1L]
  ambiguous   <- per_alias$alias[per_alias$n_symbol >  1L]
  unresolved  <- setdiff(unmatched, per_alias$alias)

  resolution <- hits %>%
    dplyr::filter(alias %in% unambiguous) %>%
    dplyr::distinct(alias, symbol)

  # Keep ambiguous and unresolved under their original symbol: they remain
  # members of the estimator, they simply cannot be checked against MitoCarta.
  out <- stats::setNames(symbols, symbols)
  if (nrow(resolution) > 0) {
    out[resolution$alias] <- resolution$symbol
  }
  final <- unique(unname(out))

  report <- tibble::tibble(
    set          = set_name,
    input_symbol = symbols,
    status       = dplyr::case_when(
      symbols %in% matched     ~ "matched",
      symbols %in% unambiguous ~ "resolved",
      symbols %in% ambiguous   ~ "ambiguous",
      TRUE                     ~ "unresolved"
    ),
    resolved_to  = unname(out[symbols])
  )

  if (length(final) < length(symbols)) {
    message(sprintf(
      "04: NOTE %s - harmonisation collapsed %d symbol(s) onto an existing member",
      set_name, length(symbols) - length(final)
    ))
  }

  list(genes = final, report = report)
}

# =============================================================================
# 3. Build the reference sets and the MYC estimators
# =============================================================================

.split_genes <- function(x) trimws(unlist(strsplit(x, ",", fixed = TRUE)))

.pathway_genes <- function(name) {
  row <- mitocarta_pathways$Genes[mitocarta_pathways$MitoPathway == name]
  if (length(row) != 1L) {
    stop("MitoPathway not found, or not unique: ", name, call. = FALSE)
  }
  .split_genes(row)
}

# --- Reference sets ----------------------------------------------------------
# MITOCARTA_ALL is the strip set. The OXPHOS sets are diagnostics: the umbrella
# includes assembly factors, the subunit set does not. Plan section 7.2 uses the
# SUBUNIT set as the primary OXPHOS measure, so that is the circularity-relevant
# number. They are reported separately and never conflated.
ref_sets <- list(
  MITOCARTA_ALL             = unique(mitocarta_inventory$Symbol),
  MITOCARTA_OXPHOS_UMBRELLA = .pathway_genes("OXPHOS"),
  MITOCARTA_OXPHOS_SUBUNITS = .pathway_genes("OXPHOS subunits"),
  MITOCARTA_MTDNA_ENCODED   = grep("^MT-", mitocarta_inventory$Symbol, value = TRUE)
)

stopifnot(length(ref_sets$MITOCARTA_MTDNA_ENCODED) == EXPECT_MITOCARTA_MTDNA)

hallmark <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
for (gs in c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")) {
  ref_sets[[gs]] <- unique(hallmark$gene_symbol[hallmark$gs_name == gs])
}

# --- MYC estimators ----------------------------------------------------------
felsher_genes <- felsher_raw_table$Symbol

myc_edges <- collectri %>% dplyr::filter(source_genesymbol == "MYC")

# is_stimulation / is_inhibition arrive as the STRINGS "True"/"False".
# Convert explicitly; do not rely on coercion. See the CollecTRI README.
myc_edges <- myc_edges %>%
  dplyr::mutate(
    stimulation = is_stimulation == "True",
    inhibition  = is_inhibition  == "True"
  )

n_myc_both <- sum(myc_edges$stimulation & myc_edges$inhibition)
if (n_myc_both > 0) {
  message(sprintf(
    "04: NOTE %d MYC edge(s) flagged as both stimulation and inhibition; these are kept in the ALL set and also counted as stimulatory in the STIM set",
    n_myc_both
  ))
}
if ("MYC" %in% myc_edges$target_genesymbol) {
  message("04: NOTE MYC appears in its own CollecTRI regulon (autoregulation); retained")
}

est_raw <- list(
  FELSHER           = felsher_genes,
  COLLECTRI_MYC_ALL = unique(myc_edges$target_genesymbol),
  COLLECTRI_MYC_STIM = unique(myc_edges$target_genesymbol[myc_edges$stimulation])
)

# --- Harmonise everything to one vocabulary ---------------------------------
harm_est <- purrr::imap(est_raw, ~ .harmonise(.x, .y))
harm_ref <- purrr::imap(ref_sets, ~ .harmonise(.x, .y))

est_sets <- purrr::map(harm_est, "genes")
ref_sets <- purrr::map(harm_ref, "genes")

harmonisation_report <- dplyr::bind_rows(
  purrr::map_dfr(harm_est, "report"),
  purrr::map_dfr(harm_ref, "report")
)

message("04: harmonisation status across all sets")
harmonisation_report %>% dplyr::count(set, status) %>% print(n = Inf)

# =============================================================================
# 4. Overlap audit
# =============================================================================

.overlap <- function(est_name, est, ref_name, ref) {
  ov <- sort(intersect(est, ref))
  tibble::tibble(
    estimator         = est_name,
    n_estimator       = length(est),
    reference         = ref_name,
    n_reference       = length(ref),
    n_overlap         = length(ov),
    frac_of_estimator = length(ov) / length(est),
    overlap_genes     = paste(ov, collapse = ";")
  )
}

audit_raw <- purrr::map_dfr(names(est_sets), function(en) {
  purrr::map_dfr(names(ref_sets), function(rn) {
    .overlap(en, est_sets[[en]], rn, ref_sets[[rn]])
  })
}) %>%
  dplyr::mutate(estimator_version = "raw", .after = estimator)

# --- Strip against the full MitoCarta inventory ------------------------------
strip_set <- ref_sets$MITOCARTA_ALL

est_stripped <- purrr::map(est_sets, ~ setdiff(.x, strip_set))

strip_summary <- tibble::tibble(
  estimator   = names(est_sets),
  n_before    = purrr::map_int(est_sets, length),
  n_removed   = purrr::map_int(est_sets, ~ length(intersect(.x, strip_set))),
  n_remaining = purrr::map_int(est_stripped, length)
) %>%
  dplyr::mutate(
    frac_removed = n_removed / n_before,
    passes_g1_count_criterion = n_remaining >= G1_MIN_SIGNATURE_SIZE
  )

# Stripped estimators against the Hallmark proliferation sets: does removing the
# mitochondrial genes change how proliferation-entangled the estimator is?
audit_stripped <- purrr::map_dfr(names(est_stripped), function(en) {
  purrr::map_dfr(
    c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT"),
    function(rn) .overlap(en, est_stripped[[en]], rn, ref_sets[[rn]])
  )
}) %>%
  dplyr::mutate(estimator_version = "mitocarta_stripped", .after = estimator)

audit <- dplyr::bind_rows(audit_raw, audit_stripped)

# =============================================================================
# 5. Report
# =============================================================================

message("\n04: === G1 OVERLAP AUDIT ===")
audit %>%
  dplyr::select(estimator, estimator_version, reference,
                n_estimator, n_reference, n_overlap, frac_of_estimator) %>%
  print(n = Inf)

message("\n04: === MITOCARTA STRIP (full 1,136-gene inventory) ===")
strip_summary %>% print(n = Inf)

message("\n04: === G1 COUNT CRITERION ===")
for (i in seq_len(nrow(strip_summary))) {
  r <- strip_summary[i, ]
  message(sprintf(
    "  %-18s %3d -> %3d genes (%d removed)  %s",
    r$estimator, r$n_before, r$n_remaining, r$n_removed,
    if (isTRUE(r$passes_g1_count_criterion)) {
      paste0("PASSES (>= ", G1_MIN_SIGNATURE_SIZE, ")")
    } else {
      paste0("FAILS (< ", G1_MIN_SIGNATURE_SIZE, ")")
    }
  ))
}
message(
  "\n  Count criterion only. The correlation-structure half of G1 (plan section 6)\n",
  "  needs expression data from script 01 and is NOT answered here. G1 is not\n",
  "  discharged until that runs."
)

# =============================================================================
# 6. Save
# =============================================================================

g1 <- list(
  audit                = audit,
  strip_summary        = strip_summary,
  estimators_raw       = est_sets,
  estimators_stripped  = est_stripped,   # script 06 consumes these
  reference_sets       = ref_sets,
  harmonisation_report = harmonisation_report,
  meta = list(
    date                  = Sys.Date(),
    strip_set             = "MITOCARTA_ALL (full 1,136-gene human inventory)",
    g1_min_signature_size = G1_MIN_SIGNATURE_SIZE,
    correlation_criterion = "NOT EVALUATED - needs expression data, script 01",
    mitocarta_md5         = "3c0bd24e362238216e142bc708e41286",
    collectri_sha256      = paste0("d72660703f7ccb8068b75994a1e74986",
                                   "b451d0c098e36d62efd7a88e631d287d"),
    felsher_blob_v1_0     = "193a18f58cc9869acae6d2dcdcd70d2bed71e861"
  )
)

saveRDS(g1, file.path(DIR_RESULTS, "g1_overlap_audit.rds"))
readr::write_csv(audit, file.path(DIR_TABLES, "g1_overlap_audit.csv"))
readr::write_csv(harmonisation_report,
                 file.path(DIR_TABLES, "g1_symbol_harmonisation.csv"))

message("\n04: written results/g1_overlap_audit.rds")
message("04: written outputs/tables/g1_overlap_audit.csv")
message("04: written outputs/tables/g1_symbol_harmonisation.csv")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  # --- The headline numbers ---------------------------------------------------
  audit %>%
    dplyr::filter(estimator == "FELSHER", estimator_version == "raw") %>%
    dplyr::select(reference, n_overlap, frac_of_estimator) %>%
    print(n = Inf)

  # Which Felsher genes are mitochondrial? Expect 6, including KARS1 - the one
  # a naive symbol join misses.
  intersect(est_sets$FELSHER, ref_sets$MITOCARTA_ALL)

  # --- Did symbol drift actually bite? ---------------------------------------
  harmonisation_report %>%
    dplyr::filter(set == "FELSHER", status != "matched") %>%
    print(n = Inf)

  # Anything the alias map could not settle, across every set. These are kept
  # under their original symbol and cannot be checked against MitoCarta.
  harmonisation_report %>%
    dplyr::filter(status %in% c("ambiguous", "unresolved")) %>%
    dplyr::count(set, status) %>%
    print(n = Inf)

  # --- OXPHOS circularity, the number that actually matters -------------------
  # Overlap with the SUBUNIT set is what would make a MYC-to-OXPHOS result
  # partly definitional. The umbrella set includes assembly factors.
  audit %>%
    dplyr::filter(reference %in% c("MITOCARTA_OXPHOS_UMBRELLA",
                                   "MITOCARTA_OXPHOS_SUBUNITS")) %>%
    dplyr::select(estimator, estimator_version, reference, n_overlap) %>%
    print(n = Inf)

  # --- Proliferation entanglement --------------------------------------------
  # Does stripping MitoCarta genes change the Hallmark overlap? If the estimator
  # is largely a proliferation readout, that is a separate problem from
  # mitochondrial circularity and stripping will not fix it.
  audit %>%
    dplyr::filter(grepl("^HALLMARK", reference)) %>%
    dplyr::select(estimator, estimator_version, reference,
                  n_overlap, frac_of_estimator) %>%
    print(n = Inf)

  # --- The two CollecTRI variants --------------------------------------------
  length(est_sets$COLLECTRI_MYC_ALL)
  length(est_sets$COLLECTRI_MYC_STIM)
  setdiff(est_sets$COLLECTRI_MYC_ALL, est_sets$COLLECTRI_MYC_STIM) %>% head() %>% print()

  # --- Sanity: reference set sizes against the READMEs ------------------------
  vapply(ref_sets, length, integer(1))

}
