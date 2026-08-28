# 07_score_mitochondrial.R
# =============================================================================
# The mitochondrial axes, the specificity battery, and its expression-matched
# null. This is the script that builds the EXPOSURE for Block C.
#
#   Instrument 1  GSVA level   on the LOG VST matrix       CO-PRIMARY
#   Instrument 2  mitoPPS      on the LINEAR matrix        CO-PRIMARY
#   Instrument 3  content      log2 sum of linear counts   descriptive only
#   Instrument 4  level_zmean  mean z of VST               descriptive only
#
# Specification: docs/2026-08-28_specificity_panel_proposal.md (status RESOLVED,
# all 8 items agreed) and plan sections 2 and 7.2 as amended 2026-08-28.
#
# =============================================================================
# THE TWO CO-PRIMARY INSTRUMENTS, AND THE RULE FOR DISAGREEMENT
# =============================================================================
# The mouse interaction (p = 0.0052) was fitted on mitoPPS, not on a level. The
# plan's original primary was a GSVA level. Those are different quantities: the
# same 87 mouse subunits on the same contrast give +0.061 unweighted, +0.201
# expression-weighted and +0.226 summed. "Both are correct. They weight
# differently."
#
#   RULE, fixed before any model is fitted: report both; claim only what BOTH
#   instruments support. An effect on one instrument alone is reported as
#   INSTRUMENT-DEPENDENT and is not a positive result.
#
# Consequence for Block F: GSVA is the instrument that travels across cohorts.
# mitoPPS is composition-dependent, so it is TCGA-internal and never compared
# numerically across cohorts or species. That asymmetry is a property of the
# instruments, not a choice.
#
# Instruments 3 and 4 exist because the mouse's own weighting question is real
# and cheap to answer here. They are DESCRIPTIVE. No claim rests on them, and
# the disagreement rule above refers to GSVA and mitoPPS only.
#
# =============================================================================
# SCALE DISCIPLINE - the most likely silent error in this repo
# =============================================================================
# GSVA wants LOG scale (VST, kcdf = "Gaussian"). mitoPPS wants LINEAR
# DESeq2-normalised counts. These are opposite requirements and the two matrices
# must not meet. Both objects assert their own `scale` field below rather than
# being trusted, and the two are checked to be on genuinely different scales
# before anything is scored.
#
# COHORT-RELATIVITY: GSVA is cohort-relative and mitoPPS is
# composition-dependent. Every TCGA sample is scored in ONE cohort-internal run
# on both. Neither set of scores is comparable to another cohort's; Block F
# meta-analyses effect estimates, never pooled scores.
#
# SPECIES: human. Human MitoCarta 3.0, human symbols. See CLAUDE.md.
#
# =============================================================================
# RUNTIME AND DISK - read before sourcing
# =============================================================================
# Section 8 (the expression-matched null) is the expensive part: 18 arms x 2,000
# matched random gene sets, GSVA-scored. Expect roughly 20-40 minutes and
# somewhere around 0.5-1 GB under results/mito_null/. The work is cached per
# arm and the section is resumable - re-sourcing skips arms whose cache file
# already exists. results/ is gitignored and regenerable; none of it needs
# backing up. For a timing pilot, set G7_NULL_NSETS low and
# G7_NULL_FORCE_REBUILD TRUE, then put both back and delete the cache before
# the real run.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
suppressPackageStartupMessages({
  library(GSVA)
})

message("\n07: mitochondrial axes, specificity battery, matched null\n",
        strrep("=", 78))

# -----------------------------------------------------------------------------
# 0. Constants and small helpers
# -----------------------------------------------------------------------------
# Pre-registered in the specificity panel proposal section 5: 2,000 draws,
# matched within 20 ventiles of mean expression. Both numbers are the mouse's
# (myc_mouse/scripts/43_substrate_specificity_and_tradeoff.R, NSET and NBIN).
G7_NULL_NSETS <- 2000L
G7_NULL_NBIN  <- 20L

# Minimum genes present in the matrix for a set to be scored at all. mitoPPS
# uses 3, matching the mouse pathway universe. GSVA uses 3 as well, because the
# per-complex sets the plan requires include Complex II (4 genes) and ATPase
# (5); scoring them is the point of having them. Anything below
# G7_SMALL_SET_FLAG is marked "small" in every output table - a 4-gene GSVA
# score is a noisy quantity and must not be read like a 100-gene one.
G7_MIN_SET_GENES  <- 3L
G7_SMALL_SET_FLAG <- 10L

G7_NULL_FORCE_REBUILD <- FALSE
DIR_MITO_NULL <- file.path(DIR_RESULTS, "mito_null")
.ensure_dir(DIR_MITO_NULL)

# The synthetic pathway holding the 13 mtDNA-encoded protein-coding genes.
# Standing convention (CLAUDE.md): never pooled with nuclear-encoded OXPHOS
# subunits, because their expression scale skews any set they sit in.
MTDNA_PATHWAY <- "mtDNA-encoded OXPHOS subunits"

.slug <- function(x) gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", x))

.rho <- function(x, y) suppressWarnings(
  stats::cor(x, y, method = "spearman", use = "complete.obs"))

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")

vst <- readRDS(file.path(DIR_RESULTS, "tcga_brca_vst.rds"))
if (!identical(vst$scale, "log_vst")) {
  stop("expected the LOG-scale VST object, got scale = '", vst$scale,
       "'. GSVA must not be given the linear matrix.", call. = FALSE)
}
E <- vst$mat

lin <- readRDS(file.path(DIR_RESULTS, "tcga_brca_linear.rds"))
if (!identical(lin$scale, "linear_deseq2_normalised")) {
  stop("expected the LINEAR object, got scale = '", lin$scale,
       "'. mitoPPS must not be given a logged matrix.", call. = FALSE)
}
L <- lin$mat

stopifnot(identical(dimnames(E), dimnames(L)))
if (!(min(L) >= 0 && max(L) > 50 * max(E))) {
  stop("the VST and LINEAR matrices are not on different scales. One of them ",
       "is not what its `scale` field claims. Stop and check script 01.",
       call. = FALSE)
}
message("   VST    : ", nrow(E), " x ", ncol(E), "  log scale, range ",
        sprintf("%.1f to %.1f", min(E), max(E)))
message("   LINEAR : ", nrow(L), " x ", ncol(L), "  linear, range ",
        sprintf("%.1f to %.0f", min(L), max(L)))

# --- Human MitoCarta 3.0, validated against the provenance README ------------
mitocarta_inventory  <- suppressWarnings(readxl::read_xls(PATH_MITOCARTA, sheet = 2))
mitocarta_background <- suppressWarnings(readxl::read_xls(PATH_MITOCARTA, sheet = 3))
mitocarta_pathways   <- suppressWarnings(readxl::read_xls(PATH_MITOCARTA, sheet = 4))

# Sheet 4's 5 all-NA padding rows are dropped at load, not defended against at
# each lookup: strsplit() turns each NA into a phantom "NA" gene and inflates
# every pathway by exactly 5, silently. See script 04.
mitocarta_pathways <- mitocarta_pathways %>% dplyr::filter(!is.na(MitoPathway))

stopifnot(
  nrow(mitocarta_inventory)  == EXPECT_MITOCARTA_GENES,
  nrow(mitocarta_background) == EXPECT_MITOCARTA_BACKGRD,
  nrow(mitocarta_pathways)   == EXPECT_MITOCARTA_PATHWAYS,
  !anyNA(mitocarta_pathways$Genes),
  !anyDuplicated(mitocarta_pathways$MitoPathway),
  all(c("Symbol", "Synonyms") %in% colnames(mitocarta_background))
)
message("   MitoCarta 3.0: ", nrow(mitocarta_inventory), " genes, ",
        nrow(mitocarta_pathways), " MitoPathways")

# --- the curated human metabolic list, for complex-level resolution ----------
# MitoCarta's "OXPHOS subunits" is one flat list. The mouse claim is that OXPHOS
# falls ACROSS ALL COMPLEXES, and testing that needs per-complex sets.
# data/genesets_metabolic_human/README.md: this is the NATIVE human sheet of a
# tracked, tag-pinned library input, not the mouse-derived outputs/gmt/human/
# tree CLAUDE.md rejects. Different objects.
gs_metab <- readr::read_csv(
  file.path(DIR_DATA, "genesets_metabolic_human", "gs_metabolic_human.csv"),
  show_col_types = FALSE, progress = FALSE
)
stopifnot(nrow(gs_metab) == 2347L,
          all(c("gene_symbol", "classification") %in% colnames(gs_metab)))
message("   GS_metabolic (human sheet): ", nrow(gs_metab), " genes, ",
        dplyr::n_distinct(gs_metab$classification), " classifications")

# --- G1 and Hallmark, for the D7 proliferation covariates --------------------
g1 <- readRDS(file.path(DIR_RESULTS, "g1_overlap_audit.rds"))
felsher_stripped <- g1$estimators_stripped$FELSHER
stopifnot(length(felsher_stripped) == 61L)

hallmark <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
hm_e2f <- unique(hallmark$gene_symbol[hallmark$gs_name == "HALLMARK_E2F_TARGETS"])
hm_g2m <- unique(hallmark$gene_symbol[hallmark$gs_name == "HALLMARK_G2M_CHECKPOINT"])
stopifnot(length(hm_e2f) == 200L, length(hm_g2m) == 200L)

# =============================================================================
# 2. Symbol harmonisation - to the MATRIX vocabulary
# =============================================================================
# Same hazard as G1, opposite direction. G1 harmonised everything TO the
# MitoCarta Sheet-3 vocabulary, because it was comparing gene sets to each
# other. Here the target is the EXPRESSION MATRIX vocabulary (GENCODE symbols
# from script 01), because sets are about to be intersected with rownames(E).
# The same rescue is needed - KARS -> KARS1 cost G1 a gene until it was caught -
# but the map runs the other way: a MitoCarta symbol absent from the matrix is
# looked up in its own Sheet-3 synonyms for a symbol that is present.
#
# Three guards, all inherited from script 04's ambiguity analysis:
#   - a synonym is rejected if it is itself the MitoCarta Symbol of a DIFFERENT
#     gene (about 1,000 aliases are);
#   - a synonym is rejected if it is already a member of the vocabulary being
#     mapped, so two annotated genes never collapse into one;
#   - if more than one synonym survives, the symbol is left alone, because there
#     is no basis for choosing.
# Anything unresolved is KEPT under its original symbol. It simply will not
# match, exactly as before, and it is reported by name rather than dropped
# silently.
#
# ONE map is built over the union of every symbol used anywhere in this script,
# and applied to each collection separately. Building it per collection would be
# wrong here: several arm names are also MitoPathway names, and the merged
# "Fatty acid oxidation" arm is NOT the "Fatty acid oxidation" pathway. Keying
# anything by set name across collections silently gives one the other's genes.
message("\n2. symbol harmonisation to the matrix vocabulary")

MATRIX_SYMBOLS <- rownames(E)
MC_SYMBOLS     <- unique(mitocarta_background$Symbol)

.syn_map <- local({
  syn <- strsplit(ifelse(is.na(mitocarta_background$Synonyms), "",
                         mitocarta_background$Synonyms), "|", fixed = TRUE)
  tibble::tibble(symbol = rep(mitocarta_background$Symbol, lengths(syn)),
                 alias  = unlist(syn)) %>%
    dplyr::filter(!is.na(alias), alias != "", alias != symbol) %>%
    dplyr::distinct()
})

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
      out[[g]]   <- cand
      status[[i]] <- "resolved"
    } else if (length(cand) > 1L) {
      status[[i]] <- "ambiguous"
    }
  }
  list(map = out,
       report = tibble::tibble(input_symbol = genes, status = status,
                               resolved_to = unname(out[genes])))
}

# =============================================================================
# 3. Gene sets
# =============================================================================
message("\n3. gene sets")

.split_genes <- function(x) trimws(unlist(strsplit(x, ",", fixed = TRUE)))

# --- 3a. the MitoCarta pathway universe, with mtDNA split out ----------------
# Faithful to myc_mouse/scripts/08_mitoPPS_analysis.R PART 2b and to Monzel et
# al.: mtDNA-encoded transcripts are orders of magnitude more abundant than
# nuclear mitochondrial ones, so any pathway containing even one of them has a
# severely inflated score. They are removed from their canonical pathways and
# held in one synthetic pathway of their own.
mito_paths_raw <- stats::setNames(
  lapply(mitocarta_pathways$Genes, .split_genes),
  mitocarta_pathways$MitoPathway
)
stopifnot(!MTDNA_PATHWAY %in% names(mito_paths_raw))

mtdna_genes <- grep("^MT-", unique(mitocarta_inventory$Symbol), value = TRUE)
stopifnot(length(mtdna_genes) == EXPECT_MITOCARTA_MTDNA)

mtdna_in_sheet4 <- grep("^MT-", unique(unlist(mito_paths_raw)), value = TRUE)
if (!all(mtdna_in_sheet4 %in% mtdna_genes)) {
  stop("Sheet 4 carries MT- symbols absent from the Sheet 2 inventory: ",
       paste(setdiff(mtdna_in_sheet4, mtdna_genes), collapse = ", "),
       call. = FALSE)
}
n_paths_hit <- sum(vapply(mito_paths_raw,
                          function(g) any(g %in% mtdna_genes), logical(1)))
message("   mtDNA-encoded protein-coding genes: ", length(mtdna_genes),
        ", present in ", n_paths_hit, " MitoPathways -> relocated to '",
        MTDNA_PATHWAY, "'")

mito_paths <- lapply(mito_paths_raw, function(g) setdiff(g, mtdna_genes))
mito_paths[[MTDNA_PATHWAY]] <- mtdna_genes

# --- 3b. the specificity battery ---------------------------------------------
# Each arm answers a NAMED alternative to H1 (plan section 2 as amended). The
# arms are the ones the mouse actually ran, taken from Human MitoCarta by
# pathway NAME - each species uses its own native MitoCarta, so nothing is
# ortholog-projected. The cross-species link is the name.
#
# `OXPHOS assembly factors` is the PRIMARY negative, not a leftover: same
# complexes, same umbrella, different function. In the mouse the subunits sit at
# percentile 0.0 of expression-matched sets while their own assembly factors sit
# at 50.2. A distant pathway cannot exclude what that pair excludes.
ARM_PATHWAYS <- list(
  "OXPHOS subunits"          = "OXPHOS subunits",
  "OXPHOS umbrella"          = "OXPHOS",
  "OXPHOS assembly factors"  = "OXPHOS assembly factors",
  "Mitochondrial ribosome"   = "Mitochondrial ribosome",
  "Nucleotide metabolism"    = "Nucleotide metabolism",
  "ROS and glutathione"      = "ROS and glutathione metabolism",
  "TCA cycle"                = "TCA cycle",
  "Amino acid metabolism"    = "Amino acid metabolism",
  "Lipid metabolism"         = "Lipid metabolism",
  "Fatty acid oxidation"     = c("Fatty acid oxidation", "Carnitine shuttle"),
  "Folate and 1-C"           = "Folate and 1-C metabolism",
  "Glycine metabolism"       = "Glycine metabolism",
  "mtDNA-encoded OXPHOS"     = MTDNA_PATHWAY
)

arms_spec <- tibble::tibble(
  arm      = names(ARM_PATHWAYS),
  pathways = vapply(ARM_PATHWAYS, paste, character(1), collapse = " + "),
  role = c("claim", "reported", "primary negative", "negative", "negative",
           "negative", "negative", "negative", "negative", "negative",
           "negative", "negative", "held separate"),
  alternative_answered = c(
    "-",
    "includes assembly factors; not the primary measure",
    "biogenesis, not respiration",
    "growth / translation capacity",
    "proliferation",
    "redox, not respiration",
    "generic mitochondrial metabolism",
    "generic mitochondrial metabolism",
    "generic mitochondrial metabolism",
    "FAO - see the Lee et al. caveat in this section",
    "one-carbon metabolism",
    "one-carbon metabolism",
    "standing convention, expression-scale skew")
)

# THE FAO ARM CARRIES A PRE-STATED CAVEAT, fixed before the result is seen.
# Lee et al. 2017 (Cell Metab 26:633) report that TNBC cancer stem cells are
# preferentially FAO-dependent. If the FAO arm fires, the honest reading is
# "consistent with Lee et al.", NOT "the OXPHOS coupling is non-specific". The
# same applies more weakly to one-carbon via MTHFD2 (G1 note section 7).
# `Carnitine shuttle` alone is 5 genes and unscoreable as a negative, so it is
# merged into `Fatty acid oxidation` and the union is reported as one arm.

missing_paths <- setdiff(unlist(ARM_PATHWAYS), names(mito_paths))
if (length(missing_paths)) {
  stop("MitoPathway name(s) not found in Human MitoCarta 3.0 Sheet 4: ",
       paste(missing_paths, collapse = ", "),
       "\nA wrong name gives a silent empty set. Fix the name, do not guess.",
       call. = FALSE)
}

# --- 3c. per-complex sets, from the curated human list -----------------------
COMPLEX_CLASSES <- c("Complex I", "Complex II", "Complex III", "Complex IV", "ATPase")
missing_cls <- setdiff(COMPLEX_CLASSES, unique(gs_metab$classification))
if (length(missing_cls)) {
  stop("classification(s) absent from gs_metabolic_human.csv: ",
       paste(missing_cls, collapse = ", "), call. = FALSE)
}
complex_raw <- stats::setNames(
  lapply(COMPLEX_CLASSES,
         function(cl) unique(gs_metab$gene_symbol[gs_metab$classification == cl])),
  paste0("Complex-level: ", COMPLEX_CLASSES))
# The mtDNA strip is not a formality here. The source curation is internally
# inconsistent about mtDNA subunits - Complex I lists all seven MT-ND genes and
# Complex III lists MT-CYB, while Complex IV lists none of MT-CO1/2/3 and ATPase
# lists neither MT-ATP6 nor MT-ATP8. Before the strip the five sets are
# therefore not comparable to each other at all; after it they are uniformly
# nuclear-encoded, which is also what the standing convention requires.
complex_sets <- lapply(complex_raw, setdiff, y = mtdna_genes)
message("   per-complex sets, before -> after the mtDNA strip:")
for (nm in names(complex_sets)) {
  message(sprintf("     %-28s %3d -> %3d", nm,
                  length(complex_raw[[nm]]), length(complex_sets[[nm]])))
}

# Cross-check against MitoCarta, which is an independent curation of the same
# biology. The five complexes should very nearly exhaust the nuclear-encoded
# OXPHOS subunits, and a large disagreement means one of the two lists is not
# what this script assumes it is.
cx_union <- unique(unlist(complex_sets, use.names = FALSE))
ox_nuc   <- setdiff(mito_paths[["OXPHOS subunits"]], mtdna_genes)
message(sprintf(
  paste("   complexes pooled: %d genes vs MitoCarta nuclear OXPHOS subunits",
        "%d;\n     %d shared, %d curated-only, %d MitoCarta-only"),
  length(cx_union), length(ox_nuc), length(intersect(cx_union, ox_nuc)),
  length(setdiff(cx_union, ox_nuc)), length(setdiff(ox_nuc, cx_union))))
if (length(intersect(cx_union, ox_nuc)) < 0.7 * length(ox_nuc)) {
  stop("the per-complex sets and MitoCarta's OXPHOS subunits agree on fewer ",
       "than 70% of genes. One of the two lists is not what it is taken to be; ",
       "do not proceed on a per-complex claim.", call. = FALSE)
}

# --- 3d. the D7 proliferation covariates -------------------------------------
# PROLIF_DISJOINT is defined against the MitoCarta-STRIPPED Felsher set and is
# specific to M-a. It must NEVER be used with M-b: a different estimator implies
# a different disjoint set, and reusing this one silently reintroduces the
# circularity it exists to remove. See docs/2026-08-28_D7_proliferation_covariate.md
PROLIF_STD      <- union(hm_e2f, hm_g2m)
PROLIF_SHARED   <- intersect(PROLIF_STD, felsher_stripped)
PROLIF_DISJOINT <- setdiff(PROLIF_STD, felsher_stripped)

D7_SHARED_NINE <- c("CTPS1", "DNMT1", "HMGA1", "NCL", "NOP56", "PRMT5",
                    "RANBP1", "TFDP1", "UCK2")
stopifnot(length(PROLIF_STD) == 327L, length(PROLIF_DISJOINT) == 318L,
          setequal(PROLIF_SHARED, D7_SHARED_NINE))
message("   PROLIF_STD ", length(PROLIF_STD), ", PROLIF_DISJOINT ",
        length(PROLIF_DISJOINT), " (D7's 9 shared genes confirmed by name)")

# --- 3e. assemble and harmonise ----------------------------------------------
arm_sets <- c(
  stats::setNames(lapply(ARM_PATHWAYS,
                         function(p) unique(unlist(mito_paths[p], use.names = FALSE))),
                  names(ARM_PATHWAYS)),
  complex_sets
)
covariate_sets <- list(PROLIF_STD = PROLIF_STD, PROLIF_DISJOINT = PROLIF_DISJOINT)

sym <- .build_symbol_map(c(unlist(arm_sets, use.names = FALSE),
                           unlist(covariate_sets, use.names = FALSE),
                           unlist(mito_paths, use.names = FALSE)))
.remap <- function(g) unique(unname(sym$map[g]))

arm_sets       <- lapply(arm_sets, .remap)
covariate_sets <- lapply(covariate_sets, .remap)
mito_paths     <- lapply(mito_paths, .remap)

harm_report <- sym$report %>%
  dplyr::mutate(in_sets = purrr::map_chr(resolved_to, function(s) {
    hit <- names(arm_sets)[vapply(arm_sets, function(g) s %in% g, logical(1))]
    paste(hit, collapse = "; ")
  }))
message("   symbol rescue over ", nrow(harm_report), " distinct symbols: ",
        sum(harm_report$status == "resolved"), " resolved via MitoCarta synonyms, ",
        sum(harm_report$status == "ambiguous"), " ambiguous, ",
        sum(harm_report$status == "unresolved"), " unresolved (all kept as-is)")

.in_matrix <- function(g) intersect(g, MATRIX_SYMBOLS)

all_scored <- c(arm_sets, covariate_sets)
coverage <- tibble::tibble(
  set       = names(all_scored),
  kind      = c(rep("arm", length(arm_sets)), rep("covariate", length(covariate_sets))),
  n_defined = vapply(all_scored, length, integer(1)),
  n_matched = vapply(all_scored, function(g) length(.in_matrix(g)), integer(1))
) %>%
  dplyr::mutate(frac  = round(n_matched / n_defined, 3),
                small = n_matched < G7_SMALL_SET_FLAG)

message("\n   set coverage in the expression matrix:")
coverage %>% as.data.frame() %>% print(row.names = FALSE)

too_small <- coverage$set[coverage$n_matched < G7_MIN_SET_GENES]
if (length(too_small)) {
  stop("set(s) with fewer than ", G7_MIN_SET_GENES,
       " genes in the matrix, cannot be scored: ",
       paste(too_small, collapse = ", "), call. = FALSE)
}
if (any(coverage$small)) {
  message("   NOTE small sets (< ", G7_SMALL_SET_FLAG, " genes): ",
          paste(coverage$set[coverage$small], collapse = ", "))
  message("        Scored because the plan requires per-complex resolution, ",
          "but a 4-gene\n        GSVA score is noisy and is flagged in every ",
          "output table.")
}

# =============================================================================
# 4. INSTRUMENT 1 (co-primary): GSVA level on the LOG VST matrix
# =============================================================================
# kcdf = "Gaussian" because the input is log-scale, per CLAUDE.md.
#
# WHY EVERY CALL CARRIES TWO PIN SETS. GSVA's enrichment statistic is a random
# walk over the genes of the expression matrix it is handed, and GSVA may
# restrict that matrix to genes appearing in the supplied gene sets. If it does,
# a call with 20 pathways and a call with 2,000 random sets walk over DIFFERENT
# gene universes and their scores are not comparable - which would silently
# invalidate the entire matched null. Two complementary half-matrix sets pin the
# universe to the full matrix in every call. (Two halves rather than one
# all-genes set: a set containing every gene leaves no genes outside it, and the
# walk's miss-penalty term is then 0/0.) Their own scores are discarded.
#
# The pin is belt; the braces is the .OBS_CHECK assertion in section 8, which
# re-scores the observed arm inside its own null batch and stops if the value
# moves at all.
message("\n4. INSTRUMENT 1: GSVA level (log VST, kcdf Gaussian)")

.PIN_A <- MATRIX_SYMBOLS[c(TRUE, FALSE)]
.PIN_B <- MATRIX_SYMBOLS[c(FALSE, TRUE)]
stopifnot(setequal(c(.PIN_A, .PIN_B), MATRIX_SYMBOLS))

.gsva_batch <- function(sets, label = "") {
  present <- lapply(sets, .in_matrix)
  n_ok    <- vapply(present, length, integer(1))
  if (any(n_ok < G7_MIN_SET_GENES)) {
    stop("GSVA batch '", label, "': set(s) below the size floor -> ",
         paste(names(present)[n_ok < G7_MIN_SET_GENES], collapse = ", "),
         call. = FALSE)
  }
  wanted <- names(present)
  present[[".PIN_A"]] <- .PIN_A
  present[[".PIN_B"]] <- .PIN_B

  par <- GSVA::gsvaParam(exprData = E, geneSets = present, kcdf = "Gaussian",
                         minSize = G7_MIN_SET_GENES, maxSize = Inf)
  s <- GSVA::gsva(par, verbose = FALSE)

  dropped <- setdiff(wanted, rownames(s))
  if (length(dropped)) {
    stop("GSVA silently dropped set(s) in batch '", label, "': ",
         paste(utils::head(dropped, 10), collapse = ", "),
         if (length(dropped) > 10) " ..." else "", call. = FALSE)
  }
  s[wanted, , drop = FALSE]
}

gsva_all  <- .gsva_batch(all_scored, "arms + covariates")
gsva_arms <- gsva_all[names(arm_sets), , drop = FALSE]
gsva_cov  <- gsva_all[names(covariate_sets), , drop = FALSE]
message("   scored ", nrow(gsva_all), " sets x ", ncol(gsva_all), " samples")

# =============================================================================
# 5. INSTRUMENT 2 (co-primary): mitoPPS on the LINEAR matrix
# =============================================================================
# The algorithm, verified against Monzel et al.'s own code
# (external/mitotyping/Code/Figure5/mitoPPS_RM.R) and the mouse's implementation
# (myc_mouse/scripts/08_mitoPPS_analysis.R PART 5):
#
#   1. pathway score  = MEAN LINEAR expression of member genes, per sample
#   2. ratio_pq_i     = score_p_i / score_q_i, for every ordered pair p != q
#   3. corrected      = ratio_pq_i / mean_over_samples(ratio_pq)
#   4. mitoPPS_p_i    = mean over q != p of corrected_pq_i
#
# mitoPPS ~ 1 means the pathway is prioritised at the dataset average; > 1
# up-prioritised; < 1 down-prioritised. It reports the SHAPE of the
# mitochondrial program and is deliberately blind to total mitochondrial
# content, which is exactly why an OXPHOS LEVEL is reported alongside it.
#
# The loops in the reference code are written here as matrix algebra - the tidy
# version materialises P x P x N rows, which at ~150 pathways and 1,095 samples
# is 24 million. The identity is exact, not an approximation, and section 5.3
# asserts it against the definitional form.
message("\n5. INSTRUMENT 2: mitoPPS composition (linear DESeq2-normalised)")

.pathway_score_matrix <- function(sets, min_genes = G7_MIN_SET_GENES) {
  present <- lapply(sets, .in_matrix)
  present <- present[vapply(present, length, integer(1)) >= min_genes]
  out <- t(vapply(present, function(g) colMeans(L[g, , drop = FALSE]),
                  numeric(ncol(L))))
  rownames(out) <- names(present)
  colnames(out) <- colnames(L)
  out
}

# --- 5.1 the pathway universe ------------------------------------------------
# All MitoCarta pathways with at least 3 genes in the matrix, mtDNA split out.
# The mouse universe exactly.
#
# >>> DOCUMENTED DEVIATION FROM THE MOUSE INSTRUMENT, NEEDS SIGN-OFF.
# The mouse added two synthetic pathways to its universe, Apoptosis-PRO and
# Apoptosis-ANTI, from a hand-curated mouse gene list. They are NOT added here,
# for two independent reasons, either of which is sufficient:
#   (1) those lists contain Bbc3 and Bcl2l1, whose human orthologs are the
#       numerator and denominator of PRIME. Putting the endpoint's own genes
#       into the construction of the exposure is the circularity gate G1 exists
#       to prevent.
#   (2) they are a MOUSE-curated gene set. Translating them by uppercasing
#       symbols is exactly the ortholog projection CLAUDE.md forbids here, and
#       it would look correct.
# MitoCarta's own human-native "Apoptosis" pathway stays in the universe, so the
# compartment is still represented. The consequence is that the human and mouse
# mitoPPS universes differ by two pathways - which is one more reason mitoPPS
# values are never compared numerically across species, only their pattern.
S_universe <- .pathway_score_matrix(mito_paths)
if (min(S_universe) <= 0) {
  stop("a MitoPathway score is zero or negative in some sample; the pairwise ",
       "ratio is undefined. Inspect before proceeding.", call. = FALSE)
}
message("   universe: ", nrow(S_universe), " MitoPathways with >= ",
        G7_MIN_SET_GENES, " genes in the matrix (of ", length(mito_paths), ")")
message("   MitoCarta's native 'Apoptosis' pathway in the universe: ",
        "Apoptosis" %in% rownames(S_universe))

# --- 5.2 the algorithm -------------------------------------------------------
.mitopps_universe <- function(S) {
  N <- ncol(S); P <- nrow(S)
  Bi <- 1 / S
  A  <- (S %*% t(Bi)) / N                # A[p,q] = mean_i S[p,i] / S[q,i]
  M  <- (1 / A) %*% Bi
  out <- S * (M - Bi) / (P - 1)          # A[p,p] == 1, so subtracting Bi drops q = p
  dimnames(out) <- dimnames(S)
  out
}

# One or many query sets against a FIXED universe that does not contain them.
# Same formula with the self-term simply absent, which is what makes an observed
# arm and its matched null exactly comparable: both go through here, against the
# same universe, with the same denominator.
.mitopps_query <- function(Sq, Su) {
  stopifnot(is.matrix(Sq), is.matrix(Su), ncol(Sq) == ncol(Su))
  if (min(Sq) <= 0 || min(Su) <= 0) {
    stop("non-positive pathway score reached .mitopps_query", call. = FALSE)
  }
  N <- ncol(Su); P <- nrow(Su)
  Bi  <- 1 / Su
  A   <- (Sq %*% t(Bi)) / N
  out <- Sq * ((1 / A) %*% Bi) / P
  dimnames(out) <- list(rownames(Sq), colnames(Sq))
  out
}

mitopps_universe <- .mitopps_universe(S_universe)
message(sprintf("   mitoPPS global mean %.4f (should be ~1), range %.3f to %.3f",
                mean(mitopps_universe), min(mitopps_universe), max(mitopps_universe)))

# --- 5.3 the identity check --------------------------------------------------
# .mitopps_universe is an algebraic rewrite of the reference implementation's
# triple loop. If the rewrite is wrong, every mitoPPS number in this project is
# wrong and nothing downstream would notice. Checked against the definitional
# form below, which is steps 2-4 transcribed literally.
.mitopps_definitional <- function(S, p) {
  q  <- setdiff(seq_len(nrow(S)), p)
  R  <- S[rep(p, length(q)), , drop = FALSE] / S[q, , drop = FALSE]  # step 2
  Rc <- R / rowMeans(R)                                              # step 3
  colMeans(Rc)                                                       # step 4
}
for (p in unique(c(1L, unname(which.max(rowMeans(S_universe))), nrow(S_universe)))) {
  d <- max(abs(mitopps_universe[p, ] - .mitopps_definitional(S_universe, p)))
  if (!is.finite(d) || d > 1e-8) {
    stop("mitoPPS matrix rewrite disagrees with the definitional form for '",
         rownames(S_universe)[p], "' (max abs diff ", signif(d, 3), ").",
         call. = FALSE)
  }
}
message("   matrix rewrite verified against the definitional form (3 pathways)")

# --- 5.4 the arms ------------------------------------------------------------
# An arm that IS a universe pathway is scored against the universe with that
# pathway held out - which is the canonical mitoPPS, since the canonical
# definition already averages over q != p. A derived arm (the merged FAO set,
# the five per-complex sets) is scored as a query against the full universe.
# Nested pathways stay in the universe in both cases, because they do in Monzel
# and in the mouse; this instrument's whole point is fidelity to the one the
# mouse interaction was fitted on.
S_arms <- .pathway_score_matrix(arm_sets)
stopifnot(setequal(rownames(S_arms), names(arm_sets)))

arm_universe_path <- vapply(names(arm_sets), function(a) {
  p <- ARM_PATHWAYS[[a]]
  if (!is.null(p) && length(p) == 1L && p %in% rownames(S_universe)) p
  else NA_character_
}, character(1))

.arm_universe <- function(a) {
  hold <- arm_universe_path[[a]]
  if (is.na(hold)) S_universe
  else S_universe[setdiff(rownames(S_universe), hold), , drop = FALSE]
}

mitopps_arms <- t(vapply(names(arm_sets), function(a)
  as.numeric(.mitopps_query(S_arms[a, , drop = FALSE], .arm_universe(a))),
  numeric(ncol(L))))
colnames(mitopps_arms) <- colnames(L)

# An arm that is a universe pathway must reproduce its canonical value exactly.
for (a in names(arm_sets)[!is.na(arm_universe_path)]) {
  d <- max(abs(mitopps_arms[a, ] - mitopps_universe[arm_universe_path[[a]], ]))
  if (d > 1e-8) {
    stop("arm '", a, "' does not reproduce its canonical mitoPPS (diff ",
         signif(d, 3), ")", call. = FALSE)
  }
}
message("   ", sum(!is.na(arm_universe_path)), " of ", length(arm_sets),
        " arms are universe pathways and reproduce their canonical values; ",
        "the rest are scored as queries")

# =============================================================================
# 6. INSTRUMENTS 3 and 4 (DESCRIPTIVE ONLY - no claim rests on these)
# =============================================================================
# The mouse's own weighting question, asked in human. figS2_oxphos_subunit_
# heatmap.R gets +0.061 unweighted, +0.201 expression-weighted and +0.226 from
# log2(sum of normalised counts) on the SAME 87 subunits and the SAME contrast.
# The sign of a headline number can depend on the weighting, so both extra
# rulers are computed and reported. They are NOT co-primary.
message("\n6. INSTRUMENTS 3 and 4 (descriptive): content and mean-z level")

content_arms <- t(vapply(arm_sets, function(g)
  log2(colSums(L[.in_matrix(g), , drop = FALSE]) + 1), numeric(ncol(L))))
colnames(content_arms) <- colnames(L)

.zmean <- function(g) {
  sub <- E[.in_matrix(g), , drop = FALSE]
  v   <- apply(sub, 1L, stats::var)
  sub <- sub[v > 0, , drop = FALSE]     # a zero-variance row makes scale() NaN
  colMeans(t(scale(t(sub))))
}
zmean_arms <- t(vapply(arm_sets, .zmean, numeric(ncol(E))))
colnames(zmean_arms) <- colnames(E)

# =============================================================================
# 7. Diagnostics, before any model exists
# =============================================================================
message("\n7. diagnostics")

myc <- readRDS(file.path(DIR_RESULTS, "tcga_brca_myc_scores.rds"))$estimators
cov <- readRDS(file.path(DIR_RESULTS, "tcga_brca_covariates.rds"))$covariates
stopifnot(identical(sort(cov$patient), sort(colnames(E))))
i_cov <- match(colnames(E), cov$patient)
i_myc <- match(colnames(E), myc$patient)

arm_summary <- tibble::tibble(
  arm       = rownames(gsva_arms),
  n_genes   = coverage$n_matched[match(rownames(gsva_arms), coverage$set)],
  small     = coverage$small[match(rownames(gsva_arms), coverage$set)],
  gsva_mean = rowMeans(gsva_arms),
  gsva_sd   = apply(gsva_arms, 1L, stats::sd),
  mitopps_mean = rowMeans(mitopps_arms),
  mitopps_sd   = apply(mitopps_arms, 1L, stats::sd),
  rho_instruments = purrr::map_dbl(rownames(gsva_arms),
                                   ~ .rho(gsva_arms[.x, ], mitopps_arms[.x, ])),
  rho_M_a    = purrr::map_dbl(rownames(gsva_arms),
                              ~ .rho(gsva_arms[.x, ], myc$M_a[i_myc])),
  rho_purity = purrr::map_dbl(rownames(gsva_arms),
                              ~ .rho(gsva_arms[.x, ], cov$purity[i_cov])),
  rho_leuko  = purrr::map_dbl(rownames(gsva_arms),
                              ~ .rho(gsva_arms[.x, ], cov$leukocyte_fraction[i_cov]))
)

message("\n   arm summary (rho_instruments is the GSVA-vs-mitoPPS agreement):")
arm_summary %>%
  dplyr::select(arm, n_genes, rho_instruments, rho_M_a, rho_purity, rho_leuko) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# The two instruments are NOT expected to agree closely - they measure level and
# composition, and the specificity panel proposal says in as many words that the
# same genes can move in opposite directions on them. A low or negative
# correlation here is information about the instruments, not a bug. It is
# printed because the "claim only what both support" rule needs it visible
# BEFORE any coefficient is seen.
message(sprintf(
  "\n   OXPHOS subunits, GSVA vs mitoPPS: rho = %.3f",
  arm_summary$rho_instruments[arm_summary$arm == "OXPHOS subunits"]))

# Breast is the worst tissue in TCGA for mitochondrial confounding - adipose is
# OXPHOS and FAO high, immune infiltrate carries its own BCL2-family profile.
# Purity and leukocyte fraction are Block C covariates for exactly this reason;
# rho_purity and rho_leuko above say how hard those covariates will have to work.

# =============================================================================
# 8. The expression-matched null - 2,000 matched random sets per arm
# =============================================================================
# WHY THIS EXISTS. Without it, "OXPHOS moved and TCA did not" is confounded by
# set size and by expression level. The mouse compares each arm against 2,000
# random gene sets drawn within baseMean ventiles, and that is what turns a
# ranking into a test: OXPHOS subunits at percentile 0.0, their own assembly
# factors at 50.2.
#
# WHAT IS DIFFERENT IN HUMAN. In the mouse the arm is an OUTCOME - a set of
# genes carrying a fold change - so its null is a cheap re-average of per-gene
# statistics. Here the arm is the EXPOSURE: it enters a model. So the null has
# to be SCORED, on both instruments, and the model refit against each draw. This
# script does the scoring. Script 09 fits the models and reports percentiles.
#
# Matching: genes are binned into 20 ventiles of mean LINEAR expression (the
# baseMean analogue, the mouse's variable) and each draw takes the same number
# of genes from each ventile as the arm has there. Drawn WITHOUT replacement
# within a ventile - the mouse drew with replacement, which is harmless for a
# mean of per-gene statistics but would silently shrink a GSVA set, since GSVA
# de-duplicates set members.
#
# One draw serves BOTH instruments, so the two percentiles are calibrated
# against the same random sets.
message("\n8. expression-matched null, ", G7_NULL_NSETS, " draws per arm")

gene_mean_lin <- rowMeans(L)
gene_mean_vst <- rowMeans(E)
rho_means <- .rho(gene_mean_lin, gene_mean_vst)
message(sprintf("   Spearman(mean linear, mean VST) across genes = %.4f", rho_means))
if (rho_means < 0.95) {
  warning("the linear and VST expression orderings disagree (rho = ",
          round(rho_means, 3), "). Ventiles matched on one are not matched on ",
          "the other, so the null is calibrated for one instrument only.",
          call. = FALSE)
}

bin_of <- stats::setNames(
  cut(rank(gene_mean_lin, ties.method = "first"),
      breaks = G7_NULL_NBIN, labels = FALSE),
  rownames(L))
by_bin <- split(names(bin_of), bin_of)

.draw_matched <- function(genes) {
  b <- bin_of[genes]
  b <- b[!is.na(b)]
  unlist(lapply(split(b, b), function(k) {
    pool <- by_bin[[as.character(k[[1]])]]
    if (length(k) > length(pool)) {
      stop("ventile ", k[[1]], " holds ", length(pool), " genes but ",
           length(k), " are needed. Reduce G7_NULL_NBIN.", call. = FALSE)
    }
    pool[sample.int(length(pool), length(k))]
  }), use.names = FALSE)
}

null_arms <- rownames(gsva_arms)

null_manifest <- dplyr::bind_rows(lapply(seq_along(null_arms), function(ai) {
  a    <- null_arms[[ai]]
  f    <- file.path(DIR_MITO_NULL, paste0("null_", .slug(a), ".rds"))
  seed <- PROJECT_SEED + 1000L + ai

  if (file.exists(f) && !G7_NULL_FORCE_REBUILD) {
    n <- readRDS(f)
    message(sprintf("   [%2d/%2d] %-26s cached (%d draws)",
                    ai, length(null_arms), a, nrow(n$gsva)))
    return(tibble::tibble(arm = a, file = f, n_sets = nrow(n$gsva),
                          n_genes = length(n$genes), seed = seed, rebuilt = FALSE))
  }

  t0    <- Sys.time()
  genes <- .in_matrix(arm_sets[[a]])

  # Seeded PER ARM, not once for the whole loop, so a resumed or partially
  # cached run reproduces exactly the same draws for a given arm regardless of
  # which arms ran before it.
  set.seed(seed)
  null_names <- sprintf("null_%04d", seq_len(G7_NULL_NSETS))
  draws <- stats::setNames(
    lapply(seq_len(G7_NULL_NSETS), function(i) .draw_matched(genes)), null_names)
  stopifnot(all(vapply(draws, length, integer(1)) == length(genes)))

  # --- GSVA, one call per arm, with the observed arm carried along as a check
  sets <- draws
  sets[[".OBS_CHECK"]] <- genes
  g <- .gsva_batch(sets, label = paste("null:", a))

  d <- max(abs(g[".OBS_CHECK", ] - gsva_arms[a, ]))
  if (d > 1e-8) {
    stop("GSVA is not independent of the set collection: the observed arm '", a,
         "' scored differently inside its own null batch (max abs diff ",
         signif(d, 3), "). The null is not comparable to the observed value. ",
         "Stop and diagnose before using any percentile.", call. = FALSE)
  }
  g <- g[null_names, , drop = FALSE]

  # --- mitoPPS, against exactly the universe the observed arm was scored on
  Sq <- t(vapply(draws, function(gg) colMeans(L[gg, , drop = FALSE]),
                 numeric(ncol(L))))
  dimnames(Sq) <- list(null_names, colnames(L))
  m <- .mitopps_query(Sq, .arm_universe(a))

  saveRDS(list(arm = a, genes = genes, seed = seed, draws = draws,
               gsva = g, mitopps = m, built = Sys.time()), f)

  message(sprintf("   [%2d/%2d] %-26s %4d draws, %3d genes, %.1f min",
                  ai, length(null_arms), a, G7_NULL_NSETS, length(genes),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  tibble::tibble(arm = a, file = f, n_sets = G7_NULL_NSETS,
                 n_genes = length(genes), seed = seed, rebuilt = TRUE)
}))

message("   null cache: ", nrow(null_manifest), " arms under results/mito_null/")

# =============================================================================
# 9. Save
# =============================================================================
# The null score matrices stay in their per-arm cache files (0.5-1 GB in total)
# rather than being bundled here, so this object stays small enough to
# load casually. Script 09 walks `null_manifest`.
message("\n9. save")

.wide <- function(m, prefix) {
  d <- as.data.frame(t(m), stringsAsFactors = FALSE, check.names = FALSE)
  colnames(d) <- paste0(prefix, .slug(colnames(d)))
  rownames(d) <- NULL
  tibble::as_tibble(d)
}

scores <- dplyr::bind_cols(
  tibble::tibble(patient = colnames(E)),
  .wide(gsva_arms,    "gsva__"),
  .wide(mitopps_arms, "mitopps__"),
  .wide(content_arms, "content__"),
  .wide(zmean_arms,   "zmean__"),
  .wide(gsva_cov,     "")
)

mito <- list(
  scores           = scores,
  gsva_arms        = gsva_arms,
  mitopps_arms     = mitopps_arms,
  content_arms     = content_arms,
  zmean_arms       = zmean_arms,
  gsva_cov         = gsva_cov,
  mitopps_universe = mitopps_universe,
  pathway_scores   = S_universe,
  arms             = arms_spec,
  arm_pathways     = ARM_PATHWAYS,
  arm_sets         = arm_sets,
  covariate_sets   = covariate_sets,
  mito_paths       = mito_paths,
  arm_universe_path = arm_universe_path,
  coverage         = coverage,
  arm_summary      = arm_summary,
  harmonisation    = harm_report,
  null_manifest    = null_manifest,
  rules = list(
    primary_exposure = "OXPHOS subunits, nuclear-encoded, on BOTH instruments",
    primary_negative = "OXPHOS assembly factors",
    instruments = c(
      gsva    = paste("CO-PRIMARY, log VST, kcdf Gaussian, cohort-relative,",
                      "travels to Block F"),
      mitopps = "CO-PRIMARY, linear DESeq2, composition, TCGA-internal only",
      content = "descriptive, log2 sum of linear counts",
      zmean   = "descriptive, mean z of VST"),
    disagreement = "report both; claim only what BOTH co-primary instruments support",
    mtdna = paste(EXPECT_MITOCARTA_MTDNA, "mtDNA-encoded genes held in", MTDNA_PATHWAY),
    null  = paste0(G7_NULL_NSETS, " sets per arm, matched within ", G7_NULL_NBIN,
                   " ventiles of mean linear expression, drawn without replacement"),
    apoptosis_arms = paste("mouse Apoptosis-PRO/ANTI NOT added to the mitoPPS",
                           "universe - see section 5.1, NEEDS SIGN-OFF")
  ),
  built = Sys.time()
)

saveRDS(mito, file.path(DIR_RESULTS, "tcga_brca_mito_scores.rds"))
readr::write_csv(scores,      file.path(DIR_TABLES, "tcga_brca_mito_scores.csv"))
readr::write_csv(arm_summary, file.path(DIR_TABLES, "mito_arm_summary.csv"))
readr::write_csv(coverage,    file.path(DIR_TABLES, "mito_set_coverage.csv"))
readr::write_csv(arms_spec,   file.path(DIR_TABLES, "mito_specificity_arms.csv"))
readr::write_csv(harm_report, file.path(DIR_TABLES, "mito_symbol_harmonisation.csv"))

message("\n07: done.")
message("    results/tcga_brca_mito_scores.rds")
message("    results/mito_null/  (", nrow(null_manifest), " arms)")
message("    outputs/tables/  5 tables")
message("\n    REMINDER for script 09: both instruments are co-primary. An ",
        "effect on\n    one alone is instrument-dependent, not a positive result.")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  m <- readRDS(file.path(DIR_RESULTS, "tcga_brca_mito_scores.rds"))

  # --- the mitoPPS rewrite, checked the slow honest way ---------------------
  # A 5-pathway toy where the triple loop is transcribed literally from Monzel
  # et al. steps 1-4. If this disagrees, section 5.3's assertion is not enough.
  set.seed(1)
  St <- matrix(runif(30, 1, 100), nrow = 5,
               dimnames = list(paste0("P", 1:5), paste0("S", 1:6)))
  brute <- St * NA
  for (p in 1:5) for (i in 1:6) {
    q   <- setdiff(1:5, p)
    num <- St[p, i] / St[q, i]
    den <- vapply(q, function(j) mean(St[p, ] / St[j, ]), numeric(1))
    brute[p, i] <- mean(num / den)
  }
  max(abs(brute - .mitopps_universe(St)))     # expect ~1e-16

  # --- do the two co-primary instruments agree, arm by arm? ----------------
  # Look at this BEFORE any model. A near-zero correlation is not a bug; it is
  # the reason the "claim only what both support" rule exists.
  m$arm_summary %>%
    dplyr::select(arm, n_genes, rho_instruments) %>%
    dplyr::arrange(rho_instruments) %>%
    as.data.frame() %>% print(row.names = FALSE)

  plot(m$gsva_arms["OXPHOS subunits", ], m$mitopps_arms["OXPHOS subunits", ],
       pch = 16, cex = 0.3, xlab = "GSVA level", ylab = "mitoPPS composition",
       main = "OXPHOS subunits on both co-primary instruments")

  # --- the primary contrast: subunits against their own assembly factors ---
  # Same complexes, same umbrella, different function. In the mouse the gap
  # between these two points IS the specificity claim.
  plot(m$gsva_arms["OXPHOS assembly factors", ], m$gsva_arms["OXPHOS subunits", ],
       pch = 16, cex = 0.3, xlab = "assembly factors", ylab = "subunits")
  .rho(m$gsva_arms["OXPHOS subunits", ], m$gsva_arms["OXPHOS assembly factors", ])

  # --- does OXPHOS move across ALL complexes, as the mouse claims? ---------
  cx <- grep("^Complex-level", rownames(m$gsva_arms), value = TRUE)
  round(cor(t(m$gsva_arms[c("OXPHOS subunits", cx), ]), method = "spearman"), 2)

  # --- mtDNA-encoded, held separate on purpose -----------------------------
  # If this tracked the nuclear subunits closely, the standing convention would
  # be costing nothing. It usually does not.
  .rho(m$gsva_arms["OXPHOS subunits", ], m$gsva_arms["mtDNA-encoded OXPHOS", ])

  # --- how confounded is the exposure by stroma? ---------------------------
  cv <- readRDS(file.path(DIR_RESULTS, "tcga_brca_covariates.rds"))$covariates
  i  <- match(colnames(m$gsva_arms), cv$patient)
  plot(cv$purity[i], m$gsva_arms["OXPHOS subunits", ], pch = 16, cex = 0.3,
       xlab = "ABSOLUTE purity", ylab = "OXPHOS subunits (GSVA)")
  boxplot(m$gsva_arms["OXPHOS subunits", ] ~ cv$PAM50[i], las = 2,
          ylab = "OXPHOS subunits (GSVA)")

  # --- the matched null, on the arm the claim rests on ---------------------
  # Script 09 turns these into percentiles of a model coefficient. Here they are
  # only scores, so all this shows is that the matching worked: the null sets
  # should straddle the observed arm rather than sitting to one side of it.
  n <- readRDS(m$null_manifest$file[m$null_manifest$arm == "OXPHOS subunits"])
  hist(rowMeans(n$gsva), breaks = 50,
       main = "2,000 expression-matched sets, mean GSVA", xlab = "")
  abline(v = mean(m$gsva_arms["OXPHOS subunits", ]), col = "red", lwd = 2)
  hist(rowMeans(n$mitopps), breaks = 50,
       main = "2,000 expression-matched sets, mean mitoPPS", xlab = "")
  abline(v = mean(m$mitopps_arms["OXPHOS subunits", ]), col = "red", lwd = 2)

  # Are the draws actually expression-matched? Member-gene mean expression
  # should overlap the arm's distribution, not the whole matrix's.
  gl <- rowMeans(readRDS(file.path(DIR_RESULTS, "tcga_brca_linear.rds"))$mat)
  boxplot(list(arm  = log10(gl[n$genes] + 1),
               null = log10(gl[unlist(n$draws[1:20])] + 1),
               all  = log10(gl + 1)), ylab = "log10 mean linear expression")

  # --- what the symbol rescue actually caught ------------------------------
  m$harmonisation %>% dplyr::count(status) %>% print(n = Inf)
  m$harmonisation %>% dplyr::filter(status == "resolved") %>% print(n = Inf)
  m$harmonisation %>%
    dplyr::filter(status != "matched", in_sets != "") %>%
    print(n = Inf)

  # --- the two proliferation covariates ------------------------------------
  # D7: they should be near-identical (9 genes of 327 differ). If they are not,
  # something is wrong with the disjoint construction, not with the biology.
  .rho(m$gsva_cov["PROLIF_STD", ], m$gsva_cov["PROLIF_DISJOINT", ])

}
