# 08_score_priming.R
# =============================================================================
# The apoptotic priming endpoints, their negative controls, and FOXO3 activity.
#
#   PRIME           log2(BBC3) - log2(BCL2L1)      the pre-specified endpoint
#   PRIME_INDEX     summed pro / summed anti, z    robust secondary
#   negatives       BID, BAX, BCL2L11, BAK1 over BCL2L1
#   FOXO3_activity  CollecTRI regulon, decoupleR ULM
#
# Specification: plan sections 7.3 and 7.4, with three points agreed 2026-08-28
# and recorded in section 0 below.
#
# =============================================================================
# SCALE DISCIPLINE - this script uses a THIRD scale, deliberately
# =============================================================================
# Script 07's rule was GSVA on VST, mitoPPS on linear. This one is different and
# the difference is the whole point:
#
#   PRIME and every ratio endpoint are computed as log2 of the LINEAR
#   DESeq2-normalised matrix. That is the plan's definition verbatim -
#   "log2(BBC3) - log2(BCL2L1)" - and it is NOT the same as the difference of
#   two VST values. VST shrinks variance gene-specifically, so it would compress
#   BBC3 (median 632) more than BCL2L1 (median 5,574), i.e. differentially
#   compress the numerator and the denominator of the endpoint. The VST
#   difference is computed as a labelled sensitivity, never as the primary.
#
#   FOXO3 activity is a decoupleR ULM score and therefore reads the LOG VST
#   matrix, exactly as M-b did in script 06.
#
# So both matrices are open in this script. They are used for different
# quantities and are never mixed within one. Every block states its scale.
#
# NO PSEUDOCOUNT. Every gene entering a ratio is asserted strictly positive in
# every sample before the log is taken. A pseudocount that is never needed is a
# silent thumb on the scale; if the assertion ever fires, that is a decision to
# make, not a default to have already taken.
#
# =============================================================================
# WHAT THIS SCRIPT DELIBERATELY DOES NOT COMPUTE
# =============================================================================
# It does not relate PRIME to MYC or to OXPHOS, marginally or otherwise. That
# relationship IS the hypothesis, and script 09 is where it is tested under the
# pre-registered model with its covariates, its D7 specifications and its
# expression-matched null. Printing a marginal correlation here would put the
# answer on screen before the test that is supposed to produce it.
#
# SPECIES: human. See CLAUDE.md.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
suppressPackageStartupMessages({
  library(decoupleR)
})

message("\n08: priming endpoints and FOXO3 activity\n", strrep("=", 78))

# =============================================================================
# 0. The three points agreed 2026-08-28, fixed here
# =============================================================================
# (1) PRIME is log2 of LINEAR normalised counts, not a VST difference. VST is a
#     sensitivity. Reason above.
#
# (2) The robust secondary index gets an explicit membership, named now rather
#     than assembled later. BCL2 is EXCLUDED from the anti-apoptotic sum, per
#     plan section 9 mitigation 1: BCL2 is the canonical estrogen-responsive
#     gene and ER+ tumours are BCL2-high, so including it makes the index partly
#     an ER readout. The with-BCL2 version is computed and reported alongside.
#
# (3) The FOXO3 regulon is MitoCarta-stripped as G1 stripped the MYC estimators,
#     AND explicitly stripped of every gene used in any priming endpoint.
#     `BBC3` and `BCL2L11` are both CollecTRI FOXO3 targets, so without this the
#     FOXO3-to-priming relationship would be partly definitional. As of the
#     2026-08-28 snapshot both are already removed by the MitoCarta strip, so
#     the priming strip is currently a no-op - it is applied anyway so the
#     guarantee survives a change in either input.

PRO_BH3_ONLY <- c("BBC3", "PMAIP1", "BID", "BIK", "BAD", "BMF", "HRK", "BCL2L11")
PRO_EFFECTOR <- c("BAX", "BAK1", "BOK")
ANTI_GUARD   <- c("BCL2L1", "MCL1", "BCL2L2", "BCL2A1", "BCL2L10")
ANTI_ER      <- "BCL2"          # reported version only, see (2)

PRIME_NUM <- "BBC3"
PRIME_DEN <- "BCL2L1"
NEG_NUM   <- c("BID", "BAX", "BCL2L11", "BAK1")   # all over PRIME_DEN

PRIMING_GENES <- unique(c(PRO_BH3_ONLY, PRO_EFFECTOR, ANTI_GUARD, ANTI_ER,
                          PRIME_NUM, PRIME_DEN, NEG_NUM))

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")

lin <- readRDS(file.path(DIR_RESULTS, "tcga_brca_linear.rds"))
if (!identical(lin$scale, "linear_deseq2_normalised")) {
  stop("expected the LINEAR object, got scale = '", lin$scale,
       "'. PRIME must not be built from a logged matrix.", call. = FALSE)
}
L <- lin$mat

vst <- readRDS(file.path(DIR_RESULTS, "tcga_brca_vst.rds"))
if (!identical(vst$scale, "log_vst")) {
  stop("expected the LOG-scale VST object, got scale = '", vst$scale, "'.",
       call. = FALSE)
}
E <- vst$mat

stopifnot(identical(dimnames(E), dimnames(L)))
if (!(min(L) >= 0 && max(L) > 50 * max(E))) {
  stop("the VST and LINEAR matrices are not on different scales. Check 01.",
       call. = FALSE)
}
message("   LINEAR : ", nrow(L), " x ", ncol(L), "  -> PRIME and the ratios")
message("   VST    : ", nrow(E), " x ", ncol(E), "  -> FOXO3 ULM, PRIME sensitivity")

missing_genes <- setdiff(PRIMING_GENES, rownames(L))
if (length(missing_genes)) {
  stop("priming gene(s) absent from the expression matrix: ",
       paste(missing_genes, collapse = ", "),
       "\nThe endpoint cannot be built as specified. Do not substitute.",
       call. = FALSE)
}
message("   all ", length(PRIMING_GENES), " priming genes present")

# =============================================================================
# 2. PRIME and the negative-control endpoints
# =============================================================================
# SCALE: log2 of LINEAR normalised counts.
#
# The negatives exist because the mouse says only PUMA/BCL-XL reverses. Each is
# the identical construction with a different pro-apoptotic numerator over the
# same BCL2L1 denominator, so a difference between PRIME and a negative cannot
# come from the denominator, the scale, or the normalisation.
message("\n2. PRIME and the negative-control endpoints (log2 linear)")

.assert_positive <- function(genes, what) {
  sub <- L[genes, , drop = FALSE]
  bad <- rownames(sub)[apply(sub, 1L, function(v) any(v <= 0))]
  if (length(bad)) {
    stop(what, ": gene(s) with a zero or negative normalised count in at least ",
         "one sample -> ", paste(bad, collapse = ", "),
         ".\nlog2 is undefined there. This needs a decision (pseudocount, or ",
         "dropping the gene), not a default.", call. = FALSE)
  }
  invisible(TRUE)
}

.assert_positive(c(PRIME_NUM, PRIME_DEN, NEG_NUM), "ratio endpoints")

.ratio <- function(num, den = PRIME_DEN) {
  as.numeric(log2(L[num, ]) - log2(L[den, ]))
}

PRIME <- .ratio(PRIME_NUM)
negatives <- vapply(NEG_NUM, .ratio, numeric(ncol(L)))
colnames(negatives) <- paste0(NEG_NUM, "_over_", PRIME_DEN)

message(sprintf(
  "   PRIME = log2(%s) - log2(%s): mean %.3f, sd %.3f, range %.2f to %.2f",
  PRIME_NUM, PRIME_DEN, mean(PRIME), stats::sd(PRIME), min(PRIME), max(PRIME)))
for (j in colnames(negatives)) {
  message(sprintf("   %-22s mean %+.3f, sd %.3f", j,
                  mean(negatives[, j]), stats::sd(negatives[, j])))
}

# --- the pre-specified sensitivity: the same endpoint on the VST -------------
# Reported, never primary. If these two disagree materially the choice of scale
# is load-bearing and belongs in Methods rather than in a comment.
PRIME_vst <- as.numeric(E[PRIME_NUM, ] - E[PRIME_DEN, ])
rho_scale <- stats::cor(PRIME, PRIME_vst, method = "spearman")
message(sprintf("   sensitivity, PRIME on VST difference: Spearman rho = %.4f",
                rho_scale))

# --- the individual limbs, saved for script 09's limb-wise fits --------------
# PRIME is close to a BBC3 readout - it correlates 0.838 with its numerator and
# -0.128 with its denominator, because log2(BBC3) varies about twice as much as
# log2(BCL2L1) - while the four negatives are genuine ratios. Script 09
# pre-specifies fits on each limb alone so that asymmetry becomes a reported
# quantity rather than a hidden one. They are saved here because 09 does no
# scoring of its own.
limbs <- t(log2(L[c(PRIME_NUM, PRIME_DEN, NEG_NUM), , drop = FALSE]))
colnames(limbs) <- paste0("log2_", c(PRIME_NUM, PRIME_DEN, NEG_NUM))
if (rho_scale < 0.95) {
  warning("PRIME on log2-linear and on the VST difference agree at rho = ",
          round(rho_scale, 3), ". The scale choice is load-bearing; report both ",
          "and say so in Methods.", call. = FALSE)
}

# =============================================================================
# 3. The robust secondary index
# =============================================================================
# SCALE: log2 of SUMMED LINEAR normalised counts, then z-scored across samples.
#
# Plan section 7.3: "summed pro-apoptotic BH3-only plus effectors over summed
# anti-apoptotic guardians, z-scored". Taken literally, which is why summing
# happens on the LINEAR scale - a sum of logs would be a product, a different
# quantity.
#
# ZEROS ARE HARMLESS HERE and that is a property of summing, not luck: HRK has
# 20 zero samples and BCL2L10 has 36, but a sum over 11 and 5 genes is strictly
# positive regardless. This is exactly why the individual-gene ratios in section
# 2 need the positivity assertion and this block does not.
#
# >>> A PROPERTY OF THE SPECIFICATION, STATED BEFORE THE RESULT.
# A sum of linear counts is ABUNDANCE-WEIGHTED by construction, so the index is
# dominated by its most abundant members rather than by equal contributions.
# MCL1 alone (median 25,503) is several times BCL2L1 (5,574), so the
# anti-apoptotic sum is largely an MCL1 readout. That is what the plan specifies
# and it is implemented as written; the per-gene contribution table below makes
# the weighting visible, and an equal-weight variant is computed alongside as a
# clearly labelled alternative. Neither replaces the other without a decision.
message("\n3. robust secondary index (log2 of summed linear, z-scored)")

PRO_ALL <- c(PRO_BH3_ONLY, PRO_EFFECTOR)

.z <- function(x) as.numeric(scale(x))
.sum_index <- function(pro, anti) {
  .z(log2(colSums(L[pro,  , drop = FALSE])) -
     log2(colSums(L[anti, , drop = FALSE])))
}

PRIME_INDEX          <- .sum_index(PRO_ALL, ANTI_GUARD)
PRIME_INDEX_withBCL2 <- .sum_index(PRO_ALL, c(ANTI_GUARD, ANTI_ER))

# Equal-weight alternative: mean of per-gene z-scores rather than a sum of
# counts. NOT the specified index. Computed so the weighting question is
# answerable from the saved object rather than by re-deriving it later.
#
# This variant DOES need a pseudocount, and the reason is instructive: it takes
# a log per gene, so HRK (zero in 20 samples) and BCL2L10 (36) would give
# log2(0) = -Inf and a whole column of NaN. The specified index escapes that by
# summing before taking the log. The pseudocount is confined to this alternative
# and never touches PRIME or the primary index.
.zrow <- function(g) t(scale(t(log2(L[g, , drop = FALSE] + 1))))
PRIME_INDEX_equalwt <- .z(
  colMeans(.zrow(PRO_ALL)) - colMeans(.zrow(ANTI_GUARD)))
stopifnot(all(is.finite(PRIME_INDEX_equalwt)))

contribution <- dplyr::bind_rows(
  tibble::tibble(side = "pro",  gene = PRO_ALL,
                 median_count = apply(L[PRO_ALL, , drop = FALSE], 1L,
                                      stats::median)),
  tibble::tibble(side = "anti", gene = ANTI_GUARD,
                 median_count = apply(L[ANTI_GUARD, , drop = FALSE], 1L,
                                      stats::median))
) %>%
  dplyr::group_by(side) %>%
  dplyr::mutate(share_of_sum = round(median_count / sum(median_count), 3)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(side, dplyr::desc(share_of_sum))

message("   per-gene share of its side's sum (this is the weighting):")
contribution %>% as.data.frame() %>% print(row.names = FALSE)

message(sprintf("\n   index vs PRIME              rho = %.3f",
                stats::cor(PRIME_INDEX, PRIME, method = "spearman")))
message(sprintf("   index vs with-BCL2 version  rho = %.3f",
                stats::cor(PRIME_INDEX, PRIME_INDEX_withBCL2, method = "spearman")))
message(sprintf("   index vs equal-weight alt   rho = %.3f",
                stats::cor(PRIME_INDEX, PRIME_INDEX_equalwt, method = "spearman")))

# =============================================================================
# 4. FOXO3 activity
# =============================================================================
# SCALE: LOG VST, kcdf-free (ULM is a linear model on the expression matrix),
# exactly as M-b in script 06.
#
# Plan section 7.4: FOXO3 is regulated by nuclear exclusion, so its own mRNA is
# a poor activity readout. The regulon score is the measure; FOXO3 mRNA is
# reported alongside as a CONTRAST, never as the measure.
message("\n4. FOXO3 activity (CollecTRI regulon, ULM on log VST)")

ct <- readr::read_tsv(PATH_COLLECTRI, show_col_types = FALSE, progress = FALSE)

# >>> THE TRAP THIS ASSERTION EXISTS FOR.
# readr parses `is_stimulation` / `is_inhibition` as LOGICAL. Script 04 reads the
# same file with utils::read.delim and gets CHARACTER, so it tests them with
# == "True". Both are correct for their own reader. Copying script 04's string
# test into a readr-based script makes every comparison FALSE, every edge take
# mor = -1, and the resulting activity score is the NEGATIVE of the truth - with
# no error and no warning. Assert the type, then assert the signs are not all
# one way.
if (!is.logical(ct$is_stimulation) || !is.logical(ct$is_inhibition)) {
  stop("CollecTRI sign columns are ", class(ct$is_stimulation),
       ", not logical. The reader has changed; fix the sign rule to match it ",
       "before going any further - a wrong rule inverts the score silently.",
       call. = FALSE)
}

# SIGN RULE, identical to script 06: mor = +1 stimulatory, -1 inhibitory, and an
# edge flagged BOTH takes +1.
foxo_edges <- ct %>%
  dplyr::filter(source_genesymbol == "FOXO3",
                !is.na(target_genesymbol), target_genesymbol != "")
n_both <- sum(foxo_edges$is_stimulation & foxo_edges$is_inhibition)
message(sprintf(
  "   FOXO3 edges: %d  (%d stimulation-only, %d inhibition-only, %d BOTH -> +1)",
  nrow(foxo_edges),
  sum(foxo_edges$is_stimulation & !foxo_edges$is_inhibition),
  sum(foxo_edges$is_inhibition & !foxo_edges$is_stimulation), n_both))

foxo_net_all <- foxo_edges %>%
  dplyr::transmute(source = "FOXO3",
                   target = target_genesymbol,
                   mor    = dplyr::if_else(is_stimulation, 1, -1),
                   likelihood = 1) %>%
  dplyr::distinct(source, target, .keep_all = TRUE)

if (all(foxo_net_all$mor > 0) || all(foxo_net_all$mor < 0)) {
  stop("every FOXO3 edge resolved to the same sign (", foxo_net_all$mor[1],
       "). That is the signature of a broken sign rule, not of biology.",
       call. = FALSE)
}
message("   after de-duplication: ", sum(foxo_net_all$mor > 0), " activating, ",
        sum(foxo_net_all$mor < 0), " repressing")

# --- the two strips ----------------------------------------------------------
g1 <- readRDS(file.path(DIR_RESULTS, "g1_overlap_audit.rds"))
mitocarta_all <- g1$reference_sets$MITOCARTA_ALL
stopifnot(length(mitocarta_all) > 1000L)

tg_raw   <- foxo_net_all$target
tg_mito  <- setdiff(tg_raw, mitocarta_all)
tg_final <- setdiff(tg_mito, PRIMING_GENES)

message("   regulon: raw ", length(tg_raw),
        " -> MitoCarta-stripped ", length(tg_mito),
        " -> priming-stripped ", length(tg_final))
removed_priming <- intersect(tg_mito, PRIMING_GENES)
message("   priming genes removed at the second strip: ",
        if (length(removed_priming)) paste(removed_priming, collapse = ", ")
        else "none (all were already MitoCarta genes)")
message("   priming genes in the RAW regulon, for the record: ",
        paste(intersect(tg_raw, PRIMING_GENES), collapse = ", "))

if (length(intersect(tg_final, PRIMING_GENES)) > 0L) {
  stop("a priming gene survived the strip. The FOXO3-to-priming relationship ",
       "would be partly definitional.", call. = FALSE)
}

.ulm_score <- function(targets, label) {
  net <- foxo_net_all %>%
    dplyr::filter(target %in% targets, target %in% rownames(E))
  message(sprintf("   %-24s %3d targets in the matrix", label, nrow(net)))
  res <- decoupleR::run_ulm(mat = E, network = net, .source = source,
                            .target = target, .mor = mor, minsize = 5L)
  v <- res %>%
    dplyr::filter(statistic == "ulm", source == "FOXO3") %>%
    dplyr::select(condition, score)
  stats::setNames(v$score, v$condition)[colnames(E)]
}

FOXO3_activity     <- .ulm_score(tg_final, "FOXO3 stripped")
FOXO3_activity_raw <- .ulm_score(tg_raw,   "FOXO3 raw")
FOXO3_mrna         <- as.numeric(E["FOXO3", ])

rho_strip <- stats::cor(FOXO3_activity, FOXO3_activity_raw, method = "spearman")
message(sprintf("   stripped vs raw regulon: rho = %.3f", rho_strip))
if (rho_strip < 0.90) {
  warning("stripping changed the FOXO3 score materially (rho = ",
          round(rho_strip, 3), "). Report both.", call. = FALSE)
}

# Plan section 7.4's premise, tested rather than assumed. A LOW correlation here
# is the expected result and the justification for using the regulon at all; a
# high one would mean the distinction is not doing any work in this cohort.
rho_mrna <- stats::cor(FOXO3_activity, FOXO3_mrna, method = "spearman")
message(sprintf("   FOXO3 activity vs FOXO3 mRNA: rho = %.3f", rho_mrna))

# =============================================================================
# 5. Diagnostics
# =============================================================================
# Deliberately NOT including anything against MYC or OXPHOS. See the header.
message("\n5. diagnostics")

cov <- readRDS(file.path(DIR_RESULTS, "tcga_brca_covariates.rds"))$covariates
stopifnot(identical(sort(cov$patient), sort(colnames(L))))
i_cov <- match(colnames(L), cov$patient)

.rho <- function(x, y) suppressWarnings(
  stats::cor(x, y, method = "spearman", use = "complete.obs"))

endpoints <- cbind(PRIME = PRIME, negatives,
                   PRIME_INDEX = PRIME_INDEX,
                   FOXO3_activity = FOXO3_activity)

diag_tbl <- tibble::tibble(
  endpoint   = colnames(endpoints),
  mean       = colMeans(endpoints),
  sd         = apply(endpoints, 2L, stats::sd),
  rho_PRIME  = apply(endpoints, 2L, function(v) .rho(v, PRIME)),
  rho_purity = apply(endpoints, 2L, function(v) .rho(v, cov$purity[i_cov])),
  rho_leuko  = apply(endpoints, 2L,
                     function(v) .rho(v, cov$leukocyte_fraction[i_cov])),
  rho_ER     = apply(endpoints, 2L, function(v)
    .rho(v, as.numeric(cov$er_call[i_cov] == "Positive")))
)
diag_tbl %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# rho_ER is here because of plan section 9: BCL2 is estrogen-responsive and ER+
# tumours are BCL2-high, which is why BCL2 is out of the primary index and why
# BCL2L1/MCL1 carry the primary measures. If PRIME itself tracks ER strongly,
# the ER adjustment in Block C is doing more work than the plan assumes and that
# should be known before the model is fitted rather than after.
message(sprintf("\n   PRIME vs ER status: rho = %.3f  (plan section 9)",
                diag_tbl$rho_ER[diag_tbl$endpoint == "PRIME"]))

# =============================================================================
# 6. Save
# =============================================================================
message("\n6. save")

priming <- tibble::tibble(
  patient                = colnames(L),
  PRIME                  = PRIME,
  PRIME_vst_sensitivity  = PRIME_vst,
  PRIME_INDEX            = PRIME_INDEX,
  PRIME_INDEX_withBCL2   = PRIME_INDEX_withBCL2,
  PRIME_INDEX_equalwt    = PRIME_INDEX_equalwt,
  FOXO3_activity         = FOXO3_activity,
  FOXO3_activity_raw     = FOXO3_activity_raw,
  FOXO3_mrna             = FOXO3_mrna
) %>%
  dplyr::bind_cols(tibble::as_tibble(negatives), tibble::as_tibble(limbs))

out <- list(
  priming      = priming,
  contribution = contribution,
  diagnostics  = diag_tbl,
  network      = foxo_net_all,
  regulon      = list(raw = tg_raw, mitocarta_stripped = tg_mito,
                      final = tg_final),
  sets = list(pro_bh3_only = PRO_BH3_ONLY, pro_effector = PRO_EFFECTOR,
              anti_guard = ANTI_GUARD, anti_er = ANTI_ER,
              prime_num = PRIME_NUM, prime_den = PRIME_DEN, neg_num = NEG_NUM),
  rules = list(
    primary_endpoint = paste("PRIME = log2(BBC3) - log2(BCL2L1), on log2 of",
                             "LINEAR normalised counts"),
    scale_sensitivity = paste("PRIME_vst_sensitivity = VST difference;",
                              "reported, never primary"),
    index = paste("summed pro over summed anti on the linear scale then",
                  "z-scored; abundance-weighted by construction, see",
                  "`contribution`"),
    bcl2 = paste("BCL2 excluded from the primary index (plan section 9);",
                 "with-BCL2 version reported alongside"),
    foxo3 = paste("CollecTRI regulon, ULM on log VST, MitoCarta-stripped then",
                  "priming-stripped; both-flagged edges take mor = +1"),
    not_computed = "nothing against MYC or OXPHOS - that is script 09's test"
  ),
  built = Sys.time()
)

saveRDS(out, file.path(DIR_RESULTS, "tcga_brca_priming.rds"))
readr::write_csv(priming,      file.path(DIR_TABLES, "tcga_brca_priming.csv"))
readr::write_csv(diag_tbl,     file.path(DIR_TABLES, "priming_diagnostics.csv"))
readr::write_csv(contribution, file.path(DIR_TABLES, "priming_index_weighting.csv"))

message("\n08: done.")
message("    results/tcga_brca_priming.rds")
message("    outputs/tables/  3 tables")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  pr <- readRDS(file.path(DIR_RESULTS, "tcga_brca_priming.rds"))
  d  <- pr$priming

  # --- does the scale choice matter? ---------------------------------------
  # If this is a straight line, the log2-linear vs VST decision is cosmetic and
  # can be reported in one sentence. If it bends at the extremes, the shrinkage
  # argument in the header is doing real work and belongs in Methods.
  plot(d$PRIME, d$PRIME_vst_sensitivity, pch = 16, cex = 0.3,
       xlab = "PRIME, log2 linear", ylab = "PRIME, VST difference")

  # --- PRIME against its own negative controls -----------------------------
  # They share a denominator, so some correlation is structural. What matters in
  # script 09 is whether the MYC:OXPHOS interaction distinguishes them, not
  # whether the endpoints themselves are independent.
  neg <- grep("_over_", names(d), value = TRUE)
  round(cor(d[, c("PRIME", neg)], method = "spearman"), 2)

  # --- how much of the endpoint is the denominator? ------------------------
  L <- readRDS(file.path(DIR_RESULTS, "tcga_brca_linear.rds"))$mat
  round(c(num = .rho(d$PRIME, log2(L["BBC3", ])),
          den = .rho(d$PRIME, log2(L["BCL2L1", ]))), 3)

  # --- the index weighting, seen rather than argued ------------------------
  pr$contribution %>% as.data.frame() %>% print(row.names = FALSE)
  plot(d$PRIME_INDEX, d$PRIME_INDEX_equalwt, pch = 16, cex = 0.3,
       xlab = "index as specified (abundance-weighted)",
       ylab = "equal-weight alternative")

  # --- FOXO3: is the regulon telling us anything its mRNA does not? --------
  # Plan section 7.4 asserts it should. Low correlation = the distinction earns
  # its place.
  plot(d$FOXO3_mrna, d$FOXO3_activity, pch = 16, cex = 0.3,
       xlab = "FOXO3 mRNA (VST)", ylab = "FOXO3 regulon activity (ULM)")

  # --- the ER confound, plan section 9 -------------------------------------
  cv <- readRDS(file.path(DIR_RESULTS, "tcga_brca_covariates.rds"))$covariates
  i  <- match(d$patient, cv$patient)
  boxplot(d$PRIME ~ cv$er_call[i], ylab = "PRIME", xlab = "ER status by IHC")
  boxplot(d$PRIME ~ cv$PAM50[i], las = 2, ylab = "PRIME")
  # and the same for the with-BCL2 index, which should track ER more strongly
  boxplot(d$PRIME_INDEX_withBCL2 ~ cv$er_call[i], ylab = "index, with BCL2")
  boxplot(d$PRIME_INDEX ~ cv$er_call[i], ylab = "index, BCL2 excluded")

  # --- the sign rule, as script 06 did for MYC -----------------------------
  # If flipping the both-flagged edges changes the ordering, the rule is
  # load-bearing and belongs in Methods rather than a comment.
  pr$network %>% dplyr::count(mor)

}
