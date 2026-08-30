# 11_build_state_variable.R
# =============================================================================
# FREEZE the composite state variable for H4 (plan section 7.5), as an
# executable, portable constructor, BEFORE any H4 outcome data is touched.
#
# Built to:
#   docs/2026-08-27_human_validation_plan.md section 7.5 (the fixed definition)
#   docs/2026-08-28_D5_cohort_selection.md section 6.1 (STATE demoted to
#     pre-specified SECONDARY; its role changed, its definition did not)
#
# =============================================================================
# WHAT THIS SCRIPT IS FOR, AND WHAT IT IS NOT
# =============================================================================
# H4 is NOT tested in TCGA. It is tested in three neoadjuvant cohorts, scored
# independently and meta-analysed (D5 section 6.2). So the deliverable here is
# not a TCGA variable - it is a FROZEN FUNCTION that script 13 calls once per
# cohort, plus the written resolution of every choice section 7.5 left open.
#
# TCGA appears in this script for exactly three reasons, none of which is a
# test:
#   1. to prove the definition runs and produces non-degenerate cells;
#   2. to record reference cell counts;
#   3. to CALIBRATE the expression-tertile BUFFER fallback against the GISTIC
#      rule, in the one cohort that has both. The neoadjuvant cohorts have no
#      CNV, so the fallback is what H4 will actually use, and whether it agrees
#      with the rule it substitutes for is something that must be known BEFORE
#      the outcome data arrives.
#
# NO OUTCOME VARIABLE IS READ BY THIS SCRIPT. TCGA has no pCR endpoint; it is
# not a neoadjuvant cohort. There is no path by which anything here could have
# been tuned to a result.
#
# =============================================================================
# THE FREEZE, AND WHY IT IS A FREEZE
# =============================================================================
# Plan section 7.5: "STATE remains frozen in script 11 and must still not be
# revised after outcome data is seen. Its ROLE changes; its DEFINITION does
# not."
#
# The definition below is section 7.5's, transcribed without amendment. In
# particular the level-3-vs-level-4 secondary contrast is kept EXACTLY as
# written, including its pre-specified direction.
#
# A TEXT NOTE, NOT A MODIFICATION: that contrast has lost much of its
# motivation. It turns on BUFFER doing something, and BUFFER's role is what H1
# clause 2 failed to support (Block B, E1) and what Block G failed to support
# functionally at n = 1,130. The correct response is to say so in the
# manuscript, NOT to revise the variable. Revising it now would be permissible
# on timing - no outcome data has been seen - but the REASON to revise would be
# Blocks B, C and G, which is exactly the "pick the cell that survives" move the
# freeze exists to prevent. Frozen as specified; the caveat goes in prose.
#
# =============================================================================
# WHAT SECTION 7.5 LEAVES OPEN, AND HOW IT IS RESOLVED HERE
# =============================================================================
# Each of these is a choice the plan does not make. Each is made NOW, before any
# H4 outcome data exists, and each is saved with the frozen object so that
# script 13 cannot quietly make it differently.
#
#   (a) WHICH OXPHOS INSTRUMENT. Both. GSVA and mitoPPS are co-primary
#       throughout this arm (CLAUDE.md), so STATE is built twice and both are
#       always reported. This is NOT a menu to choose from after the fact.
#       The constructor takes a numeric vector and does not know or care which
#       instrument produced it, so it is portable to a cohort where only one
#       instrument is computable.
#
#   (b) WHICH MYC ESTIMATOR. M-a, the Felsher signature (plan 7.1, D2). It is
#       the primary estimator, and it is the only one of the three computable in
#       an array cohort - M-c is a GISTIC call and does not exist in GEO.
#
#   (c) WHICH OXPHOS ARM. `OXPHOS subunits`, not the OXPHOS umbrella. The
#       umbrella includes assembly factors; CLAUDE.md's standing convention and
#       the plan's primary OXPHOS measure are the subunits.
#
#   (d) MEDIAN OF WHAT. Of THIS COHORT, computed on the complete cases entering
#       that cohort's H4 model, at scoring time. NOT TCGA's median. GSVA is
#       cohort-relative (CLAUDE.md) and a threshold carried across cohorts is
#       meaningless. This is enforced by construction: the constructor computes
#       its own cuts from the vector it is handed and accepts no external cut.
#
#   (e) THE TIE RULE. `> median` is high. Matches the convention already used in
#       script 10 section 3, so the two dichotomisations in this arm agree.
#
#   (f) BUFFER WHERE THERE IS NO CNV. Section 7.5 says "expression tertile" and
#       stops. Resolved as: TOP tertile of MCL1 OR top tertile of BCL2L1, both
#       within cohort, `> quantile(2/3, type = 7)`. The OR mirrors the structure
#       of the GISTIC rule in script 03 (buffer_MCL1 OR buffer_BCL2L1). BCL2 is
#       deliberately excluded - plan section 9 mitigation 1 restricts the
#       BCL-family measures to BCL2L1 and MCL1 because BCL2 is estrogen-
#       responsive.
#
#       KNOWN ASYMMETRY, recorded rather than fixed: the GISTIC rule uses
#       DIFFERENT thresholds per gene (MCL1 == +2, BCL2L1 >= +1), and the
#       tertile rule uses the same quantile for both. The two are therefore not
#       calibrated to each other and BUFFER prevalence will differ between a
#       CNV cohort and an expression cohort. Section 5 measures how much.
#
#   (g) MISSING DATA. Complete cases, matching D8. A sample missing any one of
#       MYC, OXPHOS or BUFFER gets STATE = NA and drops out of the H4 model.
#
# =============================================================================
# SCALE DISCIPLINE: no scoring happens here. Every score is READ as built by
# scripts 06-08. Median splits are invariant to the frozen z-transform used
# elsewhere in this arm, so feeding z-scores or raw scores gives the identical
# factor; z-scores are used below only for consistency with Blocks B, C and D.
#
# SPECIES: human. TCGA-BRCA. See CLAUDE.md.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

message("\n11: freeze the STATE variable (plan section 7.5)\n",
        strrep("=", 78))

# -----------------------------------------------------------------------------
# 0. The frozen specification
# -----------------------------------------------------------------------------
# Level order is section 7.5's, verbatim. Level 1 is the reference and collapses
# across OXPHOS and BUFFER by design: section 7.5's own table says MYC-low has
# "no OXPHOS dependence expected", so it is one cell, not four.
STATE_LEVELS <- c(
  "MYC_low",                          # 1  reference
  "MYC_high_OXPHOS_low",              # 2
  "MYC_high_OXPHOS_high_unbuffered",  # 3  predicted chemo-SENSITIVE
  "MYC_high_OXPHOS_high_buffered"     # 4  predicted chemo-RESISTANT
)

# The secondary contrast IS directional. Section 7.5's table predicts level 3
# sensitive and level 4 resistant, so the direction was fixed by the plan and is
# not being invented here.
STATE_CONTRAST <- c("MYC_high_OXPHOS_high_unbuffered",
                    "MYC_high_OXPHOS_high_buffered")
STATE_CONTRAST_DIRECTION <- "pCR higher in level 3 than in level 4"

BUFFER_TERTILE_Q     <- 2 / 3
BUFFER_TERTILE_GENES <- c("MCL1", "BCL2L1")

MYC_ESTIMATOR <- "M_a"
ARM_PRIMARY   <- "OXPHOS subunits"
INSTRUMENTS   <- c("gsva", "mitopps")

# =============================================================================
# 1. Inputs - consumed as built, nothing rescored
# =============================================================================
message("\n1. inputs")

mito <- readRDS(file.path(DIR_RESULTS, "tcga_brca_mito_scores.rds"))
prim <- readRDS(file.path(DIR_RESULTS, "tcga_brca_priming.rds"))$priming
myc  <- readRDS(file.path(DIR_RESULTS, "tcga_brca_myc_scores.rds"))$estimators
cov  <- readRDS(file.path(DIR_RESULTS, "tcga_brca_covariates.rds"))$covariates
bc   <- readRDS(file.path(DIR_RESULTS, "block_c_models.rds"))

pat <- colnames(mito$gsva_arms)
stopifnot(identical(sort(pat), sort(cov$patient)),
          identical(sort(pat), sort(prim$patient)),
          identical(sort(pat), sort(myc$patient)))
message("   ", length(pat), " patients, scores read as built by 06-08")

# =============================================================================
# 2. The frozen constructor
# =============================================================================
# TWO PROPERTIES THAT MATTER, both asserted below.
#
# SELF-CONTAINED. The environment is set to baseenv() so the function closes
# over NOTHING. A closure saved with environment(f) == globalenv() is rebound to
# the LOADING session's globals on readRDS, which means a cohort scored next
# month could silently pick up a different STATE_LEVELS. Setting baseenv() makes
# that impossible: every constant the function uses is one of its own arguments,
# with a literal default.
#
# NO EXTERNAL CUTS. The function computes its own median from the vector it is
# handed and offers no argument to supply one. A TCGA threshold cannot be
# applied to another cohort even by mistake - see resolution (d).
message("\n2. the frozen constructor")

build_state <- function(myc, oxphos, buffer,
                        state_levels = c(
                          "MYC_low",
                          "MYC_high_OXPHOS_low",
                          "MYC_high_OXPHOS_high_unbuffered",
                          "MYC_high_OXPHOS_high_buffered")) {
  if (!is.numeric(myc) || !is.numeric(oxphos)) {
    stop("myc and oxphos must be numeric cohort-relative scores", call. = FALSE)
  }
  if (!is.logical(buffer)) {
    stop("buffer must be logical (TRUE = buffered)", call. = FALSE)
  }
  n <- length(myc)
  if (length(oxphos) != n || length(buffer) != n) {
    stop("myc, oxphos and buffer must be the same length", call. = FALSE)
  }
  if (length(state_levels) != 4L) {
    stop("state_levels must have exactly 4 levels", call. = FALSE)
  }

  # Complete cases only (resolution g). The cuts are computed on these rows and
  # on no others, so a sample that drops out cannot move another sample's cell.
  ok <- !is.na(myc) & !is.na(oxphos) & !is.na(buffer)
  out <- rep(NA_character_, n)
  if (!any(ok)) return(factor(out, levels = state_levels))

  # Resolution (d): the cuts come from THIS cohort. There is deliberately no
  # argument by which an external threshold could be supplied.
  # Resolution (e): strictly greater than the median is "high".
  myc_high <- ok & myc    > stats::median(myc[ok])
  ox_high  <- ok & oxphos > stats::median(oxphos[ok])

  out[ok & !myc_high]                          <- state_levels[1]
  out[ok &  myc_high & !ox_high]               <- state_levels[2]
  out[ok &  myc_high &  ox_high & !buffer]     <- state_levels[3]
  out[ok &  myc_high &  ox_high &  buffer]     <- state_levels[4]

  factor(out, levels = state_levels)
}

# Resolution (f). Returns logical with NA where either gene is missing, so it
# drops into build_state()'s `buffer` argument unchanged.
buffer_from_expression <- function(mcl1, bcl2l1, q = 2 / 3) {
  if (!is.numeric(mcl1) || !is.numeric(bcl2l1)) {
    stop("mcl1 and bcl2l1 must be numeric log-scale expression", call. = FALSE)
  }
  if (length(mcl1) != length(bcl2l1)) {
    stop("mcl1 and bcl2l1 must be the same length", call. = FALSE)
  }
  if (!is.numeric(q) || length(q) != 1L || q <= 0 || q >= 1) {
    stop("q must be a single value strictly between 0 and 1", call. = FALSE)
  }
  ok  <- !is.na(mcl1) & !is.na(bcl2l1)
  out <- rep(NA, length(mcl1))
  if (!any(ok)) return(out)
  # type = 7 named explicitly: it is R's default and a silent change of default
  # would move the cut.
  c1 <- stats::quantile(mcl1[ok],   probs = q, names = FALSE, type = 7)
  c2 <- stats::quantile(bcl2l1[ok], probs = q, names = FALSE, type = 7)
  out[ok] <- (mcl1[ok] > c1) | (bcl2l1[ok] > c2)
  out
}

environment(build_state)            <- baseenv()
environment(buffer_from_expression) <- baseenv()

stopifnot(identical(environment(build_state), baseenv()),
          identical(environment(buffer_from_expression), baseenv()))

# Static check that nothing leaked. `::` is expected and is a base function;
# a non-empty $variables would mean the function reads a global and the freeze
# is not a freeze.
if (requireNamespace("codetools", quietly = TRUE)) {
  for (nm in c("build_state", "buffer_from_expression")) {
    g <- codetools::findGlobals(get(nm), merge = FALSE)
    if (length(g$variables)) {
      stop(nm, " closes over global variable(s): ",
           paste(g$variables, collapse = ", "),
           ". The freeze would not hold across sessions.", call. = FALSE)
    }
  }
  message("   both functions are self-contained (no free variables)")
} else {
  message("   codetools absent - environment asserted, globals not checked")
}

# The integrity contract. Script 13 must call THESE objects and assert the
# deparsed text still matches, rather than re-implementing section 7.5.
#
# WHAT A MISMATCH MEANS. deparse() is a formatter, not a cryptographic seal: a
# change of R version could in principle re-wrap the text and trip the check
# against an untouched function. So a mismatch means INSPECT - diff the two
# character vectors and see whether the logic moved or only the line breaks did.
# It does not by itself mean the definition was edited.
DEFINITION_SOURCE <- list(
  build_state            = deparse(build_state),
  buffer_from_expression = deparse(buffer_from_expression)
)

# =============================================================================
# 3. TCGA demonstration - reference cells, not a test
# =============================================================================
message("\n3. TCGA demonstration")

D <- tibble::tibble(patient = pat) %>%
  dplyr::left_join(dplyr::select(cov, patient, PAM50, er_call, BUFFER,
                                 complete_block_c),
                   by = "patient") %>%
  dplyr::left_join(dplyr::select(myc, patient, dplyr::all_of(MYC_ESTIMATOR)),
                   by = "patient") %>%
  dplyr::left_join(dplyr::select(prim, patient, log2_MCL1, log2_BCL2L1),
                   by = "patient")

# The frozen z-scale from Block C, for consistency with Blocks B/C/D only. A
# median split is invariant to it; see the scale-discipline note in the header.
Z <- bc$zscale
.apply_z <- function(v, key) {
  k <- Z[[key]]; stopifnot(!is.null(k)); (v - k[["mean"]]) / k[["sd"]]
}
D$MYC <- .apply_z(D[[MYC_ESTIMATOR]], MYC_ESTIMATOR)
for (ins in INSTRUMENTS) {
  m <- mito[[paste0(ins, "_arms")]]
  D[[paste0("OX_", ins)]] <- .apply_z(
    as.numeric(m[ARM_PRIMARY, D$patient]),
    paste(ins, ARM_PRIMARY, sep = "|"))
}

# BUFFER, both ways. The GISTIC form is script 03's, read not rebuilt.
D$BUFFER_gistic <- D$BUFFER
D$BUFFER_expr   <- buffer_from_expression(D$log2_MCL1, D$log2_BCL2L1,
                                          q = BUFFER_TERTILE_Q)

# Demonstrated on Block C's 938, so the cell counts are comparable with every
# other block in this arm. The full cohort is reported alongside.
idx_c <- which(D$complete_block_c)
stopifnot(length(idx_c) == 938L)

for (ins in INSTRUMENTS) {
  D[[paste0("STATE_", ins)]] <- build_state(
    myc    = D$MYC,
    oxphos = D[[paste0("OX_", ins)]],
    buffer = D$BUFFER_gistic,
    state_levels = STATE_LEVELS)
}

# STATE is NA wherever a constituent is missing, and an NA level would pivot
# into a column literally named "NA". Named explicitly instead, so the dropped
# samples are visible in the table rather than inferred from a row total.
.label_na <- function(f) {
  factor(ifelse(is.na(f), "not_computable", as.character(f)),
         levels = c(STATE_LEVELS, "not_computable"))
}

cells <- dplyr::bind_rows(lapply(INSTRUMENTS, function(ins) {
  dplyr::bind_rows(
    tibble::tibble(set = "full cohort",
                   instrument = ins,
                   STATE = .label_na(D[[paste0("STATE_", ins)]])),
    tibble::tibble(set = "Block C 938",
                   instrument = ins,
                   STATE = .label_na(D[[paste0("STATE_", ins)]][idx_c])))
})) %>%
  dplyr::count(set, instrument, STATE, .drop = FALSE) %>%
  tidyr::pivot_wider(names_from = STATE, values_from = n, values_fill = 0L)

message("   cell counts (BUFFER from GISTIC):")
cells %>% as.data.frame() %>% print(row.names = FALSE)

# Do the two instruments agree on cell assignment? This matters: "claim only
# what both instruments support" is easy for a coefficient and harder for a
# categorical, because the two STATEs are not the same partition of patients.
# Reported so script 13 knows what it is meta-analysing.
agree <- table(gsva = D$STATE_gsva[idx_c], mitopps = D$STATE_mitopps[idx_c],
               useNA = "no")
concord_state <- sum(diag(agree)) / sum(agree)
message("   instrument agreement on cell assignment (938): ",
        sprintf("%.1f%%", 100 * concord_state))

# =============================================================================
# 4. Is STATE just PAM50?
# =============================================================================
# Exposure-side sanity, no outcome. If levels 3 and 4 were purely basal-like the
# secondary contrast would be a subtype contrast wearing a mechanism's name.
message("\n4. composition")

comp <- dplyr::bind_rows(lapply(INSTRUMENTS, function(ins) {
  tibble::tibble(instrument = ins,
                 STATE   = .label_na(D[[paste0("STATE_", ins)]][idx_c]),
                 PAM50   = tidyr::replace_na(as.character(D$PAM50[idx_c]),
                                             "unknown"),
                 er_call = tidyr::replace_na(as.character(D$er_call[idx_c]),
                                             "unknown"))
}))

# Both instruments, because a claim is only made where both agree - and the two
# STATEs are different partitions of the same patients, so one table cannot
# stand in for the other.
for (ins in INSTRUMENTS) {
  message("   STATE x PAM50, ", ins, ", Block C 938:")
  comp %>% dplyr::filter(instrument == ins) %>%
    dplyr::count(STATE, PAM50, .drop = FALSE) %>%
    tidyr::pivot_wider(names_from = PAM50, values_from = n, values_fill = 0L) %>%
    as.data.frame() %>% print(row.names = FALSE)
}

# ER, because plan section 9 makes it the confound to watch in this arm.
message("   STATE x ER, GSVA, Block C 938:")
comp %>% dplyr::filter(instrument == "gsva") %>%
  dplyr::count(STATE, er_call, .drop = FALSE) %>%
  tidyr::pivot_wider(names_from = er_call, values_from = n, values_fill = 0L) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 5. Calibrating the expression fallback against the GISTIC rule
# =============================================================================
# THE POINT OF THIS SECTION. The neoadjuvant cohorts have no CNV, so H4 will use
# BUFFER_expr. TCGA is the only cohort in this arm carrying both, so it is the
# only place the substitution can be checked at all - and it has to be checked
# before the outcome data arrives, not after.
#
# If they disagree badly that is a finding about the plan's fallback, and
# revising the fallback WOULD still be permissible right now, because no outcome
# data has been seen. It would not be permissible later. That asymmetry is the
# whole reason this runs here.
#
# READ KAPPA, NOT THE PREVALENCES. The two rules are near-guaranteed to agree on
# prevalence and that agreement is arithmetic, not evidence:
#   - GISTIC BUFFER is 512 of the 938 = 54.6% (script 10 section 2).
#   - "top tertile of either of two genes" is 1 - (2/3)^2 = 55.6% whenever the
#     two genes are not strongly correlated. Checked on random input: 55.4%.
# So matching prevalence tells us nothing about whether the same PATIENTS are
# being called buffered. Concordance and kappa are the numbers that carry
# information; the prevalences are printed only so the coincidence is visible
# and cannot be mistaken for a validation.
message("\n5. BUFFER: GISTIC vs the expression-tertile fallback")

bf <- D[idx_c, ]
ok <- !is.na(bf$BUFFER_gistic) & !is.na(bf$BUFFER_expr)
tab <- table(gistic = bf$BUFFER_gistic[ok], expr = bf$BUFFER_expr[ok])

po <- sum(diag(tab)) / sum(tab)
pe <- sum(rowSums(tab) * colSums(tab)) / sum(tab)^2
kappa <- (po - pe) / (1 - pe)

message("   n with both: ", sum(ok))
message("   prevalence  GISTIC ", sprintf("%.1f%%", 100 * mean(bf$BUFFER_gistic[ok])),
        "   |  expression tertile ", sprintf("%.1f%%", 100 * mean(bf$BUFFER_expr[ok])))
message("   concordance ", sprintf("%.1f%%", 100 * po),
        "   |  Cohen's kappa ", sprintf("%.3f", kappa))
print(tab)

# STATE built on the fallback, so the two versions can be compared cell by cell.
for (ins in INSTRUMENTS) {
  D[[paste0("STATE_expr_", ins)]] <- build_state(
    myc    = D$MYC,
    oxphos = D[[paste0("OX_", ins)]],
    buffer = D$BUFFER_expr,
    state_levels = STATE_LEVELS)
}
state_agree <- table(gistic = D$STATE_gsva[idx_c],
                     expr   = D$STATE_expr_gsva[idx_c], useNA = "no")
concord_buffer <- sum(diag(state_agree)) / sum(state_agree)
message("   STATE agreement, GISTIC-BUFFER vs expression-BUFFER (GSVA): ",
        sprintf("%.1f%%", 100 * concord_buffer))

buffer_cal <- tibble::tibble(
  n              = sum(ok),
  prev_gistic    = mean(bf$BUFFER_gistic[ok]),
  prev_expr      = mean(bf$BUFFER_expr[ok]),
  concordance    = po,
  kappa          = kappa,
  state_concord  = concord_buffer)

# =============================================================================
# 6. Freeze
# =============================================================================
message("\n6. freeze")

frozen <- list(
  # --- the deliverable ---
  build_state            = build_state,
  buffer_from_expression = buffer_from_expression,
  definition_source      = DEFINITION_SOURCE,

  # --- the specification, so script 13 makes no choice of its own ---
  spec = list(
    source          = "plan section 7.5, transcribed without amendment",
    role            = paste("PRE-SPECIFIED SECONDARY (D5 section 6.1).",
                            "The primary H4 test is the CONTINUOUS three-way."),
    levels          = STATE_LEVELS,
    reference       = STATE_LEVELS[1],
    contrast        = STATE_CONTRAST,
    contrast_df     = 1L,
    contrast_dir    = STATE_CONTRAST_DIRECTION,
    myc_estimator   = MYC_ESTIMATOR,
    oxphos_arm      = ARM_PRIMARY,
    instruments     = INSTRUMENTS,
    tie_rule        = "strictly greater than the median is high",
    cut_scope       = paste("within cohort, on the complete cases of that",
                            "cohort's H4 model; NEVER a TCGA threshold"),
    buffer_cnv      = "script 03: cnv_MCL1 == 2 OR cnv_BCL2L1 >= 1",
    buffer_expr     = paste0("top tertile of MCL1 OR of BCL2L1, within cohort, ",
                             "> quantile(", signif(BUFFER_TERTILE_Q, 4),
                             ", type = 7); BCL2 excluded (plan section 9)"),
    buffer_genes    = BUFFER_TERTILE_GENES,
    missing         = "complete cases; NA in any constituent gives STATE = NA",
    caveat_for_text = paste(
      "The level-3-vs-4 contrast turns on BUFFER, and BUFFER's role is what",
      "Block B (E1) and Block G both failed to support. Frozen as specified;",
      "the lost motivation belongs in the manuscript text, not in a revision.")),

  # --- the TCGA demonstration ---
  tcga = list(
    # log2_MCL1 / log2_BCL2L1 are carried deliberately: without them the
    # BUFFER fallback cannot be diagnosed from the saved object alone.
    state = dplyr::select(D, patient, MYC, dplyr::starts_with("OX_"),
                          log2_MCL1, log2_BCL2L1,
                          BUFFER_gistic, BUFFER_expr,
                          dplyr::starts_with("STATE_")),
    cells             = cells,
    composition       = comp,
    instrument_agree  = agree,
    instrument_concord = concord_state,
    buffer_calibration = buffer_cal,
    buffer_table      = tab,
    n_full            = nrow(D),
    n_block_c         = length(idx_c)),

  no_outcome_data = paste("TCGA is not a neoadjuvant cohort and carries no pCR",
                          "endpoint. No outcome variable was read."),
  frozen_at = Sys.time())

saveRDS(frozen, file.path(DIR_RESULTS, "state_definition.rds"))
readr::write_csv(cells,      file.path(DIR_TABLES, "state_cell_counts.csv"))
readr::write_csv(buffer_cal, file.path(DIR_TABLES, "state_buffer_calibration.csv"))

message("\n11: done. STATE is frozen.")
message("    results/state_definition.rds")
message("    outputs/tables/  2 tables")
message("    Script 13 must CALL frozen$build_state and assert")
message("    identical(deparse(frozen$build_state),")
message("              frozen$definition_source$build_state)")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  s <- readRDS(file.path(DIR_RESULTS, "state_definition.rds"))

  # --- the specification, first and on its own -----------------------------
  str(s$spec, max.level = 1)
  s$spec$levels
  s$spec$contrast
  s$spec$caveat_for_text

  # --- the integrity contract, exercised the way script 13 will ------------
  identical(deparse(s$build_state), s$definition_source$build_state)
  identical(environment(s$build_state), baseenv())

  # --- cells -----------------------------------------------------------------
  s$tcga$cells %>% as.data.frame() %>% print(row.names = FALSE)
  s$tcga$instrument_agree
  s$tcga$instrument_concord

  # --- the calibration that actually decides whether H4's BUFFER is usable ---
  s$tcga$buffer_table
  s$tcga$buffer_calibration %>% as.data.frame() %>% print(row.names = FALSE)

  # Where do GISTIC and the tertile disagree? If the disagreement is one-sided
  # the fallback is biased, not merely noisy, and that changes how a null H4
  # would be read.
  st <- s$tcga$state
  d1 <- st$BUFFER_gistic & !st$BUFFER_expr    # CNV says buffered, expression not
  d2 <- !st$BUFFER_gistic & st$BUFFER_expr    # the other way
  c(gistic_only = sum(d1, na.rm = TRUE), expr_only = sum(d2, na.rm = TRUE))

  boxplot(st$log2_MCL1 ~ st$BUFFER_gistic, xlab = "GISTIC BUFFER",
          ylab = "log2 MCL1")

  # --- is the definition degenerate anywhere? ------------------------------
  # Levels 3 and 4 are the whole secondary contrast. If either is small the
  # contrast is wide, not wrong - report the CI, do not re-cut the variable.
  table(st$STATE_gsva, useNA = "ifany")
  table(st$STATE_mitopps, useNA = "ifany")

  # --- portability check, on synthetic data --------------------------------
  # The constructor must work on a cohort it has never seen, with no reference
  # to TCGA. If this errors, script 13 will error too.
  set.seed(PROJECT_SEED)
  n <- 200
  s$build_state(myc = rnorm(n), oxphos = rnorm(n),
                buffer = sample(c(TRUE, FALSE), n, TRUE)) %>% table()

  # and that missingness propagates rather than silently reclassifying
  m <- rnorm(n); m[1:5] <- NA
  table(s$build_state(m, rnorm(n), sample(c(TRUE, FALSE), n, TRUE)),
        useNA = "ifany")

  # --- the fallback on synthetic data too ----------------------------------
  s$buffer_from_expression(rnorm(n), rnorm(n)) %>% table()

}
