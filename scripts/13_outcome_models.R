# 13_outcome_models.R
# =============================================================================
# BLOCK F1 - H4, conditional chemosensitivity, in the three neoadjuvant cohorts.
# The last untested falsification criterion.
#
# Built to:
#   docs/2026-08-30_STATE_frozen_and_H4_buffer_declaration.md   section 6
#   docs/2026-08-31_H4_amendment_instruments_and_coverage.md    A1, A2, A4, A5, C3
#   docs/2026-08-28_D5_cohort_selection.md                      6.1, 6.2, 6.3
#   docs/2026-08-27_human_validation_plan.md                    Block F1
#
# =============================================================================
# READ FIRST - THE F3 GATE, AND WHY THIS SCRIPT EXISTS ANYWAY
# =============================================================================
# CLAUDE.md carries a standing instruction: "Do not build script 13 before F3
# returns. F3 is a stop gate, not a robustness check." That instruction was
# written when F3 was one undivided test and forkscale was believed to be
# METABRIC-only. Two things have happened since, and the author should confirm
# this reading before sourcing this file.
#
#   1. G3 SPLIT F3. F3-pre (exposure side, TCGA) is discharged - INTERMEDIATE,
#      the stop gate did NOT close (2026-08-30 F3-pre note). F3 proper is a
#      SURVIVAL test and is still blocked on the METABRIC identifier file.
#
#   2. F3'S GATE CANNOT BE SATISFIED FOR F1, EVER. F3 asks whether STATE adds
#      anything over MB1_forkscale. Forkscale does not exist in these three
#      cohorts and cannot be constructed in them (F3-pre note section 7). So the
#      incremental-value test is not constructible here at any point.
#
# What F3 actually protects is the PROGNOSTIC claim - Menegollo Fig 7 is a
# survival result, and F2 is the block where the redundancy bites. F1's endpoint
# is pCR under neoadjuvant chemotherapy, in cohorts with no forkscale and no
# survival. Reading the gate as blocking F1 would not defer F1; it would end it.
#
# THIS IS A JUDGEMENT CALL AND IT IS THE AUTHOR'S TO MAKE. If the answer is that
# F1 waits for METABRIC, do not source this file.
#
# =============================================================================
# WHAT IS AND IS NOT IN THIS SCRIPT
# =============================================================================
#   F1  H4 in three neoadjuvant cohorts, meta-analysed.          IN
#   F2  STATE x chemotherapy on survival, METABRIC and SCAN-B.   NOT - no data
#   F3  incremental value over forkscale, survival.              NOT - see above
#   F4  STATE vs grade/size/nodes, descriptive.                  NOT - the
#       covariates exist only in GSE25066 (A4), so a three-cohort descriptive
#       is not available, and D6 (the mouse metastasis anchor) is still open.
#
# =============================================================================
# THE STANDARD THIS RESULT IS HELD TO, WRITTEN BEFORE IT EXISTS
# =============================================================================
# H4 IS SINGLE-INSTRUMENT (A1). GSVA only; mitoPPS does not exist in any of
# these cohorts. Blocks C, B and G each required agreement between two
# instruments that demonstrably disagree. H4 cannot clear that bar and must not
# be written as though it had.
#
# POWERED FOR A LARGE EFFECT ONLY. D5 section 5: at a target-cell risk ratio of
# 2.0 the meta-analysis sits at 86% and any single cohort below 60%. Report
# confidence intervals. A null is "not powered to exclude a modest effect", NOT
# falsification - falsification requires the INFORMATIVE direction, section 8.
#
# THE TWO-WAY IS NOT A FALLBACK. Plan Block F1: dropping BUFFER buys about five
# percentage points of power and sacrifices the entire distinction from Lee et
# al. 2017. It is not fitted here as an alternative and must not be retreated to.
#
# SCALE: every cohort is log2 (script 12). GSVA with kcdf = "Gaussian".
# COHORT-RELATIVITY: each cohort is scored in its OWN single GSVA run and the
# scores are never pooled. The ESTIMATES are meta-analysed (D5 6.2, CLAUDE.md).
#
# SPECIES: human. See CLAUDE.md.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
suppressPackageStartupMessages(library(GSVA))

message("\n13: Block F1 - H4 in the neoadjuvant cohorts\n", strrep("=", 78))

# -----------------------------------------------------------------------------
# 0. Constants - the declared specification, restated as code
# -----------------------------------------------------------------------------
ARM_PRIMARY   <- "OXPHOS subunits"
ARM_NEGATIVE  <- "OXPHOS assembly factors"
MYC_SET_NAME  <- "Felsher M-a"
CI_LEVEL      <- 0.95
GSVA_MIN_SET  <- 3L        # script 07's floor

# A2: a cohort joins the specificity battery only if the assembly-factor control
# reaches this coverage. GSE25066 is 34/68 and does not.
SPECIFICITY_MIN_COV <- 0.80

# A5: missing values.
NA_GENE_MAX_FRAC <- 0.05

# The declared directions (2026-08-30 section 6.5). Two-sided tests; these are
# the PREDICTIONS, evaluated mechanically in section 8.
PRED_THREEWAY <- "negative"
PRED_TWOWAY   <- "positive"

COHORTS <- c("GSE194040", "GSE164458", "GSE25066")
PRIMARY_COHORT <- "GSE194040"

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")

neo <- readRDS(file.path(DIR_RESULTS, "neoadjuvant_cohorts.rds"))
st  <- readRDS(file.path(DIR_RESULTS, "state_definition.rds"))
gs  <- readRDS(file.path(DIR_RESULTS, "tcga_brca_mito_scores.rds"))
g1  <- readRDS(file.path(DIR_RESULTS, "g1_overlap_audit.rds"))

stopifnot(identical(names(neo$cohorts), COHORTS))

# --- the STATE integrity contract, exercised as script 11 specified ----------
# A deparse mismatch means INSPECT, not necessarily tamper: deparse is a
# formatter and an R upgrade could rewrap an untouched function. Either way it
# must not be ignored.
if (!identical(deparse(st$build_state), st$definition_source$build_state)) {
  stop("the frozen STATE constructor does not match its recorded source. ",
       "Diff deparse(st$build_state) against st$definition_source$build_state ",
       "and establish whether the logic moved or only the line breaks did, ",
       "before running anything.", call. = FALSE)
}
stopifnot(identical(environment(st$build_state), baseenv()))
message("   STATE definition verified against its frozen source")

SETS_FULL <- list(gs$arm_sets[[ARM_PRIMARY]],
                  gs$arm_sets[[ARM_NEGATIVE]],
                  g1$estimators_stripped$FELSHER)
names(SETS_FULL) <- c(ARM_PRIMARY, ARM_NEGATIVE, MYC_SET_NAME)
stopifnot(length(SETS_FULL[[MYC_SET_NAME]]) == 61L)

# A2's sensitivity universe: the genes measured in all three cohorts.
INT3 <- Reduce(intersect, lapply(neo$cohorts, function(x) rownames(x$expr)))
message("   3-cohort gene intersection: ", length(INT3))

# =============================================================================
# 2. Score each cohort - one GSVA run per cohort, never pooled
# =============================================================================
# WHY EVERYTHING GOES IN ONE CALL PER COHORT. GSVA's enrichment statistic is a
# random walk over the genes of the matrix it is handed, and GSVA may restrict
# that matrix to genes appearing in the supplied sets. Two calls with different
# sets therefore walk DIFFERENT universes and their scores are not comparable -
# which would silently invalidate the A2 primary-versus-intersection comparison,
# the one thing that comparison exists to detect.
#
# So: the per-cohort sets, the intersection sets and two half-matrix PIN sets go
# into a single call per cohort. The pins hold the universe at the full matrix.
# Two halves rather than one all-genes set, because a set containing every gene
# leaves nothing outside it and the walk's miss-penalty term becomes 0/0.
# The device is script 07 section 4's; the reasoning is repeated here because
# this script is where it would be quietly dropped.
message("\n2. scoring")

.impute <- function(M, label) {                       # A5
  na_g <- rowSums(is.na(M))
  drop <- na_g > NA_GENE_MAX_FRAC * ncol(M)
  if (any(drop)) {
    message("   ", label, ": dropped ", sum(drop), " gene(s) above ",
            100 * NA_GENE_MAX_FRAC, "% NA")
    M <- M[!drop, , drop = FALSE]
  }
  n_imp <- sum(is.na(M))
  if (n_imp) {
    med <- apply(M, 1, stats::median, na.rm = TRUE)
    idx <- which(is.na(M), arr.ind = TRUE)
    M[idx] <- med[idx[, "row"]]
    message("   ", label, ": imputed ", n_imp, " value(s) at the gene median")
  }
  stopifnot(!anyNA(M))
  M
}

.score_cohort <- function(cn) {
  M <- .impute(neo$cohorts[[cn]]$expr, cn)
  # BUFFER_c is built from these two directly, so their survival is not
  # optional: the NA rule must not silently remove an exposure limb.
  gone <- setdiff(c("MCL1", "BCL2L1"), rownames(M))
  if (length(gone)) {
    stop(cn, ": the NA rule removed ", paste(gone, collapse = ", "),
         ", which BUFFER_c is built from. Stop and reconsider A5.",
         call. = FALSE)
  }
  present <- lapply(SETS_FULL, function(g) intersect(g, rownames(M)))
  int3    <- lapply(SETS_FULL, function(g) intersect(intersect(g, INT3),
                                                     rownames(M)))
  names(int3) <- paste0(names(int3), " [int3]")
  sets <- c(present, int3)
  n_ok <- vapply(sets, length, integer(1))
  if (any(n_ok < GSVA_MIN_SET)) {
    stop(cn, ": set(s) below the GSVA floor -> ",
         paste(names(sets)[n_ok < GSVA_MIN_SET], collapse = ", "), call. = FALSE)
  }
  sym <- rownames(M)
  sets[[".PIN_A"]] <- sym[c(TRUE, FALSE)]
  sets[[".PIN_B"]] <- sym[c(FALSE, TRUE)]

  par <- GSVA::gsvaParam(exprData = M, geneSets = sets, kcdf = "Gaussian",
                         minSize = GSVA_MIN_SET, maxSize = Inf)
  s <- GSVA::gsva(par, verbose = FALSE)
  wanted <- setdiff(names(sets), c(".PIN_A", ".PIN_B"))
  dropped <- setdiff(wanted, rownames(s))
  if (length(dropped)) {
    stop(cn, ": GSVA silently dropped ", paste(dropped, collapse = ", "),
         call. = FALSE)
  }
  list(scores = s[wanted, , drop = FALSE], expr = M,
       coverage = tibble::tibble(cohort = cn, set = names(sets)[
         !names(sets) %in% c(".PIN_A", ".PIN_B")],
         n_present = n_ok[!names(n_ok) %in% c(".PIN_A", ".PIN_B")]))
}

SC <- lapply(stats::setNames(COHORTS, COHORTS), function(cn) {
  message("   ", cn, " ...")
  .score_cohort(cn)
})

coverage <- dplyr::bind_rows(lapply(SC, `[[`, "coverage")) %>%
  dplyr::mutate(n_full = vapply(sub(" \\[int3\\]$", "", set),
                                function(s) length(SETS_FULL[[s]]), integer(1)),
                frac = n_present / n_full)
coverage %>% as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 3. Exposures - built per cohort, z-scored within cohort, never pooled
# =============================================================================
# BUFFER_c is the 2026-08-30 declaration section 6.1 definition:
#   mean( z(log2 MCL1), z(log2 BCL2L1) ), within cohort. BCL2 excluded.
message("\n3. exposures")

.z <- function(v) (v - mean(v, na.rm = TRUE)) / stats::sd(v, na.rm = TRUE)

.build <- function(cn, suffix = "") {
  s <- SC[[cn]]$scores
  M <- SC[[cn]]$expr
  p <- neo$cohorts[[cn]]$pheno
  stopifnot(identical(colnames(s), p$sample_id))

  key <- function(nm) paste0(nm, suffix)
  mcl1 <- as.numeric(M["MCL1", ])
  bcl  <- as.numeric(M["BCL2L1", ])

  d <- tibble::tibble(
    cohort    = cn,
    sample_id = p$sample_id,
    pcr       = p$pcr,
    subtype   = p$subtype,
    treatment = p$treatment,
    MYC       = .z(as.numeric(s[key(MYC_SET_NAME), ])),
    OX        = .z(as.numeric(s[key(ARM_PRIMARY), ])),
    OX_neg    = .z(as.numeric(s[key(ARM_NEGATIVE), ])),
    BUFFER_c  = (.z(mcl1) + .z(bcl)) / 2,
    z_MCL1    = .z(mcl1),
    z_BCL2L1  = .z(bcl))

  # STATE, through the FROZEN constructor. The buffer argument uses the frozen
  # expression fallback, because none of these cohorts has CNV.
  d$BUFFER_bin <- st$buffer_from_expression(mcl1, bcl)
  d$STATE <- st$build_state(myc = d$MYC, oxphos = d$OX, buffer = d$BUFFER_bin,
                            state_levels = st$spec$levels)
  d
}

D  <- lapply(stats::setNames(COHORTS, COHORTS), .build)
Di <- lapply(stats::setNames(COHORTS, COHORTS), .build, suffix = " [int3]")

for (cn in COHORTS) {
  message("   ", cn, ": n = ", nrow(D[[cn]]),
          "  pCR = ", sum(D[[cn]]$pcr, na.rm = TRUE),
          "  STATE cells = ",
          paste(as.integer(table(D[[cn]]$STATE)), collapse = "/"))
}

# =============================================================================
# 4. F1 primary - the declared three-way, per cohort
# =============================================================================
# A4: the harmonised covariate set is subtype + treatment. Stage and grade do
# not exist in the primary cohort. A term constant within a cohort is dropped
# and the drop is reported - not silently absorbed.
message("\n4. F1 primary")

.usable <- function(d, v) {
  if (!v %in% names(d)) return(FALSE)
  length(unique(d[[v]][!is.na(d[[v]])])) > 1L
}

.REPORTED <- new.env(parent = emptyenv())
.covars <- function(d, cn) {
  keep <- c("subtype", "treatment")[vapply(c("subtype", "treatment"),
                                           function(v) .usable(d, v), logical(1))]
  drop <- setdiff(c("subtype", "treatment"), keep)
  # Reported once per cohort: .covars is called for every model and the same
  # message eight times reads like eight different drops.
  if (length(drop) && !exists(cn, .REPORTED, inherits = FALSE)) {
    assign(cn, TRUE, .REPORTED)
    message("   ", cn, ": ", paste(drop, collapse = " and "),
            " constant - dropped from this cohort's model")
  }
  keep
}

.frame <- function(d, vars) {
  ok <- stats::complete.cases(d[, vars, drop = FALSE])
  droplevels(as.data.frame(d[ok, , drop = FALSE]))
}

.tidy <- function(m, term, label, cohort, extra = list()) {
  co <- summary(m)$coefficients
  if (!term %in% rownames(co)) {
    stop("term '", term, "' absent from '", label, "' in ", cohort, call. = FALSE)
  }
  crit <- stats::qnorm(1 - (1 - CI_LEVEL) / 2)
  out <- tibble::tibble(
    label = label, cohort = cohort, term = term,
    n = stats::nobs(m), events = sum(m$y == 1),
    estimate = co[term, 1], se = co[term, 2],
    ci_lo = co[term, 1] - crit * co[term, 2],
    ci_hi = co[term, 1] + crit * co[term, 2],
    or = exp(co[term, 1]), p = co[term, 4])
  if (length(extra)) out <- dplyr::bind_cols(out, tibble::as_tibble(extra))
  out
}

RES <- list(); .add <- function(x) RES[[length(RES) + 1L]] <<- x

.fit_f1 <- function(d, cn, ox = "OX", label = "F1 primary", extra_cov = NULL) {
  cv <- c(.covars(d, cn), extra_cov)
  vars <- c("pcr", "MYC", ox, "BUFFER_c", cv)
  dd <- .frame(d, vars)
  f <- stats::as.formula(paste("pcr ~ MYC *", ox, "* BUFFER_c",
                               if (length(cv)) paste("+", paste(cv, collapse = " + "))
                               else ""))
  m <- stats::glm(f, data = dd, family = stats::binomial())
  list(fit = m, three = paste0("MYC:", ox, ":BUFFER_c"),
       two = paste0("MYC:", ox), label = label, cohort = cn)
}

for (cn in COHORTS) {
  z <- .fit_f1(D[[cn]], cn)
  .add(.tidy(z$fit, z$three, "F1 three-way", cn, list(scoring = "per-cohort")))
  .add(.tidy(z$fit, z$two,   "F1 two-way",   cn, list(scoring = "per-cohort")))
  # A2 sensitivity: identical model on the intersection-scored exposure
  zi <- .fit_f1(Di[[cn]], cn)
  .add(.tidy(zi$fit, zi$three, "F1 three-way", cn, list(scoring = "int3")))
  .add(.tidy(zi$fit, zi$two,   "F1 two-way",   cn, list(scoring = "int3")))
}

coefs <- dplyr::bind_rows(RES)
coefs %>% dplyr::filter(label == "F1 three-way") %>%
  dplyr::select(cohort, scoring, n, events, estimate, ci_lo, ci_hi, p) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 5. Meta-analysis - estimates pooled, scores never
# =============================================================================
# Inverse-variance fixed effect, and DerSimonian-Laird random effects, written
# out rather than taken from a package so the arithmetic is auditable:
#   w = 1/v ; theta_FE = sum(w.theta)/sum(w) ; var_FE = 1/sum(w)
#   Q = sum(w (theta - theta_FE)^2) ; df = k-1
#   C = sum(w) - sum(w^2)/sum(w) ; tau2 = max(0, (Q-df)/C)
#   w* = 1/(v+tau2) ; theta_RE = sum(w*.theta)/sum(w*) ; var_RE = 1/sum(w*)
#   I2 = max(0, (Q-df)/Q)
message("\n5. meta-analysis")

.meta <- function(est, se, label, scoring) {
  v <- se^2; w <- 1 / v; k <- length(est)
  crit0 <- stats::qnorm(1 - (1 - CI_LEVEL) / 2)
  if (k < 2L) {
    # With one study C = 0 and tau2 would be 0/0. There is no heterogeneity to
    # estimate; report the single estimate rather than a NaN.
    return(tibble::tibble(
      label = label, scoring = scoring, k = k,
      fe_estimate = est, fe_ci_lo = est - crit0 * se, fe_ci_hi = est + crit0 * se,
      fe_p = 2 * stats::pnorm(-abs(est / se)),
      re_estimate = est, re_ci_lo = est - crit0 * se, re_ci_hi = est + crit0 * se,
      re_p = 2 * stats::pnorm(-abs(est / se)),
      Q = NA_real_, Q_p = NA_real_, tau2 = NA_real_, I2 = NA_real_))
  }
  fe <- sum(w * est) / sum(w); vfe <- 1 / sum(w)
  Q  <- sum(w * (est - fe)^2); df <- k - 1
  Cc <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - df) / Cc)
  ws <- 1 / (v + tau2)
  re <- sum(ws * est) / sum(ws); vre <- 1 / sum(ws)
  crit <- stats::qnorm(1 - (1 - CI_LEVEL) / 2)
  tibble::tibble(
    label = label, scoring = scoring, k = k,
    fe_estimate = fe, fe_ci_lo = fe - crit * sqrt(vfe),
    fe_ci_hi = fe + crit * sqrt(vfe),
    fe_p = 2 * stats::pnorm(-abs(fe / sqrt(vfe))),
    re_estimate = re, re_ci_lo = re - crit * sqrt(vre),
    re_ci_hi = re + crit * sqrt(vre),
    re_p = 2 * stats::pnorm(-abs(re / sqrt(vre))),
    Q = Q, Q_p = stats::pchisq(Q, df, lower.tail = FALSE),
    tau2 = tau2, I2 = 100 * max(0, (Q - df) / Q))
}

meta <- dplyr::bind_rows(lapply(
  split(coefs, list(coefs$label, coefs$scoring), drop = TRUE),
  function(g) .meta(g$estimate, g$se, g$label[1], g$scoring[1])))
meta %>% as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 6. Specificity - I-SPY2 and BrighTNess only (A2)
# =============================================================================
# The assembly factors are this arm's primary specificity control, and a control
# at half coverage cannot do that work. GSE25066 is 34/68 and is excluded from
# this section by the rule declared in the amendment, not by inspection of its
# result.
message("\n6. specificity")

spec_ok <- coverage %>%
  dplyr::filter(set == ARM_NEGATIVE, frac >= SPECIFICITY_MIN_COV) %>%
  dplyr::pull(cohort)
message("   cohorts qualifying at coverage >= ", SPECIFICITY_MIN_COV, ": ",
        paste(spec_ok, collapse = ", "),
        "   |  excluded: ", paste(setdiff(COHORTS, spec_ok), collapse = ", "))

spec <- list()
for (cn in spec_ok) {
  z <- .fit_f1(D[[cn]], cn, ox = "OX_neg", label = "F1 three-way NEGATIVE arm")
  spec[[length(spec) + 1L]] <-
    .tidy(z$fit, z$three, "F1 three-way NEGATIVE arm", cn,
          list(scoring = "per-cohort"))
}
spec <- dplyr::bind_rows(spec)
spec %>% dplyr::select(cohort, n, events, estimate, ci_lo, ci_hi, p) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 7. F1 secondary - the frozen STATE contrast
# =============================================================================
# Level 3 (unbuffered) vs level 4 (buffered), one degree of freedom, direction
# fixed by plan section 7.5: pCR HIGHER in level 3. Reported alongside the
# continuous primary, never instead of it (D5 6.1).
message("\n7. F1 secondary - STATE")

sec <- list()
for (cn in COHORTS) {
  d <- D[[cn]]
  cv <- .covars(d, cn)
  dd <- .frame(d, c("pcr", "STATE", cv))
  # Level 3 is the reference here so the contrast is a single coefficient.
  # Re-levelled so the contrast is a single coefficient, then droplevels():
  # re-factoring against the full frozen level list would reintroduce cells that
  # are empty after complete cases, and glm returns NA coefficients for those.
  dd$STATE <- droplevels(stats::relevel(
    factor(dd$STATE, levels = st$spec$levels), ref = st$spec$contrast[1]))
  if (!st$spec$contrast[2] %in% levels(dd$STATE)) {
    message("   ", cn, ": level 4 empty after complete cases - contrast skipped")
    next
  }
  f <- stats::as.formula(paste("pcr ~ STATE",
                               if (length(cv)) paste("+", paste(cv, collapse = " + "))
                               else ""))
  m <- stats::glm(f, data = dd, family = stats::binomial())
  sec[[length(sec) + 1L]] <-
    .tidy(m, paste0("STATE", st$spec$contrast[2]),
          "F1 secondary level4 vs level3", cn,
          list(cells = paste(as.integer(table(dd$STATE)), collapse = "/")))
}
sec <- dplyr::bind_rows(sec)
if (nrow(sec)) {
  sec %>% dplyr::select(cohort, n, events, cells, estimate, ci_lo, ci_hi, or, p) %>%
    as.data.frame() %>% print(row.names = FALSE)
  meta_sec <- .meta(sec$estimate, sec$se, "F1 secondary", "per-cohort")
} else {
  meta_sec <- NULL
}

# =============================================================================
# 7b. Treatment stratification - DESCRIPTIVE, bounded by C3
# =============================================================================
# C3: I-SPY2's 179-patient control arm is 94 HRpos_HER2neg + 85 TNBC and holds
# ZERO HER2+ patients, so a control-versus-experimental contrast is estimable
# only in those two strata.
#
# NOT FITTED AS A FOUR-WAY. `MYC * OXPHOS * BUFFER_c * treatment` at 179 control
# patients with about 31 events is not an estimate, it is a number with a CI
# spanning the parameter space. The three-way is fitted SEPARATELY within
# control and within the rest and the two are reported side by side, labelled
# descriptive. No prediction attaches to this section and it is excluded from
# section 8.
message("\n7b. treatment stratification (descriptive)")

trt <- list()
d <- D[[PRIMARY_COHORT]] %>%
  dplyr::filter(subtype %in% c("HRpos_HER2neg", "TNBC"))
d$stratum <- ifelse(d$treatment == "Paclitaxel", "control", "experimental")
for (s in c("control", "experimental")) {
  ds <- d[d$stratum == s, ]
  r <- try({
    z <- .fit_f1(ds, paste0(PRIMARY_COHORT, " ", s))
    .tidy(z$fit, z$three, "F1 three-way DESCRIPTIVE by arm",
          paste0(PRIMARY_COHORT, " ", s), list(scoring = "per-cohort"))
  }, silent = TRUE)
  if (inherits(r, "try-error")) {
    message("   ", s, " (n = ", nrow(ds), "): did not fit - ",
            trimws(as.character(r))); next
  }
  trt[[length(trt) + 1L]] <- r
}
trt <- dplyr::bind_rows(trt)
if (nrow(trt)) {
  trt %>% dplyr::select(cohort, n, events, estimate, ci_lo, ci_hi, p) %>%
    as.data.frame() %>% print(row.names = FALSE)
}

# =============================================================================
# 8. Declared readings - evaluated mechanically
# =============================================================================
# Directions fixed in the 2026-08-30 declaration section 6.5 and the informative
# failure in 6.6, both before any of these numbers existed.
message("\n8. declared readings")

.sign_ok <- function(e, lo, hi, want) {
  if (want == "negative") hi < 0 else lo > 0
}

m3 <- meta %>% dplyr::filter(label == "F1 three-way", scoring == "per-cohort")
c3p <- coefs %>% dplyr::filter(label == "F1 three-way", scoring == "per-cohort",
                               cohort == PRIMARY_COHORT)
m3i <- meta %>% dplyr::filter(label == "F1 three-way", scoring == "int3")

# 6.6: the OXPHOS slope within MYC-high, evaluated at MYC = +1 SD, in the
# primary cohort. d/dOX at MYC = 1, BUFFER_c = 0  =  b_OX + b_MYC:OX
.slope_at_myc <- function(cn) {
  z <- .fit_f1(D[[cn]], cn)
  b <- stats::coef(z$fit); V <- stats::vcov(z$fit)
  k <- rep(0, length(b)); names(k) <- names(b)
  k[["OX"]] <- 1; k[[z$two]] <- 1
  est <- sum(k * b); se <- sqrt(drop(t(k) %*% V %*% k))
  crit <- stats::qnorm(1 - (1 - CI_LEVEL) / 2)
  tibble::tibble(cohort = cn, term = "d pCR / d OXPHOS at MYC = +1 SD",
                 estimate = est, se = se, ci_lo = est - crit * se,
                 ci_hi = est + crit * se,
                 p = 2 * stats::pnorm(-abs(est / se)))
}
slopes <- dplyr::bind_rows(lapply(COHORTS, .slope_at_myc))
slopes %>% as.data.frame() %>% print(row.names = FALSE)

sl <- slopes[slopes$cohort == PRIMARY_COHORT, ]
three_null <- !.sign_ok(m3$fe_estimate, m3$fe_ci_lo, m3$fe_ci_hi, "negative") &&
              !.sign_ok(m3$fe_estimate, m3$fe_ci_lo, m3$fe_ci_hi, "positive")

declared <- tibble::tribble(
  ~test, ~predicted, ~passes,
  "F1-a  meta three-way MYC:OXPHOS:BUFFER_c",
  "NEGATIVE, CI excludes 0",
  .sign_ok(m3$fe_estimate, m3$fe_ci_lo, m3$fe_ci_hi, PRED_THREEWAY),

  "F1-b  same sign in the primary cohort",
  "negative in GSE194040",
  isTRUE(c3p$estimate < 0),

  "F1-c  meta two-way MYC:OXPHOS at mean buffering",
  "POSITIVE, CI excludes 0",
  .sign_ok(meta$fe_estimate[meta$label == "F1 two-way" &
                              meta$scoring == "per-cohort"],
           meta$fe_ci_lo[meta$label == "F1 two-way" &
                           meta$scoring == "per-cohort"],
           meta$fe_ci_hi[meta$label == "F1 two-way" &
                           meta$scoring == "per-cohort"], PRED_TWOWAY),

  "F1-d  specificity: NOT the assembly factors",
  "null / weaker, qualifying cohorts",
  # NA, not TRUE, when no cohort qualifies. all(logical(0)) is TRUE, and a
  # specificity claim that passes because it was never tested is the worst
  # possible reading of this row.
  if (nrow(spec)) all(spec$ci_lo < 0 & spec$ci_hi > 0) else NA,

  "F1-e  A2 sensitivity agrees (int3 scoring)",
  "same sign as F1-a",
  isTRUE(sign(m3i$fe_estimate) == sign(m3$fe_estimate)),

  "F1-f  secondary STATE level4 vs level3",
  "level 3 higher pCR, i.e. NEGATIVE coefficient",
  if (!is.null(meta_sec))
    .sign_ok(meta_sec$fe_estimate, meta_sec$fe_ci_lo, meta_sec$fe_ci_hi,
             "negative") else NA,

  "F1-g  INFORMATIVE FAILURE (plan section 2, criterion 4)",
  "OXPHOS slope at MYC=+1SD NEGATIVE, CI excl 0, WHILE three-way null",
  isTRUE(sl$ci_hi < 0) && isTRUE(three_null))
declared %>% as.data.frame() %>% print(row.names = FALSE)

message("\n   REMINDER: H4 is SINGLE-INSTRUMENT (A1). GSVA only; mitoPPS does")
message("   not exist in these cohorts. This result does NOT meet the")
message("   two-instrument standard Blocks C, B and G were held to, and must")
message("   not be written as though it did. A null is 'not powered to exclude")
message("   a modest effect' (D5 6.3), NOT falsification - only F1-g is that.")

# =============================================================================
# 9. Save
# =============================================================================
message("\n9. save")

out <- list(
  coefficients = coefs, meta = meta, specificity = spec,
  secondary = sec, meta_secondary = meta_sec,
  treatment_descriptive = trt, slopes = slopes,
  coverage = coverage, declared = declared,
  scores = lapply(SC, `[[`, "scores"),
  frames = D, frames_int3 = Di,
  spec_cohorts = spec_ok,
  spec_rules = list(
    instrument = paste("GSVA ONLY (amendment A1). mitoPPS does not exist in",
                       "these cohorts. Single-instrument result."),
    model      = "pcr ~ MYC * OXPHOS * BUFFER_c + subtype + treatment, logistic",
    covariates = paste("A4: subtype + treatment only. Stage and grade do not",
                       "exist in the primary cohort."),
    buffer     = "mean(z(log2 MCL1), z(log2 BCL2L1)) within cohort (6.1)",
    directions = paste("three-way NEGATIVE, two-way POSITIVE, fixed 2026-08-30",
                       "section 6.5 before any of this existed"),
    scoring    = paste("one GSVA run per cohort, pinned universe; scores never",
                       "pooled, estimates meta-analysed (D5 6.2)"),
    na_rule    = "A5: drop genes >5% NA, impute rest at gene median",
    specificity = paste("A2: assembly-factor coverage >= 0.80;",
                        "GSE25066 (0.50) excluded"),
    not_done   = paste("F2, F3 proper and F4 are NOT in this script - no",
                       "METABRIC/SCAN-B, forkscale does not exist in these",
                       "cohorts, and F4's covariates exist only in GSE25066")),
  built = Sys.time())

saveRDS(out, file.path(DIR_RESULTS, "h4_outcome_models.rds"))
readr::write_csv(coefs,    file.path(DIR_TABLES, "h4_coefficients.csv"))
readr::write_csv(meta,     file.path(DIR_TABLES, "h4_meta.csv"))
readr::write_csv(declared, file.path(DIR_TABLES, "h4_declared.csv"))
readr::write_csv(coverage, file.path(DIR_TABLES, "h4_coverage.csv"))

message("\n13: done.")
message("    results/h4_outcome_models.rds")
message("    outputs/tables/  4 tables")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  h <- readRDS(file.path(DIR_RESULTS, "h4_outcome_models.rds"))

  # --- the declared readings, first and alone ------------------------------
  h$declared %>% as.data.frame() %>% print(row.names = FALSE)
  h$meta %>% as.data.frame() %>% print(row.names = FALSE)

  # --- the forest the meta-analysis is made of -----------------------------
  fo <- h$coefficients %>%
    dplyr::filter(label == "F1 three-way", scoring == "per-cohort")
  with(fo, {
    plot(estimate, seq_along(estimate), pch = 16, yaxt = "n", ylab = "",
         xlim = range(c(ci_lo, ci_hi)), xlab = "MYC:OXPHOS:BUFFER_c (log OR)")
    segments(ci_lo, seq_along(estimate), ci_hi, seq_along(estimate))
    axis(2, seq_along(cohort), cohort, las = 1, cex.axis = 0.7)
    abline(v = 0, lty = 2)
  })

  # --- does the exposure behave at all? ------------------------------------
  # Before reading any interaction, check the main effects are not degenerate.
  # A three-way built on a score with no marginal signal is not interpretable.
  for (cn in names(h$frames)) {
    d <- h$frames[[cn]]
    cat("\n", cn, "\n")
    print(summary(stats::glm(pcr ~ MYC + OX + BUFFER_c, d,
                             family = stats::binomial()))$coefficients)
  }

  # --- is BUFFER_c just one of its two genes? ------------------------------
  # The composite was chosen over max() on stated grounds; check the two limbs
  # are not simply the same variable.
  lapply(h$frames, function(d) stats::cor(d$z_MCL1, d$z_BCL2L1, method = "spearman"))

  # --- per-cohort vs int3 scoring ------------------------------------------
  # If these disagree, the exposure is sensitive to set composition and the A2
  # sensitivity is doing real work.
  h$coefficients %>%
    dplyr::filter(label == "F1 three-way") %>%
    dplyr::select(cohort, scoring, estimate, ci_lo, ci_hi) %>%
    as.data.frame() %>% print(row.names = FALSE)

  # --- subtype is the strongest predictor of pCR here; confirm it -----------
  d <- h$frames$GSE194040
  print(table(d$subtype, d$pcr))

  # --- STATE cells, and whether the contrast had anything to work with -----
  lapply(h$frames, function(d) table(d$STATE, d$pcr, useNA = "ifany"))

}
