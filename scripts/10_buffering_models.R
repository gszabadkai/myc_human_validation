# 10_buffering_models.R
# =============================================================================
# BLOCK B - the buffering test, and H1's SECOND clause.
#
# Built to docs/2026-08-29_escape_reading_declaration.md, in which every
# direction below was fixed before any of these models existed. Read that first;
# this script implements it and re-decides none of it.
#
# =============================================================================
# WHAT QUESTION THIS ANSWERS
# =============================================================================
# Plan section 1, written before any data was touched:
#
#   "If OXPHOS-high plus MYC-high is apoptotically lethal in untransformed
#    mammary epithelium, HOW DID THE HUMAN TUMOURS THAT OCCUPY THAT STATE
#    SURVIVE IT, and what did that escape cost them therapeutically?"
#
# Block C tested the coupling and it does not hold on the pre-registered
# endpoint. This script tests the escape. They are different questions and this
# one does NOT resurrect H1 - H1 is a conjunction whose first clause failed, and
# no output here may be written as if it were otherwise.
#
# =============================================================================
# WHY THIS IS NOT JUST BLOCK C WITH X AND Y SWAPPED
# =============================================================================
# Reversing a regression on expression data establishes nothing: `OXPHOS ~
# priming` and `priming ~ OXPHOS` carry identical information and identical
# confounding.
#
# `BUFFER` is not an expression level. It is SOMATIC COPY NUMBER
# (`cnv_MCL1 == +2` OR `cnv_BCL2L1 >= +1`, script 03). OXPHOS expression cannot
# cause an MCL1 amplification, so a model with BUFFER as the PREDICTOR is a
# regression on a genuinely upstream variable. That asymmetry is the entire
# methodological content of E2.
#
# And it is fragile in one specific way, so three things are mandatory rather
# than optional: BUFFER is 512 of 938 and is dominated by a BROAD 20q gain
# (`buffer_BCL2L1` 439) rather than focal amplification (`buffer_MCL1` 158). So
# every model adjusts for ANEUPLOIDY, every BUFFER-as-predictor model adjusts for
# MYC AMPLIFICATION (G2's own result is that the two co-occur), and E2 is
# calibrated against the expression-matched null - because if BUFFER predicts
# every gene set equally, the null shifts with it and OXPHOS will not stand out.
#
# SPECIES: human. See CLAUDE.md.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

message("\n10: Block B - buffering, and H1's second clause\n", strrep("=", 78))

# -----------------------------------------------------------------------------
# 0. Constants
# -----------------------------------------------------------------------------
ARM_PRIMARY  <- "OXPHOS subunits"
ARM_NEGATIVE <- "OXPHOS assembly factors"
INSTRUMENTS  <- c("gsva", "mitopps")
BLOCK_C_VARS <- c("purity", "leukocyte_fraction", "PAM50", "TP53_status", "plate")
CI_LEVEL     <- 0.95
P_EMP_ALPHA  <- 0.05

# =============================================================================
# 1. Inputs - consumed as built
# =============================================================================
message("\n1. inputs")

mito <- readRDS(file.path(DIR_RESULTS, "tcga_brca_mito_scores.rds"))
prim <- readRDS(file.path(DIR_RESULTS, "tcga_brca_priming.rds"))$priming
myc  <- readRDS(file.path(DIR_RESULTS, "tcga_brca_myc_scores.rds"))$estimators
cov  <- readRDS(file.path(DIR_RESULTS, "tcga_brca_covariates.rds"))$covariates
bc   <- readRDS(file.path(DIR_RESULTS, "block_c_models.rds"))

if (!all(c("log2_MCL1", "log2_BCL2") %in% names(prim))) {
  stop("script 08 has not been re-run since MCL1 and BCL2 were added to the ",
       "saved limbs. Re-source scripts/08_score_priming.R, then this script.",
       call. = FALSE)
}

pat <- colnames(mito$gsva_arms)
stopifnot(identical(sort(pat), sort(cov$patient)),
          identical(sort(pat), sort(prim$patient)),
          identical(sort(pat), sort(myc$patient)))

# =============================================================================
# 2. The analysis frame - identical to Block C's, by construction
# =============================================================================
# The frozen z-scaling constants and the plate map are READ FROM script 09's
# output rather than recomputed. Recomputing them would give the same numbers
# today and might not tomorrow; reading them makes every coefficient here
# directly comparable to Block C's, which is the point of having frozen them.
message("\n2. analysis frame")

D <- tibble::tibble(patient = pat) %>%
  dplyr::left_join(dplyr::select(cov, patient, purity, leukocyte_fraction,
                                 PAM50, TP53_status, plate, aneuploidy_cbio,
                                 BUFFER, buffer_MCL1, buffer_BCL2L1, MYC_amp),
                   by = "patient") %>%
  dplyr::left_join(dplyr::select(myc, patient, M_a, M_b), by = "patient") %>%
  dplyr::left_join(dplyr::select(prim, patient, log2_MCL1, log2_BCL2,
                                 log2_BCL2L1), by = "patient") %>%
  dplyr::mutate(
    PROLIF_DISJOINT = as.numeric(mito$gsva_cov["PROLIF_DISJOINT", patient]),
    PAM50           = factor(PAM50),
    TP53_status     = factor(TP53_status),
    BUFFER_amp      = as.numeric(BUFFER),
    MYC_amp_n       = as.numeric(MYC_amp)
  )

idx <- stats::complete.cases(D[, BLOCK_C_VARS]) &
  !is.na(D$BUFFER_amp) & !is.na(D$aneuploidy_cbio) & !is.na(D$MYC_amp_n)
message("   analysis set: n = ", sum(idx),
        " (Block C's 938 with BUFFER and aneuploidy present)")

# --- BUFFER's composition, reported with every result ------------------------
message(sprintf(
  paste("   BUFFER: %d of %d TRUE  |  buffer_BCL2L1 (broad 20q gain) %d,",
        "buffer_MCL1 (focal +2) %d"),
  sum(D$BUFFER_amp[idx]), sum(idx), sum(D$buffer_BCL2L1[idx]),
  sum(D$buffer_MCL1[idx])))
message("   MYC amplified: ", sum(D$MYC_amp_n[idx]),
        "  |  BUFFER within MYC-amp: ",
        sum(D$BUFFER_amp[idx] & D$MYC_amp_n[idx]),
        " of ", sum(D$MYC_amp_n[idx]))

# --- frozen scaling and plate map, read from Block C -------------------------
Z <- bc$zscale
.apply_z <- function(v, key) {
  k <- Z[[key]]
  stopifnot(!is.null(k))
  (v - k[["mean"]]) / k[["sd"]]
}
D$MYC <- .apply_z(D$M_a, "M_a")

SC <- list()
for (ins in INSTRUMENTS) {
  m <- mito[[paste0(ins, "_arms")]]
  SC[[ins]] <- t(vapply(rownames(m), function(a)
    .apply_z(as.numeric(m[a, pat]), paste(ins, a, sep = "|")), numeric(length(pat))))
  colnames(SC[[ins]]) <- pat
}
ARMS_ALL  <- rownames(SC$gsva)
ARMS_NULL <- mito$null_manifest$arm

.pool <- function(x) factor(ifelse(x %in% bc$plate_keep, x, "other"))
D$plate_f <- .pool(D$plate)
message("   plate map and ", length(Z), " scaling constants taken from Block C")

# --- an assertion that this really is Block C's frame ------------------------
# If the frame has drifted, every coefficient below is quietly on a different
# footing from the ones it will be compared with.
d_chk <- D[idx, ]
d_chk$OX <- SC$gsva[ARM_PRIMARY, idx]
d_chk$Y  <- prim$PRIME[match(pat, prim$patient)][idx]
.chk  <- stats::lm(Y ~ MYC * OX + purity + leukocyte_fraction +
                     PROLIF_DISJOINT + PAM50 + TP53_status + plate_f, d_chk)
b_chk <- stats::coef(.chk)[["MYC:OX"]]
b_ref <- bc$coefficients$estimate[bc$coefficients$label == "spine" &
                                    bc$coefficients$instrument == "gsva"]
if (!isTRUE(all.equal(b_chk, b_ref, tolerance = 1e-10))) {
  stop("this script's analysis frame does not reproduce Block C's spine ",
       "coefficient (", signif(b_chk, 6), " vs ", signif(b_ref, 6),
       "). The two are not on the same footing; stop and diagnose.",
       call. = FALSE)
}
message("   frame verified: reproduces Block C's spine coefficient exactly")

COVARS <- c("purity", "leukocyte_fraction", "PROLIF_DISJOINT", "PAM50",
            "TP53_status", "plate_f", "aneuploidy_cbio")

.tidy <- function(m, term, label, instrument, n, family = "lm") {
  co <- summary(m)$coefficients
  if (!term %in% rownames(co)) {
    stop("term '", term, "' absent from fit '", label, "'", call. = FALSE)
  }
  crit <- if (family == "lm") stats::qt(1 - (1 - CI_LEVEL) / 2, m$df.residual)
          else stats::qnorm(1 - (1 - CI_LEVEL) / 2)
  tibble::tibble(
    label = label, term = term, instrument = instrument, n = n,
    estimate = co[term, 1], se = co[term, 2],
    ci_lo = co[term, 1] - crit * co[term, 2],
    ci_hi = co[term, 1] + crit * co[term, 2],
    p = co[term, 4],
    or = if (family == "binomial") exp(co[term, 1]) else NA_real_)
}

RES <- list(); .add <- function(x) RES[[length(RES) + 1L]] <<- x

# =============================================================================
# 3. E1 - H1 clause 2, as pre-registered. PRIMARY.
# =============================================================================
# PREDICTED: the MYC_high:OXPHOS_high interaction is POSITIVE. The dangerous
# quadrant carries more buffer than its two main effects predict.
#
# `plate` is deliberately absent: the outcome is a copy-number call from SNP
# arrays and the RNA-seq plate is not a confounder of it.
message("\n3. E1 - is the MYC-high / OXPHOS-high quadrant enriched for BUFFER?")

for (ins in INSTRUMENTS) {
  d <- D[idx, ]
  d$OX <- SC[[ins]][ARM_PRIMARY, idx]
  # splits at the median of the frozen z-scores on this set
  d$MYC_high <- as.numeric(d$MYC > stats::median(d$MYC))
  d$OX_high  <- as.numeric(d$OX  > stats::median(d$OX))

  m1 <- stats::glm(BUFFER_amp ~ MYC_high * OX_high + PAM50 + purity +
                     aneuploidy_cbio, data = d, family = stats::binomial())
  .add(.tidy(m1, "MYC_high:OX_high", "E1 categorical (as written)", ins,
             nrow(d), "binomial"))

  # Continuous companion. D5: median splits cost up to 53 percentage points of
  # power against the identical truth, so the dichotomised form the plan wrote is
  # the weaker instrument. Reported alongside, never instead.
  m2 <- stats::glm(BUFFER_amp ~ MYC * OX + PAM50 + purity + aneuploidy_cbio,
                   data = d, family = stats::binomial())
  .add(.tidy(m2, "MYC:OX", "E1 continuous companion", ins, nrow(d), "binomial"))

  # The descriptive form: the quadrant against everything else, which is what
  # H1 clause 2 says in words.
  d$quadrant <- as.numeric(d$MYC_high == 1 & d$OX_high == 1)
  m3 <- stats::glm(BUFFER_amp ~ quadrant + PAM50 + purity + aneuploidy_cbio,
                   data = d, family = stats::binomial())
  .add(.tidy(m3, "quadrant", "E1 quadrant vs rest (descriptive)", ins,
             nrow(d), "binomial"))
  message(sprintf("   %-8s quadrant n = %d of %d", ins, sum(d$quadrant), nrow(d)))
}

# =============================================================================
# 4. E2 - the specificity form. THE DECLARED ADDITION.
# =============================================================================
# arm_score ~ BUFFER_amp + covariates, all 18 arms, both instruments.
# PREDICTED: POSITIVE for OXPHOS subunits and beating its matched null; NOT for
# the comparator arms, and in particular NOT for OXPHOS assembly factors.
#
# The design matrix is FIXED and only the outcome changes, so the whole null -
# 17 arms x 2,000 draws x 2 instruments - is a handful of QR solves rather than
# 68,000 model fits.
message("\n4. E2 - does BUFFER amplification predict OXPHOS specifically?")

d2 <- D[idx, ]
X  <- stats::model.matrix(
  stats::as.formula(paste("~ BUFFER_amp + MYC_amp_n +",
                          paste(COVARS, collapse = " + "))), data = d2)
QR    <- qr(X)
j_buf <- which(colnames(X) == "BUFFER_amp")
stopifnot(length(j_buf) == 1L, QR$rank == ncol(X), nrow(X) == nrow(d2))
n_obs <- nrow(X); p_obs <- ncol(X)
message("   design: ", n_obs, " x ", p_obs,
        " (BUFFER + MYC amplification + Block C covariates + aneuploidy)")

.z_own <- function(v) (v - mean(v)) / stats::sd(v)

# observed, via lm so the SE is exact; the QR path is asserted against it below
e2_obs <- dplyr::bind_rows(lapply(INSTRUMENTS, function(ins) {
  dplyr::bind_rows(lapply(ARMS_ALL, function(a) {
    dd <- d2; dd$Y <- .z_own(SC[[ins]][a, idx])
    m  <- stats::lm(stats::as.formula(paste("Y ~ BUFFER_amp + MYC_amp_n +",
                                            paste(COVARS, collapse = " + "))), dd)
    .tidy(m, "BUFFER_amp", paste("E2 arm:", a), ins, nrow(dd)) %>%
      dplyr::mutate(arm = a)
  }))
}))

b_qr <- qr.coef(QR, .z_own(SC$gsva[ARM_PRIMARY, idx]))[j_buf]
b_lm <- e2_obs$estimate[e2_obs$arm == ARM_PRIMARY & e2_obs$instrument == "gsva"]
if (!isTRUE(all.equal(unname(b_qr), b_lm, tolerance = 1e-10))) {
  stop("the QR path and lm() disagree on E2 (", b_qr, " vs ", b_lm,
       "); every null percentile would be off the observed model.", call. = FALSE)
}
message("   QR path verified against lm()")

# --- the matched null --------------------------------------------------------
e2_null <- dplyr::bind_rows(lapply(seq_along(ARMS_NULL), function(k) {
  a  <- ARMS_NULL[[k]]
  nf <- readRDS(mito$null_manifest$file[mito$null_manifest$arm == a])
  dplyr::bind_rows(lapply(INSTRUMENTS, function(ins) {
    Y <- t(nf[[ins]][, pat[idx], drop = FALSE])      # samples x draws
    Y <- apply(Y, 2L, .z_own)
    b <- qr.coef(QR, Y)[j_buf, ]
    b <- b[is.finite(b)]
    o <- e2_obs$estimate[e2_obs$arm == a & e2_obs$instrument == ins]
    tibble::tibble(
      arm = a, instrument = ins, n_draws = length(b), b_obs = o,
      null_median = stats::median(b),
      percentile = 100 * mean(b < o),
      p_emp = max(2 * min(mean(b <= o), mean(b >= o)), 1 / length(b)))
  }))
}))
message("   null: ", nrow(e2_null), " arm x instrument combinations")

# =============================================================================
# 5. E3 - Block B's expression models
# =============================================================================
# log2(MCL1)   ~ MYC * OXPHOS   PREDICTED POSITIVE, stated with low confidence:
#   buffer_MCL1 is focal high-level (+2, n = 158) where buffer_BCL2L1 is a broad
#   gain, so if any expression-level buffering survives adjustment it is MCL1's.
#
# log2(BCL2L1) ~ MYC * OXPHOS   ALREADY FITTED as script 09's `limb: log2_BCL2L1`
#   (-0.034 GSVA, -0.043 mitoPPS). NO PREDICTION IS DECLARED FOR IT. It is
#   restated below from Block C's own object rather than refitted, so nobody can
#   read it as a fresh confirmation.
message("\n5. E3 - the expression-level buffer models")

for (ins in INSTRUMENTS) {
  d <- D[idx, ]; d$OX <- SC[[ins]][ARM_PRIMARY, idx]
  for (y in c("log2_MCL1", "log2_BCL2")) {
    d$Y <- d[[y]]
    m <- stats::lm(stats::as.formula(paste("Y ~ MYC * OX +",
                                           paste(setdiff(COVARS, "aneuploidy_cbio"),
                                                 collapse = " + "))), d)
    .add(.tidy(m, "MYC:OX", paste("E3", y), ins, nrow(d)))
  }
}

e3_known <- bc$coefficients %>%
  dplyr::filter(endpoint == "log2_BCL2L1") %>%
  dplyr::transmute(label = "E3 log2_BCL2L1 (ALREADY FITTED in Block C)",
                   term = "MYC:OX", instrument, n, estimate, se, ci_lo, ci_hi,
                   p, or = NA_real_)
.add(e3_known)

# =============================================================================
# 6. E4 - the one prediction that ties Block C to the escape reading
# =============================================================================
# PREDICTED: the three-way term is POSITIVE - the MYC:OXPHOS down-regulation of
# BCL-XL found in Block C is ATTENUATED in BUFFER-amplified tumours, because
# amplification sets a floor the transcriptional effect has less room beneath.
#
# POWER CAVEAT, fixed before the fit: a three-way term on a 512/426 split is
# underpowered by construction. A NULL HERE IS UNINFORMATIVE and is reported as
# such, never as evidence against.
message("\n6. E4 - is the BCL-XL effect attenuated where BCL-XL is amplified?")

for (ins in INSTRUMENTS) {
  d <- D[idx, ]; d$OX <- SC[[ins]][ARM_PRIMARY, idx]
  m <- stats::lm(stats::as.formula(paste(
    "log2_BCL2L1 ~ MYC * OX * BUFFER_amp +",
    paste(setdiff(COVARS, "aneuploidy_cbio"), collapse = " + "))), d)
  .add(.tidy(m, "MYC:OX:BUFFER_amp", "E4 three-way", ins, nrow(d)))
}

coefs <- dplyr::bind_rows(RES)

# =============================================================================
# 7. The declared predictions, evaluated
# =============================================================================
# Directions were fixed in docs/2026-08-29_escape_reading_declaration.md before
# any of the above existed. Evaluated mechanically so the reading cannot drift.
message("\n7. the declared predictions")

.both <- function(v) length(v) == 2L && all(v)

e1 <- coefs %>% dplyr::filter(label == "E1 categorical (as written)")
p1 <- .both(e1$estimate > 0 & e1$p < 0.05)

o  <- e2_obs %>% dplyr::filter(arm == ARM_PRIMARY)
n_ <- e2_null %>% dplyr::filter(arm == ARM_PRIMARY)
p2a <- .both(o$estimate > 0 & o$p < 0.05) && .both(n_$p_emp < P_EMP_ALPHA)

oa  <- e2_obs %>% dplyr::filter(arm == ARM_NEGATIVE)
na_ <- e2_null %>% dplyr::filter(arm == ARM_NEGATIVE)
p2b <- !(.both(oa$estimate > 0 & oa$p < 0.05) && .both(na_$p_emp < P_EMP_ALPHA))

e3 <- coefs %>% dplyr::filter(label == "E3 log2_MCL1")
p3 <- .both(e3$estimate > 0 & e3$p < 0.05)

e4 <- coefs %>% dplyr::filter(label == "E4 three-way")
p4 <- .both(e4$estimate > 0 & e4$p < 0.05)

declared <- tibble::tibble(
  test = c("E1  quadrant enriched for BUFFER (H1 clause 2)",
           "E2a BUFFER predicts OXPHOS subunits, beats matched null",
           "E2b and NOT the assembly factors (specificity)",
           "E3  MCL1 expression interaction positive",
           "E4  BCL-XL effect attenuated where amplified (underpowered)"),
  predicted = c("positive", "positive", "null / weaker", "positive", "positive"),
  passes = c(p1, p2a, p2b, p3, p4))
declared %>% as.data.frame() %>% print(row.names = FALSE)

message("\n   E1, all three forms:")
coefs %>% dplyr::filter(grepl("^E1", label)) %>%
  dplyr::select(label, instrument, n, estimate, or, ci_lo, ci_hi, p) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   E2, arms ordered by coefficient (GSVA):")
e2_obs %>% dplyr::filter(instrument == "gsva") %>%
  dplyr::left_join(dplyr::select(e2_null, arm, instrument, percentile, p_emp),
                   by = c("arm", "instrument")) %>%
  dplyr::select(arm, estimate, ci_lo, ci_hi, p, percentile, p_emp) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3))) %>%
  dplyr::arrange(dplyr::desc(estimate)) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   E3 and E4:")
coefs %>% dplyr::filter(grepl("^E3|^E4", label)) %>%
  dplyr::select(label, instrument, n, estimate, ci_lo, ci_hi, p) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   REMINDER: H1 is a conjunction whose first clause failed. Nothing ",
        "here\n   resurrects it. This answers the arm's organising question, ",
        "which is a\n   different and larger thing.")

# =============================================================================
# 8. Save
# =============================================================================
message("\n8. save")

out <- list(
  coefficients = coefs,
  e2_observed  = e2_obs,
  e2_null      = e2_null,
  declared     = declared,
  buffer_composition = tibble::tibble(
    n = sum(idx), BUFFER = sum(D$BUFFER_amp[idx]),
    buffer_BCL2L1_broad_gain = sum(D$buffer_BCL2L1[idx]),
    buffer_MCL1_focal_amp = sum(D$buffer_MCL1[idx]),
    MYC_amplified = sum(D$MYC_amp_n[idx])),
  spec = list(
    document = "docs/2026-08-29_escape_reading_declaration.md",
    e3_bcl2l1 = paste("NOT a prediction - already fitted as script 09's",
                      "limb: log2_BCL2L1; restated here from Block C"),
    e4_power = paste("three-way on a 512/426 split; a null is UNINFORMATIVE,",
                     "not evidence against"),
    h1 = "H1 is a conjunction whose first clause failed; this does not resurrect it"),
  built = Sys.time()
)

saveRDS(out, file.path(DIR_RESULTS, "block_b_models.rds"))
readr::write_csv(coefs,    file.path(DIR_TABLES, "block_b_coefficients.csv"))
readr::write_csv(e2_obs,   file.path(DIR_TABLES, "block_b_e2_observed.csv"))
readr::write_csv(e2_null,  file.path(DIR_TABLES, "block_b_e2_null.csv"))
readr::write_csv(declared, file.path(DIR_TABLES, "block_b_declared.csv"))

message("\n10: done.")
message("    results/block_b_models.rds")
message("    outputs/tables/  4 tables")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  b <- readRDS(file.path(DIR_RESULTS, "block_b_models.rds"))

  # --- the declared predictions, and nothing else first --------------------
  b$declared %>% as.data.frame() %>% print(row.names = FALSE)
  b$buffer_composition %>% as.data.frame() %>% print(row.names = FALSE)

  # --- E1: does the dichotomy or the continuous form carry it? -------------
  # D5's lesson was that median splits cost power. If the continuous companion
  # is stronger than the categorical form, that is the explanation, not a
  # discrepancy.
  b$coefficients %>% dplyr::filter(grepl("^E1", label)) %>% as.data.frame()

  # --- E2 as a forest, which is the Panel b candidate ----------------------
  fo <- b$e2_observed %>% dplyr::filter(instrument == "gsva") %>%
    dplyr::arrange(estimate)
  with(fo, {
    plot(estimate, seq_along(estimate), pch = 16, yaxt = "n", ylab = "",
         xlim = range(c(ci_lo, ci_hi)), xlab = "BUFFER amplification, per SD of arm")
    segments(ci_lo, seq_along(estimate), ci_hi, seq_along(estimate))
    axis(2, seq_along(arm), arm, las = 1, cex.axis = 0.6); abline(v = 0, lty = 2)
  })

  # --- the specificity question, which is the whole point of E2 ------------
  # If the assembly factors move as much as the subunits, the association is
  # compartment-wide and is NOT an OXPHOS-specific constraint. Declaration
  # section 7 says: report and stop.
  b$e2_null %>%
    dplyr::filter(arm %in% c("OXPHOS subunits", "OXPHOS assembly factors")) %>%
    as.data.frame() %>% print(row.names = FALSE)

  # --- is E2 just aneuploidy? ----------------------------------------------
  # The model adjusts for it, but look at the raw gradient too.
  cv <- readRDS(file.path(DIR_RESULTS, "tcga_brca_covariates.rds"))$covariates
  m  <- readRDS(file.path(DIR_RESULTS, "tcga_brca_mito_scores.rds"))
  i  <- match(colnames(m$gsva_arms), cv$patient)
  boxplot(m$gsva_arms["OXPHOS subunits", ] ~ cv$BUFFER[i],
          xlab = "BUFFER amplification", ylab = "OXPHOS subunits (GSVA)")
  plot(cv$aneuploidy_cbio[i], m$gsva_arms["OXPHOS subunits", ], pch = 16, cex = 0.3,
       xlab = "aneuploidy score", ylab = "OXPHOS subunits (GSVA)")

  # --- the two BUFFER components behave differently? -----------------------
  # buffer_BCL2L1 is a broad 20q gain (n 439); buffer_MCL1 is focal +2 (n 158).
  # If the effect is carried by the broad one, say so.
  table(cv$buffer_MCL1[i], cv$buffer_BCL2L1[i], dnn = c("MCL1 +2", "BCL2L1 >=+1"))

  # --- E4, with its power caveat kept in view ------------------------------
  b$coefficients %>% dplyr::filter(grepl("^E4", label)) %>% as.data.frame()

}
