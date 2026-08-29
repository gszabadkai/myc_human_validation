# 09_interaction_models.R
# =============================================================================
# BLOCK C - the primary interaction test. This is the script the whole human arm
# exists to run.
#
# Built to docs/2026-08-29_script09_build_spec.md, signed off before any
# coefficient was computed. Read that first; this script implements it and does
# not re-decide any of it.
#
# =============================================================================
# THE DESIGN PRINCIPLE, IN ONE PARAGRAPH
# =============================================================================
# There are TWO primary fits and no others - one per co-primary instrument.
# Every other fit varies EXACTLY ONE dimension from that spine and is reported
# beside it, never selected between. That is why almost every fit below is a
# single call to .fit() with one argument changed from its default: the code
# encodes the discipline rather than relying on the reader to check it.
#
#   SPINE   PRIME ~ M_a * OXPHOS
#                 + purity + leukocyte_fraction + PROLIF_DISJOINT
#                 + PAM50 + TP53_status + plate_pooled          n = 938
#
#   OXPHOS = `OXPHOS subunits` on GSVA   -> spine G
#   OXPHOS = `OXPHOS subunits` on mitoPPS -> spine M
#
# PRIMARY TEST: the MYC:OX coefficient, two-sided, on BOTH spines. Report both;
# claim only what BOTH support. An effect on one instrument alone is
# instrument-dependent and is NOT a positive result.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
# No scoring of any kind. Exposures and endpoints are consumed exactly as scripts
# 06, 07 and 08 built them; not one gene is read here. No Block B / BUFFER
# (script 10), no STATE (11), no outcome models (13), no other cohort, no figure.
# Panel a is drawn in script 18 from this script's saved tables.
#
# No set.seed: this script draws nothing. All randomness lives in script 07 and
# is already frozen in results/mito_null/.
#
# SPECIES: human. See CLAUDE.md.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

message("\n09: Block C - the primary interaction test\n", strrep("=", 78))

# -----------------------------------------------------------------------------
# 0. Constants fixed by the spec
# -----------------------------------------------------------------------------
ARM_PRIMARY   <- "OXPHOS subunits"
ARM_NEGATIVE  <- "OXPHOS assembly factors"       # the primary negative
INSTRUMENTS   <- c("gsva", "mitopps")            # co-primary, in this order
ENDPOINT_MAIN <- "PRIME"

BLOCK_C_VARS  <- c("purity", "leukocyte_fraction", "PAM50", "TP53_status", "plate")
PLATE_MIN_N   <- 10L                             # D9
CI_LEVEL      <- 0.95
P_EMP_ALPHA   <- 0.05                            # the null's pass threshold

# Endpoint negatives and the limbs behind them (spec section 4).
NEG_ENDPOINTS <- c("BID_over_BCL2L1", "BAX_over_BCL2L1",
                   "BCL2L11_over_BCL2L1", "BAK1_over_BCL2L1")
LIMBS_PRIME   <- c("log2_BBC3", "log2_BCL2L1")
LIMBS_NEG     <- c("log2_BID", "log2_BAX", "log2_BCL2L11", "log2_BAK1")

# The three paired-null contrasts, all against the claim arm (spec 3.1).
PAIRED_CONTRASTS <- list(
  c(ARM_PRIMARY, ARM_NEGATIVE),
  c(ARM_PRIMARY, "Mitochondrial ribosome"),
  c(ARM_PRIMARY, "Nucleotide metabolism")
)

# =============================================================================
# 1. Inputs - consumed as built
# =============================================================================
message("\n1. inputs")

mito <- readRDS(file.path(DIR_RESULTS, "tcga_brca_mito_scores.rds"))
prim <- readRDS(file.path(DIR_RESULTS, "tcga_brca_priming.rds"))$priming
myc  <- readRDS(file.path(DIR_RESULTS, "tcga_brca_myc_scores.rds"))$estimators
cov  <- readRDS(file.path(DIR_RESULTS, "tcga_brca_covariates.rds"))$covariates

SCORES <- list(gsva = mito$gsva_arms, mitopps = mito$mitopps_arms)
stopifnot(ARM_PRIMARY %in% rownames(SCORES$gsva),
          ARM_PRIMARY %in% rownames(SCORES$mitopps))

ARMS_ALL  <- rownames(SCORES$gsva)
ARMS_NULL <- mito$null_manifest$arm            # 17; mtDNA is deliberately absent
message("   arms: ", length(ARMS_ALL), " scored, ", length(ARMS_NULL),
        " null-tested (", setdiff(ARMS_ALL, ARMS_NULL), " has no percentile)")

# Every table must be on the same patients, in the same order.
pat <- colnames(SCORES$gsva)
stopifnot(identical(sort(pat), sort(cov$patient)),
          identical(sort(pat), sort(prim$patient)),
          identical(sort(pat), sort(myc$patient)))

if (!all(c(LIMBS_PRIME, LIMBS_NEG) %in% names(prim))) {
  stop("script 08 has not been re-run since the limb columns were added. ",
       "Re-source scripts/08_score_priming.R, then this script.", call. = FALSE)
}

# =============================================================================
# 2. The analysis frame, the frozen scaling, and the plate map
# =============================================================================
message("\n2. analysis frame")

D <- tibble::tibble(patient = pat) %>%
  dplyr::left_join(dplyr::select(cov, patient, purity, leukocyte_fraction,
                                 PAM50, TP53_status, plate, er_call,
                                 PIK3CA_altered, PI3K_pathway_intact, MYC_amp),
                   by = "patient") %>%
  dplyr::left_join(dplyr::select(myc, patient, M_a, M_b), by = "patient") %>%
  dplyr::left_join(dplyr::select(prim, patient, dplyr::all_of(
    c(ENDPOINT_MAIN, NEG_ENDPOINTS, "PRIME_INDEX", "FOXO3_activity",
      LIMBS_PRIME, LIMBS_NEG))), by = "patient") %>%
  dplyr::mutate(
    PROLIF_DISJOINT = as.numeric(mito$gsva_cov["PROLIF_DISJOINT", patient]),
    PROLIF_STD      = as.numeric(mito$gsva_cov["PROLIF_STD", patient]),
    PAM50           = factor(PAM50),
    TP53_status     = factor(TP53_status),
    er_call         = factor(er_call)
  )
stopifnot(identical(D$patient, pat))

# --- the three analysis sets (D8) -------------------------------------------
# NOTE, and it is a refinement of the D8 note's wording rather than a change of
# decision: M2 and M3 KEEP the proliferation term. It is a GSVA score and is
# fully observed for all 1,095, so dropping it would vary the D7 dimension at the
# same time as the D8 one and break the one-at-a-time rule the whole ladder rests
# on. D8's intent - "all patients, with the covariates that are fully observed" -
# is served by keeping it.
SETS <- list(
  M1 = stats::complete.cases(D[, BLOCK_C_VARS]),
  M2 = stats::complete.cases(D[, c("leukocyte_fraction", "plate")]),
  M3 = stats::complete.cases(D[, "plate", drop = FALSE])
)
SET_COVARS <- list(
  M1 = c("purity", "leukocyte_fraction", "PROLIF", "PAM50", "TP53_status", "plate_f"),
  M2 = c("leukocyte_fraction", "PROLIF", "plate_f"),
  M3 = c("PROLIF", "plate_f")
)
for (s in names(SETS)) message(sprintf("   %s: n = %4d", s, sum(SETS[[s]])))

# --- plate pooling, one frozen map (D9) --------------------------------------
PLATE_KEEP <- local({
  t <- table(D$plate[SETS$M1])
  names(t)[t >= PLATE_MIN_N]
})
.pool_frozen <- function(x) factor(ifelse(x %in% PLATE_KEEP, x, "other"))
.pool_fresh  <- function(x) {
  t <- table(x)
  factor(ifelse(x %in% names(t)[t >= PLATE_MIN_N], x, "other"))
}
message("   plate: ", dplyr::n_distinct(D$plate[SETS$M1]), " levels in M1 -> ",
        length(PLATE_KEEP) + 1L, " after pooling at n < ", PLATE_MIN_N,
        " (", sum(!D$plate[SETS$M1] %in% PLATE_KEEP), " patients in `other`)")

# --- frozen z-scaling constants (spec 1.1) -----------------------------------
# Computed ONCE on the primary 938 and reused everywhere, so every coefficient in
# the paper is per SD of the same reference set and a stratum estimate can be
# read against the full-cohort one. Each NULL set is z-scored on its own
# constants instead, because it is a different variable - see section 5.
Z_REF <- SETS$M1
.zfix <- function(v) {
  m <- mean(v[Z_REF]); s <- stats::sd(v[Z_REF])
  if (!is.finite(s) || s == 0) stop("zero variance in a frozen-scaled variable",
                                    call. = FALSE)
  (v - m) / s
}
ZSCALE <- list()
for (nm in c("M_a", "M_b")) {
  ZSCALE[[nm]] <- c(mean = mean(D[[nm]][Z_REF]), sd = stats::sd(D[[nm]][Z_REF]))
  D[[nm]] <- .zfix(D[[nm]])
}
for (ins in INSTRUMENTS) {
  for (a in ARMS_ALL) {
    v <- as.numeric(SCORES[[ins]][a, pat])
    ZSCALE[[paste(ins, a, sep = "|")]] <- c(mean = mean(v[Z_REF]),
                                            sd = stats::sd(v[Z_REF]))
    SCORES[[ins]][a, ] <- .zfix(v)
  }
}
message("   exposures z-scored on the M1 set; ", length(ZSCALE),
        " scaling constants frozen")

# M-c is categorical and is NOT z-scored - it is a 0/1 amplification call.
D$M_c <- as.numeric(D$MYC_amp)

# =============================================================================
# 3. The fitting machinery
# =============================================================================
# Two paths that must agree:
#   .fit()      readable, uses lm(), gives estimate / SE / CI / p. ~110 fits.
#   .fast_b()   design matrix built once, two columns swapped, .lm.fit, returns
#               the interaction coefficient only. 68,000 fits.
# The agreement between them is ASSERTED for both spines in section 4, the same
# way script 07 asserted GSVA against itself inside its own null batch. Two
# implementations of one model diverging silently is the failure mode that would
# invalidate every percentile without producing an error.
message("\n3. fitting machinery")

.drop_constant <- function(df, covars) {
  keep <- vapply(covars, function(v) {
    x <- df[[v]]
    dplyr::n_distinct(x[!is.na(x)]) >= 2L
  }, logical(1))
  covars[keep]
}

.covars_for <- function(dataset, prolif, extra) {
  cvs <- SET_COVARS[[dataset]]
  if (is.null(prolif)) cvs <- setdiff(cvs, "PROLIF")
  c(cvs, extra)
}

# Assembles the modelling frame for one fit. `idx` is a logical over all 1,095.
#
# It builds ONLY the covariates this fit actually uses. That is not tidiness: the
# frame is passed through complete.cases(), so carrying an unused column with
# missing values would silently impose that column's coverage on the fit. Doing
# it the other way collapses M2 and M3 straight back to the 938 they exist to
# escape, and an all-NA PROLIF column under S2 empties the frame entirely.
.frame <- function(endpoint, myc_var, arm, instrument, prolif, idx, covars,
                   repool = FALSE) {
  cv <- D[idx, , drop = FALSE]
  df <- tibble::tibble(
    Y   = D[[endpoint]][idx],
    MYC = D[[myc_var]][idx],
    OX  = as.numeric(SCORES[[instrument]][arm, idx])
  )
  for (v in covars) {
    df[[v]] <- switch(
      v,
      plate_f     = if (repool) .pool_fresh(cv$plate) else .pool_frozen(cv$plate),
      PROLIF      = cv[[prolif]],
      PAM50       = droplevels(cv$PAM50),
      TP53_status = droplevels(cv$TP53_status),
      er_call     = droplevels(cv$er_call),
      cv[[v]])
  }
  df
}

.fit <- function(label,
                 endpoint   = ENDPOINT_MAIN,
                 myc_var    = "M_a",
                 arm        = ARM_PRIMARY,
                 instrument = INSTRUMENTS,
                 prolif     = "PROLIF_DISJOINT",
                 dataset    = "M1",
                 extra      = character(0),
                 subset     = NULL,
                 stratum    = "-",
                 repool     = FALSE,
                 random_plate = FALSE) {

  idx <- SETS[[dataset]]
  if (!is.null(subset)) idx <- idx & subset & !is.na(subset)

  dplyr::bind_rows(lapply(instrument, function(ins) {
    cvs0 <- .covars_for(dataset, prolif, extra)
    df   <- .frame(endpoint, myc_var, arm, ins, prolif, idx, cvs0, repool)
    df   <- df[stats::complete.cases(df), , drop = FALSE]
    cvs  <- .drop_constant(df, cvs0)

    if (random_plate) return(.fit_random_plate(label, df, cvs, endpoint,
                                               myc_var, arm, ins, prolif,
                                               dataset, stratum))
    fo <- stats::as.formula(paste(
      "Y ~ MYC * OX",
      if (length(cvs)) paste("+", paste(cvs, collapse = " + ")) else ""))
    m  <- stats::lm(fo, data = df)
    co <- summary(m)$coefficients
    if (!"MYC:OX" %in% rownames(co)) {
      stop("fit '", label, "' (", ins, "): no MYC:OX row. The interaction was ",
           "aliased away, which means a covariate is collinear with it.",
           call. = FALSE)
    }
    tcrit <- stats::qt(1 - (1 - CI_LEVEL) / 2, df = m$df.residual)
    tibble::tibble(
      label = label, endpoint = endpoint, myc = myc_var, arm = arm,
      instrument = ins,
      spec = if (is.null(prolif)) "S2"
             else if (identical(prolif, "PROLIF_DISJOINT")) "S1" else "S3",
      dataset = dataset, stratum = stratum, plate = "fixed",
      n = nrow(df), df_resid = m$df.residual,
      estimate = co["MYC:OX", 1], se = co["MYC:OX", 2],
      ci_lo = co["MYC:OX", 1] - tcrit * co["MYC:OX", 2],
      ci_hi = co["MYC:OX", 1] + tcrit * co["MYC:OX", 2],
      p = co["MYC:OX", 4]
    )
  }))
}

# Plate as a random intercept (D9). Sensitivity, primary arm only, never the
# battery: the null is 68,000 fits and lmer would make it hours instead of a
# minute. Skipped with a message rather than a stop if lme4 is absent.
.fit_random_plate <- function(label, df, cvs, endpoint, myc_var, arm, ins,
                              prolif, dataset, stratum) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    message("   lme4 not installed - the random-plate sensitivity is SKIPPED. ",
            "Install lme4 and re-source to obtain it.")
    return(NULL)
  }
  cvs <- setdiff(cvs, "plate_f")
  fo  <- stats::as.formula(paste(
    "Y ~ MYC * OX", if (length(cvs)) paste("+", paste(cvs, collapse = " + ")) else "",
    "+ (1 | plate_f)"))
  m  <- lme4::lmer(fo, data = df, REML = TRUE)
  co <- summary(m)$coefficients
  tibble::tibble(
    label = label, endpoint = endpoint, myc = myc_var, arm = arm,
    instrument = ins, spec = "S1", dataset = dataset, stratum = stratum,
    plate = "random", n = nrow(df), df_resid = NA_real_,
    estimate = co["MYC:OX", 1], se = co["MYC:OX", 2],
    ci_lo = co["MYC:OX", 1] - 1.96 * co["MYC:OX", 2],
    ci_hi = co["MYC:OX", 1] + 1.96 * co["MYC:OX", 2],
    p = NA_real_)      # lmerTest not assumed; the CI carries the inference
}

# --- the fast path -----------------------------------------------------------
# Only two columns change between null draws: the OXPHOS main effect and its
# interaction with MYC. They are APPENDED, so the interaction is always the last
# column and the index cannot drift - asserted anyway, because getting it wrong
# would report the wrong coefficient 68,000 times without erroring.
.fast_setup <- function(endpoint, myc_var, prolif, dataset, idx) {
  cvs0 <- .covars_for(dataset, prolif, character(0))
  df   <- .frame(endpoint, myc_var, ARM_PRIMARY, INSTRUMENTS[1], prolif, idx, cvs0)
  ok   <- stats::complete.cases(df)
  cvs  <- .drop_constant(df[ok, ], cvs0)
  fo  <- stats::as.formula(
    paste("~ MYC", if (length(cvs)) paste("+", paste(cvs, collapse = " + ")) else ""))
  X0  <- stats::model.matrix(fo, data = df[ok, , drop = FALSE])
  X   <- cbind(X0, OX = 0, `MYC:OX` = 0)
  stopifnot(colnames(X)[ncol(X)] == "MYC:OX",
            colnames(X)[ncol(X) - 1L] == "OX")
  list(X = X, y = df$Y[ok], myc = df$MYC[ok], ok = ok,
       j_ox = ncol(X) - 1L, j_int = ncol(X))
}

.fast_b <- function(S, ox_z) {
  X <- S$X
  X[, S$j_ox]  <- ox_z
  X[, S$j_int] <- S$myc * ox_z
  f <- .lm.fit(X, S$y)
  if (f$rank < ncol(X)) return(NA_real_)
  f$coefficients[S$j_int]
}

# =============================================================================
# 4. The ladder - each row varies ONE dimension from the spine
# =============================================================================
message("\n4. the ladder")

FITS <- list()
.add  <- function(x) if (!is.null(x) && nrow(x)) FITS[[length(FITS) + 1L]] <<- x

# --- row 0: the spine --------------------------------------------------------
spine <- .fit("spine")
.add(spine)
message("   [ 0] spine (2 fits)")

# --- the agreement assertion -------------------------------------------------
# Two implementations of one model. If they disagree the percentiles are void.
for (ins in INSTRUMENTS) {
  S  <- .fast_setup(ENDPOINT_MAIN, "M_a", "PROLIF_DISJOINT", "M1", SETS$M1)
  b1 <- .fast_b(S, as.numeric(SCORES[[ins]][ARM_PRIMARY, SETS$M1])[S$ok])
  b2 <- spine$estimate[spine$instrument == ins]
  if (!isTRUE(all.equal(b1, b2, tolerance = 1e-10))) {
    stop("the fast path and lm() disagree on spine ", ins, ": ", b1, " vs ", b2,
         ". Every null percentile would be computed from a different model than ",
         "the observed value. Stop here.", call. = FALSE)
  }
}
message("        fast path verified against lm() on both spines")

# --- row 1: D7 proliferation specification, M-a ------------------------------
.add(.fit("S2 no proliferation term", prolif = NULL))
.add(.fit("S3 PROLIF_STD",            prolif = "PROLIF_STD"))

# --- row 2: MYC estimator ----------------------------------------------------
# D7: PROLIF_DISJOINT is defined against the stripped Felsher set and must NEVER
# be used with M-b. For M-b, S1 and S3 are the same fit and are reported once.
.add(.fit("M-b, PROLIF_STD (S1=S3)", myc_var = "M_b", prolif = "PROLIF_STD"))
.add(.fit("M-b, no proliferation",   myc_var = "M_b", prolif = NULL))

# --- row 3: the MYC instrument ----------------------------------------------
.add(.fit("M-c, MYC 8q24 amplification", myc_var = "M_c"))

# --- row 4: missing-data ladder (D8) ----------------------------------------
.add(.fit("M2 leukocyte + prolif + plate", dataset = "M2"))
.add(.fit("M3 prolif + plate",             dataset = "M3"))

# --- row 5: plate as a random intercept (D9) --------------------------------
.add(.fit("plate random intercept", random_plate = TRUE))

# --- row 6: purity-high ------------------------------------------------------
# n = 270. Weak for a three-term interaction; a null here is uninformative
# rather than contradictory, and the spec says so in advance.
.add(.fit("purity > 0.7", subset = D$purity > 0.7, stratum = "purity>0.7",
          repool = TRUE))

# --- row 7: the specificity battery -----------------------------------------
for (a in setdiff(ARMS_ALL, ARM_PRIMARY)) .add(.fit(paste("arm:", a), arm = a))
message("   [ 7] specificity battery: ", length(ARMS_ALL) - 1L, " comparator arms")

# --- row 8: endpoint negatives ----------------------------------------------
for (e in NEG_ENDPOINTS) .add(.fit(paste("endpoint:", e), endpoint = e))

# --- row 9: the whole-family index ------------------------------------------
# NOT a robustness check on PRIME: they correlate at rho = 0.312 and were never
# measuring the same thing. Reported in its own right (spec section 5).
.add(.fit("endpoint: PRIME_INDEX", endpoint = "PRIME_INDEX"))

# --- row 10: limb-wise (spec section 4) --------------------------------------
for (e in c(LIMBS_PRIME, LIMBS_NEG)) .add(.fit(paste("limb:", e), endpoint = e))

# --- row 11: ER handling (spec section 6) ------------------------------------
for (e in c(ENDPOINT_MAIN, NEG_ENDPOINTS)) {
  .add(.fit(paste("ER-adjusted:", e), endpoint = e, extra = "er_call"))
}
for (lv in c("Positive", "Negative")) {
  .add(.fit(paste("within ER", lv), subset = D$er_call == lv,
            stratum = paste0("ER_", lv), repool = TRUE))
}

# --- row 12: strata ----------------------------------------------------------
# PTEN_altered (n = 83) and PAM50 Normal (n = 29) are NOT fitted; PTEN enters H2
# through PI3K_pathway_intact, which is the pre-specified split (spec section 7).
STRATA <- list(
  "TP53_mutant"       = D$TP53_status == "mutant",
  "TP53_wildtype"     = D$TP53_status == "wildtype",
  "PIK3CA_altered"    = D$PIK3CA_altered,
  "PIK3CA_intact"     = !D$PIK3CA_altered,
  "PI3K_intact"       = D$PI3K_pathway_intact,
  "PI3K_altered"      = !D$PI3K_pathway_intact
)
for (lv in c("BRCA_LumA", "BRCA_LumB", "BRCA_Basal", "BRCA_Her2")) {
  STRATA[[lv]] <- D$PAM50 == lv
}
for (nm in names(STRATA)) {
  .add(.fit(paste("stratum:", nm), subset = STRATA[[nm]], stratum = nm,
            repool = TRUE))
}

# --- H2's own test: FOXO3 activity as the endpoint ---------------------------
# FOXO3 activity correlates -0.508 with purity and +0.374 with leukocyte
# fraction - the strongest confound in this pipeline - so the purity-high refit
# is MANDATORY for any FOXO3 claim, not optional (spec section 8).
.add(.fit("H2 FOXO3 endpoint", endpoint = "FOXO3_activity"))
.add(.fit("H2 FOXO3, PI3K-intact", endpoint = "FOXO3_activity",
          subset = D$PI3K_pathway_intact, stratum = "PI3K_intact", repool = TRUE))
.add(.fit("H2 FOXO3, purity > 0.7 (MANDATORY)", endpoint = "FOXO3_activity",
          subset = D$purity > 0.7, stratum = "purity>0.7", repool = TRUE))

coefs <- dplyr::bind_rows(FITS)
message("   ladder complete: ", nrow(coefs), " fits over ",
        dplyr::n_distinct(coefs$label), " labels")

# =============================================================================
# 5. The expression-matched null
# =============================================================================
# Same covariates, same plate map, same 938 patients, same spine specification.
# The ONLY thing that changes between the observed value and a draw is the gene
# set (D9 consequence (a)). Each null set is z-scored on its OWN constants,
# because it is a different variable - using the arm's constants would make the
# null a test of scale rather than of the gene set.
message("\n5. expression-matched null, ", length(ARMS_NULL), " arms x 2 instruments")

S_FAST <- .fast_setup(ENDPOINT_MAIN, "M_a", "PROLIF_DISJOINT", "M1", SETS$M1)
idx_m1 <- which(SETS$M1)[S_FAST$ok]
.z_own <- function(v) (v - mean(v)) / stats::sd(v)

B_NULL <- list()      # [[instrument]][[arm]] = numeric(2000)
null_rows <- list()

for (k in seq_along(ARMS_NULL)) {
  a  <- ARMS_NULL[[k]]
  nf <- readRDS(mito$null_manifest$file[mito$null_manifest$arm == a])
  t0 <- Sys.time()

  for (ins in INSTRUMENTS) {
    M <- nf[[ins]][, pat[idx_m1], drop = FALSE]      # 2000 x n, patient-aligned
    b <- vapply(seq_len(nrow(M)),
                function(i) .fast_b(S_FAST, .z_own(as.numeric(M[i, ]))),
                numeric(1))
    if (anyNA(b)) {
      warning("arm '", a, "' (", ins, "): ", sum(is.na(b)),
              " null fit(s) were rank deficient and are dropped.", call. = FALSE)
    }
    b_obs <- coefs$estimate[coefs$label %in% c("spine", paste("arm:", a)) &
                              coefs$instrument == ins &
                              coefs$arm == a & coefs$endpoint == ENDPOINT_MAIN]
    stopifnot(length(b_obs) == 1L)

    B_NULL[[ins]][[a]] <- b
    bb <- b[!is.na(b)]
    null_rows[[length(null_rows) + 1L]] <- tibble::tibble(
      arm = a, instrument = ins, n_draws = length(bb),
      b_obs = b_obs, null_median = stats::median(bb),
      percentile = 100 * mean(bb < b_obs),
      p_emp = max(2 * min(mean(bb <= b_obs), mean(bb >= b_obs)), 1 / length(bb))
    )
  }
  message(sprintf("   [%2d/%2d] %-26s %.1f s", k, length(ARMS_NULL), a,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
null_tbl <- dplyr::bind_rows(null_rows)

# =============================================================================
# 6. The paired null - the comparative claim
# =============================================================================
# A one-set null answers "is OXPHOS extreme?". The claim is comparative, so the
# null is on the DIFFERENCE with both sets redrawn. Costs no extra fitting: the
# coefficients above already contain everything. Valid because the draws for the
# two arms come from independent seeds, so their difference is a proper null for
# the difference of two independently expression-matched sets.
message("\n6. paired null, ", length(PAIRED_CONTRASTS), " contrasts x 2 instruments")

paired_rows <- list()
for (cc in PAIRED_CONTRASTS) {
  for (ins in INSTRUMENTS) {
    ba <- B_NULL[[ins]][[cc[1]]]; bb <- B_NULL[[ins]][[cc[2]]]
    ok <- !is.na(ba) & !is.na(bb)
    d_null <- ba[ok] - bb[ok]
    d_obs  <- null_tbl$b_obs[null_tbl$arm == cc[1] & null_tbl$instrument == ins] -
              null_tbl$b_obs[null_tbl$arm == cc[2] & null_tbl$instrument == ins]
    paired_rows[[length(paired_rows) + 1L]] <- tibble::tibble(
      arm_a = cc[1], arm_b = cc[2], instrument = ins, n_draws = length(d_null),
      d_obs = d_obs, null_median = stats::median(d_null),
      percentile = 100 * mean(d_null < d_obs),
      p_emp = max(2 * min(mean(d_null <= d_obs), mean(d_null >= d_obs)),
                  1 / length(d_null))
    )
  }
}
paired_tbl <- dplyr::bind_rows(paired_rows)

# =============================================================================
# 7. The decision rule, evaluated
# =============================================================================
# Fixed in spec section 10 before any of the above existed. Evaluated
# mechanically here so the conclusion cannot drift in the writing.
message("\n7. H1 decision rule")

.on_both <- function(v) length(v) == 2L && all(v)

sp <- coefs %>% dplyr::filter(label == "spine")
c1 <- .on_both(sp$estimate > 0 & sp$ci_lo > 0)

nt <- null_tbl %>% dplyr::filter(arm == ARM_PRIMARY)
c2 <- .on_both(nt$p_emp < P_EMP_ALPHA)

pt <- paired_tbl %>% dplyr::filter(arm_a == ARM_PRIMARY, arm_b == ARM_NEGATIVE)
c3 <- .on_both(pt$p_emp < P_EMP_ALPHA)

neg  <- coefs %>% dplyr::filter(grepl("^endpoint:", label),
                                endpoint %in% NEG_ENDPOINTS)
bbc3 <- coefs %>% dplyr::filter(endpoint == "log2_BBC3")
# "CI includes zero" is AND, not OR - the OR form is a tautology and would pass
# clause 4 unconditionally.
c4 <- all(neg$ci_lo <= 0 & neg$ci_hi >= 0) &&
  .on_both(bbc3$estimate > 0 & bbc3$ci_lo > 0)

decision <- tibble::tibble(
  clause = c("1 interaction positive, CI excludes 0, BOTH instruments",
             "2 OXPHOS subunits beats its matched null, BOTH",
             "3 paired contrast vs assembly factors beats its null, BOTH",
             "4 endpoint negatives null, and BBC3 limb carries the effect"),
  passes = c(c1, c2, c3, c4)
)
decision %>% as.data.frame() %>% print(row.names = FALSE)
message("\n   H1 SUPPORTED: ", all(decision$passes),
        "  -- report the full pattern either way (plan section 15)")

message("\n   spine coefficients:")
sp %>%
  dplyr::select(instrument, n, estimate, se, ci_lo, ci_hi, p) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   battery, ordered by GSVA percentile:")
null_tbl %>%
  dplyr::select(arm, instrument, b_obs, percentile, p_emp) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3))) %>%
  dplyr::arrange(instrument, percentile) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   paired null:")
paired_tbl %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 8. Save
# =============================================================================
message("\n8. save")

out <- list(
  coefficients = coefs,
  null         = null_tbl,
  paired_null  = paired_tbl,
  decision     = decision,
  b_null       = B_NULL,
  zscale       = ZSCALE,
  plate_keep   = PLATE_KEEP,
  set_sizes    = vapply(SETS, sum, integer(1)),
  spec = list(
    document  = "docs/2026-08-29_script09_build_spec.md",
    spine     = paste("PRIME ~ M_a * OXPHOS subunits + purity +",
                      "leukocyte_fraction + PROLIF_DISJOINT + PAM50 +",
                      "TP53_status + plate_pooled, n = 938"),
    coprimary = paste("report both instruments; claim only what BOTH support.",
                      "An effect on one alone is instrument-dependent."),
    percentile = "two-sided empirical p, floor 1/n_draws",
    not_null_tested = setdiff(ARMS_ALL, ARMS_NULL),
    multiplicity = paste("primary is one pre-registered coefficient per",
                         "instrument, uncorrected; the battery is calibrated by",
                         "its matched null, not by Bonferroni")
  ),
  built = Sys.time()
)

saveRDS(out, file.path(DIR_RESULTS, "block_c_models.rds"))
readr::write_csv(coefs,      file.path(DIR_TABLES, "block_c_coefficients.csv"))
readr::write_csv(null_tbl,   file.path(DIR_TABLES, "block_c_null.csv"))
readr::write_csv(paired_tbl, file.path(DIR_TABLES, "block_c_paired_null.csv"))
readr::write_csv(decision,   file.path(DIR_TABLES, "block_c_decision.csv"))

message("\n09: done.")
message("    results/block_c_models.rds")
message("    outputs/tables/  4 tables")
message("\n    Panel a is drawn in script 18 from $coefficients and $null.")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  b <- readRDS(file.path(DIR_RESULTS, "block_c_models.rds"))

  # --- the two co-primary spines, side by side -----------------------------
  # The rule is not "either"; it is "both". Look at these together and nowhere
  # else first.
  b$coefficients %>% dplyr::filter(label == "spine") %>% as.data.frame()

  # --- the D7 three-specification panel, which is part of the result -------
  b$coefficients %>%
    dplyr::filter(endpoint == "PRIME", arm == "OXPHOS subunits",
                  myc %in% c("M_a", "M_b"), stratum == "-") %>%
    dplyr::select(label, myc, instrument, spec, n, estimate, ci_lo, ci_hi, p) %>%
    as.data.frame() %>% print(row.names = FALSE)

  # --- the specificity forest, which is Panel a ----------------------------
  fo <- b$null %>% dplyr::filter(instrument == "gsva") %>%
    dplyr::arrange(b_obs)
  with(fo, {
    plot(b_obs, seq_along(b_obs), pch = 16, yaxt = "n", ylab = "",
         xlab = "MYC:OXPHOS coefficient (per SD)", xlim = range(b_obs) * 1.3)
    axis(2, seq_along(arm), arm, las = 1, cex.axis = 0.6)
    abline(v = 0, lty = 2)
  })

  # --- is the claim arm actually distinguishable from its primary negative?
  # 07 found these two correlate at rho = 0.875 as scores. If both beat the null
  # the specificity claim is NOT supported, however good each looks alone.
  b$null %>%
    dplyr::filter(arm %in% c("OXPHOS subunits", "OXPHOS assembly factors")) %>%
    as.data.frame() %>% print(row.names = FALSE)
  b$paired_null %>% as.data.frame() %>% print(row.names = FALSE)

  # the paired null drawn, for the primary contrast
  d <- b$b_null$gsva[["OXPHOS subunits"]] -
       b$b_null$gsva[["OXPHOS assembly factors"]]
  hist(d, breaks = 60, main = "paired null: subunits - assembly factors",
       xlab = "difference in MYC:OXPHOS coefficient")
  abline(v = b$paired_null$d_obs[b$paired_null$instrument == "gsva" &
                                 b$paired_null$arm_b == "OXPHOS assembly factors"],
         col = "red", lwd = 2)

  # --- the five per-complex arms: one compartment seen five ways -----------
  # They correlate 0.81-0.97 as scores, so near-identical coefficients are
  # EXPECTED and are not five confirmations. Check they look like one signal.
  b$null %>% dplyr::filter(grepl("subunits$", arm), arm != "OXPHOS subunits") %>%
    as.data.frame() %>% print(row.names = FALSE)

  # --- limb-wise: is PRIME's effect the trigger or the guardian? -----------
  b$coefficients %>% dplyr::filter(grepl("^limb:", label)) %>%
    dplyr::select(endpoint, instrument, estimate, ci_lo, ci_hi, p) %>%
    as.data.frame() %>% print(row.names = FALSE)

  # --- strata, with their n visible ----------------------------------------
  b$coefficients %>% dplyr::filter(grepl("^stratum:", label)) %>%
    dplyr::select(stratum, instrument, n, estimate, ci_lo, ci_hi, p) %>%
    as.data.frame() %>% print(row.names = FALSE)

  # --- H2's FOXO3 test, and the purity refit that is not optional ----------
  b$coefficients %>% dplyr::filter(endpoint == "FOXO3_activity") %>%
    dplyr::select(label, instrument, n, estimate, ci_lo, ci_hi, p) %>%
    as.data.frame() %>% print(row.names = FALSE)

}
