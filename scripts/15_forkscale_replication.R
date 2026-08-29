# 15_forkscale_replication.R
# =============================================================================
# BLOCK D - forkscale replication (ED2, demoted), and the F3-pre redundancy
# check that gates Block F.
#
# Built to:
#   docs/2026-08-29_G3_result_forkscale_availability.md  - G3, and F3-pre's
#     read rules, both fixed BEFORE this script was sourced (note section 4b)
#   docs/2026-08-27_human_validation_plan.md sections 9, 10 (Block D), 12
#
# =============================================================================
# THE TWO THINGS THIS SCRIPT DOES, AND WHY THEY SHARE A FILE
# =============================================================================
# Both need the same object: the companion paper's per-sample bicluster values,
# joined to this arm's analysis set. Building that join twice would be the only
# alternative, and a forkscale assembled two ways is a forkscale that will
# eventually differ two ways.
#
#   F3-pre  How much of this arm's axis IS MB1 forkscale? An EXPOSURE-side
#           diagnostic. No outcome variable appears anywhere in it. It does not
#           pass or fail F3, which is a survival test and belongs in METABRIC
#           (plan section 3 forbids TCGA survival; data/tcga_cdr/README.md
#           quantifies why - 145 PFI events at ~2.3 years).
#
#   Block D The pre-registered three-way `MYC * OXPHOS * forkscale` on PRIME,
#           MB2 primary and MB3 as the ER-neutral control axis (plan section 9).
#
# =============================================================================
# WHAT BLOCK D NOW IS, GIVEN BLOCK C
# =============================================================================
# Block D was specified when H1 was expected to hold. It does not. The Block C
# `MYC:OXPHOS` term on PRIME is 0.020 (p 0.455) on GSVA and 0.010 (p 0.726) on
# mitoPPS. So the three-way below is a MODIFIER OF A NULL MAIN EFFECT.
#
# That is still a legitimate pre-registered stratified test - it asks whether
# the coupling exists in some fork stratum although it does not exist overall,
# exactly as script 09's TP53 and PIK3CA strata do. It is NOT a second chance at
# H1, and no output here may be written as if a positive three-way rescued it.
#
# NO DIRECTION IS PRE-SPECIFIED FOR BLOCK D. Plan section 10 fixes the model and
# not the sign, and inventing a sign now, after Block C, would be exactly the
# post-hoc move the pre-registration exists to prevent. Two-sided, reported with
# CIs. The declared readings in section 6 are about SPECIFICITY and the CONTROL
# AXIS, which the plan does fix, not about the sign of the primary term.
#
# POWER: a three-way with a continuous modifier at n = 938. Underpowered by
# construction. A null here is UNINFORMATIVE, not evidence against.
#
# SCALE DISCIPLINE: no scoring happens here. Every score is read as built by
# scripts 06-08. forkscale is a PC1 loading, not an expression value, and is
# never logged, z-scored against expression, or pooled across cohorts.
#
# SPECIES: human. TCGA-BRCA, and the companion paper's human TCGA biclusters.
# See CLAUDE.md.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

message("\n15: Block D - forkscale, and the F3-pre redundancy check\n",
        strrep("=", 78))

# -----------------------------------------------------------------------------
# 0. Constants
# -----------------------------------------------------------------------------
DIR_FORK <- file.path(DIR_DATA, "menegollo_biclusters")

ARM_PRIMARY  <- "OXPHOS subunits"
ARM_NEGATIVE <- "OXPHOS assembly factors"
INSTRUMENTS  <- c("gsva", "mitopps")
CI_LEVEL     <- 0.95

# F3-pre read rules, fixed in the G3 note section 4b before this file was run.
F3PRE_REDUNDANT   <- 0.70
F3PRE_INDEPENDENT <- 0.30
F3PRE_DECISION_VAR <- "MB1_forkscale"   # the plain form; complete, no Inf

# MB3's sign convention. Upstream is INCONSISTENT between cohorts: the METABRIC
# compile script negates MB3.pc1 before forming forkscale and the TCGA script
# does not (data/menegollo_biclusters/README.md). We keep the TCGA convention,
# as stored, because this script is TCGA-only. Anything that later puts TCGA and
# METABRIC MB3 side by side MUST flip one of them; `MB3_forkscale_mbconv` below
# is the METABRIC-convention copy, carried so nobody has to rediscover this.
MB3_TCGA_SIGN <- +1

# Expected sizes, checked at load (script 02's idiom), so a truncated or
# replaced snapshot is caught here rather than as a shifted coefficient. The
# upstream git blob SHAs each file was verified against at snapshot time are in
# data/menegollo_biclusters/README.md; they are NOT re-checked at runtime, and
# the section 2 reproduction test against the published frame is the real
# content guard.
FORK_BYTES <- c(
  "TCGA_MB1_RNASeq_data_nonorm.RData" = 242449,
  "TCGA_MB2_RNASeq_data_nonorm.RData" = 241801,
  "TCGA_MB3_RNASeq_data_nonorm.RData" = 241236,
  "TCGA.all.biclusters.RNAseq.Rdata"  = 191563
)

# =============================================================================
# 1. Inputs - consumed as built
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

for (f in names(FORK_BYTES)) {
  p <- file.path(DIR_FORK, f)
  if (!file.exists(p)) {
    stop("missing forkscale input: ", p,
         "\nSee data/menegollo_biclusters/README.md - it carries the exact\n",
         "`gh api ... git/blobs/<sha>` command to re-fetch each file.",
         call. = FALSE)
  }
  if (!identical(as.numeric(file.size(p)), as.numeric(FORK_BYTES[[f]]))) {
    stop(f, " is ", file.size(p), " bytes but should be ", FORK_BYTES[[f]],
         ". The snapshot has been replaced or truncated. Re-fetch it from the ",
         "pinned commit; see data/menegollo_biclusters/README.md.",
         call. = FALSE)
  }
}
message("   4 forkscale files present and the expected size")

# =============================================================================
# 2. Assemble forkscale, and verify it against the published frame
# =============================================================================
# The PER-BICLUSTER files are the source (1,037 patients). The assembled file is
# 849 - an upstream inner_join to two files that are not in the repository - and
# is used here ONLY to check that our arithmetic reproduces the published
# numbers. Reversing those roles would silently discard 188 patients.
message("\n2. forkscale")

.load_one <- function(mb) {
  e <- new.env()
  load(file.path(DIR_FORK, sprintf("TCGA_%s_RNASeq_data_nonorm.RData", mb)),
       envir = e)
  d <- e[[paste0("TCGA.", mb, ".RNAseq.df")]]
  stopifnot(is.data.frame(d), nrow(d) == 1037L)
  pc1 <- as.numeric(d[[paste0(mb, ".pc1")]])
  idxv <- as.numeric(d[[paste0(mb, ".index")]])
  stopifnot(setequal(idxv, seq_len(1037L)))      # a permutation, as documented
  out <- tibble::tibble(
    patient = substr(as.character(d$X), 1, 12),
    aliquot = as.character(d$X),
    pc1 = pc1, index = idxv,
    forkscale = pc1 / idxv, forkscale_log = pc1 / log(idxv))
  names(out)[3:6] <- paste0(mb, "_", names(out)[3:6])
  out
}

FK <- Reduce(function(a, b) dplyr::full_join(a, b, by = c("patient", "aliquot")),
             lapply(c("MB1", "MB2", "MB3"), .load_one))
stopifnot(nrow(FK) == 1037L, !anyDuplicated(FK$patient))

# --- the Inf guard -----------------------------------------------------------
# log(1) == 0, so the sample ranked FIRST has an infinite `_log` forkscale. The
# published assembled frame carries exactly this. Inf is not NA, so
# complete.cases() will NOT remove it: it propagates into a fit and destroys it
# silently. Converted to NA here, loudly.
for (v in grep("_forkscale_log$", names(FK), value = TRUE)) {
  bad <- !is.finite(FK[[v]])
  if (any(bad)) {
    message("   ", v, ": ", sum(bad), " non-finite (index == 1) -> NA  [",
            paste(FK$patient[bad], collapse = ", "), "]")
    FK[[v]][bad] <- NA_real_
  }
}
stopifnot(all(is.finite(FK$MB1_forkscale)),
          all(is.finite(FK$MB2_forkscale)),
          all(is.finite(FK$MB3_forkscale)))

# --- MB3, both sign conventions, both named ----------------------------------
FK$MB3_forkscale_mbconv     <- -FK$MB3_forkscale
FK$MB3_forkscale_log_mbconv <- -FK$MB3_forkscale_log
FK$MB3_forkscale     <- MB3_TCGA_SIGN * FK$MB3_forkscale
FK$MB3_forkscale_log <- MB3_TCGA_SIGN * FK$MB3_forkscale_log

# --- verification against the published assembled frame ----------------------
ea <- new.env()
load(file.path(DIR_FORK, "TCGA.all.biclusters.RNAseq.Rdata"), envir = ea)
UP <- ea$TCGA.all.biclusters.RNAseq.df2
stopifnot(is.data.frame(UP), nrow(UP) == 849L)

vchk <- dplyr::inner_join(
  dplyr::select(FK, aliquot, MB1_forkscale, MB2_forkscale),
  tibble::tibble(aliquot = as.character(UP$X),
                 up_MB1 = as.numeric(UP$MB1.forkscale),
                 up_MB2 = as.numeric(UP$MB2.forkscale)),
  by = "aliquot")
d1 <- max(abs(vchk$MB1_forkscale - vchk$up_MB1))
d2 <- max(abs(vchk$MB2_forkscale - vchk$up_MB2))
if (nrow(vchk) != 849L || d1 > 1e-12 || d2 > 1e-12) {
  stop("recomputed forkscale does not reproduce the published assembled frame ",
       "(n = ", nrow(vchk), ", max abs diff MB1 ", signif(d1, 3), ", MB2 ",
       signif(d2, 3), "). Either the snapshot changed or the definition did.",
       call. = FALSE)
}
message("   definition verified against the published frame: n = ", nrow(vchk),
        ", max abs diff ", signif(max(d1, d2), 3))
# MB3 is deliberately NOT in that check: the published TCGA frame stores it in
# the TCGA convention and we would only be checking our own sign constant.

# --- coverage ----------------------------------------------------------------
in_cohort <- sum(FK$patient %in% pat)
message("   coverage: ", in_cohort, " of 1037 forkscale patients are in the ",
        length(pat), "-patient cohort")

# =============================================================================
# 3. The analysis frame - Block C's, verified to be Block C's
# =============================================================================
# Same construction as script 10: the z-scaling constants and the plate map are
# READ from Block C rather than recomputed, and the frame is then asserted
# against Block C's own spine coefficient. Without that assertion every number
# below would be on a quietly different footing from the ones it is compared to.
message("\n3. analysis frame")

D <- tibble::tibble(patient = pat) %>%
  dplyr::left_join(dplyr::select(cov, patient, purity, leukocyte_fraction,
                                 PAM50, TP53_status, plate, er_call, BUFFER,
                                 aneuploidy_cbio, complete_block_c),
                   by = "patient") %>%
  dplyr::left_join(dplyr::select(myc, patient, M_a, M_b), by = "patient") %>%
  dplyr::left_join(dplyr::select(prim, patient, PRIME, log2_BCL2L1,
                                 log2_BCL2L11, log2_BBC3), by = "patient") %>%
  dplyr::left_join(dplyr::select(FK, patient, dplyr::ends_with("forkscale"),
                                 dplyr::ends_with("forkscale_log")),
                   by = "patient") %>%
  dplyr::mutate(
    PROLIF_DISJOINT = as.numeric(mito$gsva_cov["PROLIF_DISJOINT", patient]),
    PAM50           = factor(PAM50),
    TP53_status     = factor(TP53_status),
    er_call         = factor(er_call))

Z <- bc$zscale
.apply_z <- function(v, key) {
  k <- Z[[key]]; stopifnot(!is.null(k)); (v - k[["mean"]]) / k[["sd"]]
}
D$MYC <- .apply_z(D$M_a, "M_a")

SC <- list()
for (ins in INSTRUMENTS) {
  m <- mito[[paste0(ins, "_arms")]]
  SC[[ins]] <- t(vapply(rownames(m), function(a)
    .apply_z(as.numeric(m[a, pat]), paste(ins, a, sep = "|")),
    numeric(length(pat))))
  colnames(SC[[ins]]) <- pat
}
ARMS_ALL <- rownames(SC$gsva)

D$plate_f <- factor(ifelse(D$plate %in% bc$plate_keep, D$plate, "other"))

# Block C's 938, taken from the flag script 03 wrote rather than re-derived.
idx_c <- which(D$complete_block_c)
stopifnot(length(idx_c) == 938L)

# The forkscale set: Block C's 938 that also have forkscale.
idx <- which(D$complete_block_c & !is.na(D$MB1_forkscale))
message("   Block C analysis set: ", length(idx_c),
        "   |   with forkscale: ", length(idx),
        "  (", length(idx_c) - length(idx), " lost)")

# --- the assertion -----------------------------------------------------------
d_chk <- D[idx_c, ]
d_chk$OX <- SC$gsva[ARM_PRIMARY, idx_c]
.chk  <- stats::lm(PRIME ~ MYC * OX + purity + leukocyte_fraction +
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
            "TP53_status", "plate_f")

# `n` is taken from the FIT, never from the caller's frame. `er_call` has 46 NAs
# inside Block C's 938 (D8 is complete-case), so any Block D model that adjusts
# for ER silently fits on fewer rows than the frame it was handed. Reporting
# nrow(dd) would overstate every Block D n by 46.
.tidy <- function(m, term, label, instrument, extra = list()) {
  co <- summary(m)$coefficients
  if (!term %in% rownames(co)) {
    stop("term '", term, "' absent from fit '", label, "'", call. = FALSE)
  }
  crit <- stats::qt(1 - (1 - CI_LEVEL) / 2, m$df.residual)
  out <- tibble::tibble(
    label = label, term = term, instrument = instrument, n = stats::nobs(m),
    estimate = co[term, 1], se = co[term, 2],
    ci_lo = co[term, 1] - crit * co[term, 2],
    ci_hi = co[term, 1] + crit * co[term, 2],
    p = co[term, 4])
  if (length(extra)) out <- dplyr::bind_cols(out, tibble::as_tibble(extra))
  out
}

# =============================================================================
# 4. F3-pre - how much of this arm's axis IS MB1 forkscale?
# =============================================================================
# EXPOSURE SIDE ONLY. There is no outcome variable in this section. Read rules
# were fixed in docs/2026-08-29_G3_result_forkscale_availability.md section 4b
# before this file was sourced; they are restated as constants above and
# evaluated mechanically in section 6 so the reading cannot drift.
#
# Spearman is the primary statistic because forkscale is severely skewed by
# construction (pc1 / index, with index near 1 giving huge values). Pearson is
# reported as a companion, never as the decision.
message("\n4. F3-pre - the redundancy diagnostic (exposure side, no outcome)")

d4 <- D[idx, ]
AXIS <- list(
  `OXPHOS subunits (GSVA)`    = SC$gsva[ARM_PRIMARY, idx],
  `OXPHOS subunits (mitoPPS)` = SC$mitopps[ARM_PRIMARY, idx],
  `OXPHOS assembly (GSVA)`    = SC$gsva[ARM_NEGATIVE, idx],
  `OXPHOS assembly (mitoPPS)` = SC$mitopps[ARM_NEGATIVE, idx],
  `MYC (M_a)`                 = d4$MYC,
  `MYC (M_b)`                 = d4$M_b,
  `PRIME`                     = d4$PRIME,
  `log2 BCL2L1`               = d4$log2_BCL2L1,
  `log2 BCL2L11`              = d4$log2_BCL2L11,
  `PROLIF_DISJOINT`           = d4$PROLIF_DISJOINT,
  `aneuploidy`                = d4$aneuploidy_cbio)

FORKVARS <- c("MB1_forkscale", "MB1_forkscale_log",
              "MB2_forkscale", "MB2_forkscale_log",
              "MB3_forkscale", "MB3_forkscale_log")

f3pre_cor <- dplyr::bind_rows(lapply(FORKVARS, function(fv) {
  x <- d4[[fv]]
  dplyr::bind_rows(lapply(names(AXIS), function(a) {
    y  <- AXIS[[a]]
    ok <- is.finite(x) & is.finite(y)
    ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman"))
    tibble::tibble(
      fork = fv, axis = a, n = sum(ok),
      rho_spearman = unname(ct$estimate),
      p_spearman   = ct$p.value,
      r_pearson    = stats::cor(x[ok], y[ok]))
  }))
}))

# --- incremental variance, both directions -----------------------------------
# A correlation is symmetric; these are not. "How much of OXPHOS is forkscale"
# and "how much of forkscale is OXPHOS" answer different questions and the
# redundancy worry is about the first.
.dr2 <- function(y, add, dat) {
  base <- stats::lm(stats::reformulate(COVARS, response = "Y"),
                    data = dplyr::mutate(dat, Y = y))
  full <- stats::lm(stats::reformulate(c(COVARS, add), response = "Y"),
                    data = dplyr::mutate(dat, Y = y))
  a <- stats::anova(base, full)
  c(r2_base = summary(base)$r.squared, r2_full = summary(full)$r.squared,
    delta_r2 = summary(full)$r.squared - summary(base)$r.squared,
    p_lrt = a[["Pr(>F)"]][2])
}

f3pre_r2 <- dplyr::bind_rows(lapply(INSTRUMENTS, function(ins) {
  dplyr::bind_rows(lapply(c("MB1_forkscale", "MB2_forkscale", "MB3_forkscale"),
                          function(fv) {
    v <- .dr2(SC[[ins]][ARM_PRIMARY, idx], fv, d4)
    tibble::tibble(instrument = ins, fork = fv, direction = "OXPHOS ~ . + fork",
                   r2_base = v[["r2_base"]], r2_full = v[["r2_full"]],
                   delta_r2 = v[["delta_r2"]], p = v[["p_lrt"]])
  }))
}))

f3pre_r2_rev <- dplyr::bind_rows(lapply(INSTRUMENTS, function(ins) {
  dd <- dplyr::mutate(d4, OXv = SC[[ins]][ARM_PRIMARY, idx])
  v  <- .dr2(dd$MB1_forkscale, "OXv", dd)
  tibble::tibble(instrument = ins, fork = "MB1_forkscale",
                 direction = "forkscale ~ . + OXPHOS",
                 r2_base = v[["r2_base"]], r2_full = v[["r2_full"]],
                 delta_r2 = v[["delta_r2"]], p = v[["p_lrt"]])
}))

# --- STATE, if script 11 has run ---------------------------------------------
# STATE is FROZEN IN SCRIPT 11 and defined there and nowhere else. It is NOT
# re-derived here: two definitions of a pre-registered variable is exactly how a
# pre-registered variable stops being one.
PATH_STATE <- file.path(DIR_RESULTS, "tcga_brca_state.rds")
f3pre_state <- NULL
if (file.exists(PATH_STATE)) {
  # Written before script 11 exists, so this reads defensively rather than
  # assuming a shape. A wrong guess about 11's object must not take the whole
  # script down when the F3-pre verdict does not depend on STATE at all.
  st <- readRDS(PATH_STATE)
  sv <- if (is.data.frame(st)) st else st[["state"]]
  if (is.data.frame(sv) && all(c("patient", "STATE") %in% names(sv))) {
    f3pre_state <- tibble::tibble(
        STATE = sv$STATE[match(d4$patient, sv$patient)],
        MB1_forkscale = d4$MB1_forkscale) %>%
      dplyr::group_by(STATE) %>%
      dplyr::summarise(n = dplyr::n(),
                       median_MB1_forkscale = stats::median(MB1_forkscale,
                                                            na.rm = TRUE),
                       .groups = "drop")
    message("   STATE read from script 11 and cross-tabulated")
  } else {
    message("   ", basename(PATH_STATE), " exists but has no patient/STATE ",
            "columns - the STATE panel is skipped. Point this at script 11's ",
            "actual object.")
  }
} else {
  message("   STATE not available (script 11 has not run). The forkscale-vs-",
          "STATE\n   panel of F3-pre is DEFERRED, not skipped - re-source this ",
          "script after 11.")
}

# =============================================================================
# 5. Block D - the pre-registered three-way
# =============================================================================
# PRIME ~ MYC * OXPHOS * MB2_forkscale + covariates      PRIMARY
# PRIME ~ MYC * OXPHOS * MB3_forkscale + covariates      ER-NEUTRAL CONTROL AXIS
#
# ER mitigations, plan section 9, all three applied:
#   1. the endpoint is BBC3/BCL2L1, not BCL2 - already true of PRIME by
#      construction (script 08);
#   2. adjust for ER status, AND fit within ER status - both below;
#   3. MB3 as the control axis - the second model.
#
# On double adjustment: PAM50 is already in Block C's covariate set and largely
# encodes ER. `er_call` is added anyway because the plan says so, with a rank
# assertion; if the two were collinear the fit would fail loudly rather than
# drop a term quietly.
message("\n5. Block D - MYC x OXPHOS x forkscale on PRIME")

FORK_D <- c(MB2 = "MB2_forkscale", MB3 = "MB3_forkscale")
ARMS_D <- c(ARM_PRIMARY, ARM_NEGATIVE)

# Restrict to complete cases and droplevels() BEFORE building the design.
# Both matter and both are silent otherwise:
#   - `er_call` has 46 NAs inside Block C's 938, so every Block D fit is on 892;
#   - a factor level of `plate_f` or `PAM50` with zero rows AFTER that
#     restriction still generates an all-zero column in model.matrix(), which is
#     rank deficiency that has nothing to do with collinearity. Without
#     droplevels() the guard below would fire on a phantom.
.model_frame <- function(dd, vars) {
  ok <- stats::complete.cases(dd[, vars, drop = FALSE])
  droplevels(dd[ok, , drop = FALSE])
}

.fit_d <- function(dat, ins, arm, fv, ystr, label, rows) {
  dd <- dat
  dd$OX <- SC[[ins]][arm, rows]
  dd$FS <- dd[[fv]]
  dd$Y  <- dd[[ystr]]
  dd <- .model_frame(dd, c("Y", "MYC", "OX", "FS", "er_call", COVARS))
  f <- stats::as.formula(paste("Y ~ MYC * OX * FS + er_call +",
                               paste(COVARS, collapse = " + ")))
  X  <- stats::model.matrix(f, data = dd)
  qx <- qr(X)
  if (qx$rank < ncol(X)) {
    aliased <- colnames(X)[qx$pivot[(qx$rank + 1L):ncol(X)]]
    stop("Block D design is rank deficient for ", label, " / ", ins, " / ", arm,
         " (rank ", qx$rank, " of ", ncol(X), "). Aliased column(s): ",
         paste(aliased, collapse = ", "),
         ". Decide explicitly what to drop; do not let a term vanish into an ",
         "NA coefficient.", call. = FALSE)
  }
  m <- stats::lm(f, data = dd)
  .tidy(m, "MYC:OX:FS", label, ins,
        extra = list(arm = arm, fork = fv, endpoint = ystr))
}

RES <- list(); .add <- function(x) RES[[length(RES) + 1L]] <<- x

for (mb in names(FORK_D)) {
  for (ins in INSTRUMENTS) {
    for (arm in ARMS_D) {
      .add(.fit_d(D[idx, ], ins, arm, FORK_D[[mb]], "PRIME",
                  paste0("D ", mb, " three-way"), idx))
    }
  }
}

# --- within-ER strata, primary arm only --------------------------------------
# Mitigation 2's second half. PAM50 and plate are re-levelled inside each
# stratum by lm(); constant factors would break the fit, so strata with a
# degenerate covariate are reported as skipped rather than silently dropped.
for (mb in names(FORK_D)) {
  for (er in levels(stats::na.omit(D$er_call[idx]))) {
    rows <- idx[which(!is.na(D$er_call[idx]) & D$er_call[idx] == er)]
    if (length(rows) < 100L) {
      message("   ER stratum '", er, "' n = ", length(rows), " - skipped (< 100)")
      next
    }
    for (ins in INSTRUMENTS) {
      dd <- D[rows, ]
      dd$OX <- SC[[ins]][ARM_PRIMARY, rows]
      dd$FS <- dd[[FORK_D[[mb]]]]
      dd <- .model_frame(dd, c("PRIME", "MYC", "OX", "FS", COVARS))
      # `er_call` is deliberately NOT a covariate inside an ER stratum - it is
      # constant there by construction. Any other covariate that has become
      # constant is dropped too, and named, rather than silently aliased.
      keep <- COVARS[vapply(COVARS, function(v)
        length(unique(dd[[v]])) > 1L, logical(1))]
      if (length(keep) < length(COVARS)) {
        message("      ER=", er, " / ", ins, ": constant covariate(s) dropped - ",
                paste(setdiff(COVARS, keep), collapse = ", "))
      }
      m <- stats::lm(stats::as.formula(paste("PRIME ~ MYC * OX * FS +",
                                             paste(keep, collapse = " + "))), dd)
      .add(.tidy(m, "MYC:OX:FS", paste0("D ", mb, " within ER=", er), ins,
                 extra = list(arm = ARM_PRIMARY, fork = FORK_D[[mb]],
                              endpoint = "PRIME")))
    }
  }
}

# --- DESCRIPTIVE ONLY: the endpoints that survived Block C -------------------
# Block C's one surviving finding is `MYC x OXPHOS` -> lower BCL-XL, higher BIM.
# Asking whether that is fork-dependent is a modifier analysis on an EXISTING
# finding, not a new hypothesis - but it is NOT in plan section 10, so it is
# marked descriptive, carries NO prediction, and is EXCLUDED from every pass/fail
# in section 6. If it is ever promoted, it needs its own declaration first.
for (mb in names(FORK_D)) {
  for (ins in INSTRUMENTS) {
    for (y in c("log2_BCL2L1", "log2_BCL2L11")) {
      .add(.fit_d(D[idx, ], ins, ARM_PRIMARY, FORK_D[[mb]], y,
                  paste0("D ", mb, " DESCRIPTIVE ", y), idx))
    }
  }
}

coefs <- dplyr::bind_rows(RES)

# --- the matched null, GATED -------------------------------------------------
# The specificity requirement (plan section 2) applies to every positive OXPHOS
# result. Unlike Block B's E2, the design matrix here CHANGES with each draw -
# OX enters four columns - so this cannot be done as a QR re-solve and 17 arms x
# 2,000 draws would be 136,000 fits. It is therefore GATED on the observed
# result: run only if the primary three-way is significant on both instruments,
# and then only for the primary arm and MB2.
#
# THIS IS NOT AN OPTIONAL STEP. If the gate opens and this does not run, the
# result is not reportable.
d_prim <- coefs %>% dplyr::filter(label == "D MB2 three-way", arm == ARM_PRIMARY)
gate_open <- nrow(d_prim) == 2L && all(d_prim$p < 0.05) &&
  length(unique(sign(d_prim$estimate))) == 1L

d_null <- NULL
if (gate_open) {
  message("\n   the primary three-way is significant on both instruments - ",
          "running the matched null")
  # model.matrix() DROPS incomplete rows. The observed Block D fits lose the 46
  # patients with no ER call, so the null must be computed on the same rows or
  # cbind() below would silently recycle a 938-length column into an 892-row
  # design. Restrict first, then build.
  dn <- D[idx, ]
  dn$FS <- dn[[FORK_D[["MB2"]]]]
  dn$.row <- seq_len(nrow(dn))
  dn <- .model_frame(dn, c("PRIME", "MYC", "FS", "er_call", COVARS))
  ox_cols <- pat[idx][dn$.row]          # the null draws, on exactly these rows
  X0 <- stats::model.matrix(
    stats::as.formula(paste("~ MYC * FS + er_call +",
                            paste(COVARS, collapse = " + "))), data = dn)
  stopifnot(nrow(X0) == nrow(dn), length(ox_cols) == nrow(dn))
  yv <- dn$PRIME
  nf_path <- mito$null_manifest$file[mito$null_manifest$arm == ARM_PRIMARY]
  nf <- readRDS(nf_path)
  d_null <- dplyr::bind_rows(lapply(INSTRUMENTS, function(ins) {
    Yd <- nf[[ins]][, ox_cols, drop = FALSE]            # draws x samples
    b  <- vapply(seq_len(nrow(Yd)), function(k) {
      ox <- as.numeric(Yd[k, ])
      ox <- (ox - mean(ox)) / stats::sd(ox)
      X  <- cbind(X0, OX = ox, `MYC:OX` = dn$MYC * ox, `OX:FS` = ox * dn$FS,
                  `MYC:OX:FS` = dn$MYC * ox * dn$FS)
      fit <- .lm.fit(X, yv)
      unname(fit$coefficients[ncol(X)])
    }, numeric(1))
    b <- b[is.finite(b)]
    o <- d_prim$estimate[d_prim$instrument == ins]
    tibble::tibble(arm = ARM_PRIMARY, instrument = ins, n_draws = length(b),
                   b_obs = o, null_median = stats::median(b),
                   percentile = 100 * mean(b < o),
                   p_emp = max(2 * min(mean(b <= o), mean(b >= o)), 1 / length(b)))
  }))
} else {
  message("\n   primary three-way not significant on both instruments - ",
          "matched null not run (gate closed)")
}

# =============================================================================
# 5b. Block D2 - NOT BUILT, and why
# =============================================================================
# Plan section 11 assigns "Block D + D2" to this script. D2 is
#
#     PRIME ~ MYC * OXPHOS * (LP_score - mL_score) + covariates
#
# and needs the Pommier developmental gene sets (LP, mL), which this repo does
# NOT hold. They exist upstream as
# `Input and output data/METABRIC analysis/Pommier_dev_genesets.xlsx` in
# gszabadkai/Menegollo_Bentham.
#
# It is not implemented here because importing a GENE SET is a decision this
# repo takes deliberately and documents, not a side effect of writing a model
# (CLAUDE.md, "Gene sets - consume the snapshot, do not rebuild"). It needs a
# snapshot with a provenance README, a check that the sets are natively human,
# and an overlap audit against the MitoCarta arms in the G1 style - the LP/mL
# sets are developmental and could easily share genes with PROLIF.
#
# D2 is DEFERRED, not dropped. Do not add it to this script without that work.
message("\n5b. Block D2 is deferred - needs the Pommier sets snapshotted first")

# =============================================================================
# 6. The declared readings, evaluated
# =============================================================================
message("\n6. declared readings")

.both <- function(v) length(v) == 2L && all(v)

# --- F3-pre ------------------------------------------------------------------
rho_dec <- f3pre_cor %>%
  dplyr::filter(fork == F3PRE_DECISION_VAR,
                axis %in% c("OXPHOS subunits (GSVA)", "OXPHOS subunits (mitoPPS)"))
stopifnot(nrow(rho_dec) == 2L)
ar <- abs(rho_dec$rho_spearman)
f3pre_verdict <- if (all(ar >= F3PRE_REDUNDANT)) "REDUNDANT" else
  if (all(ar <= F3PRE_INDEPENDENT)) "INDEPENDENT" else "INTERMEDIATE"

message(sprintf("   F3-pre: |rho| = %.3f (GSVA), %.3f (mitoPPS)  ->  %s",
                ar[rho_dec$axis == "OXPHOS subunits (GSVA)"],
                ar[rho_dec$axis == "OXPHOS subunits (mitoPPS)"], f3pre_verdict))
message("   thresholds: REDUNDANT >= ", F3PRE_REDUNDANT,
        " on both, INDEPENDENT <= ", F3PRE_INDEPENDENT, " on both")

# --- Block D -----------------------------------------------------------------
# No direction is pre-specified for the primary term, so there is no pass/fail
# on its sign. What the plan DOES fix is specificity and the control axis.
dn_prim <- coefs %>% dplyr::filter(label == "D MB2 three-way", arm == ARM_PRIMARY)
dn_neg  <- coefs %>% dplyr::filter(label == "D MB2 three-way", arm == ARM_NEGATIVE)
dn_ctrl <- coefs %>% dplyr::filter(label == "D MB3 three-way", arm == ARM_PRIMARY)

d_a <- .both(dn_prim$p < 0.05) && length(unique(sign(dn_prim$estimate))) == 1L
d_b <- !.both(dn_neg$p < 0.05)
d_c <- !.both(dn_ctrl$p < 0.05)
d_d <- (!gate_open) || (!is.null(d_null) && .both(d_null$p_emp < 0.05))

declared <- tibble::tibble(
  test = c("F3-pre  MB1 forkscale vs OXPHOS, both instruments",
           "D-a  MB2 three-way non-null and same sign on both instruments",
           "D-b  and NOT the assembly factors (specificity)",
           "D-c  and NOT the MB3 control axis (fork-specific)",
           "D-e  matched null cleared, IF D-a opened the gate"),
  predicted = c(paste0("verdict, not pass/fail (", f3pre_verdict, ")"),
                "no direction pre-specified; two-sided",
                "null / weaker", "null / weaker",
                "p_emp < 0.05 on both, or gate closed"),
  passes = c(NA, d_a, d_b, d_c, d_d))
declared %>% as.data.frame() %>% print(row.names = FALSE)

message("\n   F3-pre correlations, decision variable first:")
f3pre_cor %>% dplyr::filter(fork %in% c("MB1_forkscale", "MB1_forkscale_log")) %>%
  dplyr::select(fork, axis, n, rho_spearman, p_spearman, r_pearson) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   F3-pre incremental variance:")
dplyr::bind_rows(f3pre_r2, f3pre_r2_rev) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   Block D, the pre-registered three-way terms:")
coefs %>% dplyr::filter(grepl("three-way|within ER", label)) %>%
  dplyr::select(label, arm, instrument, n, estimate, ci_lo, ci_hi, p) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   Block D, DESCRIPTIVE endpoints (no prediction, excluded above):")
coefs %>% dplyr::filter(grepl("DESCRIPTIVE", label)) %>%
  dplyr::select(label, instrument, n, estimate, ci_lo, ci_hi, p) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   REMINDER: Block D modifies a MAIN EFFECT THAT IS NULL, and a ",
        "three-way\n   with a continuous modifier at n = ", length(idx),
        " is underpowered by construction.\n   A null here is uninformative. ",
        "Nothing here resurrects H1.")

# =============================================================================
# 7. Save
# =============================================================================
message("\n7. save")

out <- list(
  forkscale    = FK,
  coefficients = coefs,
  f3pre        = list(correlations = f3pre_cor, incremental_r2 =
                        dplyr::bind_rows(f3pre_r2, f3pre_r2_rev),
                      state = f3pre_state, verdict = f3pre_verdict,
                      decision_var = F3PRE_DECISION_VAR,
                      thresholds = c(redundant = F3PRE_REDUNDANT,
                                     independent = F3PRE_INDEPENDENT)),
  d_null       = d_null,
  gate_open    = gate_open,
  declared     = declared,
  coverage     = tibble::tibble(fork_patients = nrow(FK),
                                in_cohort = in_cohort,
                                block_c_n = length(idx_c),
                                analysis_n = length(idx)),
  spec = list(
    declaration = "docs/2026-08-29_G3_result_forkscale_availability.md (4, 4b)",
    forkscale   = "pc1 / index; the _log variant is pc1 / log(index)",
    mb3_sign    = paste("TCGA convention as stored upstream;",
                        "MB3_forkscale_mbconv is the METABRIC-convention copy"),
    block_d_dir = "no direction pre-specified; two-sided",
    block_d_pow = paste("three-way with a continuous modifier at n =",
                        length(idx), "- a null is UNINFORMATIVE"),
    d2          = "deferred - needs the Pommier developmental sets snapshotted",
    f3          = "F3 proper is a survival test and belongs in METABRIC"),
  built = Sys.time())

saveRDS(out, file.path(DIR_RESULTS, "forkscale_models.rds"))
readr::write_csv(coefs,     file.path(DIR_TABLES, "block_d_coefficients.csv"))
readr::write_csv(f3pre_cor, file.path(DIR_TABLES, "f3pre_correlations.csv"))
readr::write_csv(dplyr::bind_rows(f3pre_r2, f3pre_r2_rev),
                 file.path(DIR_TABLES, "f3pre_incremental_r2.csv"))
readr::write_csv(declared,  file.path(DIR_TABLES, "block_d_declared.csv"))
if (!is.null(d_null)) {
  readr::write_csv(d_null, file.path(DIR_TABLES, "block_d_matched_null.csv"))
}

message("\n15: done.")
message("    results/forkscale_models.rds")
message("    outputs/tables/  4-5 tables")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  f <- readRDS(file.path(DIR_RESULTS, "forkscale_models.rds"))

  # --- the verdict, and nothing else first ---------------------------------
  f$f3pre$verdict
  f$declared %>% as.data.frame() %>% print(row.names = FALSE)
  f$coverage %>% as.data.frame() %>% print(row.names = FALSE)

  # --- what forkscale actually looks like ----------------------------------
  # Skewed by construction. If a model on the raw scale is being leaned on,
  # look at this first.
  hist(f$forkscale$MB1_forkscale, breaks = 60,
       main = "MB1 forkscale = pc1 / index", xlab = "")
  plot(f$forkscale$MB1_index, f$forkscale$MB1_pc1, pch = 16, cex = 0.3,
       xlab = "MB1 sort index", ylab = "MB1 PC1")

  # --- the redundancy question, drawn -------------------------------------
  mi <- readRDS(file.path(DIR_RESULTS, "tcga_brca_mito_scores.rds"))
  ok <- match(f$forkscale$patient, colnames(mi$gsva_arms))
  plot(f$forkscale$MB1_forkscale, mi$gsva_arms["OXPHOS subunits", ok],
       pch = 16, cex = 0.4, xlab = "MB1 forkscale",
       ylab = "OXPHOS subunits (GSVA)")
  # and on ranks, which is what the verdict is taken on
  plot(rank(f$forkscale$MB1_forkscale), rank(mi$gsva_arms["OXPHOS subunits", ok]),
       pch = 16, cex = 0.3, xlab = "rank MB1 forkscale", ylab = "rank OXPHOS")

  # --- do the three biclusters agree with each other? ----------------------
  # If MB1, MB2 and MB3 forkscale are near-identical, "MB3 as a control axis"
  # is not a control at all and D-c is uninformative. Check before reading it.
  fk <- f$forkscale[, c("MB1_forkscale", "MB2_forkscale", "MB3_forkscale")]
  round(stats::cor(fk, method = "spearman"), 3)
  pairs(fk, pch = 16, cex = 0.2)

  # --- the full correlation panel ------------------------------------------
  f$f3pre$correlations %>%
    dplyr::filter(fork == "MB1_forkscale") %>%
    dplyr::arrange(dplyr::desc(abs(rho_spearman))) %>%
    as.data.frame() %>% print(row.names = FALSE)

  # --- Block D three-ways as a forest --------------------------------------
  fo <- f$coefficients %>% dplyr::filter(grepl("three-way", label)) %>%
    dplyr::arrange(estimate)
  with(fo, {
    plot(estimate, seq_along(estimate), pch = 16, yaxt = "n", ylab = "",
         xlim = range(c(ci_lo, ci_hi)), xlab = "MYC:OX:forkscale")
    segments(ci_lo, seq_along(estimate), ci_hi, seq_along(estimate))
    axis(2, seq_along(label), paste(label, arm, instrument), las = 1,
         cex.axis = 0.5)
    abline(v = 0, lty = 2)
  })

  # --- the ER mitigation, checked rather than assumed ----------------------
  cv <- readRDS(file.path(DIR_RESULTS, "tcga_brca_covariates.rds"))$covariates
  i  <- match(f$forkscale$patient, cv$patient)
  boxplot(f$forkscale$MB1_forkscale ~ cv$er_call[i], xlab = "ER call",
          ylab = "MB1 forkscale")
  table(cv$PAM50[i], cv$er_call[i], dnn = c("PAM50", "ER"))

  # --- the deferred pieces, so they are not forgotten ----------------------
  f$spec$d2      # Block D2, needs the Pommier sets
  f$f3pre$state  # NULL until script 11 has run

}
