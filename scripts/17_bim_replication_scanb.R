# 17_bim_replication_scanb.R
# =============================================================================
# THE BIM REPLICATION. The only pre-registered replication in this arm, and the
# last analysis it will run.
#
# Built to:
#   docs/2026-08-29_block_c_result_H1_not_supported.md section 9  <- the declaration
#   docs/2026-08-31_scanb_bim_replication_declaration.md          <- every parameter
#   docs/2026-08-29_script09_build_spec.md sections 1.1, 3        <- z-scaling, percentile
#
# =============================================================================
# WHAT IS BEING TESTED, AND WHAT WOULD COUNT AS PASSING
# =============================================================================
# From Block C, TCGA, n = 938, spine, S1:
#
#                    GSVA               mitoPPS
#   log2(BBC3)      -0.014 (p 0.60)    -0.034 (p 0.21)    PUMA:   unmoved
#   log2(BCL2L1)    -0.034 (p 0.030)   -0.043 (p 0.009)   BCL-XL: DOWN, both
#   log2(BCL2L11)   +0.051 (p 0.005)   +0.068 (p 3e-4)    BIM:    UP, both
#
# The registered prediction: BCL2L11 POSITIVE, BCL2L1 NEGATIVE, BBC3 NULL.
# Direction is pre-stated, so the tests are ONE-SIDED IN INTERPRETATION and
# two-sided in computation.
#
# BCL2L11 replicates if and only if all four hold (declaration section 9):
#   1. the M_a:OXPHOS interaction is POSITIVE, at p < 0.05,
#   2. on BOTH co-primary instruments,
#   3. at S1,
#   4. and the arm beats its expression-matched null on this endpoint.
# ONE INSTRUMENT ALONE IS NOT A REPLICATION. That is the rule the whole arm has
# been run under and it does not soften here.
#
# If it fails: "reported as a TCGA-specific observation and dropped. No further
# variants are tried." Not another cohort, not another covariate set, not
# another estimator, not a subtype stratum. Section 6 says so in the output.
#
# =============================================================================
# SECTION 2 IS A HARD GATE AND IT RUNS BEFORE SCAN-B IS TOUCHED
# =============================================================================
# The SCAN-B covariate set is four terms short of the Block C spine - purity,
# leukocyte fraction, TP53 status and plate do not exist there. Without a
# calibration, a SCAN-B null would be uninterpretable: it could mean BIM does
# not replicate, or it could mean the covariate reduction did it.
#
# So section 2 refits the three TCGA limbs on the SAME 938 patients with the
# SAME frozen z-scaling constants, under the reduced SCAN-B covariate set, and
# stops unless BCL2L11 keeps its positive sign AND the full-spine estimate lies
# inside the reduced-set 95% CI, on BOTH instruments. The pass condition was
# fixed in the declaration before either number existed.
#
# Section 2.2 first asserts that the rebuilt TCGA frame REPRODUCES script 09's
# saved coefficients to 1e-8. If it does not, the calibration is measuring the
# rebuild rather than the covariate set, and the script stops there instead.
#
# =============================================================================
# SCALE, AND WHERE EACH NUMBER COMES FROM
# =============================================================================
#   M_a, OXPHOS subunits, PROLIF_*   GSVA on the SCAN-B VST      (log)
#   OXPHOS subunits                  mitoPPS on SCAN-B linear    (linear)
#   log2(BCL2L11 / BCL2L1 / BBC3)    log2 of SCAN-B LINEAR
#
# The endpoints are log2 of the DESeq2-normalised LINEAR counts, NOT the VST -
# because that is what script 08 did in TCGA (`limbs <- t(log2(L[LIMB_GENES,]))`)
# and a replication that changes the endpoint's scale is not a replication.
# Verified 2026-08-31: none of the three genes has a zero in either cohort, so
# log2 is defined without a pseudocount.
#
# COHORT-RELATIVITY. Every SCAN-B score is computed in ONE run over ALL 3,207
# samples and z-scored WITHIN SCAN-B. Coefficients are therefore per SD of the
# SCAN-B distribution: comparable to TCGA in sign and rough magnitude, not in
# decimal places. mitoPPS values are NEVER compared numerically between the two
# cohorts - only the sign and the CI of the interaction transfer.
#
# THERE IS NO META-ANALYSIS. k = 1. Script 13's DerSimonian-Laird machinery is
# not invoked and this is described as a single-cohort replication.
#
# FENCE (declaration section 10): no survival, no RFI, no treatment strata, no
# forkscale, no PRIME, no STATE, no BUFFER_c, no specificity battery beyond the
# single null arm. Those columns are not even in results/scanb_pheno.rds.
#
# SPECIES: human. See CLAUDE.md.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

suppressPackageStartupMessages({
  library(GSVA)
})

message("\n17: the BIM replication in SCAN-B\n", strrep("=", 78))

# =============================================================================
# 0. Constants - all from the declaration, none chosen here
# =============================================================================
ARM_PRIMARY <- "OXPHOS subunits"
CI_LEVEL    <- 0.95

# Endpoint, predicted direction, and role. Fixed in the Block C declaration.
ENDPOINTS <- tibble::tibble(
  endpoint  = c("log2_BCL2L11", "log2_BCL2L1", "log2_BBC3"),
  gene      = c("BCL2L11", "BCL2L1", "BBC3"),
  predicted = c("positive", "negative", "null"),
  role      = c("THE REPLICATION", "corroborates direction", "the control"))

# The four specifications. S1 primary; S2 and S3 are the D7 sensitivities the
# declaration requires side by side; C-alt is BEYOND the Block C covariate set
# and CANNOT become the primary (declaration section 5.1).
SPECS <- list(
  S1      = c("PROLIF_DISJOINT", "PAM50"),
  S2      = c("PAM50"),
  S3      = c("PROLIF_STD", "PAM50"),
  `C-alt` = c("PROLIF_DISJOINT", "PAM50", "age", "NHG"))
SPEC_PRIMARY <- "S1"

INSTRUMENTS <- c("gsva", "mitopps")

# TCGA calibration (section 2). Script 09's constants, restated so this script
# is readable without opening that one; they are ASSERTED against the saved
# object rather than trusted.
BLOCK_C_VARS  <- c("purity", "leukocyte_fraction", "PAM50", "TP53_status", "plate")
SPINE_COVARS  <- c("purity", "leukocyte_fraction", "PROLIF_DISJOINT", "PAM50",
                   "TP53_status", "plate_f")
REDUCED_COVARS <- SPECS$S1              # the SCAN-B mapping, exactly
PLATE_MIN_N   <- 10L
REPRO_TOL     <- 1e-8

# The expression-matched null. Script 07's constants; the scope is one arm, one
# endpoint, both instruments, S1 only (declaration section 8).
N_DRAWS       <- 2000L
N_BIN         <- 20L
MIN_SET_GENES <- 3L
MIN_SET_FRAC  <- 0.80
NULL_SEED     <- PROJECT_SEED + 1700L
# GSVA rebuilds its gene-wise kernel estimate ONCE PER CALL over the whole
# 18,153 x 3,207 matrix, independent of how many sets the call carries - so
# chunking costs a full KDE pass per chunk and buys only memory. Two chunks is
# the compromise: the observed arm rides in both and is checked against the
# 6-set call from section 3.2, which tests the pin mechanism across batch sizes
# of 6, 1,002 and 1,002 rather than merely between two similar ones.
GSVA_CHUNK    <- 1000L
DIR_SCANB_NULL <- file.path(DIR_RESULTS, "scanb_null")
NULL_FORCE_REBUILD <- FALSE

MTDNA_PATHWAY <- "mtDNA-encoded OXPHOS subunits"

.rho <- function(x, y) suppressWarnings(
  stats::cor(x, y, method = "spearman", use = "pairwise.complete.obs"))

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")

# --- TCGA, for the gate ------------------------------------------------------
mito <- readRDS(file.path(DIR_RESULTS, "tcga_brca_mito_scores.rds"))
prim <- readRDS(file.path(DIR_RESULTS, "tcga_brca_priming.rds"))$priming
myc  <- readRDS(file.path(DIR_RESULTS, "tcga_brca_myc_scores.rds"))$estimators
cov  <- readRDS(file.path(DIR_RESULTS, "tcga_brca_covariates.rds"))$covariates
b09  <- readRDS(file.path(DIR_RESULTS, "block_c_models.rds"))

# --- SCAN-B ------------------------------------------------------------------
sv <- readRDS(file.path(DIR_RESULTS, "scanb_vst.rds"))
sl <- readRDS(file.path(DIR_RESULTS, "scanb_linear.rds"))
sp <- readRDS(file.path(DIR_RESULTS, "scanb_pheno.rds"))

if (!identical(sv$scale, "log_vst")) {
  stop("expected the LOG-scale SCAN-B object, got '", sv$scale, "'.", call. = FALSE)
}
if (!identical(sl$scale, "linear_deseq2_normalised")) {
  stop("expected the LINEAR SCAN-B object, got '", sl$scale, "'.", call. = FALSE)
}
E <- sv$mat          # LOG   - GSVA only
L <- sl$mat          # LINEAR - mitoPPS and the log2 endpoints only
stopifnot(identical(dimnames(E), dimnames(L)))
message("   SCAN-B: ", nrow(E), " genes x ", ncol(E), " samples")

# The fence, re-asserted on the object this script actually reads. Script 16
# enforced it at write time; a hand-edited rds would slip past that.
if (any(grepl("surviv|relapse|chemo|endocrine|esr[12]", names(sp$pheno),
              ignore.case = TRUE))) {
  stop("FENCE BREACH: results/scanb_pheno.rds carries an outcome or treatment ",
       "column. Declaration section 10 forbids it. Re-run script 16.",
       call. = FALSE)
}
message("   fence: no outcome or treatment column in the SCAN-B phenotype")

# =============================================================================
# 2. THE HARD GATE - TCGA calibration under the reduced covariate set
# =============================================================================
# Declaration section 6. Nothing in SCAN-B is fitted until this passes.
message("\n", strrep("=", 78),
        "\n2. HARD GATE: does the TCGA result survive the SCAN-B covariate set?\n",
        strrep("=", 78))

# --- 2.1 rebuild the TCGA frame ----------------------------------------------
# Rebuilt rather than loaded, because script 09 saves coefficients and not its
# design frame. The FROZEN CONSTANTS are taken from the saved object instead of
# recomputed - the z-scaling and the plate map - so the only thing being
# reconstructed here is the assembly, and 2.2 checks that.
pat <- colnames(mito$gsva_arms)
stopifnot(identical(sort(pat), sort(cov$patient)),
          identical(sort(pat), sort(prim$patient)),
          identical(sort(pat), sort(myc$patient)))

Dt <- tibble::tibble(patient = pat) %>%
  dplyr::left_join(dplyr::select(cov, patient, purity, leukocyte_fraction,
                                 PAM50, TP53_status, plate), by = "patient") %>%
  dplyr::left_join(dplyr::select(myc, patient, M_a), by = "patient") %>%
  dplyr::left_join(dplyr::select(prim, patient,
                                 dplyr::all_of(ENDPOINTS$endpoint)),
                   by = "patient") %>%
  dplyr::mutate(
    PROLIF_DISJOINT = as.numeric(mito$gsva_cov["PROLIF_DISJOINT", patient]),
    PROLIF_STD      = as.numeric(mito$gsva_cov["PROLIF_STD", patient]),
    PAM50           = factor(PAM50),
    TP53_status     = factor(TP53_status))
stopifnot(identical(Dt$patient, pat))

zc <- b09$zscale
.apply_z <- function(v, key) {
  k <- zc[[key]]
  if (is.null(k)) stop("no frozen z-scaling constant for '", key, "'.",
                       call. = FALSE)
  (v - k[["mean"]]) / k[["sd"]]
}
Dt$M_a <- .apply_z(Dt$M_a, "M_a")

OXt <- list(
  gsva    = .apply_z(as.numeric(mito$gsva_arms[ARM_PRIMARY, pat]),
                     paste("gsva", ARM_PRIMARY, sep = "|")),
  mitopps = .apply_z(as.numeric(mito$mitopps_arms[ARM_PRIMARY, pat]),
                     paste("mitopps", ARM_PRIMARY, sep = "|")))

M1 <- stats::complete.cases(Dt[, BLOCK_C_VARS])
if (sum(M1) != b09$set_sizes[["M1"]]) {
  stop("the rebuilt M1 set has ", sum(M1), " patients, script 09 had ",
       b09$set_sizes[["M1"]], ". The covariate table has changed under this ",
       "script; the calibration cannot be trusted.", call. = FALSE)
}
Dt$plate_f <- factor(ifelse(Dt$plate %in% b09$plate_keep, Dt$plate, "other"))
message("   frame rebuilt: n = ", sum(M1), ", plate pooled to ",
        dplyr::n_distinct(Dt$plate_f[M1]), " levels (frozen map from script 09)")

.fit_tcga <- function(endpoint, instrument, covars, label) {
  df <- tibble::tibble(Y = Dt[[endpoint]][M1], MYC = Dt$M_a[M1],
                       OX = OXt[[instrument]][M1])
  for (v in covars) df[[v]] <- if (v == "PAM50") droplevels(Dt$PAM50[M1])
                               else if (v == "TP53_status") droplevels(Dt$TP53_status[M1])
                               else if (v == "plate_f") droplevels(Dt$plate_f[M1])
                               else Dt[[v]][M1]
  df <- df[stats::complete.cases(df), , drop = FALSE]
  fo <- stats::as.formula(paste("Y ~ MYC * OX +",
                                paste(covars, collapse = " + ")))
  m  <- stats::lm(fo, data = df)
  co <- summary(m)$coefficients
  if (!"MYC:OX" %in% rownames(co)) {
    stop("TCGA fit '", label, "' (", instrument, "): the interaction was ",
         "aliased away.", call. = FALSE)
  }
  tcrit <- stats::qt(1 - (1 - CI_LEVEL) / 2, df = m$df.residual)
  tibble::tibble(label = label, endpoint = endpoint, instrument = instrument,
                 n = nrow(df),
                 estimate = co["MYC:OX", 1], se = co["MYC:OX", 2],
                 ci_lo = co["MYC:OX", 1] - tcrit * co["MYC:OX", 2],
                 ci_hi = co["MYC:OX", 1] + tcrit * co["MYC:OX", 2],
                 p = co["MYC:OX", 4])
}

# --- 2.2 REPRODUCTION ASSERTION ----------------------------------------------
# If the rebuilt frame does not give script 09's numbers to 1e-8, then whatever
# section 2.3 measures is the rebuild, not the covariate reduction. This is the
# same discipline as script 07 asserting its algebraic mitoPPS rewrite against
# the definitional form, and script 09 asserting .fit against .fast_b.
message("\n2.2 reproduction check against block_c_models.rds")

spine <- dplyr::bind_rows(lapply(ENDPOINTS$endpoint, function(ep)
  dplyr::bind_rows(lapply(INSTRUMENTS, function(ins)
    .fit_tcga(ep, ins, SPINE_COVARS, "full spine")))))

saved <- b09$coefficients %>%
  dplyr::filter(endpoint %in% ENDPOINTS$endpoint, arm == ARM_PRIMARY,
                myc == "M_a", spec == "S1", dataset == "M1", stratum == "-") %>%
  dplyr::select(endpoint, instrument, saved_estimate = estimate,
                saved_se = se, saved_n = n)

chk <- dplyr::left_join(spine, saved, by = c("endpoint", "instrument")) %>%
  dplyr::mutate(d_estimate = abs(estimate - saved_estimate),
                d_se       = abs(se - saved_se),
                d_n        = n - saved_n)
if (anyNA(chk$saved_estimate)) {
  stop("script 09 has no saved limb coefficient for: ",
       paste(unique(paste(chk$endpoint, chk$instrument)[is.na(chk$saved_estimate)]),
             collapse = ", "),
       ". Re-source scripts/09_interaction_models.R.", call. = FALSE)
}
chk %>%
  dplyr::select(endpoint, instrument, n, d_n, estimate, saved_estimate,
                d_estimate, d_se) %>%
  as.data.frame() %>% print(row.names = FALSE)

if (max(chk$d_estimate, chk$d_se) > REPRO_TOL || any(chk$d_n != 0L)) {
  stop("the rebuilt TCGA frame does NOT reproduce script 09 (max |d| = ",
       signif(max(chk$d_estimate, chk$d_se), 3), ", tolerance ", REPRO_TOL,
       ").\nThe calibration below would be measuring the rebuild rather than ",
       "the covariate reduction. Fix the assembly before going further.",
       call. = FALSE)
}
message("   reproduced to within ", REPRO_TOL,
        " on all ", nrow(chk), " limb x instrument fits")

# --- 2.3 the same limbs under the SCAN-B covariate set -----------------------
message("\n2.3 refit under the SCAN-B covariate set (",
        paste(REDUCED_COVARS, collapse = " + "), ")")

reduced <- dplyr::bind_rows(lapply(ENDPOINTS$endpoint, function(ep)
  dplyr::bind_rows(lapply(INSTRUMENTS, function(ins)
    .fit_tcga(ep, ins, REDUCED_COVARS, "reduced")))))

calib <- dplyr::inner_join(
  dplyr::select(spine, endpoint, instrument,
                full_estimate = estimate, full_ci_lo = ci_lo,
                full_ci_hi = ci_hi, full_p = p),
  dplyr::select(reduced, endpoint, instrument,
                red_estimate = estimate, red_ci_lo = ci_lo,
                red_ci_hi = ci_hi, red_p = p),
  by = c("endpoint", "instrument")) %>%
  dplyr::mutate(
    sign_kept    = sign(red_estimate) == sign(full_estimate),
    full_in_red_ci = full_estimate >= red_ci_lo & full_estimate <= red_ci_hi)

calib %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- 2.4 the gate ------------------------------------------------------------
# Declaration section 6, verbatim: for BCL2L11, on BOTH instruments, the
# reduced-set estimate must (a) keep the positive sign, and (b) have the
# full-spine point estimate inside its 95% CI.
g <- calib[calib$endpoint == "log2_BCL2L11", ]
stopifnot(nrow(g) == length(INSTRUMENTS))
gate_sign <- all(g$red_estimate > 0)
gate_ci   <- all(g$full_in_red_ci)

message("\n2.4 GATE")
message("   BCL2L11 positive on both instruments under the reduced set: ",
        gate_sign)
message("   full-spine estimate inside the reduced-set 95% CI, both: ", gate_ci)

if (!(gate_sign && gate_ci)) {
  stop("\n", strrep("!", 78),
       "\nGATE FAILED. The Block C BCL2L11 result does not survive the ",
       "covariate\nset that SCAN-B can support, so a SCAN-B result could not ",
       "be read: a null\nthere would be indistinguishable from the covariate ",
       "reduction doing it.\n\n",
       "Declaration section 6: the replication is reported as NOT RUNNABLE AS ",
       "DECLARED,\nand the arm falls back to option 2 - delete the commitment ",
       "- with this\noutput as the reason. NOTHING IN SCAN-B IS FITTED.\n",
       strrep("!", 78), call. = FALSE)
}
message("   GATE PASSED. Proceeding to SCAN-B.")

# =============================================================================
# 3. Score SCAN-B
# =============================================================================
message("\n", strrep("=", 78), "\n3. scoring SCAN-B\n", strrep("=", 78))

MATRIX_SYMBOLS <- rownames(E)
.in_matrix <- function(g) intersect(g, MATRIX_SYMBOLS)

# --- 3.1 the gene sets, harmonised to the 2014 vocabulary --------------------
# Addendum A. The map is rebuilt here over the FULL union - the four declared
# sets PLUS every MitoCarta pathway gene, because mitoPPS needs the whole
# universe and script 16's map covered only the four sets. Same rule, wider
# input, and it is ASSERTED to agree with script 16's map on the overlap: if a
# wider gene universe changes a resolution, that is context-dependence and it
# has to be understood rather than absorbed.
message("\n3.1 gene sets")

mitocarta_inventory  <- suppressWarnings(readxl::read_xls(PATH_MITOCARTA, sheet = 2))
mitocarta_background <- suppressWarnings(readxl::read_xls(PATH_MITOCARTA, sheet = 3))
mitocarta_pathways   <- suppressWarnings(readxl::read_xls(PATH_MITOCARTA, sheet = 4))
stopifnot(all(c("Symbol", "Synonyms") %in% colnames(mitocarta_background)))

.split_genes <- function(x) trimws(unlist(strsplit(x, ",", fixed = TRUE)))

mito_paths_raw <- stats::setNames(
  lapply(mitocarta_pathways$Genes, .split_genes), mitocarta_pathways$MitoPathway)
mtdna_genes <- grep("^MT-", unique(mitocarta_inventory$Symbol), value = TRUE)
stopifnot(length(mtdna_genes) == EXPECT_MITOCARTA_MTDNA,
          !MTDNA_PATHWAY %in% names(mito_paths_raw))
mito_paths <- lapply(mito_paths_raw, function(g) setdiff(g, mtdna_genes))
mito_paths[[MTDNA_PATHWAY]] <- mtdna_genes

g1         <- readRDS(file.path(DIR_RESULTS, "g1_overlap_audit.rds"))
felsher_ma <- g1$estimators_stripped$FELSHER
stopifnot(length(felsher_ma) == 61L)

declared_sets <- list(
  M_a               = felsher_ma,
  `OXPHOS subunits` = mito$arm_sets[[ARM_PRIMARY]],
  PROLIF_DISJOINT   = mito$covariate_sets$PROLIF_DISJOINT,
  PROLIF_STD        = mito$covariate_sets$PROLIF_STD)

MC_SYMBOLS <- unique(mitocarta_background$Symbol)
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
    if (length(cand) == 1L) { out[[g]] <- cand; status[[i]] <- "resolved" }
    else if (length(cand) > 1L) status[[i]] <- "ambiguous"
  }
  list(map = out, report = tibble::tibble(input_symbol = genes, status = status,
                                          resolved_to = unname(out[genes])))
}

bm <- .build_symbol_map(c(unlist(declared_sets, use.names = FALSE),
                          unlist(mito_paths,    use.names = FALSE)))
symbol_map <- bm$map
print(table(bm$report$status))

mapped_in <- unname(symbol_map)[unname(symbol_map) %in% MATRIX_SYMBOLS]
if (anyDuplicated(mapped_in)) {
  stop("the symbol map sends two inputs to the same matrix row: ",
       paste(unique(mapped_in[duplicated(mapped_in)]), collapse = ", "),
       call. = FALSE)
}

# The agreement check against script 16.
common <- intersect(names(symbol_map), names(sp$symbol_map))
disagree <- common[symbol_map[common] != sp$symbol_map[common]]
if (length(disagree)) {
  stop("the symbol map disagrees with script 16's on ", length(disagree),
       " symbol(s): ", paste(utils::head(disagree, 10), collapse = ", "),
       ".\nA wider gene universe changed a resolution, which means the map is ",
       "context-dependent. Understand this before scoring anything.",
       call. = FALSE)
}
message("   agrees with script 16's map on all ", length(common),
        " shared symbols")

.remap <- function(g) {
  h <- symbol_map[g]
  unname(ifelse(is.na(h), g, h))
}

sets_scanb  <- lapply(declared_sets, function(g) .in_matrix(.remap(unique(g))))
paths_scanb <- lapply(mito_paths,    function(g) .in_matrix(.remap(unique(g))))

covg <- tibble::tibble(
  set       = names(declared_sets),
  n_defined = vapply(declared_sets, function(g) length(unique(g)), integer(1)),
  n_present = vapply(sets_scanb, length, integer(1))) %>%
  dplyr::mutate(frac = n_present / n_defined)
covg %>% as.data.frame() %>% print(row.names = FALSE)

# ENFORCED here, as script 16 promised. A set below the floor is not a
# measurement of what it is named after.
low <- covg$set[covg$frac < MIN_SET_FRAC]
if (length(low)) {
  stop("set(s) below the ", MIN_SET_FRAC, " coverage floor in SCAN-B after ",
       "harmonisation: ", paste(low, collapse = ", "),
       ". The declared model cannot be fitted as specified.", call. = FALSE)
}
message("   all declared sets at or above the ", MIN_SET_FRAC, " floor")

# --- 3.2 INSTRUMENT 1: GSVA on the VST ---------------------------------------
# The pins hold the gene universe at the full matrix in EVERY call, so a batch
# of 4 sets and a batch of 400 random ones walk the same universe. Without them
# the null is silently invalid. Two halves rather than one all-genes set,
# because a set containing every gene leaves the walk's miss-penalty at 0/0.
message("\n3.2 INSTRUMENT 1: GSVA (log VST, kcdf Gaussian)")

.PIN_A <- MATRIX_SYMBOLS[c(TRUE, FALSE)]
.PIN_B <- MATRIX_SYMBOLS[c(FALSE, TRUE)]
stopifnot(setequal(c(.PIN_A, .PIN_B), MATRIX_SYMBOLS))

.gsva_batch <- function(sets, label = "") {
  n_ok <- vapply(sets, length, integer(1))
  if (any(n_ok < MIN_SET_GENES)) {
    stop("GSVA batch '", label, "': set(s) below the size floor -> ",
         paste(names(sets)[n_ok < MIN_SET_GENES], collapse = ", "),
         call. = FALSE)
  }
  wanted <- names(sets)
  sets[[".PIN_A"]] <- .PIN_A
  sets[[".PIN_B"]] <- .PIN_B
  par <- GSVA::gsvaParam(exprData = E, geneSets = sets, kcdf = "Gaussian",
                         minSize = MIN_SET_GENES, maxSize = Inf)
  s <- GSVA::gsva(par, verbose = FALSE)
  dropped <- setdiff(wanted, rownames(s))
  if (length(dropped)) {
    stop("GSVA silently dropped set(s) in batch '", label, "': ",
         paste(utils::head(dropped, 10), collapse = ", "), call. = FALSE)
  }
  s[wanted, , drop = FALSE]
}

gsva_scores <- .gsva_batch(sets_scanb, "declared sets")
message("   scored ", nrow(gsva_scores), " sets x ", ncol(gsva_scores),
        " samples")

# --- 3.3 INSTRUMENT 2: mitoPPS on the linear matrix --------------------------
# Script 07 section 5, unchanged. mitoPPS reports the SHAPE of the
# mitochondrial program and is deliberately blind to total content, which is
# why its values never travel between cohorts - only the interaction's sign
# and CI do.
message("\n3.3 INSTRUMENT 2: mitoPPS (linear DESeq2-normalised)")

.pathway_score_matrix <- function(sets) {
  keep <- sets[vapply(sets, length, integer(1)) >= MIN_SET_GENES]
  out  <- t(vapply(keep, function(g) colMeans(L[g, , drop = FALSE]),
                   numeric(ncol(L))))
  rownames(out) <- names(keep); colnames(out) <- colnames(L)
  out
}
.mitopps_universe <- function(S) {
  N <- ncol(S); P <- nrow(S)
  Bi <- 1 / S
  A  <- (S %*% t(Bi)) / N
  M  <- (1 / A) %*% Bi
  out <- S * (M - Bi) / (P - 1)
  dimnames(out) <- dimnames(S)
  out
}
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

S_universe <- .pathway_score_matrix(paths_scanb)
if (min(S_universe) <= 0) {
  stop("a MitoPathway score is zero or negative in some SCAN-B sample; the ",
       "pairwise ratio is undefined. Inspect before proceeding.", call. = FALSE)
}
message("   universe: ", nrow(S_universe), " MitoPathways with >= ",
        MIN_SET_GENES, " genes (of ", length(paths_scanb), ")")
if (!ARM_PRIMARY %in% rownames(S_universe)) {
  stop("'", ARM_PRIMARY, "' is not in the SCAN-B mitoPPS universe.",
       call. = FALSE)
}

mitopps_universe <- .mitopps_universe(S_universe)
message(sprintf("   global mean %.4f (should be ~1), range %.3f to %.3f",
                mean(mitopps_universe), min(mitopps_universe),
                max(mitopps_universe)))

# The observed arm and every null draw must go through the SAME denominator:
# the universe with OXPHOS subunits held out. Script 07's .arm_universe rule.
U_ox   <- S_universe[setdiff(rownames(S_universe), ARM_PRIMARY), , drop = FALSE]
ox_obs <- as.numeric(.mitopps_query(S_universe[ARM_PRIMARY, , drop = FALSE], U_ox))
d_canon <- max(abs(ox_obs - mitopps_universe[ARM_PRIMARY, ]))
if (d_canon > 1e-8) {
  stop("the held-out query does not reproduce the canonical mitoPPS for '",
       ARM_PRIMARY, "' (diff ", signif(d_canon, 3), ").", call. = FALSE)
}
message("   held-out query reproduces the canonical value (diff ",
        signif(d_canon, 3), ")")

# --- 3.4 z-scaling, within SCAN-B --------------------------------------------
# Cohort-relative, so the constants are SCAN-B's own. Coefficients are per SD of
# the SCAN-B distribution and are NOT numerically comparable to TCGA's.
message("\n3.4 z-scaling within SCAN-B")

# The reference set is the S1 ANALYSIS SET, not all 3,207 - script 09's rule
# (build spec 1.1: constants computed once on the primary set and reused
# everywhere). Everything the S1 model needs is complete except PAM50, so the
# S1 set is exactly the PAM50-complete samples; asserted rather than assumed.
ZREF <- !is.na(sp$pheno$PAM50)
stopifnot(!anyNA(as.numeric(gsva_scores["M_a", ])),
          !anyNA(as.numeric(gsva_scores[ARM_PRIMARY, ])),
          !anyNA(as.numeric(gsva_scores["PROLIF_DISJOINT", ])),
          !anyNA(ox_obs))
message("   z-scaling reference set (S1 complete cases): ", sum(ZREF),
        " of ", length(ZREF))

.zfix <- function(v) {
  m <- mean(v[ZREF]); sdv <- stats::sd(v[ZREF])
  if (!is.finite(sdv) || sdv == 0) stop("zero variance in a z-scaled variable",
                                        call. = FALSE)
  (v - m) / sdv
}
ZS <- list()
for (nm in names(sets_scanb)) {
  v <- as.numeric(gsva_scores[nm, ])
  ZS[[paste("gsva", nm, sep = "|")]] <- c(mean = mean(v[ZREF]),
                                          sd = stats::sd(v[ZREF]))
}
ZS[[paste("mitopps", ARM_PRIMARY, sep = "|")]] <-
  c(mean = mean(ox_obs[ZREF]), sd = stats::sd(ox_obs[ZREF]))

M_a_s   <- .zfix(as.numeric(gsva_scores["M_a", ]))
OXs     <- list(gsva    = .zfix(as.numeric(gsva_scores[ARM_PRIMARY, ])),
                mitopps = .zfix(ox_obs))
PROLIF  <- list(PROLIF_DISJOINT = as.numeric(gsva_scores["PROLIF_DISJOINT", ]),
                PROLIF_STD      = as.numeric(gsva_scores["PROLIF_STD", ]))

message(sprintf("   GSVA vs mitoPPS on %s: rho = %.3f", ARM_PRIMARY,
                .rho(OXs$gsva, OXs$mitopps)))
message("   (the two instruments measure LEVEL and COMPOSITION and are not ",
        "expected to agree; printed before any coefficient is seen)")

# --- 3.5 the endpoints -------------------------------------------------------
# log2 of the LINEAR matrix, matching script 08's construction exactly.
message("\n3.5 endpoints (log2 of linear, as in script 08)")

absent <- setdiff(ENDPOINTS$gene, rownames(L))
if (length(absent)) {
  stop("endpoint gene(s) absent from SCAN-B: ", paste(absent, collapse = ", "),
       call. = FALSE)
}
bad <- ENDPOINTS$gene[apply(L[ENDPOINTS$gene, , drop = FALSE], 1L,
                            function(v) any(v <= 0))]
if (length(bad)) {
  stop("endpoint gene(s) with a zero or negative normalised count in at least ",
       "one sample -> ", paste(bad, collapse = ", "),
       ".\nlog2 is undefined there. This needs a decision (pseudocount, or ",
       "dropping the gene), not a default - script 08's rule.", call. = FALSE)
}
Y <- t(log2(L[ENDPOINTS$gene, , drop = FALSE]))
colnames(Y) <- ENDPOINTS$endpoint

# --- 3.6 the SCAN-B analysis frame -------------------------------------------
ph <- sp$pheno
stopifnot(identical(ph$sample_id, colnames(E)))

Ds <- tibble::tibble(
  sample_id       = ph$sample_id,
  M_a             = M_a_s,
  PROLIF_DISJOINT = PROLIF$PROLIF_DISJOINT,
  PROLIF_STD      = PROLIF$PROLIF_STD,
  PAM50           = droplevels(factor(ph$PAM50)),
  age             = ph$age,
  NHG             = droplevels(factor(ph$NHG)))
Ds <- dplyr::bind_cols(Ds, tibble::as_tibble(Y))

message("   PAM50 levels: ", paste(levels(Ds$PAM50), collapse = ", "))
message("   complete on PAM50: ", sum(!is.na(Ds$PAM50)), " of ", nrow(Ds))

# =============================================================================
# 4. The declared fits - 4 specifications x 2 instruments x 3 endpoints
# =============================================================================
# The whole model space. There are no others, and none is selected between.
message("\n", strrep("=", 78), "\n4. the 24 declared fits\n", strrep("=", 78))

.fit_scanb <- function(endpoint, instrument, spec) {
  covars <- SPECS[[spec]]
  df <- tibble::tibble(Y = Ds[[endpoint]], MYC = Ds$M_a,
                       OX = OXs[[instrument]])
  for (v in covars) df[[v]] <- Ds[[v]]
  df <- df[stats::complete.cases(df), , drop = FALSE]
  for (v in covars) if (is.factor(df[[v]])) df[[v]] <- droplevels(df[[v]])
  fo <- stats::as.formula(paste("Y ~ MYC * OX +",
                                paste(covars, collapse = " + ")))
  m  <- stats::lm(fo, data = df)
  co <- summary(m)$coefficients
  if (!"MYC:OX" %in% rownames(co)) {
    stop("SCAN-B fit ", endpoint, "/", instrument, "/", spec,
         ": the interaction was aliased away.", call. = FALSE)
  }
  tcrit <- stats::qt(1 - (1 - CI_LEVEL) / 2, df = m$df.residual)
  tibble::tibble(endpoint = endpoint, instrument = instrument, spec = spec,
                 n = nrow(df),
                 estimate = co["MYC:OX", 1], se = co["MYC:OX", 2],
                 ci_lo = co["MYC:OX", 1] - tcrit * co["MYC:OX", 2],
                 ci_hi = co["MYC:OX", 1] + tcrit * co["MYC:OX", 2],
                 p = co["MYC:OX", 4])
}

fits <- dplyr::bind_rows(lapply(ENDPOINTS$endpoint, function(ep)
  dplyr::bind_rows(lapply(INSTRUMENTS, function(ins)
    dplyr::bind_rows(lapply(names(SPECS), function(sc)
      .fit_scanb(ep, ins, sc))))))) %>%
  dplyr::left_join(dplyr::select(ENDPOINTS, endpoint, predicted, role),
                   by = "endpoint")

message("\n   ALL 24 FITS. S1 is the primary; the rest are reported, never ",
        "selected.")
fits %>%
  dplyr::mutate(dplyr::across(c(estimate, se, ci_lo, ci_hi), ~ round(.x, 4)),
                p = signif(p, 3)) %>%
  dplyr::select(endpoint, predicted, instrument, spec, n, estimate,
                ci_lo, ci_hi, p) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 5. The expression-matched null
# =============================================================================
# Declaration section 8. Scope is exactly one arm, one endpoint, both
# instruments, S1 only - not the seventeen-arm battery of Block C.
message("\n", strrep("=", 78),
        "\n5. expression-matched null: ", ARM_PRIMARY, " on log2_BCL2L11\n",
        strrep("=", 78))

.ensure_dir(DIR_SCANB_NULL)
NULL_FILE <- file.path(DIR_SCANB_NULL, "oxphos_subunits_bcl2l11.rds")

gene_mean_lin <- rowMeans(L)
rho_means <- .rho(gene_mean_lin, rowMeans(E))
message(sprintf("   Spearman(mean linear, mean VST) across genes = %.4f",
                rho_means))
if (rho_means < 0.95) {
  warning("linear and VST expression orderings disagree (rho = ",
          round(rho_means, 3), "); ventiles matched on one are not matched on ",
          "the other.", call. = FALSE)
}

bin_of <- stats::setNames(
  cut(rank(gene_mean_lin, ties.method = "first"), breaks = N_BIN,
      labels = FALSE), rownames(L))
by_bin <- split(names(bin_of), bin_of)

.draw_matched <- function(genes) {
  b <- bin_of[genes]; b <- b[!is.na(b)]
  unlist(lapply(split(b, b), function(k) {
    pool <- by_bin[[as.character(k[[1]])]]
    if (length(k) > length(pool)) {
      stop("ventile ", k[[1]], " holds ", length(pool), " genes but ",
           length(k), " are needed.", call. = FALSE)
    }
    pool[sample.int(length(pool), length(k))]
  }), use.names = FALSE)
}

ox_genes <- sets_scanb[[ARM_PRIMARY]]

if (file.exists(NULL_FILE) && !NULL_FORCE_REBUILD) {
  nullobj <- readRDS(NULL_FILE)
  if (!identical(nullobj$genes, ox_genes) || nullobj$seed != NULL_SEED ||
      nrow(nullobj$gsva) != N_DRAWS) {
    stop("the cached null was built for a different arm, seed or draw count. ",
         "Delete ", NULL_FILE, " and re-run.", call. = FALSE)
  }
  message("   cached (", nrow(nullobj$gsva), " draws)")
} else {
  t0 <- Sys.time()
  set.seed(NULL_SEED)
  draw_names <- sprintf("null_%04d", seq_len(N_DRAWS))
  draws <- stats::setNames(
    lapply(seq_len(N_DRAWS), function(i) .draw_matched(ox_genes)), draw_names)
  stopifnot(all(vapply(draws, length, integer(1)) == length(ox_genes)))
  message("   ", N_DRAWS, " matched draws, ", length(ox_genes),
          " genes each, ventile-matched on mean linear expression")

  # --- GSVA, in chunks, each carrying the observed arm as a stability check --
  # The pins hold the universe, and .OBS_CHECK is the braces: if the observed
  # arm's score moves between chunks, the universe moved and every percentile
  # would be wrong with nothing to show for it.
  message("   GSVA over ", N_DRAWS, " draws in chunks of ", GSVA_CHUNK,
          " (expect roughly 10-20 minutes; the refits that follow take ~2 s)")
  idx_chunks <- split(seq_len(N_DRAWS), ceiling(seq_len(N_DRAWS) / GSVA_CHUNK))
  gsva_null <- matrix(NA_real_, nrow = N_DRAWS, ncol = ncol(E),
                      dimnames = list(draw_names, colnames(E)))
  obs_ref <- as.numeric(gsva_scores[ARM_PRIMARY, ])
  for (ci in seq_along(idx_chunks)) {
    ii <- idx_chunks[[ci]]
    batch <- c(stats::setNames(list(ox_genes), ".OBS_CHECK"), draws[ii])
    s <- .gsva_batch(batch, paste0("null chunk ", ci))
    d <- max(abs(as.numeric(s[".OBS_CHECK", ]) - obs_ref))
    if (d > 1e-8) {
      stop("the observed arm's GSVA score moved by ", signif(d, 3),
           " inside null chunk ", ci, ". The gene universe is not being held ",
           "constant and every percentile would be invalid.", call. = FALSE)
    }
    gsva_null[ii, ] <- s[names(draws)[ii], , drop = FALSE]
    message(sprintf("     chunk %d/%d done", ci, length(idx_chunks)))
  }

  # --- mitoPPS, one matrix call ---------------------------------------------
  S_draws <- t(vapply(draws, function(g) colMeans(L[g, , drop = FALSE]),
                      numeric(ncol(L))))
  rownames(S_draws) <- draw_names
  if (min(S_draws) <= 0) {
    stop("a null draw has a non-positive pathway score; mitoPPS is undefined.",
         call. = FALSE)
  }
  mitopps_null <- .mitopps_query(S_draws, U_ox)

  nullobj <- list(genes = ox_genes, seed = NULL_SEED, gsva = gsva_null,
                  mitopps = mitopps_null, built = Sys.time())
  saveRDS(nullobj, NULL_FILE)
  message("   built in ",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
          " min -> ", basename(NULL_FILE))
}

# --- 5.1 refit the model against every draw ----------------------------------
# A design matrix built once, two columns swapped per draw. ASSERTED against
# lm() on the observed arm, because two implementations of one model diverging
# silently is the failure mode that would invalidate the percentile without
# producing an error.
message("\n5.1 refitting against each draw")

.fast_setup <- function(endpoint, spec) {
  covars <- SPECS[[spec]]
  df <- tibble::tibble(Y = Ds[[endpoint]], MYC = Ds$M_a)
  for (v in covars) df[[v]] <- Ds[[v]]
  ok <- stats::complete.cases(df)
  df <- df[ok, , drop = FALSE]
  for (v in covars) if (is.factor(df[[v]])) df[[v]] <- droplevels(df[[v]])
  X0 <- stats::model.matrix(
    stats::as.formula(paste("~ MYC +", paste(covars, collapse = " + "))),
    data = df)
  list(ok = ok, y = df$Y, myc = df$MYC, X0 = X0)
}

.fast_b <- function(S, ox) {
  X <- cbind(S$X0, OX = ox, `MYC:OX` = S$myc * ox)
  unname(.lm.fit(X, S$y)$coefficients[ncol(X)])
}

null_res <- dplyr::bind_rows(lapply(INSTRUMENTS, function(ins) {
  S  <- .fast_setup("log2_BCL2L11", SPEC_PRIMARY)
  b_obs_lm <- fits$estimate[fits$endpoint == "log2_BCL2L11" &
                            fits$instrument == ins & fits$spec == SPEC_PRIMARY]
  b_obs <- .fast_b(S, OXs[[ins]][S$ok])
  if (abs(b_obs - b_obs_lm) > 1e-9) {
    stop("the fast refit disagrees with lm() for ", ins, " (", signif(b_obs, 6),
         " vs ", signif(b_obs_lm, 6), "). Every percentile below would be ",
         "computed by a different model from the one reported.", call. = FALSE)
  }
  Nm <- if (ins == "gsva") nullobj$gsva else nullobj$mitopps
  # Each draw is z-scored on ITS OWN mean and sd - it is a different variable,
  # and not doing so would make the null a test of scale (build spec 1.1).
  b_null <- vapply(seq_len(nrow(Nm)), function(i)
    .fast_b(S, .zfix(as.numeric(Nm[i, ]))[S$ok]), numeric(1))
  tibble::tibble(
    instrument = ins, endpoint = "log2_BCL2L11", spec = SPEC_PRIMARY,
    n = length(S$y), b_obs = b_obs,
    null_mean = mean(b_null), null_sd = stats::sd(b_null),
    percentile = 100 * mean(b_null < b_obs),
    p_emp = 2 * min(mean(b_null <= b_obs), mean(b_null >= b_obs)),
    b_null = list(b_null))
}))

null_res %>% dplyr::select(-b_null) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 4))) %>%
  as.data.frame() %>% print(row.names = FALSE)
# p_emp = 2 * min(...) and either tail can be exactly 0, so the smallest
# NON-ZERO value is 2/N_DRAWS. A reported 0 means "more extreme than all
# draws" and must be written as p < 2/N_DRAWS, never as 0 (build spec 3).
P_EMP_FLOOR <- 2 / N_DRAWS
if (any(null_res$p_emp == 0)) {
  message("   NOTE: p_emp is exactly 0 for ",
          paste(null_res$instrument[null_res$p_emp == 0], collapse = ", "),
          " - the observed coefficient is more extreme than all ", N_DRAWS,
          " draws. REPORT AS p < ", P_EMP_FLOOR, ", never as 0.")
}
message("   p_emp floor with ", N_DRAWS, " draws is ", P_EMP_FLOOR)

# =============================================================================
# 6. The four replication conditions
# =============================================================================
message("\n", strrep("=", 78), "\n6. VERDICT\n", strrep("=", 78))

s1 <- fits %>% dplyr::filter(endpoint == "log2_BCL2L11",
                             spec == SPEC_PRIMARY)

# Condition 3 is STRUCTURAL, not evidential. It asserts that the two co-primary
# fits exist and are the declared pair - it cannot fail unless the script is
# broken, and it must never be read as one of three passes out of four. This is
# the F1-b error (a "condition" that only a coding mistake could fail) named
# here so it cannot be re-inherited; see the F1 note section 6.
c3 <- identical(sort(s1$instrument), sort(INSTRUMENTS)) && all(s1$n > 0)
if (!c3) {
  stop("the two co-primary S1 fits are not both present. This is a bug in ",
       "this script, not a result.", call. = FALSE)
}

c1 <- all(s1$estimate > 0)
c2 <- all(s1$p < 0.05)
# "Beats its null" = extreme in the PREDICTED direction, at the two-sided
# empirical p the build spec defines. percentile > 50 is the direction;
# p_emp < 0.05 is the test.
c4 <- all(null_res$p_emp < 0.05 & null_res$percentile > 50)

cond <- tibble::tibble(
  condition = c("1. interaction POSITIVE on both instruments",
                "2. p < 0.05 on both instruments",
                "3. both co-primary instruments fitted at S1 (STRUCTURAL)",
                "4. beats its expression-matched null on this endpoint"),
  evidential = c(TRUE, TRUE, FALSE, TRUE),
  met = c(c1, c2, c3, c4))
cond %>% as.data.frame() %>% print(row.names = FALSE)
message("   condition 3 is structural - it checks the fits exist, and is not ",
        "evidence.")

replicated <- all(cond$met)

message("\n   ", strrep("-", 70))
if (replicated) {
  message("   BCL2L11 REPLICATES in SCAN-B on both co-primary instruments.")
  message("   The one surviving finding of this arm now has an independent,")
  message("   pre-registered replication.")
} else {
  message("   BCL2L11 DOES NOT REPLICATE as declared.")
  message("   Block C section 9's failure condition applies verbatim:")
  message("     'it is reported as a TCGA-specific observation and dropped.")
  message("      No further variants are tried.'")
  message("   NOT another cohort. NOT another covariate set. NOT another")
  message("   estimator. NOT a subtype stratum. The declaration is honoured")
  message("   by reporting this, and the arm is written up as it stands.")
}
message("   ", strrep("-", 70))

# The two companions, reported with it as Block C reported them.
message("\n   the companions at S1 (interpreted as in Block C):")
fits %>% dplyr::filter(spec == SPEC_PRIMARY) %>%
  dplyr::select(endpoint, predicted, role, instrument, n, estimate,
                ci_lo, ci_hi, p) %>%
  dplyr::mutate(dplyr::across(c(estimate, ci_lo, ci_hi), ~ round(.x, 4)),
                p = signif(p, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 7. Save
# =============================================================================
message("\n7. save")

out <- list(
  gate = list(spine = spine, reduced = reduced, calibration = calib,
              reproduction = chk, passed = gate_sign && gate_ci,
              rule = paste("BCL2L11 keeps its positive sign AND the full-spine",
                           "estimate lies inside the reduced-set 95% CI, on",
                           "both instruments (declaration section 6)")),
  fits       = fits,
  null       = dplyr::select(null_res, -b_null),
  b_null     = stats::setNames(lapply(null_res$b_null, identity),
                               null_res$instrument),
  conditions = cond,
  replicated = replicated,
  scores = list(gsva = gsva_scores, mitopps_oxphos = ox_obs, zscale = ZS),
  coverage    = covg,
  symbol_report = bm$report,
  spec = list(
    declaration = "docs/2026-08-31_scanb_bim_replication_declaration.md",
    origin      = "docs/2026-08-29_block_c_result_H1_not_supported.md section 9",
    cohort      = "SCAN-B GSE202203, n = 3,207",
    model       = paste("log2(GENE) ~ M_a * OXPHOS subunits +",
                        "PROLIF_DISJOINT + PAM50"),
    coprimary   = paste("both instruments required; one alone is not a",
                        "replication"),
    meta        = "k = 1, no meta-analysis; single-cohort replication",
    endpoints   = "log2 of the LINEAR matrix, as in script 08",
    zscaling    = "within SCAN-B; not numerically comparable to TCGA",
    failure     = paste("if BCL2L11 does not replicate it is reported as a",
                        "TCGA-specific observation and dropped; no further",
                        "variants are tried"),
    fence       = paste("no survival, RFI, treatment, forkscale, PRIME, STATE",
                        "or BUFFER_c; those columns are not in the SCAN-B",
                        "phenotype at all")),
  built = Sys.time())

saveRDS(out, file.path(DIR_RESULTS, "scanb_bim_replication.rds"))
readr::write_csv(fits, file.path(DIR_TABLES, "scanb_bim_fits.csv"))
readr::write_csv(calib, file.path(DIR_TABLES, "scanb_gate_calibration.csv"))
readr::write_csv(dplyr::select(null_res, -b_null),
                 file.path(DIR_TABLES, "scanb_bim_null.csv"))

message("\n17: done.")
message("    results/scanb_bim_replication.rds")
message("    outputs/tables/  3 tables")
message("    NEXT: the result note. This arm has no further analysis.")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  r <- readRDS(file.path(DIR_RESULTS, "scanb_bim_replication.rds"))

  # --- the gate, which is the thing to read first --------------------------
  # If full_in_red_ci were FALSE anywhere for BCL2L11, the script would have
  # stopped. Look at how much the reduction moved the estimate anyway.
  r$gate$calibration %>% as.data.frame() %>% print(row.names = FALSE)

  # --- the primary, both instruments, side by side -------------------------
  # The rule is not "either". It is "both".
  r$fits %>% dplyr::filter(endpoint == "log2_BCL2L11", spec == "S1") %>%
    as.data.frame()

  # --- the four specifications, which are part of the result ---------------
  # If the interaction appears only under one, that is informative and is
  # REPORTED, never selected.
  r$fits %>% dplyr::filter(endpoint == "log2_BCL2L11") %>%
    dplyr::select(instrument, spec, n, estimate, ci_lo, ci_hi, p) %>%
    as.data.frame() %>% print(row.names = FALSE)

  # --- the null, drawn ------------------------------------------------------
  op <- par(mfrow = c(1, 2))
  for (ins in c("gsva", "mitopps")) {
    bn <- r$b_null[[ins]]
    hist(bn, breaks = 60, main = paste("null,", ins),
         xlab = "MYC:OX coefficient")
    abline(v = r$null$b_obs[r$null$instrument == ins], col = "red", lwd = 2)
  }
  par(op)

  # --- TCGA vs SCAN-B, directions only -------------------------------------
  # NOT a numerical comparison. GSVA is cohort-relative and mitoPPS is
  # composition-dependent; only the sign and the CI transfer.
  data.frame(
    endpoint   = r$fits$endpoint[r$fits$spec == "S1"],
    instrument = r$fits$instrument[r$fits$spec == "S1"],
    scanb      = round(r$fits$estimate[r$fits$spec == "S1"], 4),
    tcga_spine = round(r$gate$spine$estimate[match(
      paste(r$fits$endpoint[r$fits$spec == "S1"],
            r$fits$instrument[r$fits$spec == "S1"]),
      paste(r$gate$spine$endpoint, r$gate$spine$instrument))], 4))

  # --- the symbol map, again ------------------------------------------------
  table(r$symbol_report$status)

  # --- coverage -------------------------------------------------------------
  r$coverage %>% as.data.frame() %>% print(row.names = FALSE)

}
