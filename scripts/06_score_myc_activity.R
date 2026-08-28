# 06_score_myc_activity.R
# =============================================================================
# The three MYC estimators, and the concordance check between them.
#
#   M-a  FELSHER signature, MitoCarta-stripped, GSVA        <- PRIMARY (D2)
#   M-b  CollecTRI MYC regulon, stripped, decoupleR ULM     <- concordance
#   M-c  MYC (8q24.21) GISTIC amplification                 <- instrument
#
# =============================================================================
# SCALE DISCIPLINE
# =============================================================================
# GSVA and ULM both read the LOG-scale VST matrix from script 01
# (results/tcga_brca_vst.rds), with kcdf = "Gaussian" as CLAUDE.md requires.
# The linear matrix is for mitoPPS in script 07 and is NEVER read here. The
# object's own `scale` field is asserted below rather than trusted.
#
# COHORT-RELATIVITY: GSVA is cohort-relative. Every TCGA sample is scored in ONE
# run, which is what makes the scores internally comparable. They are NOT
# comparable to scores from METABRIC, SCAN-B or the neoadjuvant cohorts, which
# get their own runs. Meta-analyse effect estimates, never pooled scores.
#
# CIRCULARITY: the estimators used are the MitoCarta-STRIPPED versions from gate
# G1. That is the whole point of G1 - an unstripped MYC signature shares genes
# with the mitochondrial axis it is about to be related to. Raw versions are
# scored alongside as a sensitivity, never as the primary.
#
# SPECIES: human. See CLAUDE.md.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
suppressPackageStartupMessages({
  library(GSVA)
  library(decoupleR)
})

message("\n06: MYC activity estimators\n", strrep("=", 78))

# -----------------------------------------------------------------------------
# 1. Inputs
# -----------------------------------------------------------------------------
vst <- readRDS(file.path(DIR_RESULTS, "tcga_brca_vst.rds"))
if (!identical(vst$scale, "log_vst")) {
  stop("expected the LOG-scale VST object, got scale = '", vst$scale,
       "'. GSVA and ULM must not be given the linear matrix.", call. = FALSE)
}
E <- vst$mat
message("1. VST matrix: ", nrow(E), " genes x ", ncol(E), " samples (log scale)")

g1 <- readRDS(file.path(DIR_RESULTS, "g1_overlap_audit.rds"))
felsher_stripped <- g1$estimators_stripped$FELSHER
felsher_raw      <- g1$estimators_raw$FELSHER
collectri_stripped <- g1$estimators_stripped$COLLECTRI_MYC_ALL
collectri_raw      <- g1$estimators_raw$COLLECTRI_MYC_ALL
stopifnot(length(felsher_stripped) == 61L, length(felsher_raw) == 67L,
          length(collectri_stripped) == 811L, length(collectri_raw) == 886L)
message("   G1 estimators: FELSHER ", length(felsher_stripped), "/",
        length(felsher_raw), "  CollecTRI ", length(collectri_stripped), "/",
        length(collectri_raw), "  (stripped/raw)")

cov <- readRDS(file.path(DIR_RESULTS, "tcga_brca_covariates.rds"))$covariates
stopifnot(identical(sort(cov$patient), sort(colnames(E))))

# -----------------------------------------------------------------------------
# 2. M-a: Felsher signature, GSVA
# -----------------------------------------------------------------------------
# D2 (G1 note, 2026-08-28): M-a stays primary, confirmed not flipped. On the
# dimension G1 tests it is the cleaner estimator - zero overlap with MitoCarta
# OXPHOS subunits, against CollecTRI's three.
message("\n2. M-a: Felsher GSVA")

.gsva_score <- function(sets, label) {
  present <- lapply(sets, function(s) intersect(s, rownames(E)))
  for (nm in names(present)) {
    message(sprintf("   %-22s %3d of %3d genes present in the matrix",
                    nm, length(present[[nm]]), length(sets[[nm]])))
  }
  par <- GSVA::gsvaParam(exprData = E, geneSets = present,
                         kcdf = "Gaussian", minSize = 10L, maxSize = Inf)
  s <- GSVA::gsva(par, verbose = FALSE)
  message("   ", label, " scored: ", nrow(s), " set(s) x ", ncol(s), " samples")
  s
}

gsva_myc <- .gsva_score(
  list(FELSHER_stripped = felsher_stripped, FELSHER_raw = felsher_raw),
  "M-a")

M_a          <- as.numeric(gsva_myc["FELSHER_stripped", ])
M_a_raw      <- as.numeric(gsva_myc["FELSHER_raw", ])
names(M_a)     <- colnames(gsva_myc)
names(M_a_raw) <- colnames(gsva_myc)

# -----------------------------------------------------------------------------
# 3. M-b: CollecTRI MYC regulon, decoupleR ULM
# -----------------------------------------------------------------------------
# The network is read from the dated snapshot, never fetched live - see
# data/collectri_human/README.md for why (decoupleR::get_collectri() fails in
# this environment, and a live call is not reproducible anyway).
#
# SIGN RULE, fixed here: mor = +1 for a stimulatory edge, -1 for an inhibitory
# one. 92 of MYC's 891 edges are flagged BOTH stimulatory and inhibitory; they
# take +1, which is the same treatment G1 gave them when it built
# COLLECTRI_MYC_STIM. Recorded because a signed method's answer depends on it.
message("\n3. M-b: CollecTRI MYC regulon, ULM")

ct <- readr::read_tsv(PATH_COLLECTRI, show_col_types = FALSE, progress = FALSE)
myc_net_all <- ct %>%
  dplyr::filter(source_genesymbol == "MYC",
                !is.na(target_genesymbol), target_genesymbol != "") %>%
  dplyr::transmute(
    source = "MYC",
    target = target_genesymbol,
    mor    = dplyr::if_else(as.logical(is_stimulation), 1, -1),
    likelihood = 1
  ) %>%
  dplyr::distinct(source, target, .keep_all = TRUE)
message("   MYC edges in snapshot: ", nrow(myc_net_all),
        " (", sum(myc_net_all$mor > 0), " activating, ",
        sum(myc_net_all$mor < 0), " repressing)")

.ulm_score <- function(targets, label) {
  net <- myc_net_all %>%
    dplyr::filter(target %in% targets, target %in% rownames(E))
  message(sprintf("   %-22s %3d targets in the matrix", label, nrow(net)))
  res <- decoupleR::run_ulm(mat = E, network = net, .source = source,
                            .target = target, .mor = mor, minsize = 5L)
  v <- res %>%
    dplyr::filter(statistic == "ulm", source == "MYC") %>%
    dplyr::select(condition, score)
  stats::setNames(v$score, v$condition)[colnames(E)]
}

M_b     <- .ulm_score(collectri_stripped, "M-b stripped")
M_b_raw <- .ulm_score(collectri_raw,      "M-b raw")

# -----------------------------------------------------------------------------
# 4. M-c: MYC copy number
# -----------------------------------------------------------------------------
# The instrument. Not an expression measure at all, which is what makes it worth
# having: it cannot share genes with the mitochondrial axis. Taken from the ISAR
# GISTIC calls in the covariate table, at G2's primary threshold for MYC (+2).
message("\n4. M-c: MYC 8q24.21 GISTIC")
idx      <- match(colnames(E), cov$patient)
M_c_call <- cov$cnv_MYC[idx]          # -2..2, ISAR
M_c_amp  <- cov$MYC_amp[idx]          # == +2, G2's primary threshold
message("   GISTIC available for ", sum(!is.na(M_c_call)), " of ", ncol(E),
        " samples; amplified (+2): ", sum(M_c_amp, na.rm = TRUE))

# -----------------------------------------------------------------------------
# 5. Concordance
# -----------------------------------------------------------------------------
# Plan section 7.1: "if M-a and M-b disagree, report both and treat the claim as
# unsupported." That rule needs a number to be operational.
#
# >>> PRE-REGISTRATION DECISION, NEEDS SIGN-OFF.
# Proposed: M-a and M-b AGREE if their Spearman correlation across the cohort is
# positive with rho >= 0.30. Rationale: they are built from disjoint evidence -
# a 61-gene expression signature against an 811-target regulon scored by a signed
# linear model - so identical scores are not expected and would be suspicious.
# What is required is that they order the cohort similarly. 0.30 is a moderate
# bar chosen to be exceeded by genuine shared signal and failed by noise; it is
# NOT a significance threshold, since at n ~1095 almost any rho is "significant".
M_AB_MIN_RHO <- 0.30

est <- tibble::tibble(
  patient  = colnames(E),
  M_a      = M_a,      M_a_raw = M_a_raw,
  M_b      = M_b,      M_b_raw = M_b_raw,
  M_c_call = M_c_call, M_c_amp = M_c_amp
)

.rho <- function(x, y) stats::cor(x, y, method = "spearman", use = "complete.obs")
conc <- tibble::tibble(
  pair = c("M_a vs M_b", "M_a vs M_c", "M_b vs M_c",
           "M_a stripped vs raw", "M_b stripped vs raw"),
  rho  = c(.rho(est$M_a, est$M_b), .rho(est$M_a, est$M_c_call),
           .rho(est$M_b, est$M_c_call), .rho(est$M_a, est$M_a_raw),
           .rho(est$M_b, est$M_b_raw))
)
message("\n5. concordance (Spearman)")
conc %>% as.data.frame() %>% print(row.names = FALSE, digits = 3)

ab <- conc$rho[conc$pair == "M_a vs M_b"]
message(sprintf("\n   M-a vs M-b: rho = %.3f  ->  %s (threshold %.2f)",
                ab, if (ab >= M_AB_MIN_RHO) "AGREE" else "DISAGREE",
                M_AB_MIN_RHO))
if (ab < M_AB_MIN_RHO) {
  message("   Per plan section 7.1, a MYC-axis claim built on M-a alone is ",
          "NOT supported while this holds. Report both estimators.")
}

# The strip must not have changed what the estimator measures. If it did, G1's
# conclusion that the stripped signature retains its structure was wrong.
for (p in c("M_a stripped vs raw", "M_b stripped vs raw")) {
  r <- conc$rho[conc$pair == p]
  if (r < 0.90) {
    warning(p, ": rho = ", round(r, 3), ", below 0.90. MitoCarta stripping has ",
            "materially changed the estimator, which G1 did not anticipate.",
            call. = FALSE)
  }
}

# -----------------------------------------------------------------------------
# 6. Save
# -----------------------------------------------------------------------------
myc <- list(
  estimators  = est,
  concordance = conc,
  gsva_raw    = gsva_myc,
  network     = myc_net_all,
  rules       = list(primary = "M_a (FELSHER stripped, GSVA)",
                     m_ab_min_rho = M_AB_MIN_RHO,
                     mor_ambiguous = "both-flagged edges take mor = +1",
                     scale = "log VST, kcdf Gaussian, cohort-relative"),
  built       = Sys.time()
)
saveRDS(myc, file.path(DIR_RESULTS, "tcga_brca_myc_scores.rds"))
readr::write_csv(est,  file.path(DIR_TABLES, "tcga_brca_myc_estimators.csv"))
readr::write_csv(conc, file.path(DIR_TABLES, "tcga_brca_myc_concordance.csv"))

message("\n06: done. results/tcga_brca_myc_scores.rds + 2 tables")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  m <- readRDS(file.path(DIR_RESULTS, "tcga_brca_myc_scores.rds"))
  e <- m$estimators

  # --- the three estimators against each other -----------------------------
  pairs(e[, c("M_a", "M_b", "M_c_call")], pch = 16, cex = 0.3,
        main = "M-a (Felsher GSVA) / M-b (CollecTRI ULM) / M-c (GISTIC)")

  # --- does M-c separate the expression estimators? ------------------------
  # The instrument is independent of expression, so this is the sharpest check
  # that the expression scores are measuring MYC at all.
  boxplot(M_a ~ M_c_call, data = e, xlab = "MYC GISTIC call", ylab = "M-a")
  boxplot(M_b ~ M_c_call, data = e, xlab = "MYC GISTIC call", ylab = "M-b")
  kruskal.test(M_a ~ factor(M_c_call), data = e)
  kruskal.test(M_b ~ factor(M_c_call), data = e)

  # --- by subtype: MYC activity should be highest in Basal -----------------
  cv <- readRDS(file.path(DIR_RESULTS, "tcga_brca_covariates.rds"))$covariates
  d  <- merge(e, cv[, c("patient", "PAM50", "TNBC")], by = "patient")
  boxplot(M_a ~ PAM50, data = d, las = 2, ylab = "M-a")
  tapply(d$M_a, d$PAM50, median, na.rm = TRUE)

  # --- did stripping change anything? --------------------------------------
  plot(e$M_a_raw, e$M_a, pch = 16, cex = 0.3,
       xlab = "Felsher raw (67)", ylab = "Felsher stripped (61)")
  abline(0, 1, col = "red")

  # --- the 92 ambiguous edges: how much do they matter? --------------------
  # Re-score M-b with the both-flagged edges set to -1 instead of +1 and see
  # whether the ordering survives. If it does not, the sign rule is load-bearing
  # and belongs in Methods rather than a comment.
  m$network |> dplyr::count(mor)

}
