# 14_depmap_dependency.R
# =============================================================================
# BLOCK G - functional dependency in DepMap / CCLE.  -> Panel c inset
#
# Plan section 10 (Block G), section 12 ("can run in parallel from the start -
# it needs none of the TCGA work and it is the cheapest evidence in the whole
# plan"), section 3 (dataset assignment).
#
# =============================================================================
# WHY THIS MATTERS MORE NOW THAN WHEN IT WAS SPECIFIED
# =============================================================================
# Three of the four falsification criteria are met. Everything left in the TCGA
# arm is observational and cross-sectional. This block is FUNCTIONAL: CRISPR
# knockout and drug response are interventions, so a result here is a different
# KIND of evidence, not more of the same.
#
# It also speaks directly to the one finding that survived both instruments:
#
#     MYC x OXPHOS -> LOWER BCL-XL, HIGHER BIM, PUMA unmoved, no compensating
#     amplification.
#
# =============================================================================
# THE PREDICTION IS NOW TWO-SIDED, AND THAT IS THE POINT
# =============================================================================
# Plan section 10 predicts: "MYC-high / OXPHOS-high lines are selectively
# dependent on MCL1 and/or BCL2L1", i.e. a NEGATIVE `MYC:OX` coefficient on
# Chronos gene effect (more negative = more essential).
#
# Block C found the opposite thing about SUPPLY: `MYC x OXPHOS` associates with
# LESS BCL2L1 transcript. Dependency and transcript are not the same quantity -
# a cell can need a protein more precisely because it has less of it - so both
# outcomes are live and neither is a rescue of the other.
#
# THE PRE-REGISTERED PREDICTION IS STILL THE PLAN'S ONE: negative `MYC:OX` on
# MCL1 and BCL2L1 gene effect. That is what section 7 evaluates. The reverse is
# reported when it occurs and is NOT scored as a pass.
#
# =============================================================================
# SCALE DISCIPLINE - READ THIS BEFORE EDITING ANY SCORING BLOCK
# =============================================================================
# DepMap `OmicsExpressionProteinCodingGenesTPMLogp1.csv` is log2(TPM + 1).
#
#   GSVA      -> the file AS SUPPLIED. Log scale, kcdf = "Gaussian". Same as
#                scripts 06-07 on VST.
#   mitoPPS   -> LINEAR, via 2^x - 1. It cannot see the log matrix.
#
# The two must never share an object. `E_LOG` and `E_LIN` below are separate and
# are never rebuilt from one another after this point.
#
# DECLARED DEVIATION: the plan's mitoPPS input is linear DESeq2-normalised
# counts; here it is linear TPM, because that is what DepMap ships. TPM is
# additionally length-normalised. mitoPPS is a COMPOSITION measure built from
# all-pairwise pathway ratios and is deliberately robust to total content, so
# the deviation is defensible - but it is a deviation and is recorded here and
# in the saved object. It also makes the standing CLAUDE.md rule bite harder,
# not less: CCLE mitoPPS values are NEVER comparable to TCGA mitoPPS values.
# Only the PATTERN transfers.
#
# GSVA IS COHORT-RELATIVE. Every line in this script is scored in ONE run.
# CCLE scores and TCGA scores are different scales and must never be pooled.
#
# SPECIES: human. CCLE/DepMap human cell lines, human MitoCarta arms read from
# script 07. See CLAUDE.md.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
suppressPackageStartupMessages(library(data.table))

message("\n14: Block G - DepMap functional dependency\n", strrep("=", 78))

# -----------------------------------------------------------------------------
# 0. Constants
# -----------------------------------------------------------------------------
DIR_DEPMAP <- file.path(DIR_DATA, "raw", "depmap")

# Pinned release. `Model.csv`, the expression matrix and `CRISPRGeneEffect.csv`
# must come from ONE release: Model IDs are stable but the line set and the
# Chronos scaling are not, and mixing releases produces a join that looks
# complete and is not.
#
# PRISM IS THE ONE DELIBERATE EXCEPTION. Repurposing has its own release
# cadence and `Repurposing Public 24Q2` is the current one - there is no 24Q4
# Repurposing. It joins on ModelID, which is stable, so this is a different
# dataset rather than a stale copy of this one. Recorded, not silently mixed.
DEPMAP_RELEASE <- "24Q4"
PRISM_RELEASE  <- "24Q2"

PATH_MODEL  <- file.path(DIR_DEPMAP, "Model.csv")
PATH_EXPR   <- file.path(DIR_DEPMAP, "OmicsExpressionProteinCodingGenesTPMLogp1.csv")
PATH_CRISPR <- file.path(DIR_DEPMAP, "CRISPRGeneEffect.csv")
# Optional. Block G item 3. Absent -> section 6 is skipped, not faked.
PATH_PRISM  <- file.path(DIR_DEPMAP, "Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv")
PATH_PRISM_TREAT <- file.path(DIR_DEPMAP, "Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv")

LINEAGE <- "Breast"

# Block G item 2. The plan's endpoints, and the controls that make a positive
# result mean something.
DEP_PRIMARY  <- c("MCL1", "BCL2L1")
DEP_CONTROL  <- c("BCL2",      # the specificity control, as in G2 and script 08
                  "BBC3", "BCL2L11", "BAX", "BAK1",
                  "RPL3", "POLR2A")   # common essentials: a floor, not a result

# Block G item 3. Named in plan section 10, plus AZD5991 as a third MCL1 agent.
#
# COVERAGE, checked against Repurposing Public 24Q2 on 2026-08-29: 5 of these 7
# are present. **A-1331852 and A-1155463 are NOT in PRISM Repurposing** - and
# those are the two SELECTIVE BCL-XL inhibitors. What is left on the BCL-XL side
# is navitoclax, which inhibits BCL2, BCL-XL and BCL-W and therefore CANNOT
# separate BCL-XL from BCL2 on its own.
#
# That is a real limitation of item 3, not a detail. It is partly recoverable:
# venetoclax IS present and is BCL2-selective, so navitoclax and venetoclax read
# together bound the BCL-XL/BCL-W component - which is exactly the role plan
# section 10 already assigns venetoclax. Read them as a pair; do not read a
# navitoclax result as a BCL-XL result.
#
# MCL1 is well covered: S63845, AMG-176 and AZD5991 are all present.
DRUGS <- c("S63845", "AMG-176", "AZD5991",          # MCL1, all 3 present
           "A-1331852", "navitoclax", "A-1155463",  # BCL-XL; only navitoclax
           "venetoclax")                            # BCL2 specificity control

ARM_PRIMARY  <- "OXPHOS subunits"
ARM_NEGATIVE <- "OXPHOS assembly factors"
INSTRUMENTS  <- c("gsva", "mitopps")
CI_LEVEL     <- 0.95
MIN_LINES    <- 25L   # below this an interaction is not worth fitting at all
GSVA_MIN_SET <- 5L    # a GSVA score on fewer genes is noise with a name

# =============================================================================
# 1. Inputs
# =============================================================================
# Following script 02's idiom: files must already be in data/raw/, and a missing
# one stops here with the acquisition command rather than 200 lines later.
#
# NOTE ON ACQUISITION: depmap.org is behind a Cloudflare challenge, so a plain
# curl of the portal will return an HTML verification page, not data - and a
# 5 KB "CSV" that is actually HTML fails much later and confusingly. The size
# guard below catches exactly that. figshare works from the command line.
message("\n1. inputs")

# Verified 2026-08-29: figshare carries ONLY the Achilles/CRISPR half of every
# DepMap release - the omics matrices and Model.csv are portal-only, and the
# portal is behind a Cloudflare challenge. So two files need a browser and
# three do not. PRISM Repurposing, by contrast, is complete on figshare.
.acquire_help <- paste0(
  "\n\nAcquisition, into ", DIR_DEPMAP, ":\n\n",
  "  BROWSER (2 files - portal only, Cloudflare-gated):\n",
  "    https://depmap.org/portal/data_page/?tab=allData   release ",
  DEPMAP_RELEASE, "\n",
  "      Model.csv\n",
  "      OmicsExpressionProteinCodingGenesTPMLogp1.csv\n\n",
  "  SCRIPTABLE (3 files - figshare, direct):\n",
  "    curl -L -o CRISPRGeneEffect.csv \\\n",
  "      https://ndownloader.figshare.com/files/51064667      # DepMap ",
  DEPMAP_RELEASE, "\n",
  "    curl -L -o Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv \\\n",
  "      https://ndownloader.figshare.com/files/46630984      # PRISM ",
  PRISM_RELEASE, "\n",
  "    curl -L -o Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv \\\n",
  "      https://ndownloader.figshare.com/files/46630981      # PRISM ",
  PRISM_RELEASE, "\n\n",
  "Checksums and the full story are in data/depmap/README.md.")

for (p in c(PATH_MODEL, PATH_EXPR, PATH_CRISPR)) {
  if (!file.exists(p)) stop("missing input: ", p, .acquire_help, call. = FALSE)
  # An HTML challenge page is ~5 KB. Every real file here is far larger.
  if (file.size(p) < 1e5) {
    stop(basename(p), " is only ", file.size(p), " bytes. That is not the data ",
         "file - it is almost certainly the Cloudflare verification HTML page ",
         "saved under a .csv name. Re-download it.", .acquire_help,
         call. = FALSE)
  }
}
message("   3 required inputs present")

have_prism <- file.exists(PATH_PRISM) && file.exists(PATH_PRISM_TREAT)
message("   PRISM drug sensitivity: ", if (have_prism) "present" else
        "ABSENT - Block G item 3 will be skipped, not faked")

# --- the gene sets, read from where this repo already built them -------------
# NOT rebuilt. The 18 MitoCarta arms are script 07's; the MYC signature is the
# MitoCarta-stripped Felsher set that G1 produced and script 06 uses as M-a.
# Rebuilding either here would give this block a quietly different exposure from
# every TCGA result it is meant to be read beside.
mito <- readRDS(file.path(DIR_RESULTS, "tcga_brca_mito_scores.rds"))
g1   <- readRDS(file.path(DIR_RESULTS, "g1_overlap_audit.rds"))

ARM_SETS <- mito$arm_sets
MYC_SET  <- g1$estimators_stripped$FELSHER
PROLIF   <- mito$covariate_sets$PROLIF_DISJOINT
stopifnot(length(ARM_SETS) == 18L, length(MYC_SET) == 61L, length(PROLIF) > 100L)
message("   gene sets: 18 MitoCarta arms, FELSHER stripped (", length(MYC_SET),
        "), PROLIF_DISJOINT (", length(PROLIF), ")")

# =============================================================================
# 2. Cell lines
# =============================================================================
message("\n2. cell lines")

MODEL <- data.table::fread(PATH_MODEL, data.table = FALSE)
stopifnot("ModelID" %in% names(MODEL))
lin_col <- if ("OncotreeLineage" %in% names(MODEL)) "OncotreeLineage" else
  if ("lineage" %in% names(MODEL)) "lineage" else
    stop("Model.csv has no lineage column (looked for OncotreeLineage, ",
         "lineage). DepMap renamed it; fix `lin_col` rather than guessing.",
         call. = FALSE)
sub_col <- if ("OncotreeSubtype" %in% names(MODEL)) "OncotreeSubtype" else NA

MODEL$lineage <- MODEL[[lin_col]]
breast <- MODEL$ModelID[!is.na(MODEL$lineage) & MODEL$lineage == LINEAGE]
message("   ", length(breast), " ", LINEAGE, " models in Model.csv (of ",
        nrow(MODEL), ")")

# =============================================================================
# 3. Expression, on two scales that never meet
# =============================================================================
message("\n3. expression")

.read_matrix <- function(path, what) {
  m <- data.table::fread(path, data.table = FALSE)
  rn <- m[[1]]; m <- as.matrix(m[, -1, drop = FALSE])
  rownames(m) <- rn
  # DepMap column names are "SYMBOL (ENTREZ)". Take the symbol.
  colnames(m) <- trimws(sub("\\s*\\(\\d+\\)$", "", colnames(m)))
  dup <- duplicated(colnames(m))
  if (any(dup)) {
    message("   ", what, ": ", sum(dup), " duplicate symbols after stripping ",
            "Entrez ids - first occurrence kept")
    m <- m[, !dup, drop = FALSE]
  }
  m
}

EXPR <- .read_matrix(PATH_EXPR, "expression")   # models x genes, log2(TPM+1)
message("   expression: ", nrow(EXPR), " models x ", ncol(EXPR), " genes")

# A sanity check on the scale, because everything below depends on it and the
# failure is silent: log2(TPM+1) is bounded well under 20 and mostly near zero.
if (max(EXPR, na.rm = TRUE) > 40) {
  stop("expression matrix max is ", round(max(EXPR, na.rm = TRUE), 1),
       ", which is not log2(TPM+1). GSVA would be scored on the wrong scale ",
       "and mitoPPS on a doubly-exponentiated one. Check the file.",
       call. = FALSE)
}

lines <- intersect(breast, rownames(EXPR))
if (length(lines) < MIN_LINES) {
  stop("only ", length(lines), " ", LINEAGE, " lines have expression; below ",
       MIN_LINES, " an interaction model is not worth fitting.", call. = FALSE)
}
message("   ", length(lines), " ", LINEAGE, " lines with expression")

# LOG scale. GSVA only.
E_LOG <- t(EXPR[lines, , drop = FALSE])                 # genes x models
E_LOG <- E_LOG[stats::complete.cases(E_LOG), , drop = FALSE]

# LINEAR scale. mitoPPS only. Built ONCE, from E_LOG, and never mixed back.
#
# NO FLOOR IS ADDED. mitoPPS forms all-pairwise pathway ratios, so a pathway
# mean of exactly zero is a division by zero - and script 07 STOPS on that
# rather than patching it (line 615 there). The same rule holds here. A silent
# floor would be an undocumented deviation that shifts every ratio in the
# universe. The risk is real in CCLE, where small arms (CII subunits is 4
# genes) can be entirely unexpressed in a line, so the guard below names the
# offending arm and line rather than just failing.
E_LIN <- 2^E_LOG - 1
E_LIN[E_LIN < 0] <- 0                   # floating point on an exact zero only

message("   E_LOG ", nrow(E_LOG), " x ", ncol(E_LOG),
        " (log2 TPM+1, GSVA)   |   E_LIN same dims (linear TPM, mitoPPS)")

# =============================================================================
# 4. Scoring - one run over one cohort
# =============================================================================
message("\n4. scoring")

.present <- function(s) intersect(s, rownames(E_LOG))
sets_gsva <- c(list(MYC = MYC_SET, PROLIF = PROLIF), ARM_SETS)
sets_gsva <- lapply(sets_gsva, .present)
short <- names(sets_gsva)[lengths(sets_gsva) < 5L]
if (length(short)) {
  stop("gene set(s) with fewer than 5 genes present in CCLE: ",
       paste(short, collapse = ", "), ". Symbol harmonisation has failed; a ",
       "GSVA score on 2 genes is noise with a name.", call. = FALSE)
}
message("   set coverage in CCLE, min ", min(lengths(sets_gsva)), " genes")

# --- GSVA. LOG SCALE. -------------------------------------------------------
gp <- GSVA::gsvaParam(exprData = E_LOG, geneSets = sets_gsva,
                      kcdf = "Gaussian", minSize = GSVA_MIN_SET, maxSize = Inf)
GS <- GSVA::gsva(gp, verbose = FALSE)
# GSVA drops sets below minSize SILENTLY, returning a smaller matrix. A dropped
# arm would surface far downstream as a subscript error, so name it here.
missing_sets <- setdiff(names(sets_gsva), rownames(GS))
if (length(missing_sets)) {
  stop("GSVA dropped ", length(missing_sets), " set(s) below minSize: ",
       paste(missing_sets, collapse = ", "),
       ". They are silently absent from the score matrix.", call. = FALSE)
}
message("   GSVA: ", nrow(GS), " sets x ", ncol(GS), " lines")

# --- mitoPPS. LINEAR SCALE. --------------------------------------------------
# Same algorithm as script 07: pathway means, all pairwise ratios, each ratio
# corrected by its across-sample mean, mitoPPS = mean of corrected ratios. The
# universe is the 18 arms, exactly as in script 07 - mitoPPS is composition
# WITHIN a declared universe, so changing the universe changes every value.
.pathway_means <- function(E, sets) {
  t(vapply(sets, function(g) colMeans(E[g, , drop = FALSE]),
           numeric(ncol(E))))
}
.mitopps_universe <- function(S) {
  N <- ncol(S); P <- nrow(S)
  Bi <- 1 / S
  A  <- (S %*% t(Bi)) / N                # A[p,q] = mean_i S[p,i] / S[q,i]
  M  <- (1 / A) %*% Bi
  out <- S * (M - Bi) / (P - 1)          # A[p,p] == 1, so subtracting Bi drops q = p
  dimnames(out) <- dimnames(S)
  out
}
S_LIN <- .pathway_means(E_LIN, lapply(ARM_SETS, .present))
.check_positive <- function(S, what) {
  bad <- which(!is.finite(S) | S <= 0, arr.ind = TRUE)
  if (nrow(bad)) {
    lab <- utils::head(sprintf("%s / %s", rownames(S)[bad[, 1]],
                               colnames(S)[bad[, 2]]), 6L)
    stop("mitoPPS (", what, "): ", nrow(bad), " arm x line pathway mean(s) are ",
         "zero or non-finite, so the pairwise ratio is undefined. First few: ",
         paste(lab, collapse = "; "),
         ".\nDecide explicitly - drop the arm from the universe, or drop the ",
         "line - and record it. Do NOT add a floor silently: mitoPPS is a ",
         "composition measure and a floor moves every ratio in the universe.",
         call. = FALSE)
  }
}
.check_positive(S_LIN, LINEAGE)
PPS <- .mitopps_universe(S_LIN)
message("   mitoPPS: ", nrow(PPS), " arms x ", ncol(PPS), " lines")

.z <- function(v) (v - mean(v, na.rm = TRUE)) / stats::sd(v, na.rm = TRUE)

SC <- list(gsva = t(apply(GS[names(ARM_SETS), , drop = FALSE], 1, .z)),
           mitopps = t(apply(PPS, 1, .z)))
MYC_Z    <- .z(GS["MYC", ])
PROLIF_Z <- .z(GS["PROLIF", ])
stopifnot(identical(colnames(SC$gsva), colnames(SC$mitopps)),
          identical(colnames(SC$gsva), lines))

# =============================================================================
# 5. Block G item 2 - CRISPR dependency
# =============================================================================
# Chronos gene effect: 0 = no effect, -1 = median common essential. MORE
# NEGATIVE = MORE ESSENTIAL. So the plan's "selectively dependent" prediction is
# a NEGATIVE `MYC:OX` interaction.
message("\n5. CRISPR gene effect")

CRISPR <- .read_matrix(PATH_CRISPR, "CRISPR")
message("   CRISPR: ", nrow(CRISPR), " models x ", ncol(CRISPR), " genes")

dep_lines <- intersect(lines, rownames(CRISPR))
message("   ", length(dep_lines), " ", LINEAGE,
        " lines have BOTH expression and CRISPR")
if (length(dep_lines) < MIN_LINES) {
  stop("only ", length(dep_lines), " lines carry both; below ", MIN_LINES,
       " this block cannot be run as specified. Report the gap; do not ",
       "widen the lineage to manufacture n.", call. = FALSE)
}

GENES_DEP <- c(DEP_PRIMARY, DEP_CONTROL)
missing_dep <- setdiff(GENES_DEP, colnames(CRISPR))
if (length(missing_dep)) {
  message("   not screened in this release: ", paste(missing_dep, collapse = ", "))
}
GENES_DEP <- intersect(GENES_DEP, colnames(CRISPR))
stopifnot(all(DEP_PRIMARY %in% GENES_DEP))

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

.frame <- function(rows, ins, arm) {
  tibble::tibble(
    ModelID = rows,
    MYC     = MYC_Z[rows],
    OX      = SC[[ins]][arm, rows],
    PROLIF  = PROLIF_Z[rows])
}

RES <- list(); .add <- function(x) RES[[length(RES) + 1L]] <<- x

# --- 5a. breast, as the plan specifies. PRIMARY. -----------------------------
# n is small. The model is deliberately minimal: adding subtype at this n would
# spend the degrees of freedom that the interaction needs.
for (ins in INSTRUMENTS) {
  for (arm in c(ARM_PRIMARY, ARM_NEGATIVE)) {
    for (g in GENES_DEP) {
      d <- .frame(dep_lines, ins, arm)
      d$Y <- CRISPR[dep_lines, g]
      if (sum(stats::complete.cases(d)) < MIN_LINES) next
      m <- stats::lm(Y ~ MYC * OX + PROLIF, data = d)
      .add(.tidy(m, "MYC:OX", "G2 breast", ins,
                 extra = list(gene = g, arm = arm, scope = LINEAGE)))
    }
  }
}

# --- 5b. the full arm panel, primary genes only. SPECIFICITY. ----------------
# All 18 arms, so a positive OXPHOS result can be read against the others rather
# than on its own. There is NO expression-matched null in CCLE - script 07's
# nulls are TCGA-specific and are not transferable - so this is a RANK panel,
# not a calibrated p-value. If OXPHOS subunits is positive here, the matched
# null has to be BUILT in CCLE before the result is reportable. That is a named
# follow-up, not an optional extra.
for (ins in INSTRUMENTS) {
  for (arm in names(ARM_SETS)) {
    for (g in DEP_PRIMARY) {
      d <- .frame(dep_lines, ins, arm)
      d$Y <- CRISPR[dep_lines, g]
      m <- stats::lm(Y ~ MYC * OX + PROLIF, data = d)
      .add(.tidy(m, "MYC:OX", "G2 arm panel", ins,
                 extra = list(gene = g, arm = arm, scope = LINEAGE)))
    }
  }
}

# --- 5c. pan-cancer, lineage-adjusted. SECONDARY, and better powered. --------
# ~50 breast lines cannot support a two-way interaction with any confidence.
# The pan-cancer fit has 20x the n. It is SECONDARY because the plan specifies
# breast and because the MYC/OXPHOS relationship need not be lineage-invariant;
# it is reported because a breast-only null at n = 50 is uninformative and
# saying so without offering the powered comparison would be a half-report.
#
# It requires re-scoring: GSVA is cohort-relative, so a pan-cancer score is a
# different quantity from a breast-only one and the two cannot be mixed.
pan_lines <- intersect(rownames(EXPR), rownames(CRISPR))
message("   pan-cancer: ", length(pan_lines), " lines with both")
message("   re-scoring GSVA and mitoPPS over all lineages - this is the slow ",
        "step of
   this script (minutes, not seconds). It cannot be skipped ",
        "by reusing the
   breast scores: GSVA is cohort-relative.")

E_LOG_PAN <- t(EXPR[pan_lines, , drop = FALSE])
E_LOG_PAN <- E_LOG_PAN[stats::complete.cases(E_LOG_PAN), , drop = FALSE]
# From the ORIGINAL sets, not from `sets_gsva` - that one is already intersected
# with the breast matrix, and intersecting a subset would quietly shrink the
# pan-cancer sets to whatever survived in ~50 lines.
sets_pan <- lapply(c(list(MYC = MYC_SET, PROLIF = PROLIF), ARM_SETS),
                   function(s) intersect(s, rownames(E_LOG_PAN)))
GS_PAN <- GSVA::gsva(GSVA::gsvaParam(exprData = E_LOG_PAN, geneSets = sets_pan,
                                     kcdf = "Gaussian", minSize = GSVA_MIN_SET,
                                     maxSize = Inf), verbose = FALSE)

E_LIN_PAN <- 2^E_LOG_PAN - 1
E_LIN_PAN[E_LIN_PAN < 0] <- 0
S_LIN_PAN <- .pathway_means(E_LIN_PAN, lapply(ARM_SETS, function(s)
  intersect(s, rownames(E_LOG_PAN))))
.check_positive(S_LIN_PAN, "pan-cancer")
PPS_PAN <- .mitopps_universe(S_LIN_PAN)

SC_PAN <- list(gsva = t(apply(GS_PAN[names(ARM_SETS), , drop = FALSE], 1, .z)),
               mitopps = t(apply(PPS_PAN, 1, .z)))
MYC_PAN <- .z(GS_PAN["MYC", ]); PROLIF_PAN <- .z(GS_PAN["PROLIF", ])
lin_pan <- MODEL$lineage[match(pan_lines, MODEL$ModelID)]
keep_lin <- lin_pan %in% names(which(table(lin_pan) >= 10L))

for (ins in INSTRUMENTS) {
  for (arm in c(ARM_PRIMARY, ARM_NEGATIVE)) {
    for (g in DEP_PRIMARY) {
      d <- tibble::tibble(MYC = MYC_PAN[pan_lines], OX = SC_PAN[[ins]][arm, pan_lines],
                          PROLIF = PROLIF_PAN[pan_lines],
                          lineage = factor(lin_pan), Y = CRISPR[pan_lines, g])[keep_lin, ]
      d <- droplevels(d)
      m <- stats::lm(Y ~ MYC * OX + PROLIF + lineage, data = d)
      .add(.tidy(m, "MYC:OX", "G2 pan-cancer (secondary)", ins,
                 extra = list(gene = g, arm = arm, scope = "pan-cancer")))
    }
  }
}

# =============================================================================
# 6. Block G item 3 - drug sensitivity
# =============================================================================
# PRISM Repurposing log2 fold change: MORE NEGATIVE = MORE SENSITIVE, the same
# direction as Chronos, so the prediction is again a NEGATIVE `MYC:OX`.
message("\n6. drug sensitivity")

drug_res <- NULL
if (!have_prism) {
  message("   skipped - PRISM files absent. Block G item 3 is DEFERRED, not ",
          "dropped;\n   see .acquire_help at the top of this script.")
} else {
  TREAT <- data.table::fread(PATH_PRISM_TREAT, data.table = FALSE)
  nm_col <- intersect(c("Drug.Name", "name", "Name"), names(TREAT))[1]
  id_col <- intersect(c("IDs", "column_name", "broad_id"), names(TREAT))[1]
  if (is.na(nm_col) || is.na(id_col)) {
    message("   PRISM compound list has unexpected columns (",
            paste(utils::head(names(TREAT)), collapse = ", "),
            ") - item 3 skipped rather than guessed")
  } else {
    # ORIENTATION: the Repurposing matrix is COMPOUNDS x CELL LINES - rows are
    # `BRD:BRD-...` ids and columns are ModelIDs - which is the transpose of the
    # expression and CRISPR matrices. Verified 2026-08-29: 6,790 x 919. Read
    # without transposing, `intersect(lines, rownames(PR))` is empty and item 3
    # silently reports "0 lines with PRISM" rather than failing. Hence the
    # assertion, not just the t().
    PR <- data.table::fread(PATH_PRISM, data.table = FALSE)
    rn <- PR[[1]]; PR <- as.matrix(PR[, -1, drop = FALSE]); rownames(PR) <- rn
    if (!all(grepl("^BRD", utils::head(rownames(PR), 20))) ||
        !all(grepl("^ACH-", utils::head(colnames(PR), 20)))) {
      stop("PRISM matrix orientation is not compounds x ModelIDs as expected ",
           "(rows start '", substr(rownames(PR)[1], 1, 12), "', columns start '",
           substr(colnames(PR)[1], 1, 12), "'). Check the release before ",
           "transposing.", call. = FALSE)
    }
    PR <- t(PR)                                      # now ModelIDs x compounds
    # Drug.Name is UPPERCASE for some compounds and mixed for others, hence the
    # case-insensitive match on both sides.
    hit <- TREAT[tolower(TREAT[[nm_col]]) %in% tolower(DRUGS), , drop = FALSE]
    absent <- setdiff(tolower(DRUGS), tolower(hit[[nm_col]]))
    message("   ", nrow(hit), " of ", length(DRUGS),
            " named compounds found in PRISM")
    if (length(absent)) {
      message("   NOT in this PRISM release: ", paste(absent, collapse = ", "))
    }
    pr_lines <- intersect(lines, rownames(PR))
    message("   ", length(pr_lines), " ", LINEAGE, " lines with PRISM")
    if (nrow(hit) && length(pr_lines) >= MIN_LINES) {
      drug_res <- dplyr::bind_rows(lapply(seq_len(nrow(hit)), function(k) {
        col <- as.character(hit[[id_col]][k])
        if (!col %in% colnames(PR)) return(NULL)
        dplyr::bind_rows(lapply(INSTRUMENTS, function(ins) {
          d <- .frame(pr_lines, ins, ARM_PRIMARY)
          d$Y <- PR[pr_lines, col]
          if (sum(stats::complete.cases(d)) < MIN_LINES) return(NULL)
          m <- stats::lm(Y ~ MYC * OX + PROLIF, data = d)
          .tidy(m, "MYC:OX", "G3 drug", ins,
                extra = list(gene = as.character(hit[[nm_col]][k]),
                             arm = ARM_PRIMARY, scope = LINEAGE))
        }))
      }))
    }
  }
}

coefs <- dplyr::bind_rows(c(RES, list(drug_res)))

# =============================================================================
# 7. The declared predictions, evaluated
# =============================================================================
# Plan section 10's prediction, and nothing added to it. Negative `MYC:OX` on
# Chronos = MYC-high / OXPHOS-high lines are MORE dependent.
message("\n7. the declared predictions")

.both <- function(v) length(v) == 2L && all(v)
# Plain subsetting rather than dplyr::filter(): a filter whose predicate names
# `gene` and `arm` has to unquote the arguments to avoid matching the columns
# against themselves, and that is a trap worth not setting.
.g <- function(lab, gene_, arm_ = ARM_PRIMARY) {
  coefs[coefs$label == lab & coefs$gene == gene_ & coefs$arm == arm_, ]
}

g_mcl1 <- .g("G2 breast", "MCL1")
g_bclx <- .g("G2 breast", "BCL2L1")
g_neg  <- dplyr::bind_rows(.g("G2 breast", "MCL1", ARM_NEGATIVE),
                           .g("G2 breast", "BCL2L1", ARM_NEGATIVE))
g_bcl2 <- .g("G2 breast", "BCL2")

p1 <- .both(g_mcl1$estimate < 0 & g_mcl1$p < 0.05)
p2 <- .both(g_bclx$estimate < 0 & g_bclx$p < 0.05)
p3 <- !any(g_neg$estimate < 0 & g_neg$p < 0.05)
p4 <- !.both(g_bcl2$estimate < 0 & g_bcl2$p < 0.05)

declared <- tibble::tibble(
  test = c("G-a  MCL1 dependency: MYC:OX negative, both instruments",
           "G-b  BCL2L1 dependency: MYC:OX negative, both instruments",
           "G-c  and NOT the assembly factors (specificity)",
           "G-d  and NOT BCL2 (the endpoint specificity control)"),
  predicted = c("negative", "negative", "null / weaker", "null / weaker"),
  passes = c(p1, p2, p3, p4))
declared %>% as.data.frame() %>% print(row.names = FALSE)

message("\n   G2 breast, primary and control genes (", ARM_PRIMARY, "):")
coefs %>% dplyr::filter(label == "G2 breast", arm == ARM_PRIMARY) %>%
  dplyr::select(gene, instrument, n, estimate, ci_lo, ci_hi, p) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   G2 pan-cancer (secondary, better powered):")
coefs %>% dplyr::filter(label == "G2 pan-cancer (secondary)") %>%
  dplyr::select(gene, arm, instrument, n, estimate, ci_lo, ci_hi, p) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

if (!is.null(drug_res) && nrow(drug_res)) {
  message("\n   G3 drug sensitivity:")
  drug_res %>% dplyr::select(gene, instrument, n, estimate, ci_lo, ci_hi, p) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3))) %>%
    as.data.frame() %>% print(row.names = FALSE)
}

message("\n   POWER: ", length(dep_lines), " ", LINEAGE, " lines. A two-way ",
        "interaction at this n is\n   underpowered; a null in G-a or G-b is ",
        "NOT evidence against, and the plan's\n   'a negative result here ",
        "weakens H1 considerably' has to be read with that\n   in view. The ",
        "pan-cancer fit is the powered comparison.")
message("   NO MATCHED NULL EXISTS IN CCLE. The arm panel is a rank ordering, ",
        "not a\n   calibrated p-value. A positive OXPHOS result needs a CCLE ",
        "expression-matched\n   null built before it is reportable.")

# =============================================================================
# 8. Save
# =============================================================================
message("\n8. save")

out <- list(
  coefficients = coefs,
  declared     = declared,
  scores       = list(gsva = GS, mitopps = PPS, lines = lines),
  lines        = tibble::tibble(
    lineage = LINEAGE, n_model = length(breast), n_expr = length(lines),
    n_crispr = length(dep_lines), n_pan = length(pan_lines)),
  spec = list(
    release   = DEPMAP_RELEASE,
    prediction = "plan section 10: negative MYC:OX on MCL1 / BCL2L1 gene effect",
    two_sided = paste("Block C found LOWER BCL2L1 transcript under MYC x",
                      "OXPHOS; dependency and supply are different quantities",
                      "and the reverse is reported but not scored as a pass"),
    scale     = paste("GSVA on log2(TPM+1); mitoPPS on linear TPM (2^x - 1),",
                      "no floor - a non-positive pathway mean stops the script"),
    deviation = paste("mitoPPS input is linear TPM, not linear DESeq2-",
                      "normalised counts. Declared. CCLE mitoPPS values are",
                      "NEVER comparable to TCGA mitoPPS values."),
    no_null   = paste("no expression-matched null exists in CCLE; the arm",
                      "panel is a rank ordering. Build one before reporting",
                      "any positive OXPHOS result."),
    power     = paste(length(dep_lines), "breast lines; a null is not",
                      "evidence against"),
    prism     = if (have_prism)
                  paste("run; PRISM Repurposing", PRISM_RELEASE,
                        "- its own release cadence, joined on ModelID")
                else "deferred - files absent",
    prism_gap = paste("A-1331852 and A-1155463, the two SELECTIVE BCL-XL",
                      "inhibitors, are not in PRISM Repurposing. navitoclax is",
                      "BCL2/BCL-XL/BCL-W and cannot separate BCL-XL from BCL2",
                      "alone; read it against venetoclax, which is present and",
                      "BCL2-selective. MCL1 is fully covered (S63845, AMG-176,",
                      "AZD5991).")),
  built = Sys.time())

saveRDS(out, file.path(DIR_RESULTS, "depmap_dependency.rds"))
readr::write_csv(coefs,    file.path(DIR_TABLES, "block_g_coefficients.csv"))
readr::write_csv(declared, file.path(DIR_TABLES, "block_g_declared.csv"))

message("\n14: done.")
message("    results/depmap_dependency.rds")
message("    outputs/tables/  2 tables")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  g <- readRDS(file.path(DIR_RESULTS, "depmap_dependency.rds"))

  # --- the declared predictions, and the n behind them, first --------------
  g$declared %>% as.data.frame() %>% print(row.names = FALSE)
  g$lines %>% as.data.frame() %>% print(row.names = FALSE)

  # --- is there anything to detect? ---------------------------------------
  # If MCL1 and BCL2L1 are near-zero gene effect across every breast line,
  # there is no dependency to modify and G-a / G-b are null by construction
  # rather than by result. Look before interpreting.
  cr <- data.table::fread(file.path(DIR_DATA, "raw", "depmap",
                                    "CRISPRGeneEffect.csv"), data.table = FALSE)
  rn <- cr[[1]]; cr <- as.matrix(cr[, -1]); rownames(cr) <- rn
  colnames(cr) <- trimws(sub("\\s*\\(\\d+\\)$", "", colnames(cr)))
  boxplot(cr[g$scores$lines, c("MCL1", "BCL2L1", "BCL2", "RPL3")],
          ylab = "Chronos gene effect", main = "breast lines")
  abline(h = c(0, -1), lty = c(1, 2))

  # --- the interaction, drawn, because n is small --------------------------
  # At ~50 lines a single leverage point can carry the whole coefficient.
  # Plot it before believing it.
  ox  <- g$scores$gsva["OXPHOS subunits", ]
  myc <- g$scores$gsva["MYC", ]
  hi  <- myc > stats::median(myc)
  plot(ox, cr[g$scores$lines, "BCL2L1"], col = ifelse(hi, "red", "grey40"),
       pch = 16, xlab = "OXPHOS subunits (GSVA, CCLE-relative)",
       ylab = "BCL2L1 gene effect")
  legend("topright", c("MYC high", "MYC low"), col = c("red", "grey40"), pch = 16)

  # --- the arm panel as a rank ordering, which is all it is ----------------
  g$coefficients %>%
    dplyr::filter(label == "G2 arm panel", instrument == "gsva",
                  gene == "BCL2L1") %>%
    dplyr::arrange(estimate) %>%
    dplyr::select(arm, estimate, ci_lo, ci_hi, p) %>%
    as.data.frame() %>% print(row.names = FALSE)

  # --- breast vs pan-cancer, side by side ----------------------------------
  # If they disagree, the honest reading is that n = 50 could not see it, not
  # that breast is special.
  g$coefficients %>%
    dplyr::filter(gene %in% c("MCL1", "BCL2L1"), arm == "OXPHOS subunits",
                  label %in% c("G2 breast", "G2 pan-cancer (secondary)")) %>%
    dplyr::select(label, gene, instrument, n, estimate, ci_lo, ci_hi, p) %>%
    as.data.frame() %>% print(row.names = FALSE)

  # --- do the two instruments agree in CCLE at all? ------------------------
  # In TCGA they can move in opposite directions. Check here too, before any
  # "both instruments" claim is made.
  plot(g$scores$gsva["OXPHOS subunits", ],
       g$scores$mitopps["OXPHOS subunits", ], pch = 16,
       xlab = "GSVA", ylab = "mitoPPS")
  stats::cor(g$scores$gsva["OXPHOS subunits", ],
             g$scores$mitopps["OXPHOS subunits", ], method = "spearman")

  # --- the deferred pieces -------------------------------------------------
  g$spec$prism    # Block G item 3, if PRISM was absent
  g$spec$no_null  # the CCLE matched null that does not exist yet

}
