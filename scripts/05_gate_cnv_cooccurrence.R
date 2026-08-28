# 05_gate_cnv_cooccurrence.R
# =============================================================================
# GATE G2 - copy-number co-occurrence in TCGA-BRCA.  Plan section 6.
#
# Question: does MYC (8q24.21) amplification co-occur with MCL1 (1q21.3) or
# BCL2L1 (20q11.21) amplification, or with BBC3 (19q13.32) LOSS, more than
# expected?  And is the co-altered group enriched for basal / TNBC?
#
# Decision this feeds: if yes, H1 becomes Panel b and H4 gains its stratifying
# variable.  If no, H2 moves to primary and Panel b becomes the FOXO3/PI3K
# decoupling figure.  If the partners SPLIT, BCL2L1 carries the Panel b
# decision - it is the PRIME denominator and the term Lee et al. never touch.
#
# SCALE DISCIPLINE: not applicable.  This script touches NO expression data.
# Its inputs are discrete GISTIC calls in {-2,-1,0,1,2} plus clinical
# annotation.  The absence of a scale block here is deliberate, not an
# omission.  Do not add expression to this script; G2 is pure genomics.
#
# SPECIES: human.  TCGA-BRCA, human GISTIC, human symbols.  See CLAUDE.md.
#
# All design choices are pre-registered in
#   docs/2026-08-28_D4_scoping_and_G2_design.md
# and the constants live in 00_setup_packages.R.  Do not change a threshold
# here; change it there, with a dated reason, or the pre-registration record
# stops meaning anything.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
library(data.table)

message("\n05: GATE G2 - CNV co-occurrence\n", strrep("=", 78))

# -----------------------------------------------------------------------------
# 1. Sample universe
# -----------------------------------------------------------------------------
# BRCA membership comes from the GDC case list, NOT from either GISTIC file.
# Using one GISTIC run to define the sample set of the other would make the
# "independent" sensitivity analysis circular.

gdc_cases <- readr::read_tsv(PATH_GDC_BRCA_CASES, show_col_types = FALSE)
brca_patients <- unique(gdc_cases$submitter_id)

stopifnot(length(brca_patients) == EXPECT_GDC_BRCA_CASES)
message("1. GDC TCGA-BRCA cases: ", length(brca_patients))

clinical <- readr::read_tsv(PATH_TCGA_CLINICAL, show_col_types = FALSE,
                            na = "NA") %>%
  dplyr::mutate(
    ANEUPLOIDY_SCORE        = as.numeric(ANEUPLOIDY_SCORE),
    FRACTION_GENOME_ALTERED = as.numeric(FRACTION_GENOME_ALTERED)
  )
stopifnot(nrow(clinical) == EXPECT_GDC_BRCA_CASES)
stopifnot(!any(duplicated(clinical$patient_barcode)))

# --- TNBC, per the rule fixed in data/tcga_clinical/README.md ---------------
# HER2: FISH preferred where it is Positive/Negative, else IHC.  An
# IHC-equivocal case with no FISH result is UNCALLABLE, not negative.  Calling
# those negative would inflate TNBC, which is why the rule is written down
# rather than decided inline.
.clean_rec <- function(x) ifelse(x %in% c("Positive", "Negative"), x, NA_character_)

clinical <- clinical %>%
  dplyr::mutate(
    er_call   = .clean_rec(ER_STATUS_BY_IHC),
    pr_call   = .clean_rec(PR_STATUS_BY_IHC),
    her2_call = dplyr::coalesce(.clean_rec(HER2_FISH_STATUS),
                                .clean_rec(IHC_HER2)),
    receptor_callable = !is.na(er_call) & !is.na(pr_call) & !is.na(her2_call),
    tnbc = dplyr::case_when(
      !receptor_callable ~ NA,
      er_call == "Negative" & pr_call == "Negative" & her2_call == "Negative" ~ TRUE,
      TRUE ~ FALSE
    ),
    pam50 = dplyr::na_if(PAM50_SUBTYPE, "NA")
  )

stopifnot(sum(clinical$receptor_callable) == EXPECT_TNBC_CALLABLE)
stopifnot(sum(clinical$tnbc, na.rm = TRUE) == EXPECT_TNBC_N)
message("   receptor-callable: ", sum(clinical$receptor_callable),
        " | IHC-TNBC: ", sum(clinical$tnbc, na.rm = TRUE))

# -----------------------------------------------------------------------------
# 2. Read the GISTIC calls
# -----------------------------------------------------------------------------
# The ISAR file is 536 MB uncompressed and pan-cancer.  fread via a gzip pipe
# keeps it out of memory as a whole.  Do NOT read.delim this.

G2_GENES <- G2_THRESHOLDS$gene

.read_gistic <- function(path, cmd = NULL, label) {
  dt <- if (is.null(cmd)) {
    data.table::fread(path, sep = "\t", header = TRUE, showProgress = FALSE)
  } else {
    data.table::fread(cmd = cmd, sep = "\t", header = TRUE, showProgress = FALSE)
  }
  setnames(dt, 1:3, c("gene_symbol", "locus_id", "cytoband"))
  message("2. ", label, ": ", nrow(dt), " genes x ", ncol(dt) - 3L, " samples")
  dt
}

# --- Firehose (source B) is inside a tarball; extract once -------------------
if (!file.exists(PATH_GISTIC_FH_TXT)) {
  message("   extracting Firehose all_thresholded.by_genes.txt from tarball ...")
  .ensure_dir(dirname(PATH_GISTIC_FH_TXT))
  tmp <- tempfile("fh_")
  dir.create(tmp)
  utils::untar(PATH_GISTIC_FH_TAR, files = NULL, list = FALSE, exdir = tmp)
  hit <- list.files(tmp, pattern = "^all_thresholded\\.by_genes\\.txt$",
                    recursive = TRUE, full.names = TRUE)
  if (length(hit) != 1L) {
    stop("expected exactly one all_thresholded.by_genes.txt in the Firehose ",
         "tarball, found ", length(hit), call. = FALSE)
  }
  file.copy(hit, PATH_GISTIC_FH_TXT)
  unlink(tmp, recursive = TRUE)
}

gistic_isar <- .read_gistic(
  PATH_GISTIC_ISAR,
  cmd   = paste("gzip -dc", shQuote(PATH_GISTIC_ISAR)),
  label = "ISAR (A, primary)"
)
gistic_fh <- .read_gistic(PATH_GISTIC_FH_TXT, label = "Firehose (B, sensitivity)")

stopifnot(nrow(gistic_isar) == EXPECT_GISTIC_ISAR_GENES)
stopifnot(ncol(gistic_isar) - 3L == EXPECT_GISTIC_ISAR_SAMPLES)
stopifnot(nrow(gistic_fh)   == EXPECT_GISTIC_FH_GENES)
stopifnot(ncol(gistic_fh)   - 3L == EXPECT_GISTIC_FH_SAMPLES)

# -----------------------------------------------------------------------------
# 3. Subset to BRCA primary tumours, one column per patient
# -----------------------------------------------------------------------------
# Rule fixed 2026-08-28:
#   (a) patient barcode (first three fields) is in the GDC BRCA case list
#   (b) sample type is -01 (primary solid tumour)
#   (c) assert one column per patient afterwards
#
# Step (b) is what removes the 7 ISAR patients that carry both a primary and a
# metastasis (-06).  Deduplicating by "take the first column" would silently
# keep metastases for those 7.

# Split rather than regex: "(?:...)" is a PCRE construct that R's default TRE
# engine happens to accept, which is not something to depend on.
.patient_of <- function(bc) {
  vapply(strsplit(bc, "-"), function(p) paste(p[1:3], collapse = "-"), character(1))
}
.stype_of <- function(bc) {
  substr(vapply(strsplit(bc, "-"), `[`, character(1), 4), 1, 2)
}

.build_calls <- function(dt, label, expect_n) {
  bc  <- setdiff(names(dt), c("gene_symbol", "locus_id", "cytoband"))
  pat <- .patient_of(bc)
  st  <- .stype_of(bc)

  keep <- bc[pat %in% brca_patients & st == "01"]
  kpat <- .patient_of(keep)

  if (any(duplicated(kpat))) {
    stop(label, ": ", sum(duplicated(kpat)),
         " patient(s) still duplicated after the -01 filter. Inspect before ",
         "proceeding; do not dedupe arbitrarily.", call. = FALSE)
  }

  sub <- dt[gene_symbol %in% G2_GENES, c("gene_symbol", "cytoband", keep),
            with = FALSE]

  missing_genes <- setdiff(G2_GENES, sub$gene_symbol)
  if (length(missing_genes)) {
    stop(label, ": gene(s) absent from GISTIC: ",
         paste(missing_genes, collapse = ", "), call. = FALSE)
  }
  if (any(duplicated(sub$gene_symbol))) {
    stop(label, ": duplicated gene rows in GISTIC.", call. = FALSE)
  }

  m <- as.matrix(sub[, keep, with = FALSE])
  storage.mode(m) <- "integer"
  rownames(m) <- sub$gene_symbol
  colnames(m) <- kpat

  bad <- setdiff(unique(as.vector(m)), c(-2L, -1L, 0L, 1L, 2L))
  if (length(bad)) {
    stop(label, ": unexpected GISTIC value(s): ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  if (anyNA(m)) stop(label, ": NA in GISTIC calls.", call. = FALSE)

  message("3. ", label, ": ", ncol(m), " BRCA primary tumours, ",
          nrow(m), " genes of interest")
  stopifnot(ncol(m) == expect_n)

  # Cytobands must match the plan.  A silent build change here would move the
  # loci out from under the hypothesis.
  cb <- stats::setNames(sub$cytoband, sub$gene_symbol)
  expect_cb <- c(MYC = "8q24.21", MCL1 = "1q21.3", BCL2L1 = "20q11.21",
                 BBC3 = "19q13.32", BAX = "19q13.33")
  for (g in names(expect_cb)) {
    if (!identical(unname(cb[g]), unname(expect_cb[g]))) {
      stop(label, ": ", g, " cytoband is ", cb[g], ", expected ",
           expect_cb[g], call. = FALSE)
    }
  }
  m
}

calls_isar <- .build_calls(gistic_isar, "ISAR (A)", EXPECT_BRCA_N_ISAR)
calls_fh   <- .build_calls(gistic_fh,   "Firehose (B)", EXPECT_BRCA_N_FH)

# A must be a strict subset of B at patient level; the sensitivity comparison
# is only interpretable if we know how the sample sets relate.
pat_isar <- colnames(calls_isar)
pat_fh   <- colnames(calls_fh)
pat_both <- intersect(pat_isar, pat_fh)
stopifnot(length(setdiff(pat_isar, pat_fh)) == 0L)
message("   A n B = ", length(pat_both),
        " | A only = ", length(setdiff(pat_isar, pat_fh)),
        " | B only = ", length(setdiff(pat_fh, pat_isar)))

# -----------------------------------------------------------------------------
# 4. Binarise, per the pre-registered direction and threshold
# -----------------------------------------------------------------------------
# Four rules. "gain" genes get eq2 / ge1; "loss" genes get le_neg1 / le_neg2.
# The alternative rule for each gene is reported alongside its primary so the
# sensitivity of the result to that choice is visible rather than assumed.

.apply_rule <- function(v, rule) {
  switch(rule,
    eq2     = v ==  2L,                # gain: high-level amplification
    ge1     = v >=  1L,                # gain: gain or amplification
    le_neg1 = v <= -1L,                # loss: any loss
    le_neg2 = v == -2L,                # loss: homozygous deletion
    stop("unknown rule: ", rule, call. = FALSE)
  )
}

.binarise <- function(m, gene, rule) .apply_rule(m[gene, ], rule)

.primary_rule_for <- function(g) {
  r <- G2_THRESHOLDS$primary_rule[match(g, G2_THRESHOLDS$gene)]
  if (is.na(r)) stop("no pre-registered rule for gene: ", g, call. = FALSE)
  r
}
# Keyed on the primary RULE, not on direction. See 00_setup_packages.R.
.alt_rule_for <- function(g) unname(G2_ALT_RULE[.primary_rule_for(g)])

# -----------------------------------------------------------------------------
# 5. Descriptive frequency grid
# -----------------------------------------------------------------------------
.freq_grid <- function(m, source_label) {
  purrr::map_dfr(seq_len(nrow(G2_THRESHOLDS)), function(i) {
    g     <- G2_THRESHOLDS$gene[i]
    dirn  <- G2_THRESHOLDS$direction[i]
    prim  <- G2_THRESHOLDS$primary_rule[i]
    alt   <- .alt_rule_for(g)
    v     <- m[g, ]
    tibble::tibble(
      source    = source_label,
      gene      = g,
      direction = dirn,
      role      = G2_THRESHOLDS$role[i],
      n         = length(v),
      n_m2 = sum(v == -2L), n_m1 = sum(v == -1L), n_0 = sum(v == 0L),
      n_p1 = sum(v ==  1L), n_p2 = sum(v ==  2L),
      primary_rule = prim,
      n_primary    = sum(.apply_rule(v, prim)),
      pct_primary  = 100 * mean(.apply_rule(v, prim)),
      alt_rule     = alt,
      n_alt        = sum(.apply_rule(v, alt)),
      pct_alt      = 100 * mean(.apply_rule(v, alt))
    )
  })
}

freq_grid <- dplyr::bind_rows(
  .freq_grid(calls_isar, "A_ISAR"),
  .freq_grid(calls_fh,   "B_Firehose"),
  .freq_grid(calls_fh[, pat_both, drop = FALSE], "B_Firehose_intersection")
)

message("\n4. Frequency grid (primary rule per gene)")
freq_grid %>%
  dplyr::filter(source == "A_ISAR") %>%
  dplyr::select(gene, direction, role, primary_rule, n_primary, pct_primary,
                alt_rule, n_alt, pct_alt) %>%
  as.data.frame() %>%
  print()

# -----------------------------------------------------------------------------
# 6. Co-occurrence testing
# -----------------------------------------------------------------------------
# Three estimates per pair, reported side by side so the size of the aneuploidy
# confound is visible rather than assumed:
#   (i)   unadjusted Fisher            - shows the confound's magnitude
#   (ii)  logistic + ANEUPLOIDY_SCORE  - PRIMARY
#   (iii) CMH by aneuploidy tertile    - assumption-light check on (ii)
#
# 8q and 1q are among the most frequently gained arms in breast cancer, so an
# unconditioned test largely measures aneuploidy.  An unconditioned result is
# NOT reportable on its own.

.aneu_for <- function(patients, which = "ANEUPLOIDY_SCORE") {
  idx <- match(patients, clinical$patient_barcode)
  clinical[[which]][idx]
}

.test_pair <- function(m, exposure, partner, source_label,
                       exp_rule, prt_rule, tier, aneu_col = "ANEUPLOIDY_SCORE") {
  x <- .binarise(m, exposure, exp_rule)
  y <- .binarise(m, partner,  prt_rule)
  a <- .aneu_for(colnames(m), aneu_col)

  tab <- table(factor(x, c(FALSE, TRUE)), factor(y, c(FALSE, TRUE)))
  ft  <- stats::fisher.test(tab)

  # Fit from a named data.frame rather than bare vectors: glm() would otherwise
  # name the coefficient "x[ok]TRUE", and a lookup on that string breaks the
  # moment the call is edited.
  ok  <- !is.na(a)
  d   <- data.frame(y = y, x = x, aneu = a)[ok, , drop = FALSE]
  fit <- stats::glm(y ~ x + aneu, data = d, family = stats::binomial())
  cf  <- summary(fit)$coefficients
  term <- "xTRUE"
  stopifnot(term %in% rownames(cf))

  # Profile-likelihood CI where it converges, Wald otherwise. Separation is a
  # real possibility on the sparser cells, and a hard error there would take
  # down the whole grid.
  ci <- tryCatch(suppressMessages(stats::confint(fit))[term, ],
                 error = function(e) {
                   est <- cf[term, "Estimate"]; se <- cf[term, "Std. Error"]
                   c(est - 1.96 * se, est + 1.96 * se)
                 })
  adj_or <- exp(unname(cf[term, "Estimate"]))
  adj_lo <- exp(unname(ci[1]))
  adj_hi <- exp(unname(ci[2]))

  # CMH stratified by aneuploidy tertile.
  ter <- cut(a[ok], stats::quantile(a[ok], c(0, 1/3, 2/3, 1), na.rm = TRUE),
             include.lowest = TRUE, labels = c("low", "mid", "high"))
  arr <- table(factor(x[ok], c(FALSE, TRUE)),
               factor(y[ok], c(FALSE, TRUE)), ter)
  cmh <- tryCatch(stats::mantelhaen.test(arr, exact = FALSE),
                  error = function(e) NULL)

  tibble::tibble(
    source = source_label, tier = tier,
    exposure = exposure, exposure_rule = exp_rule,
    partner = partner, partner_rule = prt_rule,
    n = length(x), n_exposure = sum(x), n_partner = sum(y),
    n_both = sum(x & y),
    unadj_or = unname(ft$estimate), unadj_p = ft$p.value,
    n_adj = sum(ok), aneu_covariate = aneu_col,
    adj_or = adj_or, adj_lo = adj_lo, adj_hi = adj_hi,
    adj_p  = unname(cf[term, "Pr(>|z|)"]),
    cmh_or = if (is.null(cmh)) NA_real_ else unname(cmh$estimate),
    cmh_p  = if (is.null(cmh)) NA_real_ else cmh$p.value,
    passes = !is.na(adj_lo) & adj_lo > 1
  )
}

# --- primary tests: three pairs, one per partner -----------------------------
myc_rule <- .primary_rule_for("MYC")

primary_pairs <- c("MCL1", "BCL2L1", "BBC3")

res_primary <- purrr::map_dfr(primary_pairs, function(p) {
  .test_pair(calls_isar, "MYC", p, "A_ISAR",
             myc_rule, .primary_rule_for(p), tier = "primary")
})

# --- regional control: BAX, same rule as BBC3, 19q13 -------------------------
# If the BBC3 estimate is indistinguishable from this one, the BBC3 finding is
# arm-level 19q behaviour and not PUMA-specific.  This is a control, NOT a
# hypothesis, and is labelled as such in the output.
res_control <- .test_pair(calls_isar, "MYC", "BAX", "A_ISAR",
                          myc_rule, .primary_rule_for("BAX"),
                          tier = "regional_control")

# --- secondary: full gene x threshold grid -----------------------------------
# GUARD: ORs at different thresholds are NOT comparable across genes. Compare
# within a threshold only.
res_grid <- purrr::map_dfr(c(primary_pairs, "BAX"), function(p) {
  purrr::map_dfr(c(.primary_rule_for(p), .alt_rule_for(p)), function(r) {
    purrr::map_dfr(c("eq2", "ge1"), function(mr) {
      .test_pair(calls_isar, "MYC", p, "A_ISAR", mr, r, tier = "secondary_grid")
    })
  })
})

# --- sensitivity: source B, full and intersection ----------------------------
res_sens <- dplyr::bind_rows(
  purrr::map_dfr(primary_pairs, function(p) {
    .test_pair(calls_fh, "MYC", p, "B_Firehose",
               myc_rule, .primary_rule_for(p), tier = "sensitivity_source")
  }),
  purrr::map_dfr(primary_pairs, function(p) {
    .test_pair(calls_fh[, pat_both, drop = FALSE], "MYC", p,
               "B_Firehose_intersection",
               myc_rule, .primary_rule_for(p), tier = "sensitivity_source")
  })
)

# --- sensitivity: FGA instead of aneuploidy score ----------------------------
res_fga <- purrr::map_dfr(primary_pairs, function(p) {
  .test_pair(calls_isar, "MYC", p, "A_ISAR",
             myc_rule, .primary_rule_for(p), tier = "sensitivity_fga",
             aneu_col = "FRACTION_GENOME_ALTERED")
})

cooccurrence <- dplyr::bind_rows(res_primary, res_control, res_grid,
                                 res_sens, res_fga) %>%
  dplyr::group_by(tier) %>%
  dplyr::mutate(adj_p_bh = ifelse(tier == "primary", NA_real_,
                                  stats::p.adjust(adj_p, "BH"))) %>%
  dplyr::ungroup()

message("\n5. Primary co-occurrence tests (aneuploidy-adjusted)")
res_primary %>%
  dplyr::select(partner, partner_rule, n_both, unadj_or, adj_or, adj_lo,
                adj_hi, adj_p, passes) %>%
  as.data.frame() %>%
  print()

message("\n   19q13 regional control (BAX), compare against BBC3 above")
res_control %>%
  dplyr::select(partner, partner_rule, n_both, unadj_or, adj_or, adj_lo,
                adj_hi, adj_p) %>%
  as.data.frame() %>%
  print()

# -----------------------------------------------------------------------------
# 7. Stratified analysis: PAM50 and TNBC
# -----------------------------------------------------------------------------
# BRCA_Normal is artefact-prone and usually reflects low tumour cellularity.
# It is reported as its own stratum and must never become a silent reference
# level.

.stratified <- function(m, partner, strat_col, source_label) {
  idx  <- match(colnames(m), clinical$patient_barcode)
  lev  <- clinical[[strat_col]][idx]
  keep <- !is.na(lev)
  purrr::map_dfr(sort(unique(lev[keep])), function(l) {
    sel <- which(keep & lev == l)
    if (length(sel) < 20L) {
      return(tibble::tibble(source = source_label, strat_var = strat_col,
                            stratum = as.character(l), partner = partner,
                            n = length(sel), skipped = TRUE))
    }
    out <- .test_pair(m[, sel, drop = FALSE], "MYC", partner, source_label,
                      myc_rule, .primary_rule_for(partner),
                      tier = paste0("stratified_", strat_col))
    dplyr::mutate(out, strat_var = strat_col, stratum = as.character(l),
                  skipped = FALSE)
  })
}

stratified <- dplyr::bind_rows(
  purrr::map_dfr(primary_pairs, ~ .stratified(calls_isar, .x, "pam50", "A_ISAR")),
  purrr::map_dfr(primary_pairs, ~ .stratified(calls_isar, .x, "tnbc",  "A_ISAR"))
)

# Breslow-Day: does the association genuinely DIFFER by subtype, rather than
# eyeballing five odds ratios?
.breslow_day <- function(m, partner, strat_col) {
  idx  <- match(colnames(m), clinical$patient_barcode)
  lev  <- clinical[[strat_col]][idx]
  keep <- !is.na(lev)
  x <- .binarise(m, "MYC", myc_rule)[keep]
  y <- .binarise(m, partner, .primary_rule_for(partner))[keep]
  arr <- table(factor(x, c(FALSE, TRUE)), factor(y, c(FALSE, TRUE)),
               droplevels(factor(lev[keep])))
  # Drop strata with an empty margin; Breslow-Day is undefined there.
  ok <- apply(arr, 3, function(t) all(rowSums(t) > 0) && all(colSums(t) > 0))
  arr <- arr[, , ok, drop = FALSE]
  if (dim(arr)[3] < 2L) return(tibble::tibble(partner = partner,
                                              strat_var = strat_col,
                                              bd_p = NA_real_, n_strata = dim(arr)[3]))
  mh  <- stats::mantelhaen.test(arr, exact = FALSE)
  psi <- unname(mh$estimate)
  bd  <- sum(vapply(seq_len(dim(arr)[3]), function(k) {
    t <- arr[, , k]
    n1 <- sum(t[2, ]); n0 <- sum(t[1, ]); m1 <- sum(t[, 2]); nn <- sum(t)
    f <- function(a) a * (nn - n1 - m1 + a) - psi * (n1 - a) * (m1 - a)
    lo <- max(0, m1 + n1 - nn); hi <- min(n1, m1)
    if (f(lo) * f(hi) > 0) return(NA_real_)
    a_hat <- stats::uniroot(f, c(lo, hi))$root
    v <- 1 / sum(1 / c(a_hat, n1 - a_hat, m1 - a_hat, nn - n1 - m1 + a_hat))
    (t[2, 2] - a_hat)^2 / v
  }, numeric(1)), na.rm = TRUE)
  tibble::tibble(partner = partner, strat_var = strat_col,
                 bd_stat = bd, n_strata = dim(arr)[3],
                 bd_p = stats::pchisq(bd, dim(arr)[3] - 1L, lower.tail = FALSE))
}

homogeneity <- dplyr::bind_rows(
  purrr::map_dfr(primary_pairs, ~ .breslow_day(calls_isar, .x, "pam50")),
  purrr::map_dfr(primary_pairs, ~ .breslow_day(calls_isar, .x, "tnbc"))
)

# -----------------------------------------------------------------------------
# 8. Reverse question: is the co-altered group enriched for Basal / TNBC?
# -----------------------------------------------------------------------------
.enrichment <- function(m, partner, target_col, target_val = NULL) {
  idx <- match(colnames(m), clinical$patient_barcode)
  tgt <- if (identical(target_col, "pam50")) {
    clinical$pam50[idx] == target_val
  } else {
    clinical$tnbc[idx]
  }
  label <- if (is.null(target_val)) target_col else paste0(target_col, ":", target_val)
  co  <- .binarise(m, "MYC", myc_rule) &
         .binarise(m, partner, .primary_rule_for(partner))
  ok  <- !is.na(tgt)
  ft  <- stats::fisher.test(table(factor(co[ok], c(FALSE, TRUE)),
                                  factor(tgt[ok], c(FALSE, TRUE))))
  tibble::tibble(
    partner = partner, target = label,
    n = sum(ok), n_coaltered = sum(co[ok]),
    n_coaltered_target = sum(co[ok] & tgt[ok]),
    or = unname(ft$estimate), lo = ft$conf.int[1], hi = ft$conf.int[2],
    p = ft$p.value
  )
}

enrichment <- dplyr::bind_rows(
  purrr::map_dfr(primary_pairs, ~ .enrichment(calls_isar, .x, "pam50", "BRCA_Basal")),
  purrr::map_dfr(primary_pairs, ~ .enrichment(calls_isar, .x, "tnbc"))
)

message("\n6. Co-altered group enrichment for Basal / TNBC")
enrichment %>% as.data.frame() %>% print()

# -----------------------------------------------------------------------------
# 9. BBC3 descriptive block - the empty deep-deletion cell
# -----------------------------------------------------------------------------
# Reported because the -2 cell is what motivated moving BBC3 to <= -1. Do NOT
# "fix" this by running a Fisher test on it: the joint count with MYC
# amplification is zero, so the test is empty rather than merely underpowered.
bbc3_deep <- local({
  x <- .binarise(calls_isar, "MYC", myc_rule)
  y <- .binarise(calls_isar, "BBC3", "le_neg2")
  tibble::tibble(
    gene = "BBC3", rule = "le_neg2 (homozygous deletion)",
    n = length(y), n_deleted = sum(y),
    n_both_with_MYC_amp = sum(x & y),
    note = "descriptive only; no test is run, the joint cell is empty"
  )
})
message("\n7. BBC3 homozygous deletion (descriptive, no test)")
bbc3_deep %>% as.data.frame() %>% print()

# -----------------------------------------------------------------------------
# 10. Gate verdict
# -----------------------------------------------------------------------------
# Criterion fixed 2026-08-28 BEFORE any statistic was computed (design note
# 2.8): a partner passes if the aneuploidy-adjusted OR exceeds 1 with a 95% CI
# excluding 1, at its primary threshold and direction, in source A. Source B
# must agree in direction; a disagreement is reported, never overridden.
verdict <- res_primary %>%
  dplyr::select(partner, adj_or, adj_lo, adj_hi, adj_p, passes) %>%
  dplyr::left_join(
    res_sens %>%
      dplyr::filter(source == "B_Firehose") %>%
      dplyr::transmute(partner, b_or = adj_or),
    by = "partner"
  ) %>%
  # Direction agreement between source A and source B, on the log-OR scale.
  # Comparing A against A would always agree; the join above is what makes this
  # an actual cross-source check.
  dplyr::mutate(b_agrees_direction = sign(log(adj_or)) == sign(log(b_or)))

message("\n", strrep("=", 78))
message("G2 VERDICT (criterion: aneuploidy-adjusted OR > 1, 95% CI excluding 1)")
verdict %>% as.data.frame() %>% print()
message("\nIf the partners split, BCL2L1 carries the Panel b decision.")
message("BBC3 is NOT sufficient alone - check it against the BAX regional control.")
message(strrep("=", 78))

# -----------------------------------------------------------------------------
# 11. Save
# -----------------------------------------------------------------------------
g2 <- list(
  freq_grid    = freq_grid,
  cooccurrence = cooccurrence,
  stratified   = stratified,
  homogeneity  = homogeneity,
  enrichment   = enrichment,
  bbc3_deep    = bbc3_deep,
  verdict      = verdict,
  thresholds   = G2_THRESHOLDS,
  samples      = list(isar = pat_isar, firehose = pat_fh, intersection = pat_both),
  generated    = Sys.time()
)

saveRDS(g2, file.path(DIR_RESULTS, "g2_cnv_cooccurrence.rds"))
readr::write_csv(freq_grid,    file.path(DIR_TABLES, "g2_frequency_grid.csv"))
readr::write_csv(cooccurrence, file.path(DIR_TABLES, "g2_cooccurrence.csv"))
readr::write_csv(stratified,   file.path(DIR_TABLES, "g2_stratified.csv"))
readr::write_csv(enrichment,   file.path(DIR_TABLES, "g2_enrichment.csv"))

message("\n05: done. results/g2_cnv_cooccurrence.rds + 4 tables in outputs/tables/")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  # --- the three primary tests, side by side -------------------------------
  res_primary %>%
    dplyr::select(partner, partner_rule, n_exposure, n_partner, n_both,
                  unadj_or, adj_or, adj_lo, adj_hi, adj_p, passes) %>%
    print(width = Inf)

  # --- how big is the aneuploidy confound? ---------------------------------
  # The gap between unadj_or and adj_or IS the confound. If they are far apart,
  # say so explicitly in the write-up rather than reporting only the adjusted
  # number.
  res_primary %>%
    dplyr::transmute(partner, unadj_or, adj_or,
                     shrinkage = unadj_or - adj_or) %>%
    print()

  # --- is BBC3 distinguishable from its 19q13 regional control? ------------
  dplyr::bind_rows(
    dplyr::filter(res_primary, partner == "BBC3"),
    res_control
  ) %>%
    dplyr::select(partner, adj_or, adj_lo, adj_hi, adj_p) %>%
    print()

  # --- threshold sensitivity, within a threshold only ----------------------
  # Reminder: do NOT compare across the partner_rule column between genes.
  res_grid %>%
    dplyr::filter(exposure_rule == "eq2") %>%
    dplyr::select(partner, partner_rule, n_both, adj_or, adj_lo, adj_hi) %>%
    print(n = 20)

  # --- does the association differ by subtype? -----------------------------
  homogeneity %>% print()
  stratified %>%
    dplyr::filter(strat_var == "pam50", !skipped) %>%
    dplyr::select(partner, stratum, n, n_both, adj_or, adj_lo, adj_hi) %>%
    print(n = 30)

  # --- source A vs source B ------------------------------------------------
  cooccurrence %>%
    dplyr::filter(tier %in% c("primary", "sensitivity_source")) %>%
    dplyr::select(source, partner, n, adj_or, adj_lo, adj_hi) %>%
    print(n = 20)

  # --- aneuploidy covariate swap -------------------------------------------
  dplyr::bind_rows(res_primary, res_fga) %>%
    dplyr::select(tier, partner, aneu_covariate, adj_or, adj_lo, adj_hi) %>%
    print()

  # --- sanity: the -01 filter really did remove the metastases -------------
  bc <- setdiff(names(gistic_isar), c("gene_symbol", "locus_id", "cytoband"))
  table(.stype_of(bc[.patient_of(bc) %in% brca_patients]))

}
