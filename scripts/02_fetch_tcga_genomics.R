# 02_fetch_tcga_genomics.R
# =============================================================================
# TCGA-BRCA genomics and protein: mutation status, tumour purity, leukocyte
# fraction, arm-level aneuploidy, and RPPA. Produces tidy per-patient tables
# for script 03 to assemble into the covariate table.
#
# SCALE DISCIPLINE: not applicable. No RNA expression is touched here. RPPA is
# already level-4 normalised and is carried through as supplied - it is NOT
# re-normalised, logged, or mixed with the expression matrices from script 01.
#
# SPECIES: human. TCGA-BRCA. See CLAUDE.md.
#
# SOURCES: five PanCanAtlas files, all single downloads with recorded SHA-256,
# snapshotted under data/raw/tcga_pancanatlas/ (gitignored). Provenance in
# data/tcga_pancanatlas/README.md.
#
# Deliberately all PanCanAtlas, matching the ISAR GISTIC already used for G2.
# Mixing a per-project GDC MAF with PanCanAtlas ABSOLUTE and Thorsson would put
# the covariate set on inconsistent calling pipelines for no gain.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))
suppressPackageStartupMessages(library(data.table))

message("\n02: TCGA-BRCA genomics\n", strrep("=", 78))

DIR_PCA <- file.path(DIR_DATA, "raw", "tcga_pancanatlas")

PATH_MC3   <- file.path(DIR_PCA, "mc3.v0.2.8.PUBLIC.maf.gz")
PATH_ABS   <- file.path(DIR_PCA, "TCGA_mastercalls.abs_tables_JSedit.fixed.txt")
PATH_LEUK  <- file.path(DIR_PCA, "TCGA_all_leuk_estimate.masked.20170107.tsv")
PATH_RPPA  <- file.path(DIR_PCA, "TCGA-RPPA-pancan-clean.txt")
PATH_ARM   <- file.path(DIR_PCA, "PANCAN_ArmCallsAndAneuploidyScore_092817.txt")
PATH_MC3_CACHE <- file.path(DIR_PCA, "mc3_brca_subset.rds")

for (p in c(PATH_MC3, PATH_ABS, PATH_LEUK, PATH_RPPA, PATH_ARM)) {
  if (!file.exists(p)) {
    stop("missing input: ", p,
         "\nSee data/tcga_pancanatlas/README.md for URLs and checksums.",
         call. = FALSE)
  }
}

# Expected sizes, so a PARTIAL download is caught here rather than 100 lines
# later as a misleading gzip or fread error. A truncated .gz surfaces as
# "unexpected end of file" from gzip, which fread reports as a possible full
# disk - a message that points nowhere near the actual cause.
EXPECT_BYTES <- c(
  "mc3.v0.2.8.PUBLIC.maf.gz"                    = 753339089,
  "TCGA_mastercalls.abs_tables_JSedit.fixed.txt" =    901812,
  "TCGA_all_leuk_estimate.masked.20170107.tsv"   =    560473,
  "TCGA-RPPA-pancan-clean.txt"                   =  18901234,
  "PANCAN_ArmCallsAndAneuploidyScore_092817.txt" =   1079849
)
for (p in c(PATH_MC3, PATH_ABS, PATH_LEUK, PATH_RPPA, PATH_ARM)) {
  nm  <- basename(p)
  got <- file.size(p)
  want <- EXPECT_BYTES[[nm]]
  if (!identical(as.numeric(got), as.numeric(want))) {
    stop(nm, " is ", format(got, big.mark = ","), " bytes but should be ",
         format(want, big.mark = ","), " (",
         sprintf("%.0f%%", 100 * got / want), " complete).\n",
         "The download is unfinished or corrupt. Wait for it, or re-fetch:\n",
         "  curl -L https://api.gdc.cancer.gov/data/<UUID> -o ", p, "\n",
         "UUIDs are in data/tcga_pancanatlas/README.md.", call. = FALSE)
  }
}
message("all 5 inputs present and complete")

.patient_of <- function(bc) {
  vapply(strsplit(bc, "-"), function(p) paste(p[1:3], collapse = "-"), character(1))
}
.stype_of <- function(bc) {
  substr(vapply(strsplit(bc, "-"), `[`, character(1), 4), 1, 2)
}

# Sample universe: the patients script 01 actually kept. Anchoring to 01 rather
# than to the GDC case list means the covariate tables cannot silently include
# patients that have no expression.
vst <- readRDS(file.path(DIR_RESULTS, "tcga_brca_vst.rds"))
brca_patients <- colnames(vst$mat)
rm(vst); invisible(gc())
message("patients carried forward from script 01: ", length(brca_patients))

# -----------------------------------------------------------------------------
# 1. Mutations (MC3)
# -----------------------------------------------------------------------------
# >>> PRE-REGISTRATION DECISION, NEEDS SIGN-OFF BEFORE THIS SCRIPT IS RUN.
#
# "TP53-mutant" and "PIK3CA/PTEN-altered" are strata in H2 and H3, so what
# counts as altered has to be fixed before any model is fitted.
#
# Rule, fixed here:
#   a. Only MC3 calls with FILTER == "PASS" are used. MC3 ships non-PASS calls
#      (wga, native_wga_mix, broad_PoN_v2, ...) and including them inflates
#      mutation counts inconsistently across genes.
#   b. A gene is MUTATED in a patient if it carries at least one call in
#      MUT_NONSILENT below. Silent, UTR, flank, intron and IGR are excluded.
#   c. PTEN_ALTERED = PTEN mutated OR GISTIC deep deletion (-2). PTEN is lost by
#      deletion at least as often as by mutation, so a mutation-only definition
#      would misclassify the H2 stratum. The CNV half is joined in script 03,
#      where the GISTIC calls live; this script emits the mutation half and
#      names the dependency.
#   d. PIK3CA_ALTERED = PIK3CA mutated. It is an oncogene activated by point
#      mutation; amplification is not used.
#
# Alternative considered and rejected: hotspot-only definitions for PIK3CA
# (H1047R/E545K/E542K). Cleaner biologically, but the plan's stratum is
# "PIK3CA-wild-type / PTEN-intact" as a whole, and a hotspot-only rule would put
# non-hotspot mutants into the wild-type arm, which is the wrong error to make.

MUT_NONSILENT <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
  "In_Frame_Del", "In_Frame_Ins", "Splice_Site", "Nonstop_Mutation",
  "Translation_Start_Site"
)
STRATUM_GENES <- c("TP53", "PIK3CA", "PTEN")
# Also carried, so the CNV analyses can be checked against mutation later.
EXTRA_GENES   <- c("MYC", "MCL1", "BCL2L1", "BBC3", "BAX", "BAK1", "BCL2")

if (file.exists(PATH_MC3_CACHE)) {
  message("\n1. cached MC3 BRCA subset found")
  maf <- readRDS(PATH_MC3_CACHE)
} else {
  message("\n1. reading MC3 MAF (753 MB gz; subset and cached, runs once)")
  maf_all <- data.table::fread(
    cmd    = paste("gzip -dc", shQuote(PATH_MC3)),
    sep    = "\t", header = TRUE, showProgress = FALSE,
    select = c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode",
               "HGVSp_Short", "FILTER")
  )
  message("   MC3 total calls: ", format(nrow(maf_all), big.mark = ","))
  maf_all[, patient := .patient_of(Tumor_Sample_Barcode)]
  maf <- maf_all[patient %in% brca_patients]
  message("   calls in our BRCA patients: ", format(nrow(maf), big.mark = ","))
  saveRDS(maf, PATH_MC3_CACHE)
  rm(maf_all); invisible(gc())
}

message("   FILTER values present: ",
        paste(sprintf("%s=%d", names(sort(table(maf$FILTER), decreasing = TRUE))[1:3],
                      sort(table(maf$FILTER), decreasing = TRUE)[1:3]), collapse = " "))

maf_pass <- maf[FILTER == "PASS" & Variant_Classification %in% MUT_NONSILENT]
message("   PASS + non-silent calls: ", format(nrow(maf_pass), big.mark = ","))
message("   patients with any such call: ", length(unique(maf_pass$patient)))

mut_status <- tibble::tibble(patient = brca_patients)
for (g in c(STRATUM_GENES, EXTRA_GENES)) {
  hit <- unique(maf_pass$patient[maf_pass$Hugo_Symbol == g])
  mut_status[[paste0(g, "_mut")]] <- mut_status$patient %in% hit
}
message("\n   mutation frequency in ", length(brca_patients), " patients:")
for (g in c(STRATUM_GENES, EXTRA_GENES)) {
  n <- sum(mut_status[[paste0(g, "_mut")]])
  message(sprintf("     %-8s %4d  (%4.1f%%)", g, n, 100 * n / nrow(mut_status)))
}

# Patients with NO MC3 coverage must not be scored as wild-type. Absence of a
# call is not evidence of absence of a mutation if the patient was never
# sequenced, so they are marked NA rather than FALSE.
covered <- brca_patients %in% unique(maf$patient)
message("\n   patients with MC3 coverage: ", sum(covered), " of ", length(covered))
if (any(!covered)) {
  message("   ", sum(!covered), " patient(s) set to NA, not FALSE (no MC3 coverage)")
  for (g in c(STRATUM_GENES, EXTRA_GENES)) {
    mut_status[[paste0(g, "_mut")]][!covered] <- NA
  }
}

# -----------------------------------------------------------------------------
# 2. Tumour purity and ploidy (ABSOLUTE)
# -----------------------------------------------------------------------------
abs_tab <- data.table::fread(PATH_ABS, sep = "\t", header = TRUE,
                             showProgress = FALSE)
data.table::setnames(abs_tab, make.names(names(abs_tab)))
abs_tab[, patient := .patient_of(sample)]
abs_tab[, stype   := .stype_of(sample)]
abs_brca <- abs_tab[stype == "01" & patient %in% brca_patients]
abs_brca <- abs_brca[!duplicated(patient)]

purity <- tibble::tibble(patient = brca_patients) %>%
  dplyr::left_join(
    tibble::as_tibble(abs_brca[, .(patient, purity, ploidy,
                                   genome_doublings = Genome.doublings)]),
    by = "patient"
  )
message("\n2. ABSOLUTE purity: ", sum(!is.na(purity$purity)), " of ",
        nrow(purity), " patients")
message(sprintf("   median purity %.2f | purity > 0.7 in %d patients",
                median(purity$purity, na.rm = TRUE),
                sum(purity$purity > 0.7, na.rm = TRUE)))

# -----------------------------------------------------------------------------
# 3. Leukocyte fraction (Thorsson)
# -----------------------------------------------------------------------------
# Headerless: cancer type, aliquot barcode, leukocyte fraction.
leuk <- data.table::fread(PATH_LEUK, sep = "\t", header = FALSE,
                          col.names = c("cancer_type", "aliquot", "leuk_frac"),
                          showProgress = FALSE)
leuk[, patient := .patient_of(aliquot)]
leuk[, stype   := .stype_of(aliquot)]
leuk_brca <- leuk[stype == "01" & patient %in% brca_patients]
# Several aliquots per patient are possible; average them rather than picking
# one arbitrarily - this is a continuous estimate, not a call.
leuk_p <- leuk_brca[, .(leukocyte_fraction = mean(leuk_frac, na.rm = TRUE)),
                    by = patient]

leukocyte <- tibble::tibble(patient = brca_patients) %>%
  dplyr::left_join(tibble::as_tibble(leuk_p), by = "patient")
message("\n3. leukocyte fraction: ", sum(!is.na(leukocyte$leukocyte_fraction)),
        " of ", nrow(leukocyte), " patients")

# -----------------------------------------------------------------------------
# 4. Arm-level aneuploidy (Taylor 2018)
# -----------------------------------------------------------------------------
# The authoritative version of the aneuploidy burden used to condition G2. The
# cBioPortal ANEUPLOIDY_SCORE snapshot in data/tcga_clinical/ came from the same
# study; this file additionally carries the per-arm calls, so 8q and 1q can be
# inspected directly - which matters because G2's confound is arm-level gain.
arm <- data.table::fread(PATH_ARM, sep = "\t", header = TRUE, showProgress = FALSE)
data.table::setnames(arm, make.names(names(arm)))
arm[, patient := .patient_of(Sample)]
arm[, stype   := .stype_of(Sample)]
arm_brca <- arm[stype == "01" & patient %in% brca_patients]
arm_brca <- arm_brca[!duplicated(patient)]

keep_arms <- intersect(c("X8q", "X1q", "X20q", "X19q"), names(arm_brca))
aneuploidy <- tibble::tibble(patient = brca_patients) %>%
  dplyr::left_join(
    tibble::as_tibble(arm_brca[, c("patient", "Aneuploidy.Score", keep_arms),
                               with = FALSE]),
    by = "patient"
  )
message("\n4. aneuploidy score: ", sum(!is.na(aneuploidy$Aneuploidy.Score)),
        " of ", nrow(aneuploidy), " patients")
message("   per-arm calls carried: ", paste(keep_arms, collapse = ", "))

# -----------------------------------------------------------------------------
# 5. RPPA
# -----------------------------------------------------------------------------
# Level-4 normalised as supplied. NOT re-normalised here.
rppa <- data.table::fread(PATH_RPPA, sep = "\t", header = TRUE,
                          showProgress = FALSE)
rppa[, patient := .patient_of(SampleID)]
rppa[, stype   := .stype_of(SampleID)]
rppa_brca <- rppa[stype == "01" & patient %in% brca_patients]
rppa_brca <- rppa_brca[!duplicated(patient)]

prot_cols <- setdiff(names(rppa_brca),
                     c("SampleID", "TumorType", "patient", "stype"))
rppa_mat <- as.matrix(rppa_brca[, prot_cols, with = FALSE])
rownames(rppa_mat) <- rppa_brca$patient

message("\n5. RPPA: ", nrow(rppa_mat), " patients x ", ncol(rppa_mat), " proteins")

# The proteins Block E depends on. Reported rather than assumed, because the
# panel composition decides what Block E can actually test.
BLOCK_E_PROTEINS <- c("BAK", "BAX", "BCL2", "BCLXL", "BID", "BIM", "CMYC",
                      "FOXO3A", "FOXO3A_pS318S321", "MCL1", "PUMA")
present <- intersect(BLOCK_E_PROTEINS, colnames(rppa_mat))
absent  <- setdiff(BLOCK_E_PROTEINS, colnames(rppa_mat))
message("   Block E proteins present: ", paste(present, collapse = ", "))
if (length(absent)) {
  message("   Block E proteins ABSENT : ", paste(absent, collapse = ", "))
  message("   NOTE: PRIME's numerator (PUMA/BBC3) is not on the RPPA panel, so ",
          "Block E can confirm the BCL-XL / MYC / FOXO3 arms at protein level ",
          "but NOT the PRIME ratio itself. State this in Methods.")
}

# -----------------------------------------------------------------------------
# 6. Save
# -----------------------------------------------------------------------------
genomics <- list(
  mutations   = mut_status,
  mc3_covered = tibble::tibble(patient = brca_patients, mc3_covered = covered),
  purity      = purity,
  leukocyte   = leukocyte,
  aneuploidy  = aneuploidy,
  rppa        = rppa_mat,
  rules       = list(nonsilent = MUT_NONSILENT, filter = "PASS",
                     stratum_genes = STRATUM_GENES, extra_genes = EXTRA_GENES),
  built       = Sys.time()
)
saveRDS(genomics, file.path(DIR_RESULTS, "tcga_brca_genomics.rds"))
readr::write_csv(mut_status, file.path(DIR_TABLES, "tcga_brca_mutation_status.csv"))

message("\n02: done. results/tcga_brca_genomics.rds")
message("    PTEN_ALTERED still needs the GISTIC deep-deletion half - script 03.")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  g <- readRDS(file.path(DIR_RESULTS, "tcga_brca_genomics.rds"))

  # --- do the strata look like breast cancer? ------------------------------
  # TP53 ~30-35%, PIK3CA ~35%, PTEN ~5% is the expected shape. A big deviation
  # means the FILTER or non-silent rule is wrong.
  colMeans(g$mutations[, -1], na.rm = TRUE) * 100

  # --- purity distribution, and the >0.7 sensitivity subset ----------------
  hist(g$purity$purity, breaks = 30, main = "ABSOLUTE purity", xlab = "purity")
  abline(v = 0.7, col = "red", lty = 2)
  sum(g$purity$purity > 0.7, na.rm = TRUE)

  # --- leukocyte fraction vs purity: should be strongly negative -----------
  plot(g$purity$purity, g$leukocyte$leukocyte_fraction,
       pch = 16, cex = .4, xlab = "ABSOLUTE purity", ylab = "leukocyte fraction")
  cor(g$purity$purity, g$leukocyte$leukocyte_fraction,
      use = "complete.obs", method = "spearman")

  # --- 8q / 1q arm calls, the G2 confound, seen directly -------------------
  table(g$aneuploidy$X8q, g$aneuploidy$X1q, useNA = "ifany")

  # --- RPPA: is the BCL2 family measurable and behaved? --------------------
  summary(g$rppa[, intersect(c("BCLXL","BAX","BAK","BID","BIM","CMYC","FOXO3A"),
                             colnames(g$rppa))])

  # --- coverage overlap across every source --------------------------------
  data.frame(
    source = c("mutations", "purity", "leukocyte", "aneuploidy", "rppa"),
    n = c(sum(g$mc3_covered$mc3_covered),
          sum(!is.na(g$purity$purity)),
          sum(!is.na(g$leukocyte$leukocyte_fraction)),
          sum(!is.na(g$aneuploidy$Aneuploidy.Score)),
          nrow(g$rppa))
  )

}
