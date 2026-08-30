# Neoadjuvant cohorts for H4

Provenance for the three cohorts D5 selected. Files live under
`data/raw/neoadjuvant/`, which is **gitignored and not on origin**. This README
is the record; re-download from the commands below.

Selected in `docs/2026-08-28_D5_cohort_selection.md` by a documented search
(447 GEO series screened, 34 passing), not by assertion. Consumed by
`scripts/12_fetch_neoadjuvant_cohorts.R`.

| Cohort | Role | n | Design |
|---|---|---|---|
| GSE194040 I-SPY2-990 | **PRIMARY** | 988 | randomised, 14 arms, 179-patient paclitaxel control |
| GSE164458 BrighTNess | replication | 482 | randomised, TNBC only, 3 arms |
| GSE25066 Hatzis | third | 508 | single-arm taxane-anthracycline |

## Download

Seven files, about 180 MB total. Read on 2026-08-30; sizes are the GEO
directory listing's.

```sh
mkdir -p data/raw/neoadjuvant && cd data/raw/neoadjuvant
B=https://ftp.ncbi.nlm.nih.gov/geo

# --- GSE194040 I-SPY2-990 (PRIMARY) ---
# expression: combined gene-level, n=988, both platforms in one file (47M)
curl -L -O $B/series/GSE194nnn/GSE194040/suppl/GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt.gz
# phenotype: TWO series matrices, one per platform (51K + 30K)
curl -L -O $B/series/GSE194nnn/GSE194040/matrix/GSE194040-GPL20078_series_matrix.txt.gz
curl -L -O $B/series/GSE194nnn/GSE194040/matrix/GSE194040-GPL30493_series_matrix.txt.gz

# --- GSE164458 BrighTNess (replication) ---
curl -L -O $B/series/GSE164nnn/GSE164458/suppl/GSE164458_BrighTNess_RNAseq_log2_Processed_ASTOR.txt.gz   # 73M
curl -L -O $B/series/GSE164nnn/GSE164458/matrix/GSE164458_series_matrix.txt.gz                            # 16K

# --- GSE25066 Hatzis (third) ---
# The series matrix carries BOTH the normalised table and the phenotype.
# Do NOT fetch GSE25066_RAW.tar - it is 1.1 GB of CEL files and is not needed.
curl -L -O $B/series/GSE25nnn/GSE25066/matrix/GSE25066_series_matrix.txt.gz                               # 59M

# --- GPL96 annotation, for the GSE25066 probe collapse ---
# NOTE the path: GPL96 buckets as GPLnnn, not GPL96nnn.
curl -L -O $B/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz                                                 # 4.3M
```

The GEO HTTP gateway serves neither `HEAD` nor `Content-Range`, so exact byte
sizes are not obtainable and script 12 guards on a **size floor** plus parsed
dimensions instead of script 02's byte-exact check.

## Landmines, all verified 2026-08-30 and all handled in script 12

**1. Nothing joins on the obvious column.** Measured against the real files:

| Cohort | Key | Transform | Match |
|---|---|---|---|
| GSE194040 | `!Sample_title` | strip `ISPY2_` prefix | 654/654 on GPL20078 |
| GSE164458 | `!Sample_title` | strip `_RNAseq` suffix | 482/482 |

`patient id` looks like the I-SPY2 key and is **not**: it matches 653 of
GPL20078's 654, because one sample's `patient id` disagrees with its own title.
BrighTNess titles are `102001_RNAseq` against expression columns `102001`, so
the raw overlap is **zero**. Script 12 discovers the key by scoring every column
under three declared transformations and stops if none clears 90%.

**2. All three are log2.** GSVA is satisfied; **mitoPPS is not available** -
there is no linear DESeq2-normalised matrix for any of them, two being
microarrays and the third deposited already logged. See script 12's header and
section 6.

**3. GSE194040 is ComBat-adjusted as deposited** (`adjustment method: ComBat
Adjust`). A property of the source, not a choice. State it in Methods.

**4. GSE25066 needs a probe collapse, and loses BBC3 to it.** GPL96 HG-U133A.
The rule, fixed before any model: drop empty symbols, drop `///` multi-mapping
probes, then keep the highest-mean probe per gene. That removes `211692_s_at`,
the **only** BBC3 probe, which also maps to `MIR3190`/`MIR3191`. BBC3 is
therefore absent from this cohort - reported, not patched.

**5. GSE25066 carries four expression-derived published predictors** -
`set_class`, `chemosensitivity_prediction`, `dlda30_prediction`,
`rcb_0_i_prediction`. **Forbidden as covariates or comparators**: each is a
function of the same matrix H4's exposure is built from.

**6. Subtype is not harmonisable as PAM50.** `pam50_class` exists only in
GSE25066. The harmonised variable is **receptor subtype** (HR/HER2): four groups
in GSE194040, constant TNBC in GSE164458, and IHC-derived in GSE25066.

**7. GSE25066 is a pooled meta-cohort** - its `source` field includes `ISPY`,
meaning I-SPY**1**. That is a different trial from I-SPY2 (GSE194040), so there
is no patient overlap between the two, but `source` is carried so the pooling is
visible.

## Endpoint

`pcr = 1` / `RD = 0`, harmonised in script 12; an unrecognised value stops the
script. Marginals recorded at fetch time:

- GSE194040: **319 pCR of 988 (32.3%)** - asserted against D5, which published it
- GSE164458: 236 of 482 (49.0%)
- GSE25066: parsed from `pathologic_response_pcr_rd`

**No score-versus-pCR association is computed in script 12**, by design. See its
header, and D5 section 2, which held the same line during cohort selection.
