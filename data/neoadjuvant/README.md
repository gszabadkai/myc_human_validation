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


## Checksums

Downloaded 2026-08-30, 192 MB total. All seven pass `gzip -t`.

| File | Bytes | SHA-256 |
|---|---|---|
| `GPL96.annot.gz` | 4,522,748 | `88e0b22362bac779eb220b3b185c80faa6510a92b9358eaad159a561ab4351c4` |
| `GSE164458_BrighTNess_RNAseq_log2_Processed_ASTOR.txt.gz` | 76,442,245 | `7209cd0e5528e1fdb1fb6ffe69f91d7ff7e666df6a9453664a6a78dc6f8bd9e8` |
| `GSE164458_series_matrix.txt.gz` | 15,956 | `8e80fd0c211f8eedfe6aad4538b11f02b93a17992b64f3a3ab5d0a4f1ba57d24` |
| `GSE194040-GPL20078_series_matrix.txt.gz` | 52,000 | `8a7e5d56cb7523812a44b58daa7236b4f09a726cc028d0ab0b4384a759702a98` |
| `GSE194040-GPL30493_series_matrix.txt.gz` | 30,724 | `20aa7bd22d2352a4d96dfbd381b7100443be10dedfbab56223a786b96b80cb19` |
| `GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt.gz` | 49,349,197 | `9d8d1d35dc7f9ccee367cf37ae3403335eca3d7665a498e3037f4a4614b972c4` |
| `GSE25066_series_matrix.txt.gz` | 61,627,399 | `beef319fba17f38e05ff0994513763dfbc7c9cfa9f6c4b2b7c31dcdd036060fb` |

**Script 12 guards on a size FLOOR, not on these exact bytes, and that is
deliberate.** GEO *regenerates series matrices* - the three here carry 2026-07
and 2026-08 timestamps against source data from 2011-2022 - so a byte-exact
guard would fail on a harmless upstream rebuild. The two supplementary
expression files are static (2021 and 2022) and their checksums should not
move. If a series-matrix checksum changes, re-run script 12 and check that the
asserted dimensions and the 319 I-SPY2 pCR count still hold; if they do, the
rebuild was cosmetic and this table should be updated.

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
