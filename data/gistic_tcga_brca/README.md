# TCGA-BRCA GISTIC2 thresholded copy-number calls

Provenance for the gate G2 copy-number inputs. Snapshotted 2026-08-28.

**The data files themselves are not here.** They are large and live under
`data/raw/gistic_tcga_brca/`, which is gitignored and not on `origin`. This
README is the committed record. Both files are openly downloadable with no
account, so the URLs plus the SHA-256 sums below make them recoverable
deterministically. Re-download and verify rather than hunting for a backup.

Two runs were taken deliberately, per decision on 2026-08-28:

- **A = primary**, PanCanAtlas ISAR-corrected
- **B = pre-specified sensitivity**, Broad Firehose BRCA-only

They are not interchangeable. See "Why two" below.

---

## A. PanCanAtlas ISAR-corrected GISTIC2 (PRIMARY)

| | |
|---|---|
| File | `data/raw/gistic_tcga_brca/ISAR_GISTIC.all_thresholded.by_genes.txt.gz` |
| Source page | https://gdc.cancer.gov/about-data/publications/pancan-aneuploidy |
| Download URL | https://api.gdc.cancer.gov/data/a9dae2ab-9462-4f9a-9730-5bad520ab2d7 |
| GDC UUID | `a9dae2ab-9462-4f9a-9730-5bad520ab2d7` |
| Listed as | "ISAR-corrected GISTIC all_thresholded.by_genes" |
| Retrieved | 2026-08-28 |
| Size (gz) | 37,966,065 bytes |
| SHA-256 (gz) | `bbd82167738f5c196c56eeded12d9b2e10867858f42129b6d8bd9827b8890331` |
| Size (raw) | 536,214,210 bytes |
| SHA-256 (raw) | `dde3270dec9dc35b46ccc90d7674538f710ff2dc91895dc392ba45c121eb02a3` |
| Inner filename | `ISAR_GISTIC.all_thresholded.by_genes.txt` (gzip header, mtime 2017-12-11) |
| Dimensions | 24,203 genes x 9,991 samples, **pan-cancer (all 33 TCGA types)** |
| Publication | Taylor AM et al. 2018, Cancer Cell 33:676 (PanCanAtlas aneuploidy) |

Requires subsetting to BRCA. See "Defining the BRCA sample set".

## B. Broad Firehose BRCA-TP GISTIC2 Level 4 (SENSITIVITY)

| | |
|---|---|
| File | `data/raw/gistic_tcga_brca/gdac.broadinstitute.org_BRCA-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0.tar.gz` |
| Download URL | https://gdac.broadinstitute.org/runs/analyses__2016_01_28/data/BRCA-TP/20160128/gdac.broadinstitute.org_BRCA-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0.tar.gz |
| Firehose run | `analyses__2016_01_28`, stamp `2016012800.0.0` |
| Retrieved | 2026-08-28 |
| Size | 53,856,803 bytes |
| SHA-256 | `2ee43313d90a5c9c1efbd8f4d33103ddb85effd4014957227bed226b305ce597` |
| Member of interest | `<dir>/all_thresholded.by_genes.txt`, 60,251,425 bytes |
| Dimensions | 24,776 genes x 1,080 samples, **BRCA primary tumours only** |

The tarball also carries `all_lesions.conf_99.txt`, `amp_genes.conf_99.txt`,
`del_genes.conf_99.txt`, `focal_data_by_genes.txt`, `broad_values_by_arm.txt`
and `sample_cutoffs.txt`, none of which G2 uses. Do not silently switch to
`focal_data_by_genes.txt`; it is a different quantity.

---

## Why two, and why ISAR is primary

G2 asks whether MYC (8q24.21) amplification co-occurs with MCL1 (1q21.3) or
BCL2L1 (20q11.21) amplification **more than expected**.

8q and 1q are two of the most frequently gained arms in breast cancer. In the
Firehose calls, MCL1 is gained or amplified in 72% of tumours and MYC in 63%.
A co-occurrence test on uncorrected calls therefore largely measures aneuploidy:
tumours that gain arms gain both of these arms. ISAR (In Silico Admixture
Removal) corrects for arm-level and whole-genome ploidy background, which is
exactly the confound in play, so it is the primary source.

Firehose is retained as a **pre-specified sensitivity**, fixed before any
co-occurrence statistic was computed, because it is BRCA-native (GISTIC run on
breast alone rather than pan-cancer) and because its larger sample set is a
useful check that the ISAR result is not an artefact of the 1,043-sample subset.

ISAR alone is not sufficient. G2 must additionally condition on aneuploidy
burden (`ANEUPLOIDY_SCORE`, `FRACTION_GENOME_ALTERED`; see
`data/tcga_clinical/`). Without it, a confident positive is expected regardless
of biology.

---

## Defining the BRCA sample set

Membership is defined by the **GDC case list**, not by either GISTIC file, so
that A and B are subset by an independent authority rather than one defining the
other. The list is snapshotted at
`data/tcga_clinical/gdc_brca_cases_2026-08-28.tsv` (1,098 cases).

Rule, fixed 2026-08-28:

1. Keep columns whose patient barcode (`TCGA-XX-YYYY`, first three fields) is in
   the GDC TCGA-BRCA case list.
2. Keep **sample type `-01` only** (primary solid tumour).
3. Assert one column per patient after step 2.

Verified coverage under that rule:

| | A (ISAR) | B (Firehose) |
|---|---|---|
| Columns before filtering | 9,991 (pan-cancer) | 1,080 |
| Sample types present in BRCA subset | 1,043 x `-01`, 7 x `-06` | 1,080 x `-01` |
| Distinct BRCA patients after rule | **1,043** | **1,080** |
| GDC BRCA cases with no data | 55 | 18 |

The 7 duplicated patients in ISAR are each a primary (`-01`) plus a metastasis
(`-06`). Step 2 removes them without needing an arbitrary tie-break. **Do not
dedupe by picking the first column** - that silently keeps metastases for some
patients.

A's patients are a strict subset of B's: `A n B` = 1,043, A-only = 0,
B-only = 37. So the sensitivity comparison is well defined. Report n for each
run, and additionally run B restricted to the 1,043 intersection so that the
comparison isolates the calling method rather than the sample set.

---

## Genes of interest, verified present in both

Cytobands agree with the plan (section 6) in both files.

| Gene | Cytoband | Role |
|---|---|---|
| `MYC` | 8q24.21 | exposure, M-c CNV instrument |
| `MCL1` | 1q21.3 | BUFFER |
| `BCL2L1` | 20q11.21 | BUFFER, and the `PRIME` denominator |
| `BBC3` | 19q13.32 | `PRIME` numerator (PUMA) - descriptive only, see below |
| `TP53` | 17p13.1 | H3 stratum |
| `PIK3CA` | 3q26.32 | H2 stratum |
| `PTEN` | 10q23.31 | H2 stratum |

Also present and pulled for the endpoint negatives: `BCL2` (18q21.33),
`BID` (22q11.21), `BAX` (19q13.33), `BCL2L11` (2q13), `BAK1` (6p21.31),
`PMAIP1` (18q21.32).

### Observed distributions (ISAR, BRCA subset, n = 1,050 columns before the `-01` filter)

| Gene | -2 | -1 | 0 | +1 | +2 | amp (+2) |
|---|---|---|---|---|---|---|
| MYC | 0 | 28 | 343 | 452 | 227 | 21.6% |
| MCL1 | 0 | 23 | 233 | 618 | 176 | 16.8% |
| BCL2L1 | 0 | 63 | 505 | 454 | 28 | **2.7%** |
| BBC3 | **5** | 245 | 545 | 244 | 11 | 1.0% |

Two consequences, both recorded as decisions rather than left to the script:

- **`BCL2L1` high-level amplification is rare (2.7%).** A `BUFFER` defined as
  "MCL1 amp OR BCL2L1 amp" at threshold +2 is overwhelmingly MCL1. The amplitude
  threshold changes the result and must be pre-registered, not chosen after
  seeing the co-occurrence table.
- **`BBC3` homozygous deletion occurs in 5 of 1,050 tumours.** It is reported
  descriptively and is **not** powered as a co-occurrence test; a null on n = 5
  would be uninterpretable. Its shallow-loss frequency (23%) tracks BAX at
  19q13.33 (226), i.e. it is regional 19q behaviour and not focal to BBC3.

## Read-time traps

- The header's first three fields are `Gene Symbol`, `Locus ID`, `Cytoband`;
  sample columns start at column 4.
- Values are integers in `{-2,-1,0,1,2}` and arrive as **character** from
  `read.delim` unless coerced. Coerce explicitly and assert the value set.
- Some `Locus ID` values are negative (e.g. `PIK3CA` is `-1426`). That is a
  Firehose/GISTIC annotation quirk, not corruption. Do not filter on it.
- ISAR is pan-cancer and 536 MB uncompressed. Read it with
  `data.table::fread(cmd = "gzip -dc ...")` selecting only the needed columns,
  or stream it. Do not `read.delim` the whole thing.
- Barcodes are aliquot-level in the GISTIC files and patient-level in the
  clinical snapshot. Truncate to three fields to join, after applying the `-01`
  filter.
