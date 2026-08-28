# TCGA-BRCA STAR counts, via UCSC Xena

Provenance for the expression matrix consumed by script 01. Snapshotted
2026-08-28.

**The data file is not here.** It is 138 MB and lives at
`data/raw/tcga_expression_xena/TCGA-BRCA.star_counts.tsv.gz`, which is gitignored
and not on `origin`. This README is the committed record; the URL plus the
SHA-256 make it recoverable deterministically.

## File

| | |
|---|---|
| File | `data/raw/tcga_expression_xena/TCGA-BRCA.star_counts.tsv.gz` |
| URL | https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BRCA.star_counts.tsv.gz |
| Hub | UCSC Xena GDC hub (`gdc-hub`) |
| Retrieved | 2026-08-28 |
| Size | 138,503,699 bytes |
| SHA-256 | `058d79121460c535b73312247a55d3108d18ea6ddd0ccc1070b5dd93ceeedaa4` |
| Dimensions | 60,660 genes x **1,226 samples** |
| Row ids | versioned Ensembl (`ENSG00000000003.15`), **GENCODE v36** |
| Column ids | TCGA sample barcodes with vial (`TCGA-B6-A1KC-01B`) |
| Values | **`log2(count + 1)`**, not raw counts |

## Why Xena and not the GDC API

The plan (section 8) says "GDC harmonised STAR counts via `TCGAbiolinks` (or
recount3)". The GDC per-file route was implemented first and abandoned as
impractical: **5.2 GB across 1,231 files**, measured at 0.05 files/s, i.e. about
**7 hours**, over a transport that failed twice (R's 60-second default download
timeout, then repo-root pollution from TCGAbiolinks writing to `getwd()`).

Xena re-hosts the **same GDC STAR-Counts quantification** as a single file:
same pipeline, same GENCODE v36, same versioned ENSG ids. 138 MB, one request.

What it costs, stated rather than waved away:

1. **Provenance is one hop longer** - Xena's snapshot of a GDC release rather
   than the GDC API directly. Mitigated by this README and the SHA-256, the same
   discipline used for CollecTRI and GISTIC.
2. **Values arrive log2-transformed** and must be inverted. See below.
3. **Xena has already collapsed some aliquots** upstream, by a rule not stated
   in its documentation. See "The averaged columns", which is the important part
   of this file.

## Verification: the inversion is exact

Raw counts are recovered as `round(2^x - 1)`. This was **verified against the
GDC files themselves**, not assumed.

42 GDC per-file downloads survive from the abandoned attempt, under
`data/raw/tcga_brca_expression/`. Four were compared against the inverted Xena
matrix on a **blind match** - no barcode was supplied, the check simply asked
which Xena column, if any, equals the GDC sample exactly:

```
2f51534b-248   EXACT MATCH -> TCGA-AC-A62X-01A   (4,000 genes)
75668632-f0c   EXACT MATCH -> TCGA-E2-A10E-01A   (4,000 genes)
ad01e71e-ed1   EXACT MATCH -> TCGA-E9-A1NA-01A   (4,000 genes)
36e48370-2ec   EXACT MATCH -> TCGA-E9-A1NH-01A   (4,000 genes)
```

Each identified exactly one column, with zero total absolute difference. The
Xena values are the GDC `unstranded` raw counts, log2(x+1)-transformed.

**Keep at least a few of those GDC `.tsv` files.** Script 01 also takes its gene
annotation from one of them (`gene_id` / `gene_name` / `gene_type` at the
matching GENCODE v36), because Xena ships only ENSG ids and its published
probeMap is GENCODE v22 against v36 data.

## The averaged columns - the thing that would have gone wrong

Inverting the whole matrix does **not** give integers everywhere. Deviations
from the nearest integer are either ~0 or **exactly 0.5**, never in between:

```
values examined            4,904,000   (4,000 genes x 1,226 samples)
deviation > 1e-6               9,729   (0.198%)
deviation > 0.49               9,729   (0.198%)   <- identical count
```

That bimodality is a signature, not floating-point error. The half-integers are
confined to **5 columns of 1,226**, and in each one ~49% of genes are
half-integers - exactly what averaging two integers produces:

| Column | half-integer genes | other columns for that patient |
|---|---|---|
| `TCGA-A7-A26J-01A` | 1973 (49.3%) | `-01B` |
| `TCGA-A7-A13E-01A` | 1967 (49.2%) | `-01B`, `-11A` |
| `TCGA-A7-A26E-01A` | 1965 (49.1%) | `-01B` |
| `TCGA-A7-A13D-01A` | 1921 (48.0%) | `-01B` |
| `TCGA-A7-A0DB-01A` | 1903 (47.6%) | `-01C`, `-11A` |

**Xena has averaged replicate aliquots into those five `-01A` columns.** They are
not raw counts and must not reach DESeq2.

Every affected column is a `-01A` for a patient who also has a clean `-01B` or
`-01C`. So the obvious tie-break - "keep the lowest vial letter" - would have
selected the averaged column and discarded the clean one, **for all five
patients, silently**. That rule was written before this check and was wrong.

Script 01 therefore:

1. detects averaged columns by testing for half-integers after inversion, rather
   than hard-coding the five barcodes, so a future Xena release that averages
   different samples is caught;
2. drops an averaged column when the patient has a clean alternative;
3. keeps and **reports** an averaged column only if a patient has no clean one
   (no patient is in that position in this release);
4. applies the lowest-vial rule only to what survives;
5. hard-stops if any non-integrality is found *outside* the averaged columns,
   since that would mean the values are not `log2(raw count + 1)` at all.

## Read-time traps

- Values are `log2(count + 1)`. Do not feed them to DESeq2 directly, and do not
  log them again for GSVA - script 01 inverts to counts and then derives both a
  VST (log) and a size-factor-normalised linear matrix from them.
- Column ids are sample-level **with vial**, not patient ids and not aliquot
  barcodes. Truncate to three fields to key on patient, but only after the
  averaged-column and vial rules above.
- 1,226 columns include normal-tissue samples (`-11A`). Script 01 keeps sample
  type `-01` only.
- Row ids carry the GENCODE version suffix. Do not strip it before joining to
  the GDC annotation; both sides carry it and the join is exact.
