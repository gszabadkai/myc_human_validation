# TCGA-BRCA clinical snapshot (subtype, receptor status, aneuploidy burden)

Snapshotted 2026-08-28 for gate G2. Small enough to commit, so unlike the GISTIC
files these are on `origin`.

Everything here comes from **live web services**. A live call at analysis time
would silently drift, which is not acceptable under this arm's pre-registration
discipline. Same reasoning as the CollecTRI snapshot in `data/collectri_human/`:
take a dated copy, record the exact query, consume the copy.

## Files

| File | Rows | SHA-256 |
|---|---|---|
| `gdc_brca_cases_2026-08-28.tsv` | 1,098 | `15aa46890064814a18d74d9950694ad3ae75df81d33c706b44f389ed12c3d9d9` |
| `tcga_brca_clinical_snapshot_2026-08-28.tsv` | 1,098 | `40fd740362f5aef02920a7eae9abe393e6669f8deea88d5fd99de9b5385ec181` |

### `gdc_brca_cases_2026-08-28.tsv`

The authoritative TCGA-BRCA case list. Two columns, `id` (GDC UUID) and
`submitter_id` (patient barcode). Used to define BRCA membership when subsetting
the pan-cancer ISAR GISTIC file, so that membership does not depend on either
GISTIC run. See `data/gistic_tcga_brca/README.md`.

Query, verbatim:

```
GET https://api.gdc.cancer.gov/cases
  filters={"op":"in","content":{"field":"project.project_id","value":["TCGA-BRCA"]}}
  fields=submitter_id
  size=2000
  format=TSV
```

### `tcga_brca_clinical_snapshot_2026-08-28.tsv`

One row per GDC BRCA case (1,098), left-joined to the attributes below. Missing
values are the literal string `NA`.

| Column | Source study | Non-missing |
|---|---|---|
| `patient_barcode` | GDC | 1,098 |
| `PAM50_SUBTYPE` | `brca_tcga_pan_can_atlas_2018` (patient) | 981 |
| `ER_STATUS_BY_IHC` | `brca_tcga` (patient) | 1,048 |
| `PR_STATUS_BY_IHC` | `brca_tcga` (patient) | 1,047 |
| `IHC_HER2` | `brca_tcga` (patient) | 919 |
| `HER2_FISH_STATUS` | `brca_tcga` (patient) | 421 |
| `ANEUPLOIDY_SCORE` | `brca_tcga_pan_can_atlas_2018` (**sample**) | 1,041 |
| `FRACTION_GENOME_ALTERED` | `brca_tcga_pan_can_atlas_2018` (**sample**) | 1,068 |

Queries, verbatim (one call per attribute):

```
GET https://www.cbioportal.org/api/studies/<STUDY>/clinical-data
  ?clinicalDataType=<PATIENT|SAMPLE>
  &attributeId=<ATTRIBUTE>
  &projection=SUMMARY
  &pageSize=3000
```

Note the two studies. PAM50 and the aneuploidy metrics come from the PanCanAtlas
study; **receptor IHC is not in it** (its 60 attributes contain no ER/PR/HER2),
so those come from the legacy `brca_tcga` study, which carries 138 attributes
including the full IHC and FISH panel.

`ANEUPLOIDY_SCORE` and `FRACTION_GENOME_ALTERED` are **sample-level**, not
patient-level; querying them with `clinicalDataType=PATIENT` returns zero rows
without erroring. They were re-queried at sample level and keyed back to patient.
All BRCA sample ids in the response end `-01`, so the collapse is unambiguous.

## Two definitions must be pre-specified before use

### PAM50

`PAM50_SUBTYPE` values are `BRCA_LumA` (499), `BRCA_LumB` (197), `BRCA_Basal`
(171), `BRCA_Her2` (78), `BRCA_Normal` (36). This is the PanCanAtlas subtype
call, the same table `TCGAbiolinks::PanCancerAtlas_subtypes()` returns; taking it
from cBioPortal avoids installing TCGAbiolinks for G2.

`BRCA_Normal` (normal-like) is an artefact-prone class usually reflecting low
tumour cellularity. Decide explicitly whether it is a stratum or excluded. Do not
let it fall silently into a reference level.

A free cross-check exists: the G3 file
`TCGA.all.biclusters.RNAseq.Rdata` carries PAM50 for its 849 samples
(173 Basal / 123 Her2 / 237 LumA / 250 LumB / 66 Normal). Concordance on the
overlap is worth reporting once.

### TNBC

There is no TNBC column. It must be derived, and the derivation rule is a
decision, not an implementation detail. Rule fixed 2026-08-28:

- ER = `ER_STATUS_BY_IHC`, PR = `PR_STATUS_BY_IHC`; `Positive` / `Negative`
  only, `Indeterminate` is uncallable.
- HER2 = `HER2_FISH_STATUS` where it is `Positive` or `Negative`; otherwise
  `IHC_HER2`.
- **IHC-equivocal HER2 with no FISH result is uncallable, not negative.**
- TNBC = ER `Negative` AND PR `Negative` AND HER2 `Negative`.

Under that rule: **951 patients callable on all three, 161 TNBC.**

The attrition from 1,049 to 951 is driven by the 179 IHC-equivocal cases lacking
FISH. Calling those negative would inflate TNBC and is the reason the rule is
written down here rather than decided inside a script.

161 IHC-TNBC against 171 PAM50-Basal. The two sets overlap heavily but are not
identical, which is expected; report the cross-tabulation rather than treating
them as interchangeable. Among the 951 callable patients:

| PAM50 | TNBC | non-TNBC |
|---|---|---|
| `BRCA_Basal` | 117 | 32 |
| `BRCA_Her2` | 13 | 53 |
| `BRCA_LumA` | 4 | 429 |
| `BRCA_LumB` | 1 | 166 |
| `BRCA_Normal` | 9 | 21 |
| no PAM50 call | 17 | 89 |

So 44 of 161 TNBC tumours are not Basal, and 32 of 149 callable Basal tumours
are not TNBC. The divergence is large enough that "basal / TNBC" in plan section
6 cannot be treated as one variable. G2 reports both, separately.

## Read-time traps

- Values arrive as strings throughout, including `ANEUPLOIDY_SCORE` (integer-
  like) and `FRACTION_GENOME_ALTERED` (double-like). Coerce explicitly.
- `HER2_FISH_STATUS` includes the literal values `[Not Evaluated]`,
  `Indeterminate` and `Equivocal`. Only `Positive` and `Negative` are usable.
- Missing values are the literal `NA` string, not empty. Read with
  `na.strings = "NA"` or coerce deliberately.
- These are **patient** barcodes. The GISTIC files are aliquot-level. Join only
  after applying the `-01` sample-type filter described in
  `data/gistic_tcga_brca/README.md`.

---

## `gdc_brca_rnaseq_aliquots_2026-08-28.tsv`

Added 2026-08-28 to recover **sequencing plate**, a pre-specified Block C
covariate (plan section 8). Plate is field 6 of the ALIQUOT barcode, and the
Xena expression matrix carries only sample-level ids, so it cannot be recovered
from the expression data alone.

| | |
|---|---|
| Rows | 1,231 (1,111 primary tumour, 113 normal, 7 metastatic) |
| Columns | `sample_id`, `patient`, `aliquot_barcode`, `plate`, `sample_type` |
| Distinct plates, primary | 43 |
| SHA-256 | `703ac03b6bb997dad0b67c28f842b99598d3fbc1d57ec6c1a78bb66e686128ad` |

Query, verbatim:

```
GET https://api.gdc.cancer.gov/files
  filters={"op":"and","content":[
    {"op":"in","content":{"field":"cases.project.project_id","value":["TCGA-BRCA"]}},
    {"op":"in","content":{"field":"data_type","value":["Gene Expression Quantification"]}},
    {"op":"in","content":{"field":"analysis.workflow_type","value":["STAR - Counts"]}}]}
  fields=cases.samples.portions.analytes.aliquots.submitter_id,
         cases.samples.submitter_id,cases.samples.sample_type
  size=2000  format=TSV
```

`cases.samples.submitter_id` is sample-level (`TCGA-A8-A09E-01A`), which is
exactly the Xena column format, so the join to script 01's `sample_map` is direct.

### It independently confirms the Xena averaging

Five sample ids carry **two** RNA-seq aliquots each, the second re-sequenced on
plate `A277`:

```
TCGA-A7-A0DB-01A   A00Z, A277
TCGA-A7-A13D-01A   A12P, A277
TCGA-A7-A13E-01A   A12P, A277
TCGA-A7-A26E-01A   A169, A277
TCGA-A7-A26J-01A   A169, A277
```

These are **exactly** the five columns script 01 flags as Xena aliquot averages,
identified there by an entirely different route - the half-integer signature in
the inverted log2 values. Two independent lines of evidence, same five samples.
See `data/tcga_expression_xena/README.md`.
