# TCGA Clinical Data Resource (CDR)

Curated, quality-controlled clinical outcome endpoints for all 33 TCGA cohorts.
Snapshotted here so that any TCGA clinical-outcome variable used anywhere in this repo
has a single provenanced source.

**Read the scope note below before using this for survival modelling.** The plan
forbids the TCGA survival analysis; this file is here for descriptive work (F4) and
for stage/grade/histology covariates, not to reopen that decision.

**Do not edit in place.** Re-snapshot if a newer CDR release is adopted.

## Provenance

- File: `TCGA-CDR-SupplementalTableS1.xlsx`
- Upstream: NCI GDC PanCanAtlas publication page,
  https://gdc.cancer.gov/about-data/publications/pancanatlas
- Direct URL: `https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81`
- Citation: Liu J, Lichtenberg T, Hoadley KA, et al. "An Integrated TCGA Pan-Cancer
  Clinical Data Resource to Drive High-Quality Survival Outcome Analytics."
  *Cell* 2018;173(2):400-416.e11.
- Snapshot date: 2026-08-29
- Size: 2,945,129 bytes
- SHA-256: `ea594c0fbb6731477c7ac511fab449ca9c38b0d42d269591ed9f5c4090e75a5a`
- MD5: `a4591b2dcee39591f59e5e25a6ce75fa`

## Contents

Eight sheets. Sheet 1 `TCGA-CDR` is the table: **11,160 x 34**, one row per patient
across all cohorts, keyed on `bcr_patient_barcode` (12-character patient barcode) with
`type` as the cohort code.

Columns of interest:

- Endpoints, as event/time pairs: `OS`/`OS.time`, `DSS`/`DSS.time`, `DFI`/`DFI.time`,
  `PFI`/`PFI.time`. Times are in **days**.
- Covariates: `age_at_initial_pathologic_diagnosis`, `gender`, `race`,
  `ajcc_pathologic_tumor_stage`, `clinical_stage`, `histological_type`,
  `histological_grade`, `menopause_status`, `margin_status`, `residual_tumor`,
  `tumor_status`, `vital_status`, `new_tumor_event_*`.
- `Redaction` flags patients the TCGA project redacted. Filter on it.

Sheet 2 `TCGA-CDR_Notes` carries the authors' usage guidance; `ExtraEndpoints`,
`Table4/5_PHAssumptionTests`, `TSS_Info` and `Fig2EFG_AdditionalInfo` are supporting.

`readxl::read_excel()` renames the unnamed first column to `...1` - it is a row index,
not data. Address columns by name.

## BRCA scope note - why this is not a licence to run survival in TCGA

Plan section 3 states, in bold: **"Do not run the survival analysis in TCGA."** The
numbers in this file are the reason. Counted from sheet 1, `type == "BRCA"`, 1,097 rows:

| Endpoint | Events | n with both event and time | Median follow-up (days) |
|---|---|---|---|
| OS  | 151 | 1,096 | 843 |
| DSS |  83 | 1,077 | 843 |
| DFI |  84 |   952 | 764 |
| PFI | 145 | 1,096 | 773 |

Roughly 2.3 years of median follow-up and 145-151 events. That is far too thin for the
three-way interaction the plan specifies, and thinner still once stratified. **F2 and F3
belong in METABRIC and SCAN-B.**

The CDR authors' own recommendation (sheet 2, verbatim): *"For clinical outcome
endpoints, we recommend the use of PFI for progression-free interval, and OS for overall
survival. Both endpoints are relatively accurate. Given the relatively short follow-up
time, PFI is preferred over OS."* If a TCGA sensitivity analysis is ever run, **PFI is
the endpoint**, and it is a sensitivity analysis, not a result.

Note also, from sheet 2: BRCA has neither `treatment_outcome_first_course` nor
`residual_tumor`, so its `DFI` was derived from `margin_status` alone.

## What this is legitimately used for here

- **F4** (plan section 10) - descriptive association of `STATE` with stage, grade,
  histological type and nodal status, following the Menegollo Fig 7A-C template.
  Descriptive and supplementary; no survival model.
- Stage and grade covariates, where the existing `data/tcga_clinical/` snapshot (8
  columns: PAM50, ER/PR/HER2, aneuploidy, FGA) does not reach.
