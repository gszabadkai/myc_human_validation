# TCGA PanCanAtlas genomics and protein

Provenance for the script 02 inputs. Snapshotted 2026-08-28.

**The data files are not here.** They live under
`data/raw/tcga_pancanatlas/`, which is gitignored and not on `origin`. This
README is the committed record; the URLs plus SHA-256 sums make them recoverable
deterministically.

All five are **single-file downloads from the GDC API**, open access, no token.
This is deliberate after the script 01 experience: the GDC *per-file* route
(1,231 files) was abandoned as impractical, but single large files over the same
API are fine as long as R's 60-second default download timeout is raised.

## Files

| File | Size | SHA-256 | GDC UUID |
|---|---|---|---|
| `mc3.v0.2.8.PUBLIC.maf.gz` | 753 MB | *see below* | `1c8cfe5f-e52d-41ba-94da-f15ea1337efc` |
| `TCGA_mastercalls.abs_tables_JSedit.fixed.txt` | 0.9 MB | `f430a975433d82e0098d7405619d4f12a0c765fcd97e7d63cc9b1de7f2d763cd` | `4f277128-f793-4354-a13d-30cc7fe9f6b5` |
| `TCGA_all_leuk_estimate.masked.20170107.tsv` | 0.6 MB | `5a8268caedbf8dc98a75be0528d583238d7355761d9fc746e42002f223a982d9` | `6f75c9d7-5134-4ed1-b8f3-72856c98a4e8` |
| `TCGA-RPPA-pancan-clean.txt` | 18.9 MB | `06246573836865589134bd9424189f81b0d9fb436fcbf5e72024225442c400de` | `fcbb373e-28d4-4818-92f3-601ede3da5e1` |
| `PANCAN_ArmCallsAndAneuploidyScore_092817.txt` | 1.1 MB | `c8a7f3dd059a5dc3b66b85ba53c0a6d0cd6c20be71e72b843b0eff7318e8a23c` | `4c35f34f-b0f3-4891-8794-4840dd748aad` |

Download pattern: `https://api.gdc.cancer.gov/data/<UUID>`

Source pages:
- MC3, ABSOLUTE, arm calls: https://gdc.cancer.gov/about-data/publications/pancanatlas
  and https://gdc.cancer.gov/about-data/publications/pancan-aneuploidy
- Leukocyte fraction, RPPA: https://gdc.cancer.gov/about-data/publications/panimmune

## Why PanCanAtlas throughout

Xena hosts a BRCA-only somatic mutation file at 2.5 MB, against MC3's 753 MB,
and it was considered. It was **rejected**: it is a per-project GDC call set,
whereas MC3 is the PanCanAtlas consensus across seven callers. The copy-number
(ISAR GISTIC), purity (ABSOLUTE), immune (Thorsson) and aneuploidy inputs are
all PanCanAtlas already, and mixing a differently-called MAF into that set buys
nothing but inconsistency. Size is a one-off cost; the subset is cached.

## Contents and read-time traps

### `mc3.v0.2.8.PUBLIC.maf.gz`

Pan-cancer, ~3.6M calls. Script 02 reads five columns, subsets to the BRCA
patients script 01 kept, and caches the result as `mc3_brca_subset.rds`, so the
753 MB file is parsed once.

- **Not every call is usable.** MC3 ships non-PASS calls (`wga`,
  `native_wga_mix`, `broad_PoN_v2`, ...). Script 02 uses `FILTER == "PASS"` only.
- Barcodes are aliquot-level; truncate to three fields for patient.
- **Absence of a call is not wild-type** if the patient was never sequenced.
  Script 02 marks patients with no MC3 coverage as `NA`, not `FALSE`.

### `TCGA_mastercalls.abs_tables_JSedit.fixed.txt`

ABSOLUTE purity/ploidy. Columns contain **spaces** (`call status`,
`Genome doublings`); script 02 runs `make.names()` on load. Keyed on `sample`,
aliquot-level.

### `TCGA_all_leuk_estimate.masked.20170107.tsv`

**Headerless**, three columns: cancer type, aliquot barcode, leukocyte fraction.
Multiple aliquots per patient occur; script 02 averages them, since this is a
continuous estimate rather than a call.

### `TCGA-RPPA-pancan-clean.txt`

Level-4 normalised, ~200 proteins. **Carried through as supplied - not
re-normalised, not logged, and never mixed with the expression matrices from
script 01.**

Relevant to Block E, confirmed present in the panel:

```
BAK  BAX  BCL2  BCLXL  BID  BIM  CMYC  FOXO3A  FOXO3A_pS318S321
```

`FOXO3A_pS318S321` is the AKT phosphorylation site, so **Block E can test H2's
FOXO3 nuclear-exclusion mechanism at protein level directly** rather than
inferring it from regulon activity.

**`PUMA`/`BBC3` is NOT on the panel.** So Block E can confirm the BCL-XL, MYC and
FOXO3 arms at protein level but **cannot confirm the `PRIME` ratio itself**,
whose numerator is PUMA. State this in Methods rather than letting the reader
assume protein-level confirmation of the endpoint.

### `PANCAN_ArmCallsAndAneuploidyScore_092817.txt`

Taylor et al. 2018, the authoritative version of the aneuploidy burden used to
condition G2. The `ANEUPLOIDY_SCORE` already snapshotted in
`data/tcga_clinical/` came via cBioPortal from the same study; this file
additionally carries **per-arm calls**, so 8q and 1q can be inspected directly -
which matters because G2's confound is arm-level gain at exactly those arms.

Column names contain spaces and arm names start with digits (`1p`, `8q`);
`make.names()` turns these into `X1p`, `X8q`.

## Pre-registration note

Script 02 fixes the definition of "mutated" before any model is fitted:
`FILTER == "PASS"` plus a non-silent `Variant_Classification`. `PTEN_ALTERED`
additionally requires the GISTIC deep-deletion half, joined in script 03,
because PTEN is lost by deletion at least as often as by mutation. Full
reasoning is in the script header. **That rule needs sign-off before script 02
is run.**
