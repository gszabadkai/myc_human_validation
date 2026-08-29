# DepMap / CCLE - Block G inputs

Provenance for the DepMap files consumed by `scripts/14_depmap_dependency.R`
(Block G, plan section 10). The files themselves live in
`data/raw/depmap/`, which is **gitignored** - like every other large download in
this repo, they are re-fetchable from the URLs below and are not on `origin`.

**Release: DepMap Public 26Q1** (release notes snapshotted at
`docs/README_depmap-public-26Q1.txt`). `Model.csv`, the expression matrix and
`CRISPRGeneEffect.csv` must all come from ONE release: Model IDs are stable
across releases but the line set and the Chronos scaling are not, and mixing
releases produces a join that looks complete and is not. If you move to a newer
release, change `DEPMAP_RELEASE` in script 14, re-fetch **all three**, and record
the new checksums here.

## Files

| File | Release | Required | Used for |
|---|---|---|---|
| `Model.csv` | 26Q1 | yes | lineage assignment (`OncotreeLineage == "Breast"`) |
| `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | 26Q1 | yes | MYC / OXPHOS scoring. **log2(TPM+1)** |
| `CRISPRGeneEffect.csv` | 26Q1 | yes | Chronos gene effect, Block G item 2 |
| `Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv` | 24Q2 | optional | Block G item 3, drug sensitivity |
| `Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv` | 24Q2 | optional | compound name -> column id |

Without the two PRISM files, script 14 **skips** item 3 and says so. It does not
substitute anything for it.

## Acquisition

Established by enumerating the figshare API and reading the 26Q1 release notes
on 2026-08-30, not assumed:

- **figshare's latest DepMap mirror is 24Q4 (December 2024).** It never mirrored
  25Qx or 26Q1, and even for the releases it does carry it holds only the
  Achilles/CRISPR half - checked across 24Q4 (`27993248`), 24Q2 (`25880521`) and
  earlier. **So for the current release the portal is the only source.**
- **The DepMap portal is behind a Cloudflare challenge.** A `curl` of any portal
  download URL returns a ~5 KB HTML verification page, not data. Saved under a
  `.csv` name it fails much later and very confusingly, so script 14 size-checks
  every required input on load. The matching `storage.googleapis.com` paths
  return 403. Use a browser; that is the intended route.
- **PRISM Repurposing is a separate release and is complete on figshare**
  (`25917643`), including both files script 14 wants.

### Browser - the three DepMap files

https://depmap.org/portal/data_page/?tab=allData , release **26Q1**:

- `Model.csv`
- `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`
- `CRISPRGeneEffect.csv`

### Scriptable - the two PRISM files

```
cd data/raw/depmap
curl -L -o Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv \
  https://ndownloader.figshare.com/files/46630984
curl -L -o Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv \
  https://ndownloader.figshare.com/files/46630981
```

The numeric file ids are pinned deliberately: resolving them through
`/v2/articles/25917643/files` works too, but a pinned id is a pinned file, and a
figshare version bump would otherwise change the data under a stable command.
DOI for citation: `10.25452/figshare.plus.25917643`.

## What changed between 24Q4 and 26Q1 - three breaking changes

Worth recording, because script 14 was first written against 24Q4 and two of
these would have failed silently rather than erroring.

1. **The expression file was renamed.**
   `OmicsExpressionProteinCodingGenesTPMLogp1.csv` (24Q4) is now
   `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` (26Q1).
2. **Its shape changed: it is now ONE ROW PER SEQUENCING RUN, not per model.**
   In 24Q4 it was one row per model with the ModelID as an unnamed first column.
   Read the old way, the metadata columns would be taken for genes, and any
   model with more than one run would appear more than once - silently
   duplicating lines and inflating n. Script 14 filters on the model-level
   default flag and asserts one row per ModelID afterwards.
3. **There are now stranded and unstranded variants.** Script 14 uses the
   **unstranded** file (`EXPR_STRANDED <- FALSE`), because it is continuous with
   every earlier DepMap release and with the published CCLE analyses this block
   is read beside. The script reports the line count, so if the stranded file
   turns out to cover materially more lines that is checkable rather than
   assumed. Do not mix them; whichever is used is recorded in the saved object.

## The 26Q1 release notes misdescribe their own expression file

Found on 2026-08-30, when script 14's column guard fired. The notes
(`docs/README_depmap-public-26Q1.txt`) say the columns are *"ProfileID,
is_default_entry, ModelID and then a column per gene"*. The shipped file begins:

```
<unnamed 0-based index>, SequencingID, ModelConditionID, ModelID,
IsDefaultEntryForMC, IsDefaultEntryForModel, then "SYMBOL (ENTREZ)" columns
```

**There is no `ProfileID` and no `is_default_entry`, and the flags are the
strings `"Yes"`/`"No"`, not `TRUE`/`FALSE`.** Three consequences, all handled:

- Script 14 locates columns **by name with fallbacks** and decides what is a
  gene by testing that the column is **numeric**, not by matching a name
  pattern. A newly added metadata column stops the script rather than becoming a
  gene.
- The flag reader accepts logical, numeric, `Yes/No`, `True/False`, `1/0`, and
  **refuses anything else rather than guessing**. This repo has already been
  bitten once by testing `== "True"` against a column parsed as logical, where
  the result was a silent all-FALSE - and here a `"Yes"`/`"No"` file would have
  produced exactly that.
- It filters on **`IsDefaultEntryForModel`, not `IsDefaultEntryForMC`**. The
  latter is the model-*condition* default, a different and larger set; in 26Q1
  exactly one row differs between them.

Verified on the real file: 1,775 rows collapse to **1,719 models**, 19,215 gene
columns, values 0.000-17.361 (consistent with log2(TPM+1)), no NAs.

## Line counts in 26Q1, measured

| | n |
|---|---|
| Breast models in `Model.csv` | 96 |
| Breast with expression | 71 |
| Breast with CRISPR | 53 |
| **Breast with both - the Block G primary** | **51** |
| **Pan-cancer with both - the powered secondary** | **1,140** |

51 lines is above script 14's `MIN_LINES` floor of 25 but is genuinely small for
a two-way interaction, which is exactly why the lineage-adjusted pan-cancer fit
is reported alongside.

## Is the PRISM *secondary* (dose-response) screen needed? No.

Asked and settled 2026-08-30. **Repurposing Public 24Q2 has no secondary screen.**
It contains two single-dose screens, REP1M and REP300, both at 2.5 uM with a
5-day treatment (its release notes are at `docs/README_prism.txt`). The
dose-response screen with AUC and IC50 exists only in the **older** PRISM
Repurposing releases, 19Q4 (`9393293`) and 20Q2 (`20564034`), as
`secondary-screen-dose-response-curve-parameters.csv`.

Taking it would be a bad trade, and not on general grounds - on a specific one.
Checked against the 20Q2 treatment table
(`prism-repurposing-20q2-primary-screen-replicate-collapsed-treatment-info.csv`),
of the seven compounds plan section 10 names, **that release contains only
navitoclax and venetoclax**. None of the three MCL1 inhibitors is in it -
S63845, AMG-176 and AZD5991 are all absent, and the secondary screen is a subset
of the primary, so they cannot be in it either.

So the secondary screen would cost **the entire MCL1 arm**, which is the arm the
plan's prediction actually rests on ("selectively dependent on MCL1 and/or
BCL2L1"), in exchange for dose-response on two compounds that 24Q2 already
covers. It would also bring a different, older cell-line panel and 290 MB.

**Limitation to state, not to fix:** the 24Q2 readout is single-dose LFC at
2.5 uM, which is coarser than an AUC - a null could in principle be a dose
artefact rather than an absence of effect. There is no release that gives both
dose-response and MCL1 coverage, so this is a limitation of the available data.

Also **not needed**: `Repurposing_Public_24Q2_Cell_Line_Meta_Data.csv`. The
Extended Primary Data Matrix already has `depmap_id` (`ACH-...`) as its column
names, verified directly, so no cell-line join is required.

## The release gap in item 3

`Model.csv`, the expression matrix and `CRISPRGeneEffect.csv` are **26Q1**.
**PRISM is 24Q2**, and that gap is now about twenty months. It is not an
oversight: Repurposing has its own release cadence, `Repurposing Public 24Q2`
(May 2024) is still the latest, and **26Q1 contains no drug-sensitivity data at
all** - verified against its release notes. PRISM joins on `ModelID`, which is
stable across releases, so this is a separate dataset rather than a stale copy of
this one. The consequence is that models added since 24Q2 have no PRISM row:
**item 3 loses n, it does not gain bias.** Script 14 records this under
`spec$prism`. Do not "fix" it by downgrading the rest to 24Q2.

## A superseded file, deliberately renamed

`data/raw/depmap/CRISPRGeneEffect_24Q4_figshare_SUPERSEDED.csv` is the 24Q4
CRISPR matrix fetched from figshare on 2026-08-29, before the current portal
release had been checked. It is kept, renamed, rather than deleted - but it must
not be used, because it is a different release from everything else here.
Nothing inside the file records its release, so script 14 guards on the only
signal that distinguishes it: its exact byte length, 428,678,699. If a file of
exactly that size appears as `CRISPRGeneEffect.csv`, the script stops.

## Checksums

```
shasum -a 256 data/raw/depmap/*.csv
```

| File | Release | Bytes | SHA-256 |
|---|---|---|---|
| `Model.csv` | 26Q1 | 697,455 | `ea4e0b2a3bc806f81df62689a5ae75f1a100135727a3d7b8a4c7ccc8815183f8` |
| `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | 26Q1 | 305,007,605 | `0377be80c525fde98cbd2c6e8b06bdf2a4014a9683eb70182c1f8649d711021a` |
| `CRISPRGeneEffect.csv` | 26Q1 | 440,646,050 | `e610a4cefb13a82b5b256b47eb08b63ff14843f8dbd0fb164bc0a32688e5b89e` |
| `Repurposing_...Extended_Primary_Data_Matrix.csv` | 24Q2 | 72,456,953 | `3b6554cfc6c765af53088a676edc7bce00ee7d84fe808b93bbfa892de607bc3d` |
| `Repurposing_...Extended_Primary_Compound_List.csv` | 24Q2 | 719,567 | `7e78f5901c1a97d2baab0789ab89832e716388da4eacaa9f094e7d2f2f5a3463` |
| *superseded* `CRISPRGeneEffect_24Q4_figshare_SUPERSEDED.csv` | 24Q4 | 428,678,699 | `3d8f3ec6dbf2db7ff834b79b508622ec0b226f3518003fe96ecf5a4fcf167e3b` |

PRISM downloaded and byte-length verified 2026-08-29.

## Read-time findings, verified against the actual files

Checked on 2026-08-29 rather than assumed, because two of them would have failed
silently:

- **`CRISPRGeneEffect.csv` is ModelIDs x genes**, first column unnamed
  (`ACH-000001`), gene columns named `SYMBOL (ENTREZ)`.
- **The PRISM data matrix is the TRANSPOSE of that: compounds x cell lines**,
  6,790 x 919. Rows are `BRD:BRD-...` ids, columns are ModelIDs. Read without
  transposing, the line intersection is empty and Block G item 3 reports
  "0 lines with PRISM" instead of failing. Script 14 transposes and asserts the
  orientation.
- **`Repurposing_..._Compound_List.csv` keys on `IDs`** (`BRD:BRD-...`, matching
  the matrix rownames exactly) with the human name in `Drug.Name`. `Drug.Name`
  is UPPERCASE for some compounds (`VENETOCLAX`, `NAVITOCLAX`) and mixed-case
  for others (`AZD5991`, `S63845`), so match case-insensitively.

### The BCL-XL gap in item 3

**5 of the 7 compounds plan section 10 names are in this release. The two that
are missing are `A-1331852` and `A-1155463` - both SELECTIVE BCL-XL
inhibitors.**

What remains on the BCL-XL side is **navitoclax**, which inhibits BCL2, BCL-XL
and BCL-W and therefore cannot attribute a result to BCL-XL on its own. It is
partly recoverable: **venetoclax is present and is BCL2-selective**, which is
exactly the specificity-control role the plan already assigns it, so navitoclax
and venetoclax read as a pair bound the BCL-XL/BCL-W component. Read them
together; do not read a navitoclax result as a BCL-XL result.

MCL1 is fully covered: `S63845`, `AMG-176` and `AZD5991` are all present, all at
2.5 uM. Note they sit in different screens (`REP.1M` vs `REP.PRIMARY`), so each
drug is modelled separately and no cross-drug arithmetic is done.

## Read-time notes

- **Column names are `SYMBOL (ENTREZ)`** in both the expression and CRISPR
  matrices. Script 14 strips the Entrez suffix and keeps the first occurrence of
  any duplicated symbol.
- **First column is the unnamed `ModelID`.** `data.table::fread` names it `V1`.
- **Chronos sign:** 0 = no effect, -1 = median common essential. **More negative
  = more essential.** So the plan's "selectively dependent" prediction is a
  *negative* `MYC:OX` interaction.
- **PRISM log2 fold change:** more negative = more sensitive. Same direction.
- **`OmicsExpressionTPMLogp1HumanProteinCodingGenes` is log2(TPM+1).** GSVA takes it
  as supplied; mitoPPS takes `2^x - 1`. Script 14 keeps these in two objects
  that never meet, and asserts the scale on load.

## Two things script 14 will not do quietly

- **There is no expression-matched null in CCLE.** Script 07's nulls are
  TCGA-specific and are not transferable. The 18-arm panel in script 14 is a
  **rank ordering, not a calibrated p-value**. A positive OXPHOS result needs a
  CCLE matched null built before it is reportable.
- **mitoPPS here runs on linear TPM, not linear DESeq2-normalised counts.** A
  declared deviation, recorded in the script header and in the saved object.
  mitoPPS is a composition measure and is robust to total content, which is why
  it is defensible - but CCLE mitoPPS values are **never** numerically
  comparable to TCGA mitoPPS values. Only the pattern transfers.
