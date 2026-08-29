# DepMap / CCLE - Block G inputs

Provenance for the DepMap files consumed by `scripts/14_depmap_dependency.R`
(Block G, plan section 10). The files themselves live in
`data/raw/depmap/`, which is **gitignored** - like every other large download in
this repo, they are re-fetchable from the URLs below and are not on `origin`.

**Release: DepMap Public 24Q4.** Everything must come from ONE release. Model IDs
are stable across releases but the line set and the Chronos scaling are not, and
mixing releases produces a join that looks complete and is not. If you move to a
newer release, change `DEPMAP_RELEASE` in script 14, re-fetch **all** files, and
record the new checksums here.

## Files

| File | Required | Used for |
|---|---|---|
| `Model.csv` | yes | lineage assignment (`OncotreeLineage == "Breast"`) |
| `OmicsExpressionProteinCodingGenesTPMLogp1.csv` | yes | MYC / OXPHOS scoring. **log2(TPM+1)** |
| `CRISPRGeneEffect.csv` | yes | Chronos gene effect, Block G item 2 |
| `Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv` | optional | Block G item 3, drug sensitivity |
| `Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv` | optional | compound name -> column id |

Without the two PRISM files, script 14 **skips** item 3 and says so. It does not
substitute anything for it.

## Acquisition - three files are scriptable, two are not

Established by enumerating the figshare API on 2026-08-29, not assumed:

- **figshare carries only the Achilles/CRISPR half of every DepMap release.**
  Checked across 24Q4 (`27993248`), 24Q2 (`25880521`) and earlier: each article
  holds `AchillesCommonEssentialControls`, `AchillesHighVarianceGeneControls`,
  `AchillesNonessentialControls`, the two QC reports, `AvanaGuideMap`,
  `AvanaLogfoldChange`, `AvanaRawReadcounts`, `CRISPRGeneDependency` and
  `CRISPRGeneEffect` - and **nothing else**. `Model.csv` and the omics matrices
  are portal-only.
- **The DepMap portal is behind a Cloudflare challenge.** A `curl` of any portal
  download URL returns a ~5 KB HTML verification page, not data. Saved under a
  `.csv` name that fails much later and very confusingly, so script 14
  size-checks every required input on load. `storage.googleapis.com` paths for
  the same bucket return 403. There is no scriptable route; use a browser, which
  is the intended one.
- **PRISM Repurposing is complete on figshare** (`25917643`), including the two
  files script 14 wants.

### Scriptable - 3 files

```
cd data/raw/depmap
curl -L -o CRISPRGeneEffect.csv \
  https://ndownloader.figshare.com/files/51064667                        # DepMap 24Q4
curl -L -o Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv \
  https://ndownloader.figshare.com/files/46630984                        # PRISM 24Q2
curl -L -o Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv \
  https://ndownloader.figshare.com/files/46630981                        # PRISM 24Q2
```

The numeric file ids are pinned deliberately: resolving them through
`/v2/articles/<id>/files` works too, but a pinned id is a pinned file, and a
figshare version bump would otherwise change the data under a stable command.

### Browser - 2 files

https://depmap.org/portal/data_page/?tab=allData , release **24Q4**:

- `Model.csv`
- `OmicsExpressionProteinCodingGenesTPMLogp1.csv`

Move both into `data/raw/depmap/`.

## The one deliberate release mismatch

`Model.csv`, the expression matrix and `CRISPRGeneEffect.csv` are all **24Q4**.
**PRISM is 24Q2**, because Repurposing has its own release cadence and there is
no 24Q4 Repurposing - `Repurposing Public 24Q2` (2024-05-31) is the current one.
It joins on `ModelID`, which is stable across releases, so this is a separate
dataset rather than a stale copy of this one. Script 14 records it in the saved
object under `spec$prism`. Do not "fix" it by downgrading the rest to 24Q2.

## Checksums

```
shasum -a 256 data/raw/depmap/*.csv
```

| File | Release | Bytes | SHA-256 |
|---|---|---|---|
| `Model.csv` | 24Q4 | | *fill in after browser download* |
| `OmicsExpressionProteinCodingGenesTPMLogp1.csv` | 24Q4 | | *fill in after browser download* |
| `CRISPRGeneEffect.csv` | 24Q4 | 428,678,699 | `3d8f3ec6dbf2db7ff834b79b508622ec0b226f3518003fe96ecf5a4fcf167e3b` |
| `Repurposing_...Extended_Primary_Data_Matrix.csv` | 24Q2 | 72,456,953 | `3b6554cfc6c765af53088a676edc7bce00ee7d84fe808b93bbfa892de607bc3d` |
| `Repurposing_...Extended_Primary_Compound_List.csv` | 24Q2 | 719,567 | `7e78f5901c1a97d2baab0789ab89832e716388da4eacaa9f094e7d2f2f5a3463` |

Downloaded and byte-length verified 2026-08-29.

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
- **`OmicsExpressionProteinCodingGenesTPMLogp1` is log2(TPM+1).** GSVA takes it
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
