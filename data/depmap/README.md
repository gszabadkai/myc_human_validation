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

## Acquisition

**The DepMap portal is behind a Cloudflare challenge.** A plain `curl` of a
portal download URL returns an HTML verification page, not data - and a 5 KB
"CSV" that is actually HTML fails much later and very confusingly. Script 14
guards against exactly this with a file-size check. Two routes:

### `CRISPRGeneEffect.csv` - figshare, scriptable

DepMap mirrors the Achilles/CRISPR half of each release on figshare. 24Q4 is
article `27993248`.

```
cd data/raw/depmap
curl -L -o CRISPRGeneEffect.csv \
  "$(curl -s https://api.figshare.com/v2/articles/27993248/files \
     | python3 -c 'import json,sys;print([f["download_url"] for f in json.load(sys.stdin) if f["name"]=="CRISPRGeneEffect.csv"][0])')"
```

That article carries only the CRISPR files - `AchillesCommonEssentialControls`,
`AvanaGuideMap`, `AvanaLogfoldChange`, `AvanaRawReadcounts`,
`CRISPRGeneDependency`, `CRISPRGeneEffect`, and the two QC reports. **It does not
carry `Model.csv` or the omics matrices.**

### `Model.csv`, the expression matrix, PRISM - browser

https://depmap.org/portal/data_page/?tab=allData , release 24Q4. Download in a
browser and move the files into `data/raw/depmap/`.

## Checksums - FILL IN AFTER DOWNLOAD

Record these before running script 14, so that a re-download or a release change
is detectable rather than silent.

```
shasum -a 256 data/raw/depmap/*.csv
```

| File | Bytes | SHA-256 |
|---|---|---|
| `Model.csv` | | |
| `OmicsExpressionProteinCodingGenesTPMLogp1.csv` | | |
| `CRISPRGeneEffect.csv` | ~428,700,000 | |
| `Repurposing_...Data_Matrix.csv` | | |
| `Repurposing_...Compound_List.csv` | | |

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
