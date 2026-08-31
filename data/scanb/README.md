# SCAN-B, for the BIM replication

Provenance for the cohort that honours the replication declared in
`docs/2026-08-29_block_c_result_H1_not_supported.md` section 9. Files live under
`data/raw/scanb/`, which is **gitignored and not on origin**. This README is the
record; re-download from the commands below.

Every parameter of the analysis is fixed, pre-data, in
`docs/2026-08-31_scanb_bim_replication_declaration.md`. **Read that before
touching `scripts/16_fetch_scanb.R` or `17`.** Consumed by script 16.

| Accession | n | Platform | Role |
|---|---|---|---|
| **GSE202203** | 3,207 | Illumina HiSeq 2000 (GPL11154) + NextSeq 500 (GPL18573) | the BIM replication, single cohort |

Sweden Cancerome Analysis Network - Breast (SCAN-B), ClinicalTrials.gov
NCT02306096. The deposit accompanies the ESR2 paper from the Staaf group; the
cohort is population-based primary breast cancer.

## THE ACCESSION IS GSE202203, NOT GSE96058

This is the single most important line in this file, because the plan named the
wrong one and it looks correct.

Checked 2026-08-31 against the GEO FTP listings and the series-matrix
`!Sample_data_processing` fields:

| Accession | n | Deposited quantity | Counts? |
|---|---|---|---|
| GSE96058 (Brueffer 2018) | 3,273 | cufflinks FPKM, deposited as `log2(FPKM + 0.1)` | **no** |
| GSE81538 (Brueffer 2016) | 405 | FPKM | **no** |
| **GSE202203** | **3,207** | **raw gene counts** and TPM | **yes** |

`mitoPPS` needs linear DESeq2-normalised counts (CLAUDE.md). A deposit that
arrives already logged cannot carry it - that is exactly what amendment A1
established for the three neoadjuvant cohorts, and it applies to GSE96058 for
the same reason. **GSE202203 is the only SCAN-B deposit on which the
declaration's co-primary clause can be honoured.**

GSE202203 and GSE96058 are overlapping draws from the same study. Nothing in
this repo reads GSE96058, so there is no double-counting to guard against.

## Download

Three files, about 205 MB. Sizes are the GEO directory listing's, read
2026-08-31.

```sh
mkdir -p data/raw/scanb && cd data/raw/scanb
B=https://ftp.ncbi.nlm.nih.gov/geo

# --- expression: RAW COUNTS. Not the TPM file. ---
curl -L -O $B/series/GSE202nnn/GSE202203/suppl/GSE202203_RawCounts_gene_3207.tsv.gz   # 204M

# --- phenotype: TWO series matrices, one per sequencer. Both are needed. ---
curl -L -O $B/series/GSE202nnn/GSE202203/matrix/GSE202203-GPL11154_series_matrix.txt.gz  # 264K, 2,913 samples
curl -L -O $B/series/GSE202nnn/GSE202203/matrix/GSE202203-GPL18573_series_matrix.txt.gz  #  46K,   294 samples
```

**Do not fetch `GSE202203_TPM_Raw_gene_3207.tsv.gz` (402 MB).** It is the same
data on the wrong scale, it is not read by any script here, and having it in the
directory is an invitation to point script 16 at it. Script 16's section 4.1
would catch that - a TPM matrix has a maximum far below the counts floor - but
the guard is the backstop, not the plan.

## Checksums

Downloaded and verified 2026-08-31. All three pass `gzip -t`. Script 16 section
1 recomputes `md5` and `sha256` on every run, so a silently changed
re-download is visible.

```
GSE202203_RawCounts_gene_3207.tsv.gz     214,227,990 bytes
  md5     c4620ee5c20d82ff697ed3beccbd78ec
  sha256  db38c17fad4b744a2a1ad93dbf957ea5bec3d645400fe1d3172bafe09f5c3709

GSE202203-GPL11154_series_matrix.txt.gz      270,141 bytes   2,913 samples
  md5     6758c933a90f1898352d357bcd21d87c
  sha256  59613566126997b912878e9da78ba96760b376d150823c0d9e52e115436f1350

GSE202203-GPL18573_series_matrix.txt.gz       47,209 bytes     294 samples
  md5     d359544a1c90f59b782e0ee9e647656d
  sha256  762b437c00935c91e3b1635a1acd07040c8a8d6d7f4d3ad6e56a568fd4442079
```

Shape verified on the files themselves, not assumed:

```
counts       3,208 header fields, 3,208 data fields (plain layout), 19,644 genes
symbols      19,644 distinct - no duplicate collapse needed
samples      2,913 + 294 = 3,207
endpoints    BCL2L11, BCL2L1, BBC3, MYC all present (also MCL1, BAX, BID, BAK1)
mtDNA        all 13 MT- protein-coding genes present
```

## Landmines

Five, four of them already met in the neoadjuvant cohorts. Script 16 carries a
guard for each; none of these is visible in the output if the guard is removed.

1. **The counts are not integers.** `A1BG 412.35820192`. They are StringTie `-e`
   estimated counts, the prepDE-style output, and DESeq2 requires integers. They
   are **rounded**, and script 16 reports what fraction of values that changed
   and the maximum shift. Gene-level counts are all that is deposited, so
   tximport with transcript lengths is not an option. Documented deviation from
   TCGA, where GDC counts are already integers.

2. **The samples are split across two platforms.** Both series matrices are
   needed; 2,913 + 294 = 3,207. Script 16 asserts the NextSeq count and the
   total, so a missing platform file cannot pass unnoticed.

3. **The join key.** Expression columns are `S000001`. The phenotype's
   `scanb external id` is `Q009012.C009079.S000008.l.r.m.c.lib.g.k2.a`.
   `!Sample_title` happens to be exactly `S000008`, so `title`/`raw` matches
   100%, but that is luck and not a contract. Script 16 declares four
   transformations - including one that extracts the `S`-token from the external
   id - scores every phenotype column under each, and stops below 90%. In
   BrighTNess the obvious key had **zero** raw overlap.

4. **The vroom connection buffer.** GEO series matrices put one field per sample
   on a single line, so metadata lines run to hundreds of kilobytes; GSE25066's
   was 185,443 bytes against vroom's 131,072 default, which stops readr with a
   message that reads like a corrupt download and is not one. Script 16 raises
   it to 8 MB per call and restores it afterwards.

5. **Ragged characteristic lines.** GEO does not align `key: value` fields
   across samples. In GSE25066 that mislabelled 198 of 508 endpoints. Script 16
   reads the key **per sample** and pivots on it. GSE202203's fields parse
   cleanly (verified 2026-08-31 on the GPL18573 matrix: 294 samples, 25
   columns, no raggedness), but the parser is the same one and stays.

   Note one cosmetic consequence: the ESR1/ESR2 lines carry a second colon
   inside the value (`esr1 expression: log2(tpm+0.1): 5.97`), so the parsed
   value keeps the `log2(tpm+0.1): ` prefix. Both columns are forbidden and
   dropped, so this is never cleaned up.

6. **The gene symbols are a 2014 vintage.** SCAN-B is annotated against UCSC
   knownGenes downloaded 22 September 2014; MitoCarta 3.0 uses current HGNC
   symbols. The ATP synthase subunits were renamed wholesale in 2018, so
   **nineteen of the eighty-nine genes in `OXPHOS subunits` - the exposure - do
   not match by name**, and coverage lands at 0.787 against a 0.80 floor.

   The genes are all there under their 2014 names. Script 16 section 7.3
   harmonises them through MitoCarta's own `Synonyms` column, forward direction
   only, refusing ambiguous and colliding candidates - script 07 section 2's
   map, applied to this matrix instead of TCGA's. Verified 2026-08-31: all 19
   resolve, none ambiguous, none colliding, coverage 0.787 -> 1.000.

   This is the landmine with the worst failure mode of the six. Accepting 70 of
   89 without noticing would have left a Complex V with no F1 head and no
   c-ring, still labelled "OXPHOS subunits", and nothing downstream would have
   looked wrong.

   The Felsher and Hallmark sets are not MitoCarta genes so their few renames
   (`H2AX`, `POLR1G`, `VARS1` and a handful more) are **not** rescued - they sit
   at 0.95-0.98, above the floor, and are reported by name. A general alias
   source was deliberately not introduced: choosing one after seeing which genes
   were missing is what the coverage floor exists to make unnecessary.

## The fence

`GSE202203` carries `overall survival days/years/event`, `relapse free interval
days/years/event`, `endocrine treated`, `chemo treated`, and expression-derived
`esr1`/`esr2` values.

**None of them reaches the saved object.** They are named in `SCANB_FORBIDDEN`
in script 16, dropped in section 3.1, and section 8 asserts they are absent
before saving. SCAN-B was fetched for the BIM replication and nothing else
(declaration section 10).

The ESR1/ESR2 columns are forbidden for a second, independent reason: they are
functions of the same matrix the exposure is built from, so using them as
covariates would be circular. Same rule as `HATZIS_FORBIDDEN` in script 12.

They are dropped rather than hidden behind a flag on purpose. A toggle gets
flipped; a named list has to be edited in a commit with a dated note beside it.
**If F2 (treatment-stratified survival) is ever decided, that decision amends
`SCANB_FORBIDDEN` and says so in its own note.**

## Source metadata, as deposited

```
!Series_title      Clinical associations of ESR2 (estrogen receptor beta) expression
                   across thousands of primary breast tumors [SCAN-B, n=3207]
Assembly           GRCh38 / hg38
Alignment          HISAT2 v2.1.0, UCSC knownGenes (22 Sep 2014), GENCODE 27 model
Quantification     StringTie v1.3.3b, -e --rf, protein-coding transcripts only,
                   novel transcripts discarded, collapsed on gene symbols
PAM50              Parker method, fixed reference cohort (depositors' own call)
```

Because the deposit is already restricted to protein-coding transcripts and
collapsed on symbols, there is no `gene_type` column to filter on and none is
needed. Script 16 checks instead that the 13 mtDNA-encoded protein-coding genes
survived the collapse, since MitoCarta's mtDNA pathway depends on them.
