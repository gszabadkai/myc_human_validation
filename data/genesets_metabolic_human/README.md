# Curated human metabolic gene list (`GS_metabolic`)

Snapshotted 2026-08-28 for the script 07 specificity work. Small enough to
commit, so this one is on `origin`.

| | |
|---|---|
| File | `gs_metabolic_human.csv` |
| Rows | 2,347 genes, one per row |
| Columns | `gene_name`, `gene_symbol`, `entrez_id`, `classification` |
| Classifications | 74 |
| SHA-256 | `891427debea8403a24e15a5d43eb309c4b32481ef18309bd21bd5d6606c50bb8` |

## Provenance

Extracted from the **`human` sheet** of
`data/raw/user_curated/GS_metabolic_genes_list.xlsx` in the geneset library:

```
repo   /Users/gs/G/data/MK_myc_2022/mammary_geneset_library
tag    v1.0
path   data/raw/user_curated/GS_metabolic_genes_list.xlsx
blob   21f79d9a04ae88dbeff8794fd5a4f34a51809047
```

Retrieved read-only with `git -C <repo> show v1.0:<path>`, per CLAUDE.md. Rows
with a missing symbol or classification were dropped (2,348 -> 2,347);
whitespace trimmed. Nothing else was changed.

## This is human-native, and that is the point

The xlsx has two sheets, `human` and `mouse orthologs`. The library reads the
**human** sheet and maps it *to* mouse (`R/04_load_user_curated_sets.R`, which
requires columns `Gene Symbol` and `Classification` from the `human` sheet).

So the human version is the **original**, not a round trip. This is not the
mouse-derived `outputs/gmt/human/` tree that CLAUDE.md rejects - that tree is
mouse-native sets pushed through `mouse_to_human()`, gitignored and unpinned by
the tag. This is a tracked raw input carrying its own human data. Different
objects; see the amended wording in CLAUDE.md.

## What it is for, and what it is NOT for

**NOT the specificity battery.** That is MitoCarta, matching the mouse arms
named in `myc_mouse/scripts/43_substrate_specificity_and_tradeoff.R`. See
`docs/2026-08-28_specificity_panel_proposal.md`.

**It is here for complex-level resolution**, which MitoCarta does not provide:

```
Complex I 49 | Complex II 4 | Complex III 11 | Complex IV 27
ATPase 5     | Proton Transport 49 | Ubiquinone 13
```

MitoCarta's `OXPHOS subunits` is one flat list of 102 genes. The mouse claim is
that OXPHOS subunit LFC and mitoPPS drop **across all complexes**
(`fig2_wt_mito_contraction.R`), and testing that in human needs per-complex
sets. This file has them.

It also carries the wider metabolic classifications (`Fatty Acid` 139,
`Krebs` 26, `Redox` 68, `Glutathione` 34, `Folate` 12, `Glycine` 7 ...). Those
are **not** used for the primary battery - MitoCarta is - but they are available
as an independent second opinion if a specificity result needs one, since they
are a different curation of the same biology.

## Read-time traps

- `classification` is a flat character column, not a hierarchy. `Complex I` and
  `Proton Transport` are siblings, and a gene appears **once**: the 74 classes
  partition the 2,347 genes rather than overlapping. Do not assume set overlap
  semantics from MSigDB or MitoCarta.
- Symbols are human HGNC as curated, and have not been re-harmonised against the
  MitoCarta Sheet-3 vocabulary the way G1 harmonised the Felsher and CollecTRI
  sets. **Run them through the same alias resolution before intersecting with
  anything**, or a deprecated symbol will silently drop (`KARS` -> `KARS1` cost
  G1 a gene until it was caught).
- Every gene appears exactly once, so class sizes sum to 2,347. A gene present
  in the expression matrix but absent here is not "not metabolic" - the list is a
  curation, not an exhaustive annotation.
