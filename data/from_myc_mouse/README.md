# Frozen mouse result tables

Result tables carried across from the `myc_mouse` repo so the human arm can be compared
against the mouse findings it is validating.

**These are frozen mouse RESULT tables for cross-species comparison only. They are never
analysis inputs and they are never gene sets. Do not edit them here. If the mouse
analysis changes, re-snapshot from a new commit rather than patching files in place.**

Nothing in this directory may be loaded as a scoring input, a gene set, an annotation
table, or an ortholog map. They exist to be read alongside a human result, not to enter
a human model. See the "THIS IS A HUMAN REPO" section of `CLAUDE.md`.

## Provenance

- Source repo: `myc_mouse`, at `/Users/gs/G/data/MK_myc_2022/myc_mouse`
- Branch: `paper-final`
- Commit (dereferenced): `6a9c7dd513800a2a433934314a87d161ce98caa2`
- Commit date: 2026-08-18
- Snapshot date: 2026-08-27
- Method: one-off `--add-dir` session, treated as read-only. Nothing in `myc_mouse` was
  created, modified, or deleted. Files were copied out; nothing is sourced across repos
  at runtime.

**The source tree was dirty at snapshot time.** One unrelated tracked file
(`figures/panels/rebuild_panels.R`, 3 insertions / 4 deletions) plus session logs and
untracked docs were uncommitted. None of them are snapshot paths; the copied files were
verified unmodified against the commit above.

**The commit SHA does not pin the CSV tables.** `outputs/` and `results/` are untracked
in `myc_mouse` - these are generated artefacts. The SHA pins the *producing scripts*, not
the files. Source mtimes are recorded per file below as the closest available handle.

## Contents

```
data/from_myc_mouse/
  README.md                                  # this file
  interaction_sig_genes_p05.csv              # DESeq2 genotype x timepoint interaction genes
  interaction_three_sets_summary.csv         # set-level enrichment summary for the above
  mitopps_anova_sig.csv                      # mitoPPS pathway ordering (ANOVA)
  raw_pathway_anova_sig.csv                  # same ANOVA on raw, un-normalised pathway scores
  mouse_oxphos_puma_coupling_estimates.csv   # transcribed MYC x OXPHOS -> PUMA/Bcl-xL estimates
  mouse_oxphos_puma_coupling_estimates.md    # caveats that must travel with the above
```

| file | source path in `myc_mouse` | producing script | source mtime |
|---|---|---|---|
| `interaction_sig_genes_p05.csv` | `outputs/interaction_analysis/` | `11_interaction_gene_characterisation.R` | 2026-02-28 00:11 |
| `interaction_three_sets_summary.csv` | `outputs/interaction_analysis/` | `11_interaction_gene_characterisation.R` | 2026-02-28 00:11 |
| `mitopps_anova_sig.csv` | `outputs/mitopps/` | `08_mitoPPS_analysis.R:1246` | 2026-07-24 12:35 |
| `raw_pathway_anova_sig.csv` | `outputs/mitopps/` | `08_mitoPPS_analysis.R:1241` | 2026-07-24 12:35 |
| `mouse_oxphos_puma_coupling_estimates.csv` | (transcribed, see below) | `43_substrate_specificity_and_tradeoff.R`, `46_axis_ruler_test.R` | n/a |

### interaction_sig_genes_p05.csv

357 genes plus header. Columns: `ensembl_id, baseMean, log2FoldChange, lfcSE, stat,
pvalue, padj, weight, gene_symbol`.

**This is nominal p < 0.05, not FDR < 0.05.** The `p05` in the filename means raw
p-value. No gene in the table clears FDR 0.05 - the smallest `padj` is about 0.16. Do
not describe it as an FDR-significant gene set.

Mouse Ensembl IDs and mouse symbols. Any cross-species comparison needs an ortholog map,
which lives outside this repo and is not snapshotted here.

### interaction_three_sets_summary.csv

Three rows, set-level context for the table above: genes unique to Myc+ (1662 total, 357
at p<0.05, 4.30x expectation), unique to Myc- (906, 207, 4.60x), and overlapping (961,
21, 0.44x).

### mitopps_anova_sig.csv and raw_pathway_anova_sig.csv

125 and 122 rows. Columns: `pathway, effect, F_value, p_value, padj`.

Kept as a pair deliberately. `mitopps_anova_sig.csv` is the mitoPPS pathway ordering;
`raw_pathway_anova_sig.csv` is the same ANOVA before composition normalisation. The
difference between them is the composition adjustment doing its work.

**mitoPPS values are not comparable across cohorts or species.** mitoPPS reports the
*shape* of the mitochondrial program, not its level, and its baseline is
composition-dependent. Only the pathway *ordering* transfers - never the numbers. Report
an OXPHOS level metric separately. See the scale-discipline section of `CLAUDE.md`.

### mouse_oxphos_puma_coupling_estimates.csv (+ .md)

The MYC x OXPHOS coupling on the PUMA/Bcl-xL ratio - the p = 0.0052 interaction.

**Transcribed by hand, not copied.** The producing scripts write only `.rds` result
objects (`substrate_specificity_tradeoff.rds`, 5.3 KB; `axis_ruler_test.rds`), never a
CSV. Those objects were deliberately not copied: they are result objects rather than
summary tables. The values were transcribed from the two mouse-repo markdown tables that
report them, and every row carries its `source_doc` and `source_line`.

**Read the companion `.md` before quoting p = 0.0052.** The empirical permutation p is
0.083, the effect does not survive adding a timepoint term (p 0.088), and the
specificity is ruler-dependent. The mouse repo's own reading is that this is a lead, not
an established interaction.

## Not included, and why

Recorded so their absence is not mistaken for an error.

- **No `.rds` result objects.** Summary tables only. No raw count matrices, no DESeq2
  objects.
- **No gene sets of any kind.** No mouse GMTs, no mouse MitoCarta sheet, no ortholog
  table. Human sets come from `data/genesets_from_library_human/` and nowhere else.
- **No mouse figures.** The panel PDFs stay in `myc_mouse`.

## Reading a further mouse file later

Without granting directory access again:

    git -C /Users/gs/G/data/MK_myc_2022/myc_mouse show <ref>:<path>

Note this reaches tracked files only. The CSVs above are untracked in `myc_mouse` and
are not retrievable this way.
