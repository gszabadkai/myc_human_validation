# CLAUDE.md — Human validation of the MYC / OXPHOS / apoptotic-priming axis

Rules and context for Claude Code sessions in this repo. Read this first, then the
plan: `docs/2026-08-27_human_validation_plan.md`.

## What this repo is

The human arm of the MMTV-Myc mouse study. Bulk human breast cancer transcriptomes
(TCGA-BRCA primary; METABRIC, SCAN-B, neoadjuvant GEO cohorts, DepMap/CCLE), used to
test whether the mouse MYC x OXPHOS coupling to PUMA/BCL-XL priming holds in human
tumours, how surviving tumours escaped it, and what that escape costs therapeutically.

Feeds the manuscript "Mitochondria integrate oncogenic and metabolic transcriptional
programs to shape breast cancer progression" (Nature Metabolism Letter), as 2-3 main
panels plus Extended Data. Companion paper: Menegollo, Bentham et al., Cancer Res 2024
(CAN-23-3172). Mouse repo: `myc_mouse` (separate, see "Cross-repo rules").

## THIS IS A HUMAN REPO

The single most important rule. `myc_mouse` and this repo were deliberately split to
prevent wrong-species errors.

- Human MitoCarta 3.0, human GMTs, human gene symbols only.
- Never load a mouse GMT, a mouse MitoCarta sheet, or the mouse ortholog table here.
- Never write into `myc_mouse` from this repo.
- If a file in `data/` has `mouse` in its name and is not under `data/from_myc_mouse/`,
  stop and ask. Everything under `data/from_myc_mouse/` is a frozen mouse *result*
  table for cross-species comparison, never a gene set and never an analysis input.

## Current phase

**All three gates are discharged. Scripts 00–11, 14 and 15 are written and run.
`STATE` is frozen.** Updated 2026-08-30.

| Gate | Outcome | Note |
|---|---|---|
| **G1** MYC signature overlap | passed, with D2 | `docs/2026-08-28_G1_result_and_decisions.md` |
| **G2** MYC/MCL1/BCL2L1 CNV co-occurrence | **passed**; BBC3 fails | `docs/2026-08-28_G2_result.md` |
| **G3** forkscale availability | **passed**, wider than assumed — TCGA too, not METABRIC-only | `docs/2026-08-29_G3_result_forkscale_availability.md` |
| **F3-pre** MB1 forkscale redundancy | **INTERMEDIATE** — rho 0.529 / 0.418; the stop gate does not close | `docs/2026-08-30_F3pre_and_block_D_result.md` |

Where the science stands, so that nothing is re-derived by accident:

- **H1 is falsified on both clauses** (Block C and Block B). **H2 and H3 are not
  supported.** The reverse "escape" reading is not supported either. **Three of the four
  falsification criteria are met; only H4 is untested.** Criterion 1 is now discharged on
  its **complete** stratum list — Block D added the forkscale stratum, which Block C could
  not fit because G3 was not yet discharged.
- **Block G (DepMap) failed at adequate power.** G-a (MCL1) and G-b (BCL2L1) both fail at
  n = 1,130 with CIs of about ± 0.04. This was the orthogonal functional test, not more
  observational TCGA correlation.
- **Block D is null** on both instruments, gate closed, matched null correctly not run.
  Underpowered by construction (n = 844, continuous modifier); a null there is
  uninformative and no direction was pre-specified for it.
- **One finding survives both instruments:** `MYC x OXPHOS` associates with **lower
  BCL-XL and higher BIM**, PUMA unmoved, with no compensating amplification. BIM is
  pre-specified in the mouse arm (`myc_mouse/scripts/44`), not a fishing result. In
  DepMap it is **UNRESOLVED** — read the section 5.1 correction in the Block G note
  before touching it. **No claim of breast-specificity is licensed.**
- Two results were reported and **deliberately not pursued** (`Glycine metabolism`,
  `Folate and 1-C`). Chasing either is the fifth post-hoc hypothesis section 2 forbids.
  A third joins them: `D MB2 DESCRIPTIVE log2_BCL2L1` (p 0.056, one instrument).
- **One live commitment:** the BIM replication declared in the Block C note section 9.
  Honour it or delete it; do not amend it. Block D's descriptive `log2_BCL2L11` rows do
  **not** discharge it.

Next: script 12 — fetch GSE194040, GSE164458, GSE25066, where H4 and the BIM
replication share their downloads. **Nothing in the H4 declaration may be revised once
any of those outcome columns has been read.**

**`STATE` is frozen** (script 11, `results/state_definition.rds`) with the definition
unamended, as a portable constructor script 13 must **call** rather than re-implement.
Two findings from the freeze that constrain H4 and belong in the text:

- **H4's BUFFER is an expression construct, not the copy-number BUFFER that G2 passed
  on and Blocks B and G failed to support** — the tertile fallback agrees with the
  GISTIC rule at only kappa 0.221. Kept, not tuned. A null H4 does not re-test Block B
  or G, and a positive H4 would not rescue them.
- **The level-3-vs-4 contrast is subtype-confounded against its own prediction.** PAM50
  adjustment is mandatory, and a positive in the predicted direction is conservative.

The continuous `BUFFER_c` that D5's primary H4 test requires is declared pre-data in
`docs/2026-08-30_STATE_frozen_and_H4_buffer_declaration.md` section 6, with
`MYC:OXPHOS:BUFFER_c` predicted **negative**.

**A limitation registered before the H4 data arrives:** forkscale does not exist for the
neoadjuvant cohorts and cannot be constructed there. If H4 returns a positive, the
fork-redundancy question must be argued in text from the TCGA rho of 0.53 — it cannot be
tested where the claim is made. See the F3-pre note section 7.

Do not build script 13 before F3 returns. **F3 is a stop gate, not a robustness check.**
F3-pre (exposure side) is discharged; **F3 proper is a survival test, belongs in
METABRIC, and is still blocked on the METABRIC sample-identifier file.**

## Workflow — "Option A" (do not deviate)

- Claude Code **writes and edits** the numbered pipeline scripts. It does **not run
  them.** The author sources them in Positron interactively.
- Infrastructure tasks (git, snapshotting files, provenance READMEs, editing this file,
  planning docs) Claude Code may execute directly. The line: numbered analysis scripts
  are hand-run by the author; plumbing is not.
- Every numbered script ends with an `if (FALSE) { ... }` sandbox block — skipped by
  `source()`/`Rscript`, run line-by-line in Positron for inspection.
- Plan first: in a planning session produce the build spec and wait for approval before
  writing scripts. Commit per verified phase — git is the safety net.
- When in doubt, ask. The cost of clarifying is small; the cost of destroying ambiguous
  state is large.

## Pre-registration discipline

This arm was specified before any model was fitted, because the obvious analyses here
can only confirm. Treat the plan's section 2 as binding.

- The hypotheses (H1–H4), the primary endpoint `PRIME = log2(BBC3) - log2(BCL2L1)`, the
  specificity battery, and the falsification criteria are **fixed**. Do not add a fifth
  post-hoc hypothesis if the first four fail.
- The `STATE` variable (plan section 7.5) is frozen in script 11 and **must not be
  revised after outcome data has been seen**. A three-way conjunction on median splits
  gives eight cells; picking the worst one post hoc is p-hacking.
- Every positive OXPHOS result requires its accompanying negatives (pathway negatives
  and endpoint negatives) reported alongside. A positive without the negatives is not
  reportable.
- The F3 incremental-value check against MB1 forkscale is a **stop gate**, not a
  robustness check. If `STATE` adds nothing over `MB1_forkscale`, there is no Panel c.

## Scale discipline — the most likely silent error

- **GSVA / ssGSEA** want log-scale input: VST, `kcdf = "Gaussian"`.
- **mitoPPS** wants linear DESeq2-normalised counts.

These are opposite requirements. They must not share an input object. State the scale in
a comment at the top of every scoring block.

Two further traps:
- **GSVA is cohort-relative.** Score all samples of a cohort in one run. Scores from
  separately-scored cohorts are not comparable and must never be pooled. In Block F, five
  cohorts are scored independently — meta-analyse the effect estimates, not the scores.
- **mitoPPS baseline is composition-dependent.** It reports the *shape* of the
  mitochondrial program, not its level, and is deliberately robust to total content.
  Never compare mitoPPS values numerically across cohorts or species; only the pattern
  transfers. Report an OXPHOS *level* metric separately.

## Gene sets — consume the snapshot, do not rebuild

Three provenanced inputs, each with its own README. Consume as-is: do not rebuild them
here, do not edit them in place. If a source changes, re-snapshot rather than patch.

| Input | Location | Source |
|---|---|---|
| Human MitoCarta 3.0 | `data/mitocarta_human/` | Broad, `Human.MitoCarta3.0.xls`, 1,136 genes / 149 pathways |
| CollecTRI regulons | `data/collectri_human/` | OmniPath web service, dated snapshot, 1,201 TFs |
| Felsher MYC signature | `data/genesets_from_library_human/` | `mammary_geneset_library` tag `v1.0` (`cbd8f16d…`), 67 genes |
| Menegollo bicluster forkscale | `data/menegollo_biclusters/` | `gszabadkai/Menegollo_Bentham` @ `8fdbb343…`, TCGA n=1,037 + METABRIC PC1 |
| TCGA Clinical Data Resource | `data/tcga_cdr/` | GDC PanCanAtlas, Liu et al. 2018 Cell. **F4 and covariates only — not TCGA survival** |

### The library's human GMT tree is mouse-derived — do not use it

`mammary_geneset_library` ships `outputs/gmt/human/`. **It must not be loaded here.** The
library is mouse-native; that tree is built by running `mouse_to_human()` over the mouse
sets (`R/14_export_gmts.R`). The files carry human symbols and mouse provenance, which is
exactly the failure mode this repo split exists to prevent — and it looks correct.

Two concrete instances, both already avoided:

- `01_mitocarta_human.gmt` is Mouse MitoCarta ortholog-projected: 1,083 genes against the
  real 1,136. Use `data/mitocarta_human/` instead.
- its Felsher set is a human→mouse→human round trip: 62 genes against the native 67,
  dropping NPM1 and RPLP0 and adding the paralog artefact EIF5AL1.

`outputs/` is also gitignored in the library, so nothing there is reachable by
`git show` and nothing there is pinned by the tag. Tracked library files are, and that is
the only supported route:
`git -C /Users/gs/G/data/MK_myc_2022/mammary_geneset_library show v1.0:<path>`

**The rejection is of `outputs/gmt/human/`, not of the library.** Amended
2026-08-28. That tree is mouse-native sets pushed through `mouse_to_human()`,
gitignored, and unpinned by the tag - all three reasons apply to it and to
nothing else. Tracked *raw inputs* carrying their own native human data are a
different object:
`data/raw/user_curated/GS_metabolic_genes_list.xlsx` has a `human` sheet (2,347
genes, 74 classifications) which the library reads and maps *to* mouse. Taking
the human sheet is not a round trip. **Check the sheet, not the repository** -
and snapshot what you take, with the blob SHA. See
`data/genesets_metabolic_human/`.

**RESOLVED 2026-08-28.** Script 07's specificity panel comes from **Human
MitoCarta 3.0**, using the same MitoPathway names as the mouse arms in
`myc_mouse/scripts/43_substrate_specificity_and_tradeoff.R`. The cross-species
link is the pathway *name*; each species uses its own native MitoCarta, so no
projection is involved. See `docs/2026-08-28_specificity_panel_proposal.md`.

### Standing convention

- mtDNA-encoded protein-coding genes (13, `MT-` prefix) sit in a separate synthetic
  "mtDNA-encoded OXPHOS subunits" pathway and are never pooled with nuclear-encoded
  OXPHOS subunits — expression-scale skew.
- MitoCarta's `OXPHOS` umbrella includes assembly factors. The plan's primary OXPHOS
  measure is `OXPHOS subunits`. These are different sets; pick deliberately. The
  **assembly factors are the primary specificity control**, not a leftover: in the
  mouse they sit at percentile 50.2 of expression-matched sets while the subunits
  sit at 0.0, and being the same complexes they exclude what a distant pathway
  cannot (D-2026-08-28, specificity panel proposal section 2).
- **OXPHOS is measured on two instruments, both primary.** GSVA level (portable
  across cohorts) and mitoPPS (composition; the instrument the mouse interaction
  was fitted on). They answer different questions and the same genes can move in
  opposite directions on them. Report both; claim only what both support.

## R coding rules

Copied from `docs/R_CODING_INSTRUCTIONS.md` (snapshotted from `myc_mouse`). The full
file governs; these three cause the most damage.

1. **Never `print(n = X)` after `head()`.** `head()` may coerce a tibble to a
   data.frame, causing `n` to be read as `na.print`. Use `head(X) %>% print()`, or
   `tibble %>% print(n = X)` as a separate call.
2. **Always `dplyr::count()`**, never bare `count()` — namespace conflicts.
3. **ASCII-only strings in scripts.** Handle latin1/cp1252 at read time with
   `fileEncoding = "latin1"` and `iconv()`.

No `renv`; packages are installed system-wide.

## Project structure

```
scripts/    numbered R pipeline (00-18), see plan section 11
docs/       the plan, R_CODING_INSTRUCTIONS.md, dated decision notes
data/
  raw/                          large downloads (gitignored, NOT on origin)
  from_myc_mouse/               frozen mouse result tables + provenance README
  mitocarta_human/              Human MitoCarta 3.0 xls + provenance README
  collectri_human/              CollecTRI snapshot (tsv.gz) + provenance README
  genesets_from_library_human/  library v1.0 human assets + provenance README
  menegollo_biclusters/         companion-paper MCbiclust pc1/index + provenance README
  tcga_cdr/                     TCGA Clinical Data Resource + provenance README
functions/  shared utilities
external/   vendored reference code (mitotyping / Monzel et al.) — reference, not to edit
results/    intermediate .rds (gitignored, generated at runtime)
outputs/    figures and tables (gitignored, generated at runtime)
```

`results/` and `outputs/` are regenerable by re-running scripts. **`data/raw/` is not**
— it holds slow TCGA/GEO downloads, is gitignored, and needs its own backup outside git.

## Git discipline

- `main` is the trunk. Feature branches off `main` as needed.
- Do **not** reproduce the `myc_mouse` branch structure (`new-analysis` trunk, `main`
  historical). That is an artifact of a two-pipeline consolidation in that repo, not a
  model.
- No worktrees until there is a reason for one.
- Read-only git ops (`git show`, `git diff`, `git log`, `git ls-tree`) are always fine.
  Stop-and-check before any destructive action; never force-push a shared branch.
- Backup is `origin`, not Google Drive. This repo lives outside the Drive sync tree
  deliberately — see `docs/2026-08-27_human_validation_plan.md` section 5. Do not move
  it into `/Users/gs/G/`.

## Cross-repo rules

- The mouse repo is at `/Users/gs/G/data/MK_myc_2022/myc_mouse`. It is **not** attached
  to normal sessions and must not be added to `additionalDirectories`.
- `--add-dir` grants read **and write** access; there is no read-only mode. It was used
  once, for the initial snapshot, and then dropped.
- To read a mouse file later without granting directory access:
  `git -C /Users/gs/G/data/MK_myc_2022/myc_mouse show <ref>:<path>`
- Anything that becomes an *input* to this repo gets snapshotted with a provenance
  README recording the source commit SHA. Nothing is sourced across repos at runtime.

## Reference

- Plan: `docs/2026-08-27_human_validation_plan.md` — hypotheses, gates, models, panel
  budget, falsification criteria.
- Companion paper: Menegollo, Bentham et al., Cancer Res 2024 (CAN-23-3172).
- Format target: Nature Metabolism Letter — 200-word introductory paragraph, 2,500-word
  main text, 2-4 display items total, up to 10 Extended Data figures.
