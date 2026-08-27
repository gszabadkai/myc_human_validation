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

Setup and gates. Nothing downstream of the gates should be built yet.

Run first, in order (see plan sections 6 and 12):
- **G1** — MYC signature overlap audit (script 04). Decides whether the Felsher
  signature survives as a primary MYC estimator.
- **G2** — MYC/MCL1/BCL2L1 CNV co-occurrence (script 05). Pure GISTIC, no expression.
  Decides whether H1 or H2 becomes the Panel b figure.
- **G3** — forkscale availability from `github.com/gszabadkai/Menegollo_Bentham`.

Do not build script 09 before G1 and G2 return. Do not build script 13 before the F3
incremental-value check returns.

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

- Human sets live in `data/genesets_from_library_human/`, a snapshot of
  `mammary_geneset_library` at tag `v1.0`, dereferenced commit
  `cbd8f16d2b0f95c5d4e86bed6aa112e42538a34b`.
- Consume as-is. Do not rebuild them here, do not edit them in place. If the library
  changes, re-snapshot from a new tag rather than patching files.
- The tag is the point of truth. See `data/genesets_from_library_human/README.md`.
- mtDNA-encoded protein-coding genes (13, `MT-` prefix) sit in a separate synthetic
  "mtDNA-encoded OXPHOS subunits" pathway and are never pooled with nuclear-encoded
  OXPHOS subunits — expression-scale skew.

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
  genesets_from_library_human/  library v1.0 snapshot + provenance README
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
