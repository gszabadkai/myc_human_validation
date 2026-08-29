---
date: 2026-08-29
status: G3 DISCHARGED
relates-to:
  - 2026-08-27_human_validation_plan.md (sections 6 G3, 10 Block D and F3, 12, 13)
  - 2026-08-29_handoff.md (section 3, blocker a - now resolved)
  - data/menegollo_biclusters/README.md
  - data/tcga_cdr/README.md
decides:
  - G3 - PASSES. Stored forkscale exists and has been snapshotted.
  - It exists for TCGA, not METABRIC-only. The plan and Block D assumed the opposite.
  - No MCbiclust re-derivation is needed. The fidelity risk in plan section 13 is void.
  - F3 splits - an exposure-side redundancy pre-check in TCGA, then F3 proper in METABRIC.
  - One upstream file is still outstanding, for METABRIC only. See section 5.
next-action: script 15 (F3-pre + Block D), script 14 (DepMap) in parallel
---

# G3 result - forkscale availability

G3 was recorded as "passed" in a single line of the 2026-08-28 handoff with **no note
behind it and nothing snapshotted**. That was the state the 2026-08-29 handoff found and
flagged as blocker (a). This note discharges G3 properly: the check was redone, the data
located, verified and snapshotted, and the definition confirmed by recomputation.

## 1. What the gate asked

Plan section 6:

> Check `github.com/gszabadkai/Menegollo_Bentham` for stored TCGA bicluster assignments
> and forkscale values.
>
> **Decision:** if present, ED2 is cheap and high-fidelity. If absent, re-deriving via
> MCbiclust is a genuine fidelity risk (different run, seed, sample QC) and must be
> documented as a limitation.

## 2. Result - PASSES, and wider than the plan assumed

The repository is private; `gh` is authenticated as the repo owner, so it is readable
without cloning. Stored per-sample bicluster values are present for **both** cohorts,
and for TCGA they are complete.

**TCGA.** `TCGA_MB{1,2,3}_RNASeq_data_nonorm.RData` each hold a 1,037 x 6 frame with the
aliquot barcode, `pc1` and `index` for that bicluster. That is everything forkscale is
built from. A pre-assembled frame, `TCGA.all.biclusters.RNAseq.Rdata`, carries forkscale
already computed, but only for 849 samples (an upstream `inner_join` to two files that
are not in the repo). **The per-bicluster files are the primary source; the assembled
frame is a cross-check.**

**METABRIC.** The sort orders (`METABRIC_sort_data.RData`) and PC1 vectors
(`METABRIC_PC1_GSEA.RData`, n = 1,981) are both in the repo. The sample identifier
vector is not - see section 5.

All six files are snapshotted under `data/menegollo_biclusters/`, pinned to upstream
commit `8fdbb3437ae5537055d5d5429411bdb3b333c04a`, each verified byte-for-byte against
its upstream git blob SHA.

**This voids a plan section 13 landmine.** "Re-derived forkscale: if G3 fails and
MCbiclust is re-run, the result is a new result, not a replication." G3 did not fail.
Nothing is re-derived; the published values are used as published.

## 3. The definition, confirmed by recomputation

```
forkscale     = pc1 / index
forkscale.log = pc1 / log(index)     # the companion paper's preferred variant
```

Recomputing `MB1.pc1 / MB1.index` from the per-bicluster file reproduces the assembled
frame's `MB1.forkscale` with **max absolute difference 0** across all 849 shared rows.
The definition is verified, not inferred from the variable name.

Three traps, all recorded in full in `data/menegollo_biclusters/README.md`:

- `forkscale.log` is **`Inf`** for the sample at `index == 1`. The published frame
  carries one such value in MB1 and one in MB3. `Inf` is not `NA`, so it will not be
  removed by `complete.cases()` and will silently destroy a fit.
- **MB3's sign convention differs between cohorts** - METABRIC negates `MB3.pc1` before
  forming forkscale, TCGA does not. The plan uses MB3 as the ER-neutral control axis in
  Block D, so this must be fixed explicitly and stated.
- `forkscale` is severely skewed by construction (observed MB1 range `[-42.6, 8.8]`).

## 4. What this changes downstream

**Block D (script 15) gets better.** Plan section 10 wrote Block D around METABRIC
forkscale and noted that the TCGA PARADIGM analysis "ran on 250 TCGA samples; that is
thin for an interaction with this covariate set". That constraint is gone. Continuous
forkscale is available for 1,037 TCGA patients, **all 1,037 of whom are inside our
1,095-patient covariate table**. Overlap with the 938-patient Block C/B analysis set is
therefore at least 880, to be reported exactly by script 15.

**F3 splits into two steps.** F3 asks whether `STATE` adds anything over
`MB1_forkscale`, i.e. whether this arm's outcome claim is redundant with Menegollo
Fig 7. That is a survival question, and plan section 3 forbids TCGA survival - the
snapshotted CDR gives BRCA 151 OS and 145 PFI events at ~2.3 years median follow-up
(`data/tcga_cdr/README.md`). So:

- **F3-pre, TCGA, exposure side only, runnable now.** With forkscale joined to the
  analysis set, measure directly how much of this arm's axis *is* MB1 forkscale:
  correlation of `MB1.forkscale` with OXPHOS on both instruments, with M-a/M-b/M-c,
  with PRIME and with STATE; and whether the Block C interaction estimate survives
  adding forkscale as a covariate. **No outcome variable is touched.** This neither
  consumes the pre-registered outcome test nor requires `STATE` to have been frozen
  against outcome data. If the axis turns out to be forkscale under another name, F3
  proper is a formality and Block F stops early and cheaply.
- **F3 proper, METABRIC.** `m0/m1/m2` on OS and RFS as plan section 10 specifies, once
  METABRIC and the outstanding file (section 5) are in.

F3-pre is a **diagnostic on the exposure**, not a fifth hypothesis. It is declared here,
before it is run, and it has no outcome variable in it.

## 5. What is still outstanding - METABRIC only

METABRIC forkscale is one file short. The upstream compile script builds it as

```r
sort.df <- data.frame(sample = Complete_METABRIC_Clinical_Features_Data[, 1], ...)
```

The sort and PC1 objects are **positional indices into that clinical frame**, so without
the identifier vector the 1,981 METABRIC forkscale values cannot be attached to patients.
`Complete_METABRIC_Clinical_Features_Data` lives in `METABRIC_DATA.RData`, which is not
in the upstream repository (the upstream README links it on Google Drive).

Any **one** of the following closes the gap, in descending order of preference:

1. `METABRIC_starting_data_final_corr_groups.Rdata` - the compile script's output. Best:
   `all.clinical.df` already carries `sample` plus `MB{1,2,3}.forkscale` and the fork
   calls, so nothing needs reconstructing.
2. `METABRIC_data_w_forkscales` or `essential_METABRIC_MCbiclust_data.Rdata` - later
   outputs of the same script, same content.
3. `METABRIC_DATA.RData` - the raw input. Sufficient, but requires re-running the
   `sort.df` construction here, which we would then have to verify.

TCGA requires nothing further.

## 6. Status

**G3 discharged, PASSES.** Blocker (a) of the 2026-08-29 handoff is closed for TCGA and
reduced to a single named file for METABRIC.
