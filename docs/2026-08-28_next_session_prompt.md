---
date: 2026-08-28
status: disposable - delete once used
purpose: resume prompt for the next Claude Code session
---

# Resume prompt

Paste the block below into a fresh session started plainly in
`/Users/gs/code/myc_human_validation`:

    claude --model opus

---

```
Read CLAUDE.md first, then docs/2026-08-28_G1_result_and_decisions.md.
The plan is docs/2026-08-27_human_validation_plan.md - read sections 2, 6, 10
and 14; you do not need the whole thing.

Context, so you do not re-derive it:

- Setup and the G1/G3 gates are DONE. Do not re-run them, do not re-audit the
  repo state, do not re-check what was snapshotted. It is all recorded in
  docs/2026-08-28_G1_result_and_decisions.md and in git log.
- G1 passed its count criterion: stripped Felsher 61 genes, zero overlap with
  MitoCarta OXPHOS subunits. D2 resolved, M-a stays primary.
- G1's correlation-structure criterion is still OPEN and needs script 01.
  G1 is not fully discharged. Do not describe it as closed.
- G3 passed: continuous forkscale for 849 TCGA samples exists in
  github.com/gszabadkai/Menegollo_Bentham, complete, no NAs.
- The mouse repo is deliberately NOT attached. Do not ask for it. Do not add it
  to additionalDirectories. Read-only access if genuinely needed:
  git -C /Users/gs/G/data/MK_myc_2022/myc_mouse show <ref>:<path>
- CLAUDE.md carries the standing conventions - Option A, if (FALSE) blocks,
  ASCII, scale discipline, species hygiene. I am not repeating them. Follow them.

This session is D4 and then G2. Nothing downstream of G2.

--- Task 1. D4 - Lee et al. 2017 scoping. This BLOCKS G2, do it first. ---

Lee et al. 2017, Cell Metabolism 26:633. Read it properly and answer one
question: does its MYC/MCL1/OXPHOS finding PRE-EMPT hypothesis H1, or SUPPORT
it?

Report back:
  - what exactly it claims, in what cohort or model, with what evidence
  - whether the MYC + MCL1 co-amplification-to-OXPHOS link is already
    established there, or only partially
  - whether our H1 as written in plan section 2 would be novel after it
  - specifically: does it test an INTERACTION, does it test SPECIFICITY against
    other mitochondrial axes, and does it touch apoptotic priming or PUMA/BCL-XL

The consequence, from plan section 14: if it pre-empts, H1 becomes a citation,
H2 becomes the novel mechanism and H4 becomes the novel consequence - which the
plan itself notes is arguably a better paper. If it supports, H1 stands and G2
decides whether it earns Panel b.

Do not skim an abstract and guess. If you cannot get the full text, say so
plainly rather than inferring from the abstract.

--- Task 2. G2 - CNV co-occurrence (plan section 6, script 05) ---

Only after D4 returns. Pure GISTIC, no expression modelling.

In TCGA-BRCA thresholded GISTIC calls:
  - does MYC (8q24.21) amplification co-occur with MCL1 (1q21.3) or BCL2L1
    (20q11.21) amplification more than expected? Fisher / log-odds, stratified
    by PAM50
  - is the co-amplified group enriched for basal / TNBC?

Decision this feeds: if yes, H1 becomes Panel b and H4 gains its stratifying
variable. If no, H2 moves to primary and Panel b becomes the FOXO3/PI3K
decoupling figure.

Before writing script 05, tell me what data it needs and whether each piece is
available - same discipline as the G1 input check. Known issue: TCGAbiolinks is
NOT installed, and PAM50 comes from TCGAbiolinks::PanCancerAtlas_subtypes().
GISTIC all_thresholded.by_genes comes from PanCanAtlas. Neither is in data/ yet.
Stop and tell me what needs downloading or installing before you write anything.

Note G2 needs NO expression data. Do not write or trigger
01_fetch_tcga_expression.R - the STAR counts stay unfetched until the gates
have decided what gets built.

--- Also flag, do not act on ---

  - D7 (proliferation covariate overlaps the MYC exposure, section 6 of the G1
    note) is PROPOSED and must be confirmed before script 09 exists.
  - D5 (primary neoadjuvant cohort) is unconfirmed.
  - Script 07's specificity panel has no agreed source now that the library's
    human GMT tree is rejected.

Remind me of these at the end of the session.
```

---

## State at end of 2026-08-28

Commits on `main`, all pushed to `origin`:

```
53604d0  Fix MitoCarta Sheet-4 padding rows corrupting pathway lookups
6382e77  Add 00_setup_packages.R and 04 (gate G1 overlap audit)
f63f14d  Snapshot the three G1 inputs with provenance
e92e9f1  Snapshot mouse artefacts from myc_mouse 6a9c7dd
```

Generated but gitignored, so NOT on origin - regenerate by re-running script 04
if lost:

```
results/g1_overlap_audit.rds
outputs/tables/g1_overlap_audit.csv
outputs/tables/g1_symbol_harmonisation.csv
```

`data/` is fully committed and on origin, including the 10 MB MitoCarta workbook
and the CollecTRI snapshot. `data/raw/` is still empty - no TCGA download has
been made.

Packages installed this session: `OmnipathR` 4.0.0 (and its dependencies).
Still missing and needed later: `TCGAbiolinks` (scripts 01, 02, and G2's PAM50).
