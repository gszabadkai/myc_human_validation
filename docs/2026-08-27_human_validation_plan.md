---
date: 2026-08-27
version: 2 (supersedes v1 of same date)
tags: [project/myc_mouse, project/human_validation, tcga, metabric, scanb, depmap, mcbiclust, apoptotic-priming, planning]
status: draft-for-review
relates-to:
  - myc_mouse_finalisation_plan.md (AP7 MB-fork validation; this doc is the human-side counterpart)
  - docs/library_reference/2026-08-22_consensus_myc_double_hit_thread.md
  - Menegollo et al. 2024, Cancer Research (companion paper)
next-action: script 01 (discharges G1's correlation criterion)
resolved: D1, D2, D3, D4, D5, D7 (dated notes in docs/)
open: D6, G1's correlation criterion, script 07's specificity panel source
changes-in-v2:
  - added H4 (conditional chemosensitivity) and Block F (clinical outcome)
  - added Block G (DepMap dependency), an orthogonal functional test
  - added section 4 (display-item budget) after D3 resolved: Nature Metabolism Letter,
    2-3 main panels plus Extended Data
  - forkscale replication demoted from centrepiece to Extended Data as a consequence
---

# Human validation of the MYC / OXPHOS / apoptotic-priming axis

Build spec for the human arm. Written for a Claude Code session under the standing
Option A workflow: Claude Code writes and edits scripts, Gyorgy sources them
interactively in Positron. Claude Code does not auto-run numbered pipeline scripts.

---

## 1. What this arm is for, and what it is not for

The mouse result being verified is **not** "MYC drives OXPHOS" and **not** "OXPHOS
correlates with PUMA". It is a three-part conditional claim:

1. The OXPHOS-to-PUMA/BCL-XL coupling is **conditional on MYC** (mouse interaction
   p = 0.0052).
2. The coupling is **specific to OXPHOS subunits** among mitochondrial and redox axes
   ("no such link exists for other redox and metabolic axes").
3. It runs through **FOXO3 and is p53-independent** (iMMEC Tre-MYC cells are p53-null;
   Nutlin fails to induce p21/PUMA/BAX).

Any human analysis that does not test an interaction, does not test specificity, and
does not test p53-independence is not verifying this paper's claim. It is
re-establishing MYC-driven mitochondrial biogenesis, which is already established.

### The survivor-bias problem, stated plainly

The mouse model says the permissive window **closes before tumours form**. TCGA
contains only tumours. Every TCGA sample has already passed the bottleneck the paper
describes. The naive prediction (MYC-high tumours should be OXPHOS-low, having
selected against mitochondria) is contradicted by our own Menegollo data: MB2_UF is
both MYC/miRNA-driven and mitochondria-high.

That contradiction is the asset, not the obstacle. This arm is therefore built around
a different question:

> If OXPHOS-high plus MYC-high is apoptotically lethal in untransformed mammary
> epithelium, how did the human tumours that occupy that state survive it, and what
> did that escape cost them therapeutically?

**Non-goal:** the permissive-window claim itself cannot be made in TCGA. Do not try.
See section 3.

---

## 2. Pre-specified hypotheses and falsification criteria

Written before any model is fitted. Both approaches discussed in chat were, as
proposed, structured so that they could only confirm. This section is the fix.

Let:

- `MYC` = MYC transcriptional activity (three estimators, section 7.1)
- `OXPHOS` = nuclear-encoded OXPHOS subunit level (section 7.2)
- `PRIME` = log2(BBC3) - log2(BCL2L1), the pre-specified confirmatory endpoint
- `BUFFER` = MCL1 / BCL2L1 amplification and/or expression

### H1 - coupling preserved, neutralised downstream

`MYC:OXPHOS` interaction on `PRIME` is positive and significant, **and** the
MYC-high / OXPHOS-high quadrant is enriched for `BUFFER`.

Reading: the trigger relationship is intact in human tumours; survival was bought by
anti-apoptotic buffering. Second hit = BCL2-family amplification.

Prior support: MYC and MCL1 are frequently co-amplified in TNBC, and the pair raise
mitochondrial OXPHOS (Lee et al. 2017, Cell Metab 26:633). This is either a strong
prior in our favour or a partial scooping risk. **Establish which before building the
figure** - see D4.

### H2 - coupling broken by upstream decoupling (PI3K/AKT)

`MYC:OXPHOS` interaction is null or attenuated overall, but is **restored in the
PIK3CA-wild-type / PTEN-intact stratum**, and FOXO3 regulon activity is low in the
MYC-high / OXPHOS-high quadrant.

Reading: PI3K/AKT excludes FOXO3 from the nucleus, switching off the PUMA trigger
without touching OXPHOS. Mechanistically identical to the mouse axis, achieved
genetically. Given PIK3CA is the most commonly mutated gene in breast cancer, this is
arguably the most likely outcome and it is the strongest possible result for the
paper, because it validates the *mechanism* rather than the correlation.

### H3 - p53-independence

The `MYC:OXPHOS` interaction on `PRIME` persists in the **TP53-mutant** stratum.

This is the sharpest available test and TCGA is uniquely suited to it (~35% TP53-mutant,
MC3 calls). If the coupling vanishes in TP53-mutant tumours, the mouse mechanism as
written does not transfer.

### H4 - conditional chemosensitivity (NEW in v2)

**This is the counter-intuitive prediction and the one worth a main panel.**

The naive expectation is that MYC-high / OXPHOS-high tumours are aggressive. The model
predicts the opposite, conditionally:

| State | Trigger present? | Predicted response to apoptotic chemotherapy |
|---|---|---|
| MYC-high, OXPHOS-high, unbuffered | yes | **sensitive** (high pCR, chemo benefit) |
| MYC-high, OXPHOS-high, buffered (MCL1/BCL2L1 amp) | neutralised | resistant |
| MYC-high, OXPHOS-low | absent | resistant |
| MYC-low | n/a | no OXPHOS dependence expected |

So: **OXPHOS-high is protective unless buffered.** Formally, a three-way
`MYC x OXPHOS x BUFFER` structure on a treatment-responsive endpoint, and a
`score x treatment` interaction where randomised or treatment-stratified data exist.

**Tested continuously as primary** (section 7.5, amended 2026-08-28), with the
categorical `STATE` contrast as pre-specified secondary. Powered for a large
conditional effect only - see Block F1 and the D5 note.

Why this framing and not the simpler one:

- **Reduced priming predicts chemoresistance trivially.** That is the entire
  BH3-profiling literature. A marginal priming-to-resistance association proves nothing
  new. The novel content is the *conditionality on OXPHOS*, and that it is MYC-specific.
- **Prognostic association is confounded beyond repair.** MYC-high + OXPHOS-high is
  approximately basal-like, and basal-like is aggressive. Any survival main effect will
  be swamped.
- **Redundancy risk with our own companion paper.** Menegollo Fig 7 already reports that
  MB1_UF (mitochondria-high, non-MYC) predicts worse survival, and that MB1_LF
  OXPHOS/glutamine patterns identify a near-100% 5-year survival group under
  endocrine+chemotherapy. A new mitochondria-high prognostic score that does not add
  information beyond MB1_UF forkscale is not a new result. **Incremental value must be
  tested explicitly** (nested Cox, likelihood ratio, C-index delta). See Block F.

### Falsification

The human arm does **not** support the mouse model if all of the following hold:

- `MYC:OXPHOS` on `PRIME` is null in the full cohort **and** in every pre-specified
  stratum (TP53, PIK3CA, PAM50, forkscale);
- no `BUFFER` enrichment in the MYC-high / OXPHOS-high quadrant at either CNV or
  expression level;
- FOXO3 regulon activity shows no relationship to `OXPHOS` in any stratum;
- **H4 fails in the informative direction**: OXPHOS-high predicts chemo*resistance*
  regardless of buffering status, which would invert the model rather than merely fail
  to support it.

If that is the outcome, report it, or drop the human arm. Do not reach for a fifth
post-hoc hypothesis.

### Specificity requirement (applies to every positive result)

A positive OXPHOS result is only reportable alongside the accompanying negatives:

- **Pathway negatives:** the identical interaction fitted with FAO / carnitine shuttle,
  one-carbon and glycine cleavage, mitoribosome, TCA, ROS defence, and mtDNA-encoded
  OXPHOS held separately. Expected: OXPHOS-subunit-specific.
- **Endpoint negatives:** BID/BCL2L1, BAX/BCL2L1, BCL2L11/BCL2L1, BAK1/BCL2L1.
  The mouse says only PUMA/BCL-XL reverses. Human should show the same.

---

## 3. Dataset assignment

Do not force everything onto TCGA. These are different questions needing different
cohorts.

| Question | Cohort | Why |
|---|---|---|
| Survivor adaptation / second hit (H1-H3) | **TCGA-BRCA** | CNV, MC3 mutation calls, ABSOLUTE purity, RPPA. Nothing else has all four. |
| Replication of the coupling | **METABRIC** | n ~2000, biclusters native, long follow-up. |
| Independent third | **SCAN-B (Brueffer)** | RNA-seq, large n, treatment-stratified. Already used in Menegollo Fig 7H-L. |
| **Chemoresistance (H4), primary** | **I-SPY2-990 (GSE194040), n=988** | Largest public neoadjuvant cohort with pre-treatment expression and complete pCR (319 events, 32.3%). Randomised, 14 arms with a 179-patient paclitaxel control, so `score x treatment` is testable. RPPA companion GSE196093 (736 of the same patients). **D5 RESOLVED 2026-08-28 - see that note; selected by documented search, not assertion.** |
| **Chemoresistance (H4), replication** | **BrighTNess (GSE164458), n=482** | Randomised, TNBC-only, RNA-seq. Tests whether BCL2L1's subtype-breadth (G2) holds in the matched stratum. |
| **Chemoresistance (H4), third** | **GSE25066 (Hatzis), n=508** | Demoted from primary: single-arm, so it cannot test treatment interaction, and its 2011 HG-U133A array carries one multi-mapping `BBC3` probe. Retained as a genuinely independent replication (different platform generation and era). |
| **Chemoresistance (H4), secondary** | **SCAN-B endocrine-only vs endocrine+chemo strata** | Template already exists in the companion paper. Not randomised - see caveat below. |
| **Functional orthogonal test** | **DepMap / CCLE / PRISM / GDSC** | MCL1 and BCL2L1 dependency, and MCL1i/BCL-XLi sensitivity, against MYC and OXPHOS status. Not observational. |
| Permissive window / pre-malignant trajectory | **Hannon (Rebbeck, Xian et al. 2022)** | The only human dataset with the mouse's axis. Count matrix requested. **Out of scope until it lands.** |

### Survival-endpoint caveat, important

**Do not run the survival analysis in TCGA.** TCGA-BRCA is event-poor with short
follow-up; the PanCanAtlas Clinical Data Resource (Liu et al. 2018, Cell) recommends
PFI over OS for BRCA for exactly this reason. Verify the recommended endpoint against
the CDR paper before use. Survival work belongs in METABRIC and SCAN-B.

### Confounding-by-indication caveat

Chemotherapy assignment in METABRIC and SCAN-B is not randomised; sicker patients get
chemotherapy. A `score x chemo` interaction in these cohorts is therefore suggestive,
not causal. This is the reason the neoadjuvant pCR cohorts are the **primary** H4 test
and the treatment-stratified survival analysis is secondary.

### Power caveat

Treatment-interaction tests need roughly four times the sample size of a main-effect
test. Budget accordingly and report CIs rather than leaning on p-values.

---

## 4. Display-item budget (D3 resolved)

**Decision:** the human arm stays in the Nature Metabolism Letter, as 2-3 main panels
plus Extended Data.

Format constraints (nature.com/natmetab/content, checked 2026-08-27):

- Introductory paragraph up to 200 words, referenced, replacing the abstract
- Main text up to 2,500 words, excluding intro, Methods, references, figure legends
- **2-4 display items total for the whole paper**
- Up to 10 Extended Data figures
- ~40 references
- No headings except Methods

### Flag before proceeding

The current draft carries in vivo mouse RNA-seq, iMMEC MYC induction, Bcl-xL rescue,
selected survivors, PGC-1a overexpression across two cell models plus Hs578t and MYAZ,
and an in vivo mitochondria-lowering model. That is already at or past four display
items before the human arm exists. Adding a human figure may force the choice between
compressing the mouse and cell work substantially, or moving to Article format
(8 display items, subheadings allowed, more mechanistic credit). Worth deciding
deliberately rather than discovering at submission.

### Panel allocation for the human arm (3 panels, one display item)

The arc is: the coupling exists, tumours escaped it, escape has a therapeutic cost.

| Panel | Content | Block |
|---|---|---|
| **a** | `MYC x OXPHOS` interaction on PRIME, shown as a coefficient forest alongside the pathway and endpoint negatives. One panel carries both the result and its specificity. | C |
| **b** | The escape route - whichever G2 selects: BUFFER co-amplification (H1) or FOXO3/PI3K decoupling (H2). | B or C-strata |
| **c** | H4 - conditional chemosensitivity. pCR rate by MYC/OXPHOS/BUFFER state in the neoadjuvant cohort, with the DepMap dependency inset. | F + G |

**If only two panels survive:** keep **a** and **c**. They are the load-bearing pair
(the axis is real; it matters clinically). Panel **b** is the most intellectually
interesting and also the most exposed to Lee et al. 2017, so it is the right one to
demote.

### Extended Data allocation

- ED1: MYC estimator concordance (M-a/M-b/M-c), overlap audit from G1
- ED2: **forkscale replication (Block D and D2)** - *demoted from centrepiece*
- ED3: RPPA protein confirmation (Block E)
- ED4: METABRIC and SCAN-B replication of Panel a
- ED5: purity, immune and subtype sensitivity analyses

Note the demotion in ED2 explicitly. Under a 3-panel budget the MB-fork projection
cannot hold a main panel. It is still worth doing, and it is the natural place to cite
the companion paper in the main text ("consistent with the multistate switch described
previously"), but it does not earn main-figure space against the coupling and the
clinical consequence.

---

## 5. Repository setup and species hygiene (D1 RESOLVED)

**Decision:** new sibling repo `myc_human_validation`, at
`/Users/gs/code/myc_human_validation` - **outside the Google Drive sync tree**.

Three reasons, in order of importance:

1. **Species hygiene.** The library README records an explicit prior decision: the human
   GMT tree was deliberately not snapshotted into `myc_mouse` to remove an easy
   wrong-species mistake. The finalisation plan separately flags the Mouse-MitoCarta
   name-clash as one of the two things most likely to bite a Claude Code session that
   globs the data directory. Adding human MitoCarta and human GMTs to that tree makes
   both worse.
2. **Google Drive.** `git_icon_gdrive_issue.md` documents `Icon\r` litter inside `.git`
   at `/Users/gs/G/data/MK_myc_2022/myc_mouse`, recurring on every Drive sync, and calls
   a repo inside a live Drive sync a known corruption risk needing a proper fix. That
   fix was deferred because moving an existing repo has cost. A new repo has none. Take
   it now.
3. **Data volume.** TCGA-BRCA STAR counts, GISTIC, MC3 and RPPA run to several GB.
   Syncing that to Drive burns quota and slows every git operation.

Backup comes from the `origin` remote, not from Drive. Drive mirrors `.git` internals
mid-write and propagates corruption outward; it cannot reconstruct a coherent repo state
from per-file version history. **Caveat:** `results/` and `outputs/` are gitignored and
therefore not on origin. They are regenerable by re-running scripts, so that is
acceptable. `data/raw/` is *not* cheaply regenerable (slow TCGA downloads) and needs its
own backup - Time Machine or an external drive, not a folder holding the `.git`.

### 5.1 Bootstrap

Run once, in Terminal.

```bash
mkdir -p /Users/gs/code/myc_human_validation
cd /Users/gs/code/myc_human_validation
git init

mkdir -p scripts docs functions external \
         data/raw data/from_myc_mouse data/genesets_from_library_human \
         results outputs
touch data/raw/.gitkeep results/.gitkeep outputs/.gitkeep

cat > .gitignore <<'EOF'
# Large downloads and generated artefacts - not on origin, back up data/raw/ separately
data/raw/*
!data/raw/.gitkeep
results/*
!results/.gitkeep
outputs/*
!outputs/.gitkeep

# R session cruft
.Rhistory
.RData
.Rproj.user/

# macOS / Drive litter (insurance; should not occur outside Drive)
.DS_Store
Icon?
EOF

git add .
git commit -m "Initial repo structure for human validation arm"
git branch -M main
```

Then create the remote and push. With the GitHub CLI:

```bash
gh repo create myc_human_validation --private --source=. --remote=origin --push
```

Or, if the empty remote repo was made in the browser first:

```bash
git remote add origin <url>
git push -u origin main
```

Verify the backup actually exists before relying on it:

```bash
git remote -v
git log origin/main -1
```

### 5.2 CLAUDE.md

Save the drafted `CLAUDE.md` to the repo root:
`/Users/gs/code/myc_human_validation/CLAUDE.md`, then commit it.

**Do not import the mouse `CLAUDE.md`** via
`CLAUDE_CODE_ADDITIONAL_DIRECTORIES_CLAUDE_MD=1`. The R coding rules are copied into the
new file deliberately; the mouse branch model, worktree rules and data conventions would
be actively wrong here and would mislead a session.

### 5.3 One-off snapshot session from myc_mouse

`--add-dir` grants **read and write** access, not read-only. Adding `myc_mouse` to every
session hands Claude Code write access to the mouse repo, which is exactly the
contamination this repo split exists to prevent.

So: attach it once, snapshot, detach. Never put it in `additionalDirectories`.

```bash
cd /Users/gs/code/myc_human_validation
claude --model opus --add-dir /Users/gs/G/data/MK_myc_2022/myc_mouse
```

In that session, ask Claude Code to copy into this repo and write provenance:

| From `myc_mouse` | To here | Note |
|---|---|---|
| `scripts/R_CODING_INSTRUCTIONS.md` | `docs/` | Copy, do not import |
| `functions/generate_heatmap.R` | `functions/` | |
| `external/mitotyping/` | `external/mitotyping/` | Monzel et al. mitoPPS reference code. Upstream of `myc_mouse`, not part of it - re-vendoring from the original source is equally valid |
| mouse interaction gene table (`interaction_sig_genes_p05.csv`) | `data/from_myc_mouse/` | Cross-species comparison |
| mouse MYC x OXPHOS effect estimate and mitoPPS pathway ordering | `data/from_myc_mouse/` | Summary tables only, not raw objects |

Then `data/from_myc_mouse/README.md`, following the `genesets_from_library` provenance
pattern exactly: source repo, branch, **commit SHA (dereferenced)**, snapshot date,
snapshot method, and a "do not edit, re-snapshot instead" line.

This is infrastructure work, so under Option A Claude Code may execute it directly.

Afterwards, launch plainly and never re-attach:

```bash
claude --model opus
```

For occasional later reference, this reads a mouse file without granting any directory
access:

```bash
git -C /Users/gs/G/data/MK_myc_2022/myc_mouse show <ref>:<path>
```

### 5.4 Gene sets - two upstreams, only one is myc_mouse

Human GMTs do **not** come from `myc_mouse`; that repo deliberately has none. Snapshot
the human tree from `mammary_geneset_library` at tag `v1.0` (dereferenced commit
`cbd8f16d2b0f95c5d4e86bed6aa112e42538a34b`) into `data/genesets_from_library_human/`,
with a README recording the same tag and commit. Do not rebuild, do not edit in place;
if the library changes, re-snapshot from a new tag.

### 5.5 Git discipline

- `main` is the trunk. Feature branches off `main` as needed.
- **Do not reproduce the mouse branch structure.** `new-analysis` as trunk with `main`
  as historical reference is an artifact of the two-pipeline consolidation in that repo,
  not a model to copy.
- Not a submodule of `myc_mouse`. The two repos move independently and submodule
  coupling buys nothing here.
- No worktrees until there is a reason for one.
- Commit per verified phase, as in `myc_mouse`. Git is the safety net.

---

## 6. Gates (run first, cheap, they decide what gets built)

### G1 - overlap audit (~1 hour)

Compute and report:

- `|Felsher MYC signature INTERSECT MitoCarta 3.0 human|`, as a fraction of the signature
- `|Felsher INTERSECT Hallmark E2F_TARGETS / G2M_CHECKPOINT|`
- the same for the CollecTRI MYC regulon

Then produce MitoCarta-stripped versions of both MYC estimators, recording how many
genes were removed.

**Decision:** if the stripped Felsher signature falls below ~50 genes or loses its
correlation structure, the signature-ranking approach is not salvageable as a primary
estimator and the CollecTRI regulon plus the 8q24 CNV instrument carry the MYC axis.

Without this step, any MYC-to-OXPHOS result is partly definitional.

### G2 - CNV co-occurrence (~half day)

Pure GISTIC, no expression modelling. In TCGA-BRCA thresholded calls:

- Does `MYC` (8q24.21) amplification co-occur with `MCL1` (1q21.3) or `BCL2L1`
  (20q11.21) amplification more than expected? Fisher / log-odds, stratified by PAM50.
- Is the co-amplified group enriched for basal / TNBC?

**Decision:** if yes, H1 becomes Panel b and H4 gains its stratifying variable. If no,
H2 moves to primary and Panel b becomes the FOXO3/PI3K decoupling figure. Either way
this is one afternoon and it reshapes everything downstream.

### G3 - forkscale availability (~1 hour)

Check `github.com/gszabadkai/Menegollo_Bentham` for stored TCGA bicluster assignments
and forkscale values.

**Decision:** if present, ED2 is cheap and high-fidelity. If absent, re-deriving via
MCbiclust is a genuine fidelity risk (different run, seed, sample QC) and must be
documented as a limitation. Under the 3-panel budget, a failed G3 is a reason to drop
ED2 rather than to spend a week on it.

---

## 7. Measurement definitions

### 7.1 MYC activity - three estimators, concordance required

| ID | Estimator | Role |
|---|---|---|
| M-a | Felsher signature, MitoCarta-stripped, GSVA on VST | **Primary.** Cross-species continuity with the mouse arm. |
| M-b | CollecTRI / DoRothEA MYC regulon via `decoupleR` (VIPER or ULM), MitoCarta-stripped | Less circular. Required concordance check. |
| M-c | MYC 8q24.21 GISTIC amplification (categorical) | Quasi-instrument, assigned independently of the transcriptome. |

**Rule:** the headline claim requires directional concordance across all three. If M-a
and M-b disagree, report both and treat the claim as unsupported.

### 7.2 Mitochondrial axes

- **Level (primary):** mean z-score of nuclear-encoded MitoCarta 3.0 OXPHOS subunits.
  Per standing convention the 13 mtDNA-encoded protein-coding genes sit in a separate
  synthetic "mtDNA-encoded OXPHOS subunits" pathway and are never pooled with the
  nuclear set (expression-scale skew).
- **Shape (secondary):** mitoPPS. It normalises each pairwise pathway ratio by the
  global average across samples, so it reports the *shape* of the mitochondrial program
  and is deliberately robust to total content. Part of the mouse claim is about OXPHOS
  *level*. Report both. **Never compare mitoPPS values numerically across cohorts or
  species** - the baseline is composition-dependent.
- **Specificity panel:** FAO / carnitine shuttle, one-carbon and glycine cleavage,
  mitoribosome, TCA, ROS defence. Size-matched where possible.

### 7.3 Priming endpoints

- **Confirmatory (pre-specified):** `PRIME = log2(BBC3) - log2(BCL2L1)`. Single, named
  in advance, mirroring the mouse pre-specification argument. This is what protects the
  result from a multiple-comparison criticism.
- **Robust secondary:** summed pro-apoptotic BH3-only plus effectors over summed
  anti-apoptotic guardians, z-scored.
- **Negative-control endpoints:** as listed in section 2.

**Known weakness, state it in Methods:** BBC3 mRNA is a poor proxy for PUMA protein.
The internal handoff already documents RNA/protein discordance. This is why ED3 (RPPA)
is not optional.

### 7.4 FOXO3 - activity, not expression

FOXO3 is regulated by nuclear exclusion; its own mRNA is a bad activity readout. Use a
**FOXO3 target-gene regulon score** (CollecTRI / DoRothEA via `decoupleR`), with
mitochondrial genes excluded. Report FOXO3 mRNA alongside as a contrast, not as the
measure.

### 7.5 The composite state variable for H4

Pre-specify **before** looking at outcome data. A three-way conjunction on median splits
gives eight cells; picking the worst post hoc is p-hacking.

> **AMENDED 2026-08-28 (D5).** The **primary** H4 test is now the **continuous**
> three-way `MYC x OXPHOS x BUFFER` interaction. `STATE` below and its
> level-3-vs-level-4 contrast are the **pre-specified secondary**, reported alongside.
>
> Reason: simulation at n=988 shows median splits cost up to **53 percentage points**
> of power against the identical truth (D5 note section 5.2). This amendment is made
> before script 11 exists, before any H4 outcome data has been touched, and from a
> simulation containing no project data - it is a power decision, not a response to a
> result. Both specifications are fixed now and both will be reported.
>
> `STATE` remains frozen in script 11 and must still not be revised after outcome data
> is seen. Its **role** changes; its **definition** does not.

Fixed definition:

```
STATE = factor, four levels, in this order:
  1. MYC_low                                       (reference)
  2. MYC_high & OXPHOS_low
  3. MYC_high & OXPHOS_high & BUFFER_absent        (predicted chemo-SENSITIVE)
  4. MYC_high & OXPHOS_high & BUFFER_present       (predicted chemo-RESISTANT)
```

Splits at cohort median for MYC and OXPHOS; BUFFER by GISTIC amplification where CNV
exists, else by expression tertile.

Secondary H4 contrast: **level 3 vs level 4**. That single contrast is the whole
prediction and it is one degree of freedom. It is reported alongside the continuous
primary, never instead of it, and never selected between after the fact.

---

## 8. Covariates and their sources

Breast is the worst tissue in TCGA for mitochondrial confounding. Adipose is OXPHOS
and FAO high; immune infiltrate carries its own BCL2-family profile.

| Covariate | Source | Why |
|---|---|---|
| Tumour purity | ABSOLUTE (Aran 2015 / PanCanAtlas), via `TCGAbiolinks` | Non-negotiable for any mito or apoptosis score in breast. |
| Leukocyte / stromal fraction | Thorsson 2018 immune landscape (PanCanAtlas) | Already computed for all TCGA. Cheaper than CIBERSORTx. |
| PAM50 / ER status | `TCGAbiolinks::PanCancerAtlas_subtypes()` | **BCL2 is estrogen-responsive.** See section 9. |
| Proliferation index | Hallmark E2F_TARGETS + G2M_CHECKPOINT GSVA. **Two versions, see D7 note.** `PROLIF_STD` = the union, 327 genes. `PROLIF_DISJOINT` = 318 genes, the 9 shared with the stripped Felsher estimator removed. | Nuclear mito genes track proliferation; so does any MYC score. **The covariate shares 14.8% of M-a's genes**, so adjusting for it partially adjusts away the exposure - overlapping measurement, not confounding. |
| TP53, PIK3CA, PTEN status | MC3 MAF (PanCanAtlas) | H2 and H3 strata. |
| GISTIC calls (MYC, MCL1, BCL2L1) | PanCanAtlas GISTIC2 `all_thresholded.by_genes` | G2, BUFFER, H4. |
| Stage, grade, nodal status, size | TCGA clinical / METABRIC | Block F covariates. |
| Treatment (chemo, endocrine, radiotherapy) | METABRIC / SCAN-B annotation | Block F treatment interaction. |
| Plate / batch | TCGA barcode | Standard. |

Sensitivity: repeat the primary model in the purity-high subset (ABSOLUTE > 0.7). If
the interaction exists only at low purity, it is stroma.

Expression handling: GDC harmonised STAR counts via `TCGAbiolinks` (or recount3).
Raw counts -> VST for GSVA/ssGSEA. DESeq2-normalised **linear** counts for mitoPPS.
Opposite scale requirements; the two must not share an input object.

---

## 9. The ER confound (read before touching Block D)

MB2_UF is ER-negative by construction in Menegollo. BCL2 is a canonical
estrogen-responsive gene and ER+ tumours are BCL2-high. **Any BCL2-family ratio
compared across MB2_UF and MB1_UF is partly an ER readout, not a MYC readout.**

Mitigations, all three applied:

1. Restrict primary BCL-family measures to **BCL2L1 and MCL1**, not BCL2. This also
   matches the mouse, which is a BCL-XL story.
2. Adjust for ER status, and additionally fit within-ER-status.
3. Use **MB3 as the control axis**. MB3 is an independent switch (a switch in
   mitochondrial *function*, where MB1/MB2 are switches in *biogenesis*, Menegollo
   Fig 1D), so it dissociates mitochondrial content from the MYC/ER axis. Gyorgy's
   instinct to include MB3 was right; this is the reason to state in Methods.

---

## 10. Models

### Block C - primary interaction (TCGA) -> Panel a

```
PRIME ~ MYC * OXPHOS
        + purity + leukocyte_fraction + <proliferation, see below>
        + PAM50 + TP53_status + plate
```

Test: the `MYC:OXPHOS` coefficient, two-sided. Report effect size and CI; with n ~1000
almost anything reaches significance.

> **D7 RESOLVED 2026-08-28 - the proliferation term is estimator-specific.** The
> covariate shares **14.8%** of M-a's genes, so adjusting for it partially adjusts away
> the exposure and biases `MYC:OXPHOS` toward the null. The fix is asymmetric: removing
> the shared genes costs **2.8%** of the covariate for M-a but **24.5%** for M-b, and
> for M-b it strips exactly the canonical proliferation machinery (`AURKA`, `CDK1`,
> `PCNA`, `MYBL2`, ...). So:
>
> | Fit | Proliferation term | Applies to |
> |---|---|---|
> | **S1 primary** | `PROLIF_DISJOINT` (318 genes) | **M-a** |
> | **S1 primary** | `PROLIF_STD` (327 genes) | **M-b** |
> | **S2 sensitivity** | none - term dropped | both |
> | **S3 sensitivity** | `PROLIF_STD` (327 genes) | both |
>
> **Report all three side by side; this is part of the result, not an appendix.**
> Survives all three = robust to construction. Appears only in S2 = report that
> plainly, never select it. `PROLIF_DISJOINT` is defined against the *stripped*
> Felsher set and must never be reused for M-b. See
> `docs/2026-08-28_D7_proliferation_covariate.md`.

Then in order:

1. **Specificity battery** - refit with each pathway negative in place of OXPHOS, and
   each endpoint negative in place of PRIME. Coefficient forest, which is Panel a.
2. **Stratified refits** - TP53-mutant vs wild-type (H3); PIK3CA/PTEN-altered vs intact
   (H2); within PAM50.
3. **Instrument version** - replace `MYC` with M-c (8q24 amplification).
4. **Purity-high sensitivity.**

### Block B - buffering test -> Panel b (if G2 positive)

```
BUFFER_amp  ~ MYC_high * OXPHOS_high + PAM50 + purity      # logistic, CNV
MCL1_expr   ~ MYC * OXPHOS + covariates
BCL2L1_expr ~ MYC * OXPHOS + covariates
```

### Block F - clinical outcome (NEW in v2) -> Panel c

**F1. Chemoresistance, primary. Neoadjuvant cohorts, pCR endpoint.**

```
# PRIMARY (continuous, D5 2026-08-28)
pCR ~ MYC * OXPHOS * BUFFER + stage + grade + ER_status + arm    # logistic

# SECONDARY (categorical, section 7.5)
pCR ~ STATE + stage + grade + ER_status + arm                    # logistic
```

Primary test: the three-way `MYC:OXPHOS:BUFFER` coefficient. Secondary contrast:
`STATE` level 3 (unbuffered) vs level 4 (buffered). Prediction: level 3 has **higher**
pCR than level 4, and higher than level 2 (OXPHOS-low). Report both, always.

`arm` is not optional in I-SPY2: 14 arms spanning 17.3% to 59.1% pCR, several
HER2-targeted and so confounded with subtype. Adjust or stratify.

**Three cohorts, meta-analysed - not one cohort with two afterthoughts.** Score each
independently (GSVA is cohort-relative) and meta-analyse the three interaction
estimates. Pooling scores across cohorts remains forbidden. Report per-cohort
estimates, the pooled estimate, and heterogeneity.

TNBC-only sensitivity in BrighTNess. Note this now tests whether BCL2L1's
subtype-breadth (G2) holds in the matched stratum, rather than assuming a subtype
match - see the D5 note section 1.

**Powered for a large effect only.** Pre-registered: at a target-cell risk ratio of
2.0 the meta-analysis sits at 86% power and any single cohort below 60% (D5 note
section 5). Report confidence intervals, not p-values alone. A null H4 is
"not powered to exclude a modest effect", **not** falsification.

**The two-way is not a fallback.** Dropping `BUFFER` buys about five percentage
points of power while sacrificing the entire distinction from Lee et al. 2017. Do not
retreat to it.

**F2. Treatment interaction, secondary. METABRIC and SCAN-B.**

```
Surv(time, event) ~ STATE * chemotherapy + age + stage + grade + nodes + PAM50
```

Test the `STATE:chemotherapy` interaction, not the main effect. Prediction: the
survival benefit of chemotherapy is largest in level 3 and smallest in level 4.

State the confounding-by-indication limitation explicitly. This is supporting evidence
for F1, not a substitute.

**F3. Incremental value over the companion paper. Mandatory.**

```
m0: Surv(...) ~ clinical covariates
m1: m0 + MB1_forkscale
m2: m1 + STATE
```

Likelihood-ratio test m2 vs m1, plus C-index delta with bootstrap CI. If `STATE` adds
nothing over `MB1_forkscale`, the outcome claim is redundant with Menegollo Fig 7 and
should not be made. **Run this before building Panel c.**

**F4. Aggressiveness, descriptive only.**

Association of `STATE` with grade, size, nodal status, and nuclear pleomorphism,
following the Menegollo Fig 7A-C template. Supplementary, not a panel. Prognostic
association alone is confounded by subtype and is not the claim.

### Block G - DepMap functional dependency (NEW in v2) -> Panel c inset

Orthogonal to everything above because it is functional rather than observational, and
it is immediately available.

In CCLE/DepMap breast lines:

1. Score `MYC` and `OXPHOS` as in section 7 (CCLE expression). Menegollo already built
   CCLE GSVA scoring infrastructure for fork assignment - reuse it.
2. Regress **CRISPR gene-effect (Chronos)** for `MCL1` and `BCL2L1` on `MYC * OXPHOS`.
3. Repeat with **drug sensitivity** (PRISM / GDSC): MCL1 inhibitors (S63845, AMG-176),
   BCL-XL inhibitors (A-1331852, navitoclax), venetoclax as a BCL2 specificity control.

Prediction: MYC-high / OXPHOS-high lines are selectively dependent on MCL1 and/or
BCL2L1. A negative result here weakens H1 considerably, which is what makes it worth
running.

This also supplies the translational hook, which is what a Nature Metabolism editor
will look for.

### Block D - forkscale replication -> ED2 (demoted)

Use **continuous forkscale**, not fork membership. The PARADIGM analysis in Menegollo
Fig 3A ran on 250 TCGA samples; that is thin for an interaction with this covariate set.

```
PRIME ~ MYC * OXPHOS * MB2_forkscale + covariates
PRIME ~ MYC * OXPHOS * MB3_forkscale + covariates    # ER-neutral control axis
```

### Block D2 - developmental analogue via cell of origin -> ED2

```
PRIME ~ MYC * OXPHOS * (LP_score - mL_score) + covariates
```

LP maps to MB2_UF, mL to MB1_UF. The closest available structural translation of the
mouse stage effect. Frame as an analogue, not an equivalence.

### Block E - RPPA confirmation -> ED3

TCPA / PanCanAtlas RPPA level 4 for TCGA-BRCA carries BCL-XL, BAX, BAK and
caspase-cleavage readouts. Repeat the primary interaction with the protein endpoint.
This is the answer to the documented RNA/protein discordance.

---

## 11. Script plan

New repo `myc_human_validation` (pending D1). Same conventions as `myc_mouse`.

```
00_setup_packages.R
01_fetch_tcga_expression.R        # STAR counts, VST, DESeq2-normalised linear
02_fetch_tcga_genomics.R          # MC3, GISTIC, RPPA, ABSOLUTE, Thorsson
03_build_covariate_table.R        # one sample x covariate table + QC report
04_snapshot_human_genesets.R      # library v1.0 human GMTs + overlap audit   <- G1
05_gate_cnv_cooccurrence.R        # pure GISTIC                               <- G2
06_score_myc_activity.R           # M-a, M-b, M-c
07_score_mitochondrial.R          # OXPHOS level, mitoPPS, specificity panel,
                                  #   PROLIF_STD + PROLIF_DISJOINT (D7)
08_score_priming.R                # PRIME, robust index, negatives, FOXO3 regulon
09_interaction_models.R           # Block C S1/S2/S3 (D7) + specificity
                                  #   battery + strata                      -> Panel a
10_buffering_models.R             # Block B                                  -> Panel b
11_build_state_variable.R         # section 7.5, frozen before Block F
12_fetch_neoadjuvant_cohorts.R    # GSE194040 primary, GSE164458, GSE25066
13_outcome_models.R               # Block F1-F4                              -> Panel c
14_depmap_dependency.R            # Block G                                  -> Panel c inset
15_forkscale_replication.R        # Block D + D2                             -> ED2
16_rppa_confirmation.R            # Block E                                  -> ED3
17_metabric_scanb_replication.R   # ED4
18_figures.R
```

Results as `.rds` in `results/`. Figures to named subdirectories under `outputs/`.
Every script carries an `if (FALSE)` sandbox block. ASCII-only strings; handle any
latin1/cp1252 at read time with `fileEncoding = "latin1"` and `iconv()`.

Standing R rules apply: never `print(n = X)` after `head()`; always `dplyr::count()`.

---

## 12. Execution order

```
G1 (04) -> G2 (05) -> G3 (external check)
   |
   +-- if G2 positive: H1 -> Panel b. Build 06,07,08 -> 10 -> 09
   +-- if G2 negative: H2 -> Panel b. Build 06,07,08 -> 09 with PIK3CA strata,
                                       FOXO3 regulon foregrounded
   |
   +-- H3 (TP53 strata) runs inside 09 either way. Cheap, most mechanistically
       specific test available. Do not defer it.

Then, for Panel c:
   11 (freeze STATE) -> F3 incremental-value check FIRST
       |
       +-- if STATE adds nothing over MB1_forkscale: STOP. No Panel c.
       +-- else: 12 -> 13 (F1 primary, F2 secondary) -> 14 (DepMap)

14 (DepMap) can run in parallel from the start - it needs none of the TCGA work
   and it is the cheapest evidence in the whole plan.

15 (forkscale) and 16 (RPPA) after 09, as Extended Data.
17 (METABRIC/SCAN-B) last, as replication, not discovery.
```

Do not build 09 before G1 and G2 return. Do not build 13 before F3 returns.

---

## 13. Known landmines

- **Species contamination.** Human MitoCarta, human GMTs, human ortholog tables. This is
  why the sibling repo is recommended.
- **Scale confusion.** GSVA wants log-scale (VST, `kcdf = "Gaussian"`). mitoPPS wants
  linear DESeq2-normalised counts. Opposite requirements. Enforce with an explicit
  comment in 07.
- **GSVA cohort-relativity.** Score all samples of a cohort in one run. Not portable
  across separately-scored cohorts. This matters most in Block F, where five cohorts are
  scored independently and must never be pooled at the score level - meta-analyse the
  effect estimates instead.
- **mitoPPS baseline drift.** Dataset composition sets the "normal". Only the pattern
  transfers across cohorts.
- **Bulk composition.** A bulk score cannot separate "more cells running the program"
  from "same cells upregulating it". Same constraint as the mouse arm. State findings as
  *associated*, not *driven*.
- **Three-way conjunction as p-hacking.** Section 7.5 exists to prevent this. Freeze
  STATE in script 11 and never revise it after seeing outcome data.
- **Redundancy with Menegollo Fig 7.** F3 is the guard. Run it first.
- **Re-derived forkscale.** If G3 fails and MCbiclust is re-run, the result is a new
  bicluster solution, not the published one.

---

## 14. Open decisions

- **D1 - RESOLVED.** New sibling repo `myc_human_validation` at
  `/Users/gs/code/myc_human_validation`, outside the Google Drive sync tree. Bootstrap,
  CLAUDE.md placement and the one-off `myc_mouse` snapshot procedure are in section 5.
- **D2 - primary MYC estimator.** M-a Felsher-stripped (cross-species continuity) vs
  M-b CollecTRI regulon (less circular). Spec says M-a primary with mandatory M-b/M-c
  concordance. Confirm or flip.
- **D3 - RESOLVED.** Human arm stays in the Nature Metabolism Letter as 2-3 main panels
  plus Extended Data. See section 4. Consequences: forkscale work demoted to ED2; the
  overall display-item budget for the paper needs a decision (section 4, flag).
- **D4 - Lee et al. 2017 scoping.** Read properly and establish whether the
  MYC/MCL1/OXPHOS finding pre-empts H1 or supports it. If it pre-empts, H1 becomes a
  citation, H2 becomes the novel mechanism, and H4 becomes the novel *consequence* -
  which is arguably a better paper. **Do this before G2.**
- **D5 - RESOLVED 2026-08-28.** The original candidate set was undocumented, so D5 was
  reopened rather than confirmed. A pre-specified search (GEO E-utilities 447 series ->
  34 passing; PubMed; Consensus) gives: **I-SPY2-990 (GSE194040, n=988) primary**,
  BrighTNess (GSE164458, n=482) replication, GSE25066 (n=508) third, I-SPY1 (GSE22226,
  n=150) dropped. GSE25066 is demoted because it is single-arm and so cannot test
  treatment interaction. Two consequent decisions: the H4 primary test becomes
  **continuous** (section 7.5 amended), and H4 uses **all three cohorts meta-analysed**.
  See `docs/2026-08-28_D5_cohort_selection.md`.
- **D7 - RESOLVED 2026-08-28.** Raised by G1: the proliferation covariate shares 14.8%
  of M-a's genes, so adjusting for it partially adjusts away the exposure. Overlapping
  measurement, not confounding. Handling is **estimator-specific** because the fix is
  asymmetric - a disjoint covariate costs 2.8% of the covariate for M-a but 24.5% for
  M-b. M-a's primary model takes `PROLIF_DISJOINT`; M-b keeps `PROLIF_STD`; all three
  specifications reported for both. See
  `docs/2026-08-28_D7_proliferation_covariate.md`.
- **D6 - NEW. Does the mouse arm have a metastasis phenotype to anchor F4 to?** The
  manuscript's framing mentions metastatic capacity. If there is a mouse readout, F4
  gains a cross-species anchor and might justify supplementary space. If not, drop F4.

---

## 15. What would make this arm fail well

If the coupling does not survive in human tumours, the defensible framing is already
available: established tumours are survivors of the bottleneck, so the *absence* of the
coupling in TCGA is consistent with the window having closed. That framing is only
credible if the falsification criteria in section 2 were written first, and if the
Hannon pre-malignant arm is named as the proper test.

If H4 fails in the informative direction - OXPHOS-high predicts chemoresistance
regardless of buffering - the model is inverted rather than unsupported, and that is a
result worth reporting honestly rather than burying. It would mean OXPHOS in established
human tumours does something other than set the apoptotic threshold.

That is the reason for this document.
