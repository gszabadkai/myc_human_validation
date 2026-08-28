---
date: 2026-08-28
status: D4 RESOLVED (option 1); G2 design decisions FIXED
relates-to:
  - 2026-08-27_human_validation_plan.md (sections 2, 3, 6, 8, 14)
  - 2026-08-28_G1_result_and_decisions.md
decides:
  - D4 - RESOLVED. Supports, not pre-empts. H1 stays a hypothesis, novelty
    relocated onto PRIME
  - G2 amplitude thresholds - FIXED, gene-specific and direction-specific
  - G2 sample definition and TNBC rule - FIXED
  - G2 pass criterion - FIXED
  - BCL2L1 carries the Panel b decision if the two partners split
  - BBC3 - added to G2 by amendment, tested in the LOSS direction
next-action: run script 05 (G2)
---

# D4 scoping, and the G2 design decisions that follow

---

## Part 1. D4 - Lee et al. 2017

Lee KM, Giltnane JM, Balko JM, ... Arteaga CL. "MYC and MCL1 Cooperatively
Promote Chemotherapy-Resistant Breast Cancer Stem Cells via Regulation of
Mitochondrial Oxidative Phosphorylation." Cell Metab 2017;26(4):633-647.e7.
PMID 28978427. PMC5650077. doi:10.1016/j.cmet.2017.09.009

Full text read from PubMed Central, not the abstract. Passages where the text
extractor dropped italicised gene symbols were re-verified against the PMC HTML.

### 1.1 What the paper establishes

Claim chain: MYC + MCL1 co-amplification -> increased mtOXPHOS -> increased ROS
-> HIF-1alpha accumulation -> cancer stem cell (CSC) expansion -> chemoresistance.

| Layer | Evidence | Strength |
|---|---|---|
| MYC/MCL1 co-amplification in resistant TNBC | 18 patients, serial biopsies, Foundation Medicine 236-gene capture. 8/18 (44%) co-amplified in at least one biopsy; 17/18 (94%) MYC and/or MCL1. Prior Balko cohort: post-NAC residual MYC 54%, MCL1 35%, and 83% of MYC-amplified co-amplified MCL1 | descriptive, small n |
| MYC -> mtOXPHOS | DOX-inducible MYC (MDA-MB-468), siMYC (MDA-MB-436, SUM159PT). Seahorse OCR, TEM, mtDNA content and mass | causal, cell line |
| MCL1 -> mtOXPHOS | MCL1-dMTS mutant; BH3-mimetic VU0659158; LC-MS/MS TCA intermediates | causal, cell line; the paper's strongest result |
| Cooperativity | double knockdown exceeds either single; MYC induction on MCL1 background exceeds GFP background | causal, cell line |
| Human tumours | GSVA of TCGA-BRCA n=816: hallmark OXPHOS, ROS and hypoxia higher in MYC- and/or MCL1-altered tumours, and higher in both-altered than either alone. Also METABRIC TNBC n=320, CCLE n=20, NanoString n=34 | correlative, and **supplementary** (Fig. S6 / Table S1) |
| Reversal | digoxin, N-acetylcysteine, S3QEL2, HIF1A siRNA reduce CSCs and re-sensitise to paclitaxel; xenograft ELDA | causal, in vivo |

All 18 tumours were TP53-mutant. Incidental to their argument, but it means the
axis is demonstrated in a p53-null background - a weak prior favouring H3's
premise.

### 1.2 The three tests that matter to us

**Interaction: NO.** Statistics are Student's t throughout, with one-way ANOVA
for the multi-group GSVA comparison. There is no interaction term and no
continuous model anywhere in the paper. The "both exceed either alone" claim is a
four-group ANOVA on a GSVA score - a super-additivity comparison, not an
interaction test.

**Specificity against other mitochondrial axes: NO, and one result cuts against
us.** The pathway work is hallmark GSVA (OXPHOS, ROS, hypoxia, EMT). Hallmark
OXPHOS is not MitoCarta OXPHOS subunits, and nothing is held against FAO,
one-carbon, mitoribosome, TCA or mtDNA-encoded sets. Worse: their XF Mito Fuel
Flex result shows TNBC CSCs are preferentially **FAO-dependent** (etomoxir
suppresses OCR more than BPTES or UK5099). FAO / carnitine shuttle is one of our
pre-specified *pathway negatives*. **Carry this into the specificity battery
alongside the MTHFD2 caveat** in the G1 note - if the FAO negative fires, Lee et
al. supplies a ready prior explanation and we must argue our way out of it.

**Apoptotic priming, PUMA, BCL-XL: NO - and deliberately routed around.** This is
the most important finding for this arm. The paper's central mechanistic move is
to *dissociate* MCL1's respiratory function from its anti-apoptotic function:

- MCL1-dMTS reduces mtOXPHOS, ROS and CSCs **without affecting anti-apoptotic
  activity**;
- the BH3-mimetic VU0659158 **induces apoptosis but does not reduce OCR or
  mammosphere formation**.

`PUMA` appears twice, both in background prose listing MCL1's binding partners;
never measured. `BBC3` never appears. `BCL2L1` appears exactly once in the whole
paper - as an anti-Bcl-xL antibody used as a mitoplast marker in submitochondrial
fractionation. No BH3 profiling, no priming assay, no PRIME-like ratio.

### 1.3 Verdict

**SUPPORTS, with H1 clause 2 partially pre-empted for MCL1 only.**

- **H1 clause 1** (`MYC:OXPHOS` interaction on `PRIME`) is **not pre-empted**.
  The endpoint does not exist in Lee et al.
- **H1 clause 2** (BUFFER enrichment in the MYC-high / OXPHOS-high quadrant) is
  **substantially pre-empted for MCL1**, in TCGA, by their Fig. S6. Same 2x2,
  opposite direction of conditioning. **Not pre-empted for BCL2L1.**

So the novelty of H1 rests entirely on `PRIME`. Framed as "MYC/MCL1
co-amplification associates with OXPHOS", it is their supplementary figure.
Framed as "MYC and OXPHOS jointly raise apoptotic priming, and BCL2-family
amplification neutralises it", it is new.

Lee et al. establishes MYC/MCL1 -> OXPHOS and then argues it is *not* about
apoptosis. This arm claims it *is*. Those are competing mechanisms for the same
lesion, and ours is testable against theirs. That is the opening.

### 1.4 What they did NOT do that G2 does

**Lee et al. never test co-occurrence against a null.** They report conditional
frequencies (83% of MYC-amplified are MCL1 co-amplified) and compare their
cohort's marginal frequencies against TCGA's (MYC 77.7% vs 44%, MCL1 66.6% vs
24%, both p<0.0001 Fisher). Neither is an excess-over-expected test. Whether MYC
and MCL1 amplification co-occur more than chance in TCGA-BRCA is unanswered in
the paper.

G2 therefore adds, in order of importance:

1. **BCL2L1 / 20q11.21**, untouched by Lee et al. and the `PRIME` denominator;
2. formal co-occurrence testing against a null, conditioned on aneuploidy;
3. PAM50 and TNBC stratification.

G2 is not redundant. But its selling point is no longer "are MYC and MCL1
co-amplified" - that is accepted - and the script and any figure must lead with
BCL2L1 and with excess-over-expected.

### 1.5 Consequence for H4 - sharpened, not threatened

Lee et al. predicts OXPHOS-high -> chemo*resistance*. H4 predicts OXPHOS-high is
protective *unless buffered*. These are not opposed: Lee's OXPHOS-high state is
by construction the MCL1-co-amplified state, i.e. precisely the **buffered row**
of the H4 table, which they populate with prior evidence and an independent
mechanism (CSC / HIF-1alpha rather than buffering, same direction).

H4's novel content is the **unbuffered row** - MYC-high / OXPHOS-high *without*
MCL1/BCL2L1 amplification being chemo-*sensitive*. Lee et al. never test that
stratum, because their cohort is selected on the co-amplification.

Carry forward honestly: plan section 2's falsification criterion says H4 failing
"in the informative direction" - OXPHOS-high resistant regardless of buffering -
would invert the model. **Lee et al. is substantive prior evidence in that
inverting direction.** H4's buffering stratification is therefore load-bearing,
not a robustness check, and D5 (cohort choice) must be settled with that in mind.

### 1.6 D4 RESOLVED - Option 1, confirmed 2026-08-28

**H1 stays a hypothesis. Its novelty is relocated onto `PRIME` and must be stated
that way in the text. Lee et al. is cited in the introduction as prior support
for the upstream half - MYC/MCL1 to OXPHOS - which we do not claim to discover.**

The alternative offered by plan section 14 (H1 demoted to a citation, H2 promoted
to the novel mechanism, H4 to the novel consequence) is **not** taken. Lee et al.
pre-empts the association but not the endpoint, and the endpoint is where this
arm's claim lives.

Consequences that bind downstream writing:

- No sentence may present "MYC/MCL1 co-amplification associates with OXPHOS" as a
  finding of this work.
- Any H1 result is stated as an apoptotic-priming claim, with the OXPHOS
  association attributed to Lee et al.
- The MCL1 arm of G2 is framed as **replication**; the BCL2L1 arm is the novel
  test. See 2.8.

---

## Part 2. G2 design decisions - FIXED 2026-08-28

All fixed before any co-occurrence statistic was computed, except as disclosed in
2.6.

### 2.1 GISTIC source - option C, both runs

| | A - PanCanAtlas ISAR (PRIMARY) | B - Firehose BRCA-TP (SENSITIVITY) |
|---|---|---|
| Dimensions | 24,203 x 9,991 pan-cancer | 24,776 x 1,080 BRCA-only |
| BRCA patients after rule | 1,043 | 1,080 |

ISAR is primary because 8q and 1q are among the most frequently gained arms in
breast cancer - MCL1 is gained-or-amplified in 76% of tumours - so an uncorrected
co-occurrence test largely measures aneuploidy. ISAR corrects arm-level and
ploidy background, which is the confound in play.

Full provenance, SHA-256 sums and read-time traps: `data/gistic_tcga_brca/README.md`.

### 2.2 Sample definition

BRCA membership from the **GDC case list** (1,098, snapshotted), independently of
either GISTIC run, then sample type `-01` only, then assert one column per
patient.

The 7 duplicated ISAR patients are each a primary plus a metastasis (`-06`), so
the `-01` filter removes them without an arbitrary tie-break. A is a strict
subset of B (A n B = 1,043, A-only = 0, B-only = 37), so B must additionally be
run on the 1,043 intersection to isolate calling method from sample set.

### 2.3 Amplitude thresholds - GENE-SPECIFIC

**Decision: gene-specific, author's call, overriding my recommendation.**

| Gene | Direction | Primary threshold | Base rate (ISAR, n=1,043) |
|---|---|---|---|
| `MYC` | gain | `== +2` | 225 (21.6%) |
| `MCL1` | gain | `== +2` | 174 (16.7%) |
| `BCL2L1` | gain | `>= +1` | 479 (45.9%) |
| `BBC3` | **loss** | `<= -1` | 247 (23.7%) |
| `BAX` (control) | **loss** | `<= -1` | 226 (21.7%) |

Rationale, as given: at `+2`, `BCL2L1` is n=28 and uninterpretable, and BCL2L1 is
this arm's *a priori* primary focus - it is the `PRIME` denominator, and the term
Lee et al. never touch. Keeping MYC and MCL1 at the standard `+2` preserves the
published axis on its published definition, so the contrast between the
replication arm and the novel arm stays legible.

This is a structural power argument plus an a-priori-focus argument, not a choice
fitted to an observed effect size. Recorded as such.

**Mandatory guard.** Odds ratios computed at different thresholds are **not
comparable across genes**. Script 05 must therefore report the **full gene x
threshold grid** (`{MYC, MCL1, BCL2L1} x {+2, >=+1}`), with the primary cell per
gene marked. Cross-gene comparison is licensed only within a threshold, never
across. Any figure obeys the same rule.

My original recommendation was `+2` primary with `>=+1` as sensitivity, on the
grounds that `+2` is the standard GISTIC definition of high-level amplification.
It was overridden with reasons. Both positions are on the record.

### 2.4 Aneuploidy conditioning - REQUIRED, not optional

A naive Fisher test on these calls will return a confident positive largely
because aneuploid tumours gain both arms. G2 must condition on aneuploidy burden
(`ANEUPLOIDY_SCORE`, `FRACTION_GENOME_ALTERED`, snapshotted in
`data/tcga_clinical/`) in addition to using ISAR. This matters *more* at the
`>= +1` threshold chosen for BCL2L1, where the base rate is 45.9%.

An unconditioned co-occurrence result is not reportable.

### 2.5 PAM50 and TNBC - both, separately

TNBC rule (fixed): ER and PR from `*_STATUS_BY_IHC`; HER2 from
`HER2_FISH_STATUS` where Positive/Negative, else `IHC_HER2`; **IHC-equivocal
without FISH is uncallable, not negative**. Gives 951 callable, 161 TNBC.

They are not interchangeable, contrary to the plan's "basal / TNBC" phrasing: 44
of 161 TNBC are non-Basal and 32 of 149 callable Basal are non-TNBC. Report both,
separately, plus the cross-tabulation. See `data/tcga_clinical/README.md`.

`BRCA_Normal` (n=36) is artefact-prone and usually reflects low cellularity.
Report as its own stratum; never allow it to become a silent reference level.

### 2.6 BBC3 - AMENDMENT to G2, tested in the LOSS direction

**This broadens a pre-specified gate and is recorded as an amendment, dated,
with its reason. It is legitimate because G2 is descriptive and no outcome data
has been seen.**

Reason for adding: `PRIME = log2(BBC3) - log2(BCL2L1)` has two terms, and testing
CNV determinants of only the denominator is asymmetric.

**Direction: loss, not amplification.** `BBC3` (PUMA) is the pro-apoptotic
numerator, so the `PRIME`-lowering event is deletion. `BCL2L1` and `MCL1` are
tested for gain; `BBC3` is tested for loss. The tests are not the same shape and
must not be described as "co-amplification" collectively.

Observed, ISAR, n = 1,043:

```
BBC3  19q13.32   -2:   5    -1: 242    0: 541   +1: 244   +2:  11
                 loss (<= -1): 247  (23.7%)
```

**Threshold for BBC3: `<= -1`.** Homozygous deletion is 5 tumours, and its
crude joint count with `MYC` amplification is **zero** - that test is not
underpowered, it is empty. `<= -1` gives 23.7%, which is testable and is the
direct symmetric counterpart of the `>= +1` threshold chosen for `BCL2L1`. The
`-2` counts are reported descriptively alongside.

**Regional control, REQUIRED.** BBC3's shallow-loss rate is close to background
for its neighbourhood and for other BCL2-family loci:

| Gene | Cytoband | loss `<= -1` |
|---|---|---|
| `BBC3` | 19q13.32 | 247 (23.7%) |
| `BAX` | 19q13.33 | 226 (21.7%) |
| `BCL2L11` | 2q13 | 241 (23.1%) |
| `BID` | 22q11.21 | 489 (46.9%) |

So a `MYC` x `BBC3`-loss association could be arm-level 19q behaviour rather than
PUMA-specific. ISAR plus aneuploidy conditioning is intended to absorb that, but
it must be tested rather than trusted: **`MYC` x `BAX`-loss at `<= -1` is run as
a regional control for 19q13.** If the `BBC3` estimate is indistinguishable from
the `BAX` estimate, the finding is regional and is reported as such. This is the
arm's standing "every positive needs its negatives" rule applied at CNV level.

Carry the BBC3 call to script 08 as a covariate, to show `log2(BBC3)` variation
is not CNV-driven.

**Correction, 2026-08-28.** This section previously read "homozygous BBC3
deletion occurs in 11 of 1,043 tumours" and concluded BBC3 should be descriptive
only. **11 is the amplification (`+2`) count; homozygous deletion is 5.** The
error was caught by the author asking whether the loss direction had been
considered. The direction had been analysed correctly throughout, but an
inconsistent threshold standard had been applied - `BCL2L1` was moved to `>= +1`
for interpretability while `BBC3` was left at `-2`, where it is empty. Correcting
the standard is what promotes BBC3 from descriptive to tested. Recorded rather
than amended away, because this is a pre-registration document.

### 2.7 Disclosure - threshold decision was not made blind

To make the threshold question concrete, crude unadjusted joint counts were
computed **before** the threshold was fixed:

```
[+2 only]  MYC x MCL1    both=53    OR=1.77
[+2 only]  MYC x BCL2L1  both=16    OR=5.14
[>= +1  ]  MYC x MCL1    both=515   OR=1.10
[>= +1  ]  MYC x BCL2L1  both=391   OR=4.38
```

These are unadjusted for aneuploidy and are **not** the gate result. The
threshold was then chosen with them visible. What preserves the record is that
the competing recommendations, and their reasons, were stated before the numbers
existed and are time-ordered in the session transcript. Recorded here rather than
omitted.

### 2.8 Pass criterion and the Panel b decision - FIXED

The plan gives G2 no numeric criterion, only "more than expected". Fixed here,
before the statistic was computed:

> **G2 passes for a partner if the aneuploidy-adjusted odds ratio exceeds 1 with
> a 95% CI excluding 1, at that gene's primary threshold and direction, in
> source A (ISAR).** Source B must agree in direction. A B disagreement is
> reported, never overridden.

The partners can split, and plausibly will, since `MCL1` is Lee et al.'s axis and
`BCL2L1` is ours. Plan section 6 assumes they move together. They may not, so:

**If the partners split, `BCL2L1` carries the Panel b decision.** It is the
`PRIME` denominator and the term Lee et al. never touch. `MCL1` is reported as
replication of Lee et al. either way.

`BBC3` loss is tested on the same criterion but is **not** sufficient on its own
to carry Panel b, because of the 19q13 regional-control caveat in 2.6.

### 2.9 Statistical design

For each pair, three estimates reported side by side, so the size of the
aneuploidy confound is visible rather than assumed:

| Estimate | Role |
|---|---|
| unadjusted Fisher, OR + 95% CI | shows the confound's magnitude |
| logistic `partner ~ MYC + ANEUPLOIDY_SCORE`, adjusted OR | **primary** |
| Cochran-Mantel-Haenszel, stratified by aneuploidy tertile | assumption-light check |

Breslow-Day across PAM50 strata, to test whether the association genuinely
differs by subtype rather than eyeballing five ORs.

Multiplicity: the primary tests are `MYC` x `MCL1` (+2), `MYC` x `BCL2L1` (>= +1)
and `MYC` x `BBC3` (<= -1), reported with exact p and no correction. The full
gene x threshold grid, the strata and the `BAX` regional control are secondary,
BH-adjusted within family and labelled exploratory in the output table.

`FRACTION_GENOME_ALTERED` is the pre-specified alternative to `ANEUPLOIDY_SCORE`;
the primary uses the latter and confirms with the former. n is reported per
model, since aneuploidy is missing for roughly 57 cases.

---

## 3. Still open

- ~~D4 framing consequence~~ - RESOLVED, Option 1. See 1.6.
- **D5** - primary neoadjuvant cohort. Now more urgent: section 1.5 makes H4's
  buffered/unbuffered split load-bearing against Lee et al., so the cohort must
  support it.
- **D7** - proliferation covariate overlaps the MYC exposure. PROPOSED. Confirm
  before script 09.
- **G1 correlation criterion** - needs script 01. G1 is not fully discharged.
- **Script 07 specificity panel source** - undecided, and now additionally needs
  an interpretive note on the FAO negative (section 1.2).
- **`TCGAbiolinks`** - not installed. Not needed for G2; needed for scripts
  01/02/03.
