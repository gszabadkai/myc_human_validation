---
date: 2026-08-30
status: RECORDED + DECLARATION - section 6 is a pre-data declaration
relates-to:
  - 2026-08-27_human_validation_plan.md (sections 7.5, 9, 11, 12; H4)
  - 2026-08-28_D5_cohort_selection.md (6.1 continuous primary, 6.2 meta-analysis,
    6.3 what is detectable)
  - 2026-08-29_block_b_result_escape_not_supported.md (BUFFER's role, E1 and E4)
  - 2026-08-30_block_g_result_and_bcl2l11_declaration.md (BUFFER's role, functionally)
decides:
  - STATE is FROZEN as plan section 7.5 specifies. Definition unamended.
  - The expression-tertile BUFFER fallback is KEPT, not tuned. kappa = 0.221.
  - H4's BUFFER is an EXPRESSION construct and is NOT the copy-number BUFFER
    that G2 passed on and that Blocks B and G failed to support.
  - The continuous BUFFER for H4's primary test is declared here, pre-data,
    with its direction fixed (section 6).
next-action: script 12 - fetch GSE194040, GSE164458, GSE25066. Nothing in this
  note may be revised after any of those outcome columns has been read.
---

# STATE is frozen, and H4's BUFFER is not the BUFFER that failed

Script 11 run 2026-08-30. Results in `results/state_definition.rds`.

**No outcome variable was read.** TCGA is not a neoadjuvant cohort and carries
no pCR endpoint, so nothing below could have been tuned to a result. The
declaration in section 6 is written after the exposure-side calibration in
section 3 and before any H4 outcome data exists. That ordering is the point of
the note and it is auditable from the git history.

---

## 1. STATE is frozen

Plan section 7.5's definition, transcribed without amendment, as a portable
constructor rather than a TCGA column - H4 is tested in three neoadjuvant
cohorts scored independently (D5 section 6.2), so the deliverable is a function
script 13 calls once per cohort.

Reference cells, BUFFER from GISTIC:

```
                                    MYC_low   OX_low   unbuffered   buffered   NA
Block C 938   GSVA                      460      170          129        179    0
Block C 938   mitoPPS                   460      186          110        182    0
full cohort   GSVA                      520      188          140        192   55
full cohort   mitoPPS                   520      203          122        195   55
```

Two properties asserted at runtime, both of which held:

- **The constructor closes over nothing.** `environment()` is `baseenv()` and
  `codetools::findGlobals` reports no free variables. A closure saved with
  `environment == globalenv()` is rebound to the *loading* session's globals on
  `readRDS`, which would let a cohort scored later silently use a different
  definition.
- **There is no argument for supplying a cut.** The function computes its own
  median from the vector it is handed. GSVA is cohort-relative (CLAUDE.md), so
  a TCGA threshold applied elsewhere is meaningless; now it cannot happen by
  accident.

The integrity contract returns `TRUE` on both checks
(`identical(deparse(...))`, `identical(environment(...), baseenv())`).

**Instrument agreement on cell assignment is 91.7%** (938). All 78 disagreements
sit at the level-2 boundary - whether OXPHOS is high - and **zero** patients
move between levels 3 and 4, which is as it must be, since BUFFER does not
depend on the instrument.

### 1.1 What was resolved, that section 7.5 left open

Seven choices, all recorded in `frozen$spec` and in the script header. Six are
forced; one (f) is a genuine new specification and is the subject of section 3.

| | Choice | Resolution |
|---|---|---|
| a | OXPHOS instrument | **both**, co-primary. Not a menu to pick from later. |
| b | MYC estimator | **M-a**. The only one of the three computable in an array cohort - M-c is a GISTIC call. |
| c | OXPHOS arm | **`OXPHOS subunits`**, not the umbrella (which carries assembly factors). |
| d | median of what | **this cohort**, at scoring time, on that model's complete cases. Never a TCGA threshold. |
| e | tie rule | **`> median` is high**, matching script 10 section 3. |
| f | BUFFER without CNV | **top tertile of MCL1 OR of BCL2L1**, within cohort, `> quantile(2/3, type = 7)`. BCL2 excluded (plan section 9). |
| g | missing data | complete cases (D8). NA in any constituent gives `STATE = NA`. |

## 2. The level-3-vs-4 contrast is frozen as specified, and its motivation is not

Kept exactly as written, direction included: pCR higher in level 3 than level 4.

It has nevertheless lost most of its motivation. It turns on BUFFER doing
something, and BUFFER's role is what Block B (E1) failed to support at the
quadrant level and what Block G failed to support functionally at n = 1,130.

**That belongs in the manuscript text, not in a revision.** Revising now would
be permissible on timing - no outcome data has been seen - but the *reason* to
revise would be Blocks B, C and G, which is precisely the "keep the cell that
survives" move the freeze exists to prevent. Recorded in
`frozen$spec$caveat_for_text`.

## 3. The calibration: the fallback does not reproduce the rule it replaces

The neoadjuvant cohorts have no CNV, so resolution (f) is what H4 will actually
use. TCGA is the only cohort in this arm carrying both, so it is the only place
the substitution can be checked - and it had to be checked now, because
revising the fallback is permissible before outcome data and not after.

```
n = 938

                 expr FALSE   expr TRUE
GISTIC FALSE          223         203      specificity  223/426 = 52.3%
GISTIC TRUE           156         356      sensitivity  356/512 = 69.5%

prevalence   GISTIC 54.6%   expression tertile 59.6%
concordance  61.7%          Cohen's kappa 0.221
```

Nearly half the tumours GISTIC calls unbuffered are called buffered by
expression, and 30% of the buffered ones are missed. Specificity is barely
above a coin flip. The disagreement is mildly one-sided - the fallback
over-calls.

**The matching prevalences are arithmetic, not validation.** GISTIC BUFFER is
54.6%; "top tertile of either of two genes" is `1 - (2/3)^2 = 55.6%` for any two
genes that are not strongly correlated, measured at 55.4% on random input. They
were printed only so the coincidence could not be mistaken for agreement. Kappa
is the number that carries information.

### 3.1 Why this is mechanically expected

Script 10's header already says what BUFFER is made of: it is
`cnv_MCL1 == +2 OR cnv_BCL2L1 >= +1`, and it is **dominated by a broad 20q gain**
(`buffer_BCL2L1` 439) rather than focal amplification (`buffer_MCL1` 158).

An arm-level gain moves any single gene's expression weakly. Low concordance is
therefore the predicted result of substituting expression for copy number here,
not the symptom of a badly chosen quantile.

### 3.2 The number that matters for the secondary contrast is 36.7%, not 88.0%

STATE agreement between the GISTIC-BUFFER and expression-BUFFER versions is
88.0%. That is the reassuring-looking number and it is the wrong one.

BUFFER cannot move a patient in or out of levels 1 or 2; it only decides level 3
versus level 4. So all `938 - 825 = 113` disagreements fall inside the 308
patients in levels 3+4. **36.7% of the patients the secondary contrast is made
of change sides** depending on which BUFFER definition is used.

### 3.3 DECISION: keep the fallback, do not tune it

Recorded as a decision, with reasons, so it is not revisited as a matter of
taste:

1. **No quantile fixes it.** The signal being asked for - broad arm-level gain
   raising one gene's expression into a top tertile - is genuinely weak
   (section 3.1). Threshold tuning would move the prevalence, not the kappa.
2. **It would optimise against one platform.** A cut chosen to maximise
   agreement with TCGA GISTIC on RNA-seq is then applied to two array cohorts.
3. **GISTIC is not ground truth for this quantity.** Amplification and
   expression are different constructs, and script 10 makes that distinction
   the entire methodological content of E2. Fitting one to the other treats a
   real difference as measurement error.

So the fallback stands as frozen in script 11, and the 0.221 is carried as a
measured property of the substitution rather than a defect to be repaired.

## 4. Consequence: H4's BUFFER is not the BUFFER that failed

This is the interpretive load the 0.221 carries, and it must be in the text.

- **G2** passed on copy number: MYC amplification co-occurs with MCL1 and
  BCL2L1 amplification (adj OR 1.69 and 1.77).
- **Block B (E1)** showed that co-occurrence is not conditional on OXPHOS -
  copy-number BUFFER again.
- **Block G** tested the functional consequence of MCL1 and BCL2L1 dependency.
- **H4, in the neoadjuvant cohorts, will test high MCL1 / BCL2L1 EXPRESSION.**

At kappa 0.221 these are not interchangeable labels for one variable.
**A null H4 does not re-test what Block B and Block G already found, and a
positive H4 would not rescue them.** The two sets of results are about
different constructs and must be reported as such.

One thing that follows in H4's favour, and is worth stating because it is not
obvious: since none of the three cohorts has CNV, all three use the *same*
expression construct. The meta-analysis is therefore internally uniform - it
does not mix a copy-number BUFFER in one cohort with an expression BUFFER in
another.

## 5. The secondary contrast is subtype-confounded, against the prediction

STATE is not a subtype variable in disguise - both cells are genuinely mixed:

```
GSVA, Block C 938        Basal   Her2   LumA   LumB   Normal      ER-  ER+
level 3  unbuffered         27     18     54     27        3       34   88   (27.9% ER-)
level 4  buffered           59     22     40     58        0       68   97   (41.2% ER-)
```

But level 4 - the cell predicted **resistant** - is the basal-enriched, more
ER-negative one, and basal tumours have markedly higher pCR under neoadjuvant
chemotherapy.

**Subtype therefore pushes level 4's pCR up, against H4's prediction.** Two
consequences, both fixed here rather than after the fact:

1. **PAM50 adjustment is mandatory**, not a robustness check, in every H4 model
   and in the secondary contrast.
2. **A positive result in the predicted direction would be conservative** -
   obtained despite the confounder pulling the other way. A null is
   correspondingly harder to read, and per D5 section 6.3 must be reported as
   "not powered to exclude a modest effect", never as evidence of absence.

---

# 6. DECLARATION - the continuous BUFFER for H4's primary test

**Everything in this section is fixed before any H4 outcome data exists.**

Plan section 7.5 defines a categorical BUFFER. D5 section 6.1 made the
**primary** H4 test the continuous three-way, for power. Neither says what
BUFFER becomes when it goes continuous, and nothing else in the plan resolves
it. This declaration fills that gap.

**This is not a revision of the freeze.** STATE stays exactly as script 11 froze
it, tertile BUFFER included, and remains the pre-specified secondary. Section 6
specifies a quantity the plan requires for the primary test and never defined.

### 6.1 Definition

```
BUFFER_c = mean( z(log2 MCL1), z(log2 BCL2L1) )
```

z-scored **within cohort**, on the complete cases entering that cohort's H4
model. Log scale. The two genes are MCL1 and BCL2L1 - the same pair as the
categorical rule. **BCL2 is excluded**, per plan section 9 mitigation 1: it is
estrogen-responsive and would make the measure partly an ER readout.

### 6.2 Why the mean and not the maximum

The categorical rule combines the two genes with an **OR**, which argues on its
face for `max(z1, z2)`. Rejected, before any fit:

- The OR is an artefact of **dichotomisation** - some rule is needed to combine
  two binary calls - not a claim that buffering capacity is a maximum.
- Biologically it is wrong. Two anti-apoptotic proteins both being high is more
  apoptotic buffering than one being high; `max()` cannot see that.
- `max()` of two z-scores is skewed and unstable, and it discards the second
  gene entirely.

The mean is the standard composite construction and it is the one declared.

### 6.3 Companion, mandatory

The **two genes are additionally reported separately**, each in its own
three-way, alongside the composite. Not as an alternative to select between
after the fact - as the same "report the components" discipline script 10
applies to `buffer_MCL1` and `buffer_BCL2L1`. If the composite fires and neither
component does, that is reported.

### 6.4 Scale and cohort rules

- **Within-cohort z, always.** GSVA is cohort-relative and expression scales
  differ by platform; a z computed in one cohort and applied in another is not a
  quantity. This is the same rule as resolution (d).
- **Never pooled.** Three cohorts, three estimates, meta-analysed (D5 section
  6.2 and CLAUDE.md). Scores are never concatenated across cohorts.
- Per-SD coefficients, so the three estimates are on a common scale for the
  meta-analysis.

### 6.5 Directions, fixed now

Primary model, per cohort, logistic on pCR, PAM50-adjusted (section 5):

```
pCR ~ MYC * OXPHOS * BUFFER_c + covariates
```

with MYC = M-a, OXPHOS = `OXPHOS subunits`, both z-scored within cohort, on
**both instruments**.

| Term | H4 predicts | Reasoning |
|---|---|---|
| `MYC:OXPHOS:BUFFER_c` | **NEGATIVE** | The `MYC x OXPHOS` effect on pCR is positive when buffering is low and is attenuated as buffering rises. That is "OXPHOS-high is protective unless buffered", written as a coefficient. |
| `MYC:OXPHOS` | **POSITIVE** | The two-way at `BUFFER_c = 0`, i.e. at mean buffering. |

**Claim only what both instruments support**, as everywhere else in this arm.
Report confidence intervals, not p-values alone (D5 section 6.3).

### 6.6 The informative failure, in coefficient terms

Plan section 2's fourth falsification criterion is "H4 fails in the informative
direction: OXPHOS-high predicts chemo*resistance* regardless of buffering."
Written so it cannot be reinterpreted after the numbers are in:

> The OXPHOS effect on pCR **within MYC-high** - the model's OXPHOS slope
> evaluated at `MYC = +1 SD` - is **negative with a CI excluding zero**, on both
> instruments, **while `MYC:OXPHOS:BUFFER_c` is null**.

That combination inverts the model rather than failing to support it, and per
plan section 15 it is reportable as a result rather than a non-result. It is
**not** the same as a null three-way, and a null three-way must not be written
up as if it were.

### 6.7 What this declaration does NOT do

- It does not revise STATE, the tertile BUFFER, or the level-3-vs-4 contrast.
- It does not specify the rest of the H4 model. Covariates beyond PAM50, the
  `score x treatment` interaction against I-SPY2's 179-patient paclitaxel
  control arm, the meta-analysis machinery and the handling of the three
  cohorts' platform differences are **script 13's build spec** and are not
  settled here.
- It does not add a hypothesis. H4 is one of the original four; this fixes how
  an already-declared hypothesis is measured.

## 7. Ordering, for the audit trail

1. Script 11 written and committed (`b533eee`) - the tertile fallback fixed
   before it was calibrated.
2. Script 11 run. Section 3's kappa observed. **Exposure side only; TCGA has no
   pCR endpoint.**
3. This note, including the section 6 declaration, written and committed.
4. Only then: script 12 fetches GSE194040, GSE164458, GSE25066.

The kappa result motivated *specifying* the continuous form now rather than
discovering the gap inside script 13. It carries no information about any
outcome, so it cannot have biased the direction fixed in section 6.5.
