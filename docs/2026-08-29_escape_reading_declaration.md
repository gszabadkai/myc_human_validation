---
date: 2026-08-29
status: DECLARED - directions fixed before any of these models were fitted
relates-to:
  - 2026-08-27_human_validation_plan.md (sections 1, 2 H1 clause 2, 7.5 H4, 10 Block B)
  - 2026-08-28_G2_result.md
  - 2026-08-29_block_c_result_H1_not_supported.md
decides:
  - Nothing new is hypothesised. This fixes the DIRECTION of tests already
    pre-registered, plus one declared addition (E2), before script 10 exists.
next-action: scripts/10_buffering_models.R
---

# The escape reading - directions declared before Block B

## 0. What this is, and what it is not

**It is not a fifth hypothesis.** Plan section 2 forbids one and the prohibition
is honoured here. Every test below is either (a) already pre-registered, or (b)
one declared addition, E2, whose direction is fixed in this document before any
model exists.

**It is not a reversal of Block C.** H1 is a conjunction and its first clause
failed. Nothing here rescues H1. What follows answers a *different and larger*
question, which the plan states as the reason the arm exists at all.

## 1. The organising question was always this

Plan section 1, written before any data was touched:

> *"If OXPHOS-high plus MYC-high is apoptotically lethal in untransformed
> mammary epithelium, **how did the human tumours that occupy that state survive
> it**, and what did that escape cost them therapeutically?"*

That is the escape reading, verbatim, as the arm's organising question. The
`MYC:OXPHOS` interaction tested in Block C was **one clause of one hypothesis**,
not the whole arm.

Why the framing suits the data: TCGA contains only survivors of the bottleneck
the mouse describes. A *causal* claim about what happened during transformation
is unavailable here - plan section 1 says so and calls it a non-goal. A
*constraint* claim is not: "you only find the dangerous state where the danger
has been neutralised" is a statement about survivors, and a survivor cohort is
the right design for it.

## 2. Where it already lives in the pre-registration

| Element | Status |
|---|---|
| **H1 clause 2** - "the MYC-high / OXPHOS-high quadrant is enriched for `BUFFER`" | pre-registered, **NOT RUN** (Block B) |
| **H4** - "OXPHOS-high is protective *unless buffered*" | pre-registered, **NOT RUN** (Block F) |
| **G2** - MCL1 / BCL2L1 amplification co-occurs with MYC amplification | **PASSED**: adj OR 1.69 (1.17-2.44) and 1.77 (1.26-2.47) |

Two of the three are untested. That is the whole reason this can be pursued
without inventing anything.

## 3. What carries directional weight, and what does not

**Reversing a regression establishes nothing.** `OXPHOS ~ priming` and
`priming ~ OXPHOS` on expression data carry identical information and identical
confounding. If the escape reading were only an algebraic rearrangement it would
add exactly zero, and that must be said plainly rather than glossed.

**`BUFFER` is not an expression level. It is DNA.** Script 03 defines it as
`cnv_MCL1 == +2` OR `cnv_BCL2L1 >= +1` - somatic copy number. OXPHOS expression
cannot cause an MCL1 amplification, so a model with `BUFFER_amp` as predictor is
a regression on a genuinely upstream variable, not a reversed one. **That
asymmetry is the entire methodological content of this document.**

**The composition of `BUFFER` is load-bearing and is reported with every result.**
Within the 938: `BUFFER` TRUE 512 / FALSE 426, of which `buffer_BCL2L1` = 439 and
`buffer_MCL1` = 158. It is dominated by a **broad 20q gain**, not by focal
amplification. So:

- **aneuploidy adjustment is mandatory**, not optional (G2 built and validated it);
- **MYC amplification must be adjusted**, because G2's own positive result is that
  the two co-occur;
- **the expression-matched null in E2 is load-bearing**, because a buffer-to-OXPHOS
  association could be aneuploidy wearing a mechanism's clothes. If `BUFFER`
  predicts every gene set equally, the null distribution shifts with it and
  OXPHOS will not stand out. That is what the null is for here.

## 4. The other half of the idea is not testable, and G2 already settled it

The escape could be buffer amplification **or** loss of the BH3-only trigger.
The second has no genetic instrument in TCGA:

```
BBC3 loss   adjusted OR 1.36   CI 0.966-1.91   FAIL
            indistinguishable from its BAX 19q13 regional control (1.41)
            homozygous deletion n = 5; joint count with MYC amplification = 0
```

So "reduced BH3" can only be examined at expression level, where it is symmetric
with OXPHOS and carries no directional weight whatsoever. **The escape reading is
declared for the buffer side only.** Whether `BCL2L11` copy loss is recurrent
enough to serve as an instrument has not been checked and is not promised here.

## 5. The tests, with directions fixed

### E1 - H1 clause 2, as pre-registered. PRIMARY.

```
BUFFER_amp ~ MYC_high * OXPHOS_high + PAM50 + purity + aneuploidy
             logistic, n = 938, both co-primary instruments
```

**PREDICTED: the `MYC_high:OXPHOS_high` interaction is POSITIVE.** The dangerous
quadrant carries more buffer than its two main effects predict.

Splits at the median of the frozen z-scores on the 938. `plate` is deliberately
**excluded**: the outcome is a copy-number call from SNP arrays and the RNA-seq
plate is not a confounder of it.

Two companions, reported alongside, never instead:
- the **continuous** form `BUFFER_amp ~ MYC * OXPHOS`, same predicted direction.
  D5 established that median splits cost up to 53 percentage points of power, so
  the dichotomised form is the weaker instrument even though it is the one the
  plan wrote.
- the **descriptive** quadrant-versus-rest contrast, which is what the H1 clause
  says in words.

### E2 - the specificity form. DECLARED ADDITION.

```
arm_score ~ BUFFER_amp + purity + leukocyte_fraction + PROLIF_DISJOINT
          + PAM50 + TP53_status + plate_pooled + aneuploidy
          all 18 arms, both instruments, n = 938
```

**PREDICTED: POSITIVE for `OXPHOS subunits`, and beating its 2,000-set
expression-matched null; NOT positive for the comparator arms - in particular
NOT for `OXPHOS assembly factors`.**

This is the declared addition, and it is the test that turns "buffered tumours
tolerate more OXPHOS" from a correlation into a specificity claim. It reuses the
null already cached in `results/mito_null/`; the design matrix is fixed and only
the outcome changes, so all 68,000 fits are one QR decomposition.

`plate` **is** included here, because the outcome is an RNA-derived score.

### E3 - Block B's expression models. ONE IS ALREADY ANSWERED.

```
log2(MCL1)   ~ MYC * OXPHOS + <Block C covariates>     NOT YET FITTED
log2(BCL2L1) ~ MYC * OXPHOS + <Block C covariates>     ALREADY FITTED
```

**`log2(BCL2L1)` has already been run**, as script 09's `limb: log2_BCL2L1`. It
gave **-0.034 (p 0.030)** on GSVA and **-0.043 (p 0.009)** on mitoPPS. **No
prediction is declared for it, because the result is known.** Declaring one now
would be writing a prediction backwards, which is the exact failure this whole
document exists to avoid.

What that known result means for the escape reading, stated now: it **runs
against** the naive expression-level version of it. `MYC x OXPHOS` associates
with *less* BCL-XL transcript, not more.

**PREDICTED for `log2(MCL1)`, which has not been fitted: POSITIVE.** Stated with
low confidence and with its rationale exposed - `buffer_MCL1` is a focal
high-level amplification (`== +2`, n = 158) where `buffer_BCL2L1` is a broad gain,
so if any expression-level buffering survives adjustment it should be MCL1's.

**If MCL1 is also null or negative, the escape reading is a DNA-level claim
only.** That is a coherent and reportable outcome - gene dosage can set a floor
that transcriptional regulation modulates around - but it must be stated as a
narrowing, not presented as if the expression level had supported it.

### E4 - the one prediction that ties Block C to the escape reading

```
log2(BCL2L1) ~ MYC * OXPHOS * BUFFER_amp + <Block C covariates>
```

**PREDICTED: the three-way term is POSITIVE** - i.e. the `MYC:OXPHOS`
down-regulation of BCL-XL found in Block C is **attenuated in BUFFER-amplified
tumours**, because amplification sets a floor the transcriptional effect has less
room beneath.

**Power caveat, stated before the fit:** a three-way term on a 512/426 split is
underpowered by construction. **A null here is uninformative and will be reported
as such, not as evidence against.** Said now so it cannot later sound like an
excuse.

### E5 - H4. Not in script 10.

Already pre-specified in plan section 7.5 as amended by D5. It is the sharpest
form of the escape reading, because a treatment-response endpoint is measured
*after* the tumour state and so breaks the symmetry cross-sectional data can
never break. It belongs to Block F and its own stop gate.

## 6. Rules that bind this document

- **Nothing here changes H1's status.** H1 is a conjunction whose first clause
  failed. Block B answers the arm's organising question; it does not resurrect a
  failed conjunction, and no text may imply otherwise.
- **Script 09 discarded the OXPHOS and MYC main effects.** They were fitted and
  thrown away; only `MYC:OX` was kept. Those main effects are relevant to a
  constraint reading, and retrieving them is a one-line change. **They must not
  be inspected before this note is committed**, and when retrieved they inform
  E2's interpretation without redefining it.
- **Every result reports `BUFFER`'s composition alongside it** (512 of 938, of
  which 439 BCL2L1 gain and 158 MCL1 amplification), so no reader mistakes a
  broad 20q event for focal buffering.
- **No further variants.** If E1 and E2 fail, the escape reading is reported as
  unsupported and the arm rests on G2's DNA-level co-occurrence alone.

## 7. Failure conditions, fixed now

| If | Then |
|---|---|
| E1's interaction is null | H1 clause 2 fails. The quadrant carries no excess buffer and the escape reading has no quadrant evidence. |
| E2 is positive for OXPHOS subunits **and** equally for the comparator arms | the association is compartment-wide, not an OXPHOS-specific constraint. Report and stop. |
| E2's OXPHOS coefficient does not beat its matched null | it is a set-size or expression-level artefact. Report and stop. |
| E2 is positive for `OXPHOS assembly factors` as strongly as for the subunits | the specificity claim fails exactly as it did in Block C, and for the same reason. |
| E3 MCL1 is null | the escape reading narrows to DNA only. Reportable; not to be dressed up. |
| E4 is null | uninformative, per the power caveat above. Not evidence against. |

If E1 and E2 both fail, that is the answer, and plan section 15 already describes
how to report it.
